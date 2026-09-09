"""Running the program under test.

Two things here matter more than they look:

* the command is built as an argv **list**, never a shell string, so filenames containing
  spaces or shell metacharacters cannot be misinterpreted;
* a timed-out run is killed by **process group**, not just the direct child.  Three tests
  abort on purpose, and a hung ``mpirun`` otherwise leaves orphaned ranks behind that go on
  to poison later tests.
"""

import os
import shlex
import signal
import subprocess
from dataclasses import dataclass
from pathlib import Path

#: Template for launching a parallel run.  ``{nprocs}`` is substituted; the rest is split
#: with shlex, so Intel MPI, srun and mpiexec can all be configured from the command line.
DEFAULT_MPI_LAUNCHER = "mpirun -np {nprocs}"

STDOUT_LOG = "stdout.log"
STDERR_LOG = "stderr.log"


class RunnerError(Exception):
    """Raised when a run could not be started at all."""


@dataclass
class RunResult:
    """Outcome of one invocation."""

    argv: list[str]
    returncode: int
    work_dir: Path
    timed_out: bool = False
    timeout: int | None = None

    @property
    def stdout_log(self) -> Path:
        return self.work_dir / STDOUT_LOG

    @property
    def stderr_log(self) -> Path:
        return self.work_dir / STDERR_LOG

    @property
    def command(self) -> str:
        """The command as you would type it, for error messages."""
        return " ".join(shlex.quote(a) for a in self.argv)

    def stderr_tail(self, lines: int = 20) -> str:
        """Last few lines of stderr, for failure messages."""
        try:
            content = self.stderr_log.read_text(errors="replace").strip()
        except OSError:
            return ""
        if not content:
            return ""
        tail = content.splitlines()[-lines:]
        return "\n".join(tail)


def build_argv(
    exe: Path,
    input_file: str,
    args=(),
    *,
    nprocs: int | None = None,
    launcher: str = DEFAULT_MPI_LAUNCHER,
) -> list[str]:
    """Assemble the command line.

    ``nprocs`` of ``None`` or 0 means a serial run with no launcher at all -- not
    ``mpirun -np 1``, which is a different thing and is what the old default category did.
    """
    argv: list[str] = []
    if nprocs:
        argv.extend(shlex.split(launcher.format(nprocs=nprocs)))
    argv.append(str(exe))
    argv.extend(str(a) for a in args)
    argv.append(str(input_file))
    return argv


def _kill_process_group(process: subprocess.Popen) -> None:
    """Kill the whole process group, then reap the child."""
    try:
        os.killpg(os.getpgid(process.pid), signal.SIGKILL)
    except (ProcessLookupError, PermissionError):
        # Already gone, or we cannot signal the group; fall back to the direct child.
        process.kill()
    try:
        process.wait(timeout=10)
    except subprocess.TimeoutExpired:
        pass


def run(
    argv: list[str],
    work_dir: Path,
    *,
    timeout: int | None = None,
    env: dict[str, str] | None = None,
) -> RunResult:
    """Run ``argv`` in ``work_dir``, capturing stdout and stderr to files."""
    if not work_dir.is_dir():
        raise RunnerError(f"work directory does not exist: {work_dir}")

    stdout_path = work_dir / STDOUT_LOG
    stderr_path = work_dir / STDERR_LOG
    timed_out = False

    with open(stdout_path, "wb") as out, open(stderr_path, "wb") as err:
        try:
            process = subprocess.Popen(
                argv, cwd=work_dir, stdout=out, stderr=err,
                env=env,
                # Put the child in its own session so we can signal the entire group,
                # including MPI ranks it spawned, if it hangs.
                start_new_session=True,
            )
        except FileNotFoundError as exc:
            raise RunnerError(f"cannot execute {argv[0]!r}: {exc}") from exc
        except OSError as exc:
            raise RunnerError(f"failed to start {argv[0]!r}: {exc}") from exc

        try:
            returncode = process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            _kill_process_group(process)
            timed_out = True
            returncode = -signal.SIGKILL

    return RunResult(
        argv=argv, returncode=returncode, work_dir=work_dir,
        timed_out=timed_out, timeout=timeout,
    )
