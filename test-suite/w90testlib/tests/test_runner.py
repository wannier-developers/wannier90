"""Unit tests for command construction and subprocess handling."""

import sys

import pytest

from w90testlib.runner import RunnerError, build_argv, run


def test_serial_runs_have_no_launcher_prefix():
    argv = build_argv("/bin/w90.x", "gaas.win", nprocs=None)
    assert argv == ["/bin/w90.x", "gaas.win"]


def test_nprocs_zero_is_serial_not_mpirun_np_1():
    """The old 'default' category ran the bare binary; -np 1 is a different thing."""
    assert build_argv("/bin/w90.x", "gaas.win", nprocs=0) == ["/bin/w90.x", "gaas.win"]


def test_parallel_runs_are_prefixed_with_the_launcher():
    argv = build_argv("/bin/w90.x", "gaas.win", nprocs=2)
    assert argv == ["mpirun", "-np", "2", "/bin/w90.x", "gaas.win"]


def test_the_launcher_is_configurable():
    argv = build_argv("/bin/w90.x", "gaas.win", nprocs=4, launcher="srun -n {nprocs} --exclusive")
    assert argv == ["srun", "-n", "4", "--exclusive", "/bin/w90.x", "gaas.win"]


def test_args_go_between_the_executable_and_the_input():
    assert build_argv("/bin/w90.x", "wannier.win", ["-pp"]) == \
        ["/bin/w90.x", "-pp", "wannier.win"]


def test_the_command_is_a_list_so_metacharacters_are_never_interpreted():
    argv = build_argv("/bin/w90.x", "odd name; rm -rf /.win")
    assert argv[-1] == "odd name; rm -rf /.win"


# --- actually running things ------------------------------------------------------------

def test_a_successful_run_captures_stdout_and_stderr(tmp_path):
    result = run([sys.executable, "-c",
                  "import sys; print('out'); print('err', file=sys.stderr)"], tmp_path)
    assert result.returncode == 0
    assert not result.timed_out
    assert result.stdout_log.read_text().strip() == "out"
    assert result.stderr_log.read_text().strip() == "err"


def test_a_failing_run_reports_its_exit_code(tmp_path):
    result = run([sys.executable, "-c", "raise SystemExit(3)"], tmp_path)
    assert result.returncode == 3


def test_stderr_tail_is_available_for_failure_messages(tmp_path):
    result = run([sys.executable, "-c",
                  "import sys; print('boom', file=sys.stderr)"], tmp_path)
    assert "boom" in result.stderr_tail()


def test_a_hanging_run_is_timed_out_and_killed(tmp_path):
    result = run([sys.executable, "-c", "import time; time.sleep(30)"], tmp_path, timeout=1)
    assert result.timed_out
    assert result.returncode != 0


def test_a_timeout_kills_the_whole_process_group_not_just_the_child(tmp_path):
    """A hung mpirun otherwise leaves orphaned ranks that poison later tests."""
    marker = tmp_path / "grandchild_alive"
    script = (
        "import subprocess, sys, time\n"
        f"subprocess.Popen([sys.executable, '-c', "
        f"\"import time,pathlib; time.sleep(20); pathlib.Path(r'{marker}').write_text('x')\"])\n"
        "time.sleep(20)\n"
    )
    result = run([sys.executable, "-c", script], tmp_path, timeout=1)
    assert result.timed_out
    # If the group was killed, the grandchild never gets to write its marker.
    import time
    time.sleep(3)
    assert not marker.exists(), "a grandchild process survived the timeout kill"


def test_a_missing_executable_is_reported_clearly(tmp_path):
    with pytest.raises(RunnerError, match="cannot execute"):
        run([str(tmp_path / "does_not_exist")], tmp_path)


def test_running_in_a_missing_directory_is_reported_clearly(tmp_path):
    with pytest.raises(RunnerError, match="work directory does not exist"):
        run([sys.executable, "-c", ""], tmp_path / "nope")
