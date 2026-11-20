#!/usr/bin/env python3
# Copyright (C) 2026 Thomas D. Kuehne
# SPDX-License-Identifier: LGPL-2.1-or-later
#
# This file is part of Wannier90 and is distributed under the GNU Lesser
# General Public License, version 2.1 or (at your option) any later version.
# It is provided without any warranty. See the LICENSE file for details.

"""Select mesh-wide projection columns and compare independent Wannier90 runs.

This optional utility uses NumPy and SciPy. It does not evaluate orbitals or
replace Wannier90's disentanglement and spread minimization.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shlex
import shutil
import signal
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

try:
    import numpy as np
    from scipy.linalg import qr
except ImportError as error:
    raise SystemExit("Install the optional utility dependencies with: pip install numpy scipy") from error


def number(text: str) -> float:
    value = float(text.replace("D", "e").replace("d", "e"))
    if not math.isfinite(value):
        raise ValueError("Non-finite numerical input")
    return value


@dataclass
class WinInput:
    """Preserve input blocks while reading and replacing scalar settings."""

    lines: list[str]
    settings: dict[str, str]
    scalar_lines: dict[int, str]
    blocks: dict[str, tuple[int, int]]

    @classmethod
    def read(cls, path: Path) -> WinInput:
        lines = path.read_text().splitlines()
        settings, scalar_lines, blocks = {}, {}, {}
        current = None
        start = 0
        for index, raw in enumerate(lines):
            line = re.split(r"[!#]", raw, maxsplit=1)[0].strip()
            match = re.fullmatch(r"(begin|end)\s*:?\s+(\w+)", line, re.IGNORECASE)
            if match:
                action, name = (part.lower() for part in match.groups())
                if action == "begin":
                    if current is not None or name in blocks:
                        raise ValueError(f"Nested or repeated input block {name}")
                    current, start = name, index
                else:
                    if current != name:
                        raise ValueError(f"Unmatched input block {name}")
                    blocks[name] = (start, index)
                    current = None
            elif line and current is None:
                match = re.fullmatch(r"(\w+)\s*(?:[:=]\s*|\s+)(.+)", line)
                if match is None:
                    raise ValueError(f"Cannot parse input line {index + 1}")
                key, value = match.groups()
                key = key.lower()
                if key in settings:
                    raise ValueError(f"Repeated input keyword {key}")
                settings[key] = value.strip()
                scalar_lines[index] = key
        if current is not None:
            raise ValueError(f"Unclosed input block {current}")
        return cls(lines, settings, scalar_lines, blocks)

    def real(self, key: str, default: float) -> float:
        return number(self.settings[key]) if key in self.settings else default

    def integer(self, key: str, default: int) -> int:
        return int(self.settings[key]) if key in self.settings else default

    def enabled(self, key: str) -> bool:
        return self.settings.get(key, "false").lower().strip(".") in ("true", "t")

    def render(self, updates: dict[str, object]) -> str:
        # The supplied AMN is authoritative. A projection block would otherwise
        # make Wannier90 expect the original (possibly larger) projector count.
        removed = set()
        if "projections" in self.blocks:
            first, last = self.blocks["projections"]
            removed.update(range(first, last + 1))
        result = [line for i, line in enumerate(self.lines)
                  if i not in removed and self.scalar_lines.get(i) not in updates]
        result.extend(f"{key} = {value}" for key, value in updates.items())
        return "\n".join(result) + "\n"


def read_amn(path: Path) -> np.ndarray:
    """Read band/projector/k-point records, rejecting missing or duplicate data."""
    with path.open() as stream:
        next(stream)
        nb, nk, nc = map(int, next(stream).split())
        if min(nb, nk, nc) < 1:
            raise ValueError(f"Invalid AMN dimensions in {path}")
        values = np.empty((nk, nb, nc), dtype=complex)
        seen = np.zeros(values.shape, dtype=bool)
        for line in stream:
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) != 5:
                raise ValueError(f"Invalid AMN record in {path}")
            band, column, point = (int(x) - 1 for x in fields[:3])
            if not (0 <= band < nb and 0 <= point < nk and 0 <= column < nc):
                raise ValueError(f"AMN index out of range in {path}")
            if seen[point, band, column]:
                raise ValueError(f"Duplicate AMN record in {path}")
            values[point, band, column] = complex(number(fields[3]), number(fields[4]))
            seen[point, band, column] = True
        if not seen.all():
            raise ValueError(f"Missing AMN records in {path}")
    return values


def write_amn(path: Path, values: np.ndarray) -> None:
    nk, nb, nc = values.shape
    with path.open("x") as stream:
        stream.write("Selected projection columns\n")
        stream.write(f"{nb} {nk} {nc}\n")
        for point in range(nk):
            for column in range(nc):
                for band in range(nb):
                    value = values[point, band, column]
                    stream.write(f"{band + 1} {column + 1} {point + 1} "
                                 f"{value.real:.17e} {value.imag:.17e}\n")


def read_eig(path: Path, nk: int, nb: int) -> np.ndarray:
    values = np.empty((nk, nb))
    seen = np.zeros((nk, nb), dtype=bool)
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        band, point, value = line.split()
        band, point = int(band) - 1, int(point) - 1
        if not (0 <= band < nb and 0 <= point < nk) or seen[point, band]:
            raise ValueError(f"Invalid or duplicate eigenvalue index in {path}")
        values[point, band] = number(value)
        seen[point, band] = True
    if not seen.all():
        raise ValueError(f"Missing eigenvalues in {path}")
    return values


def window_mask(win: WinInput, nk: int, nb: int, nw: int, eig: Path) -> np.ndarray:
    """Use the outer window, or a completely filled frozen target space."""
    if nw > nb or nw < 1:
        raise ValueError("num_wann must lie between one and the AMN band count")
    if not any(key.startswith(("dis_win_", "dis_froz_")) for key in win.settings):
        return np.ones((nk, nb), dtype=bool)
    energies = read_eig(eig, nk, nb)
    low, high = win.real("dis_win_min", -math.inf), win.real("dis_win_max", math.inf)
    allowed = (energies >= low) & (energies <= high)
    if low > high or np.any(allowed.sum(axis=1) < nw):
        raise ValueError("The outer window contains too few states")
    if "dis_froz_min" in win.settings and "dis_froz_max" not in win.settings:
        raise ValueError("dis_froz_min requires dis_froz_max")
    if "dis_froz_max" in win.settings:
        fmin = win.real("dis_froz_min", low)
        fmax = win.real("dis_froz_max", high)
        if fmin < low or fmax > high or fmin > fmax:
            raise ValueError("The frozen window must lie inside the outer window")
        frozen = (energies >= fmin) & (energies <= fmax)
        counts = frozen.sum(axis=1)
        if np.any(counts > nw):
            raise ValueError("The frozen window contains more than num_wann states")
        allowed[counts == nw] = frozen[counts == nw]
    return allowed


def projection_quality(values: np.ndarray, tolerance: float) -> tuple[int, float]:
    singular = np.linalg.svd(values, compute_uv=False)
    scales = singular[:, :1]
    relative = np.divide(singular, scales, out=np.zeros_like(singular), where=scales > 0)
    defect = int(np.count_nonzero(relative <= tolerance))
    return defect, float(relative[:, -1].min())


def select_columns(values: np.ndarray, nw: int, allowed: np.ndarray,
                   tolerance: float = 1e-10) -> tuple[np.ndarray, float]:
    """Select one common trial set, checking its rank at every k-point.

    QR acts on the stacked mesh, not independently on each point. Rows are
    masked for selection only; the exported AMN entries remain unmodified.
    Rank repair accepts exchanges only if the total rank defect decreases.
    """
    nk, nb, nc = values.shape
    if not 0 < tolerance < 1 or min(nk, nb, nc) < 1 or nw < 1 or nw > min(nb, nc):
        raise ValueError("Invalid projection count or rank tolerance")
    if allowed.shape != (nk, nb) or not np.isfinite(values).all():
        raise ValueError("Invalid candidate bank or window mask")
    masked = np.where(allowed[:, :, None], values, 0)
    best, best_defect = None, nk * nw + 1
    # A single anchor can cure a mesh-wide QR selection that loses local rank.
    for anchor in [masked.reshape(nk * nb, nc), *masked]:
        _, order = qr(anchor, mode="r", pivoting=True, check_finite=False)
        columns = np.sort(order[:nw])
        defect, minimum = projection_quality(masked[:, :, columns], tolerance)
        if defect < best_defect:
            best, best_defect = columns.copy(), defect
        if defect == 0:
            return columns, minimum
    assert best is not None
    columns = best
    while best_defect > 0:
        improved = False
        for position in range(nw):
            for candidate in range(nc):
                if candidate in columns:
                    continue
                trial = columns.copy()
                trial[position] = candidate
                defect, minimum = projection_quality(masked[:, :, trial], tolerance)
                if defect < best_defect:
                    columns, best_defect, improved = trial, defect, True
                    break
            if improved:
                break
        if not improved:
            raise ValueError(f"No common full-rank trial set found ({best_defect} missing directions)")
    return np.sort(columns), minimum


def localization_result(text: str, nw: int, disentangle: bool) -> tuple[bool, float | None]:
    """Require the native convergence messages, not merely a Final State block."""
    converged = "Wannierisation convergence criteria satisfied" in text
    if disentangle:
        converged = converged and "Disentanglement convergence criteria satisfied" in text
        converged = converged and "Disentanglement convergence criteria not satisfied" not in text
    final = text.rsplit("Final State", 1)
    if len(final) != 2:
        return False, None
    states = re.findall(r"WF centre and spread\s+(\d+)\s+\(([^)]*)\)\s+(\S+)", final[1])
    sums = re.findall(r"Sum of centres and spreads\s+\([^)]*\)\s+(\S+)", final[1])
    if sorted(int(state[0]) for state in states) != list(range(1, nw + 1)) or not sums:
        return False, None
    try:
        values = [number(state[2]) for state in states]
        for _, centre, _ in states:
            coords = centre.split(",")
            if len(coords) != 3:
                return False, None
            for coordinate in coords:
                number(coordinate)
        total = number(sums[0])
    except ValueError:
        return False, None
    valid = min(values) >= -1e-7 and abs(sum(values) - total) <= 2e-7 * max(1, nw)
    return converged and valid, total if valid else None


def link_input(source: Path, destination: Path) -> None:
    # Symlinks avoid copying the often large MMN and optional plotting grids.
    destination.symlink_to(source.resolve())


def execute(command: list[str], directory: Path, timeout: float) -> int | None:
    """Stop the entire MPI process group on timeout, not only the launcher."""
    with (
        (directory / "run.log").open("x") as log,
        subprocess.Popen(command, cwd=directory, stdout=log, stderr=subprocess.STDOUT,
                         start_new_session=os.name == "posix") as process,
    ):
        try:
            return process.wait(timeout=timeout)
        except (subprocess.TimeoutExpired, KeyboardInterrupt) as error:
            if os.name == "posix":
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
            else:
                process.kill()
            process.wait()
            if isinstance(error, KeyboardInterrupt):
                raise
            return None


def run_trials(win_path: Path, banks: list[Path], output: Path, command: list[str],
               max_iter: int = 1000, cg_steps: list[int] | None = None,
               tolerance: float = 1e-10, timeout: float = 600,
               prepare_only: bool = False) -> dict:
    win_path = win_path.resolve()
    win = WinInput.read(win_path)
    seed = win_path.stem
    if not re.fullmatch(r"[\w][\w.-]*", seed):
        raise ValueError("Use a seed name containing letters, digits, underscores, dots or hyphens")
    for key in ("restart", "exclude_bands", "select_projections"):
        if key in win.settings:
            raise ValueError(f"{key} is not supported for fresh matrix-only trials")
    for key in ("guiding_centres", "site_symmetry", "slwf_constrain", "use_bloch_phases",
                "postproc_setup", "read_ibz"):
        if win.enabled(key):
            raise ValueError(f"{key} is not supported for full-mesh matrix-only trials")
    if win.integer("slwf_num", 0) > 0:
        raise ValueError("Selective localization uses a different objective")
    if max_iter < 1 or timeout <= 0 or not math.isfinite(timeout) or not 0 < tolerance < 1:
        raise ValueError("Invalid iteration limit, timeout or rank tolerance")
    schedules = [win.integer("num_cg_steps", 5), 0] if cg_steps is None else cg_steps
    schedules = list(dict.fromkeys(schedules))
    if not schedules or min(schedules) < 0:
        raise ValueError("CG step counts must be non-negative")
    if not banks:
        raise ValueError("Supply at least one projection bank")
    arrays = [read_amn(path) for path in banks]
    nk, nb, _ = arrays[0].shape
    if any(array.shape[:2] != (nk, nb) for array in arrays):
        raise ValueError("All projection banks must share the same bands and ordered k-points")
    if win.integer("num_bands", nb) != nb:
        raise ValueError("num_bands does not match the supplied projection bank")
    nw = win.integer("num_wann", 0)
    mmn, eig = win_path.with_suffix(".mmn"), win_path.with_suffix(".eig")
    with mmn.open() as stream:
        next(stream)
        mb, mk, neighbours = map(int, next(stream).split())
    if (mb, mk) != (nb, nk) or neighbours < 1:
        raise ValueError("MMN dimensions do not match the projection banks")
    allowed = window_mask(win, nk, nb, nw, eig)
    if nb > nw or eig.exists():
        read_eig(eig, nk, nb)
    if not prepare_only:
        if not command or shutil.which(command[0]) is None:
            raise ValueError("Wannier90 command not found")
        command = [str(Path(shutil.which(command[0])).resolve()), *command[1:]]
    # Never reuse a work directory or an old checkpoint from a previous trial.
    output.mkdir(parents=True, exist_ok=False)
    summary = {"winner": None, "trials": [], "skipped": [], "prepared_only": prepare_only}
    for bank, values in zip(banks, arrays):
        try:
            columns, minimum = select_columns(values, nw, allowed, tolerance)
        except ValueError as error:
            summary["skipped"].append({"bank": str(bank.resolve()), "reason": str(error)})
            continue
        for schedule in schedules:
            directory = output / f"trial-{len(summary['trials']) + 1:03d}"
            directory.mkdir()
            settings = {"num_iter": max_iter, "num_cg_steps": schedule, "iprint": 1,
                        "conv_window": win.integer("conv_window", 5)
                        if win.integer("conv_window", 5) > 1 else 5,
                        "auto_projections": "false", "write_u_matrices": "true"}
            (directory / win_path.name).write_text(win.render(settings))
            write_amn(directory / f"{seed}.amn", values[:, :, columns])
            link_input(mmn, directory / mmn.name)
            if eig.exists():
                link_input(eig, directory / eig.name)
            if win.enabled("wannier_plot"):
                for grid in win_path.parent.glob("UNK*.*"):
                    if grid.is_file():
                        link_input(grid, directory / grid.name)
            trial = {"directory": directory.name, "bank": str(bank.resolve()),
                     "columns": (columns + 1).tolist(), "min_relative_singular_value": minimum,
                     "num_cg_steps": schedule, "converged": False, "spread": None}
            if not prepare_only:
                trial["returncode"] = execute([*command, seed], directory, timeout)
                if trial["returncode"] is None:
                    trial["reason"] = "Timed out"
                report = directory / f"{seed}.wout"
                if trial["returncode"] == 0 and report.exists():
                    trial["converged"], trial["spread"] = localization_result(
                        report.read_text(), nw, nb > nw)
                    # A winning result must retain its optimized gauge, not only AMN.
                    trial["converged"] = trial["converged"] and all(
                        (directory / name).is_file() for name in
                        (f"{seed}.chk", f"{seed}_u.mat"))
                    if nb > nw:
                        trial["converged"] = trial["converged"] and (
                            directory / f"{seed}_u_dis.mat").is_file()
            summary["trials"].append(trial)
    successful = [trial for trial in summary["trials"] if trial["converged"]]
    if successful:
        summary["winner"] = min(successful, key=lambda trial: trial["spread"])["directory"]
    (output / "selection.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("win", type=Path, help="Base input with full-mesh MMN and, if needed, EIG")
    parser.add_argument("--amn", action="append", required=True, type=Path,
                        help="Projection bank in the SAME Bloch gauge as MMN; repeat for alternatives")
    parser.add_argument("--output", required=True, type=Path, help="New directory for independent trials")
    parser.add_argument("--command", default="wannier90.x", help="Command, e.g. 'mpiexec -n 2 wannier90.x'")
    parser.add_argument("--max-iter", type=int, default=1000)
    parser.add_argument("--cg-steps", type=int, nargs="+", help="CG schedules (default: input value and zero)")
    parser.add_argument("--rank-tol", type=float, default=1e-10)
    parser.add_argument("--timeout", type=float, default=600, help="Seconds per trial")
    parser.add_argument("--prepare-only", action="store_true", help="Write inputs without running Wannier90")
    args = parser.parse_args(argv)
    try:
        summary = run_trials(args.win, args.amn, args.output, shlex.split(args.command),
                             args.max_iter, args.cg_steps, args.rank_tol, args.timeout, args.prepare_only)
    except (OSError, ValueError, StopIteration) as error:
        parser.error(str(error))
    if summary["winner"] is not None:
        print(f"Selected converged result: {args.output / summary['winner']}")
        return 0
    if args.prepare_only and summary["trials"]:
        print(f"Prepared trials: {args.output}")
        return 0
    print(f"No admissible converged result; inspect {args.output / 'selection.json'}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
