"""Projection selection and small end-to-end tests using existing material data."""

import importlib.util
import itertools
import shlex
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "wannier_multistart", ROOT / "utility" / "wannier_multistart.py")
ms = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = ms
spec.loader.exec_module(ms)


def test_selection_repairs_local_rank_and_preserves_gauge():
    values = np.array([[[100, 1, 0], [0, 0, 1]],
                       [[0, 1, 0], [100, 0, 1]]], dtype=complex)
    mask = np.ones((2, 2), dtype=bool)
    columns, minimum = ms.select_columns(values, 2, mask)
    np.testing.assert_array_equal(columns, [1, 2])
    assert minimum == pytest.approx(1)
    rotations = np.array([[[1, 1j], [1j, 1]], [[1j, 1], [-1, -1j]]]) / np.sqrt(2)
    rotated, _ = ms.select_columns(rotations @ values, 2, mask)
    np.testing.assert_array_equal(rotated, columns)
    with pytest.raises(ValueError, match="full-rank"):
        ms.select_columns(np.zeros_like(values), 2, mask)


def test_window_selection_and_unmodified_amn(tmp_path):
    win = tmp_path / "model.win"
    win.write_text("num_wann=2\ndis_win_max=2\ndis_froz_max=2\n")
    eig = tmp_path / "model.eig"
    eig.write_text("1 1 0\n2 1 1\n3 1 3\n")
    allowed = ms.window_mask(ms.WinInput.read(win), 1, 3, 2, eig)
    values = np.array([[[1, 0, 0], [0, 1, 0], [4, 5, 1000]]], dtype=complex)
    columns, _ = ms.select_columns(values, 2, allowed)
    np.testing.assert_array_equal(columns, [0, 1])
    bank = tmp_path / "selected.amn"
    ms.write_amn(bank, values[:, :, columns])
    np.testing.assert_array_equal(ms.read_amn(bank), values[:, :, columns])
    win.write_text("num_wann=2\ndis_win_max=0\n")
    with pytest.raises(ValueError, match="too few states"):
        ms.window_mask(ms.WinInput.read(win), 1, 3, 2, eig)


@pytest.mark.parametrize("records,reason", [
    ("1 1 1 1D0 0\n1 1 1 1 0\n", "Duplicate"),
    ("", "Missing"),
    ("1 1 1 NaN 0\n", "Non-finite"),
    ("2 1 1 0 0\n", "out of range"),
])
def test_invalid_amn(tmp_path, records, reason):
    path = tmp_path / "bad.amn"
    path.write_text("test\n1 1 1\n" + records)
    with pytest.raises(ValueError, match=reason):
        ms.read_amn(path)


def test_only_converged_finite_results_are_accepted():
    final = ("Final State\n WF centre and spread 1 (0, 0, 0) 1.0\n"
             "Sum of centres and spreads (0, 0, 0) 1.0\n")
    success = "<<< Wannierisation convergence criteria satisfied >>>\n"
    assert ms.localization_result(final, 1, False) == (False, 1.0)
    assert ms.localization_result(success + final, 1, False) == (True, 1.0)
    assert not ms.localization_result(success + final, 1, True)[0]
    assert not ms.localization_result(success + final.replace("0, 0, 0", "NaN, 0, 0"), 1, False)[0]
    assert not ms.localization_result(success + final.replace("spread 1", "spread 2"), 1, False)[0]


def test_prepare_preserves_sources_and_rejects_old_output(tmp_path):
    win = tmp_path / "model.win"
    win.write_text("num_wann=1\nbegin projections\nH:s\nend projections\n")
    (tmp_path / "model.mmn").write_text("test\n1 1 1\n")
    bank = tmp_path / "bank.amn"
    ms.write_amn(bank, np.ones((1, 1, 2), dtype=complex))
    before = {path: path.read_bytes() for path in (win, bank)}
    output = tmp_path / "trials"
    result = ms.run_trials(win, [bank], output, [], prepare_only=True)
    assert result["winner"] is None and len(result["trials"]) == 2
    assert "begin projections" not in (output / "trial-001/model.win").read_text()
    assert all(path.read_bytes() == data for path, data in before.items())
    with pytest.raises(FileExistsError):
        ms.run_trials(win, [bank], output, [], prepare_only=True)
    win.write_text("num_wann=1\nread_ibz=true\n")
    with pytest.raises(ValueError, match="read_ibz"):
        ms.run_trials(win, [bank], tmp_path / "ibz", [], prepare_only=True)


def test_timeout_stops_command(tmp_path):
    assert ms.execute([sys.executable, "-c", "import time; time.sleep(30)"],
                      tmp_path, 0.05) is None


@pytest.mark.parametrize("dataset,seed", [
    ("testw90_example01", "gaas"),
    ("testw90_guidingcentre_selectproj", "silicon"),
])
def test_native_selection_matches_compact_export(tmp_path, executables, run_options, dataset, seed):
    """Compare to native select_projections, including a frozen-window case."""
    source = ROOT / "test-suite/tests" / dataset
    original = ms.WinInput.read(source / f"{seed}.win")
    win = tmp_path / f"{seed}.win"
    without_selection = ms.WinInput(
        [line for i, line in enumerate(original.lines)
         if original.scalar_lines.get(i) != "select_projections"], {}, {}, {})
    win.write_text(without_selection.render({}))
    # Disable metadata-dependent guiding centres for the matrix-only comparison.
    parsed = ms.WinInput.read(win)
    win.write_text(parsed.render({"guiding_centres": "false", "conv_tol": "1e-8",
                                 "conv_window": 5}))
    for suffix in ("amn", "mmn", "eig"):
        path = source / f"{seed}.{suffix}"
        if path.exists():
            ms.link_input(path, tmp_path / path.name)
    snapshots = {path: path.read_bytes() for path in (win, source / f"{seed}.amn")}
    command = [str(executables.require("wannier90"))]
    if run_options.nprocs:
        command = shlex.split(run_options.launcher.format(nprocs=run_options.nprocs)) + command
    result = ms.run_trials(win, [source / f"{seed}.amn"], tmp_path / "trials", command,
                           max_iter=400, cg_steps=[5, 0])
    assert result["winner"] is not None, result
    winner = next(t for t in result["trials"] if t["directory"] == result["winner"])
    assert winner["spread"] == min(t["spread"] for t in result["trials"] if t["converged"])
    selected = tmp_path / "trials" / result["winner"]
    control = tmp_path / "control"
    control.mkdir()
    settings = {"select_projections": " ".join(map(str, winner["columns"])),
                "guiding_centres": "false", "conv_tol": "1e-8", "conv_window": 5,
                "num_iter": 400, "num_cg_steps": winner["num_cg_steps"], "iprint": 1}
    # Keep the original projections block and full bank for the independent run.
    lines = [line for i, line in enumerate(original.lines)
             if original.scalar_lines.get(i) not in settings]
    lines.extend(f"{key}={value}" for key, value in settings.items())
    (control / win.name).write_text("\n".join(lines) + "\n")
    for suffix in ("amn", "mmn", "eig"):
        path = source / f"{seed}.{suffix}"
        if path.exists():
            ms.link_input(path, control / path.name)
    assert ms.execute([*command, seed], control, 120) == 0
    nw = original.integer("num_wann", 0)
    reference = (control / f"{seed}.wout").read_text()
    converged, spread = ms.localization_result(reference, nw, seed == "silicon")
    assert converged
    assert winner["spread"] == pytest.approx(spread, abs=2e-7)
    # The centres must agree up to a Wannier-function permutation, not just the sum.
    def centres(text):
        final = text.rsplit("Final State", 1)[1]
        return np.array([[float(x) for x in c.split(",")] for c in
                         ms.re.findall(r"WF centre and spread\s+\d+\s+\(([^)]*)\)", final)])
    a, b = centres((selected / f"{seed}.wout").read_text()), centres(reference)
    assert min(np.max(np.abs(a - b[list(p)])) for p in itertools.permutations(range(nw))) < 1e-6
    if seed == "gaas":
        unfinished = ms.run_trials(win, [source / f"{seed}.amn"], tmp_path / "unfinished",
                                   command, max_iter=1, cg_steps=[5])
        assert unfinished["winner"] is None
        assert not unfinished["trials"][0]["converged"]
    assert all(path.read_bytes() == data for path, data in snapshots.items())
