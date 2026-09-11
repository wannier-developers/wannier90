"""Unit tests for the declarative prepare steps that replaced the per-test Makefiles."""

import bz2

import pytest

from w90testlib.prepare import PrepareError, run_prepare


def make_bz2(directory, name, content=b"payload\n"):
    path = directory / name
    path.write_bytes(bz2.compress(content))
    return path


def test_bunzip2_decompresses_and_drops_the_suffix(tmp_path):
    make_bz2(tmp_path, "Fe.mmn.bz2", b"mmn data\n")
    produced = run_prepare([{"bunzip2": "*.mmn.bz2"}], tmp_path)
    assert produced == [tmp_path / "Fe.mmn"]
    assert (tmp_path / "Fe.mmn").read_bytes() == b"mmn data\n"


def test_bunzip2_keeps_the_compressed_original(tmp_path):
    """The .bz2 is what is committed; decompressing must not consume it."""
    make_bz2(tmp_path, "Fe.mmn.bz2")
    run_prepare([{"bunzip2": "*.mmn.bz2"}], tmp_path)
    assert (tmp_path / "Fe.mmn.bz2").is_file()


def test_bunzip2_handles_every_compressed_family_used_by_the_suite(tmp_path):
    for ext in ("mmn", "uHu", "amn", "spn", "sHu", "sIu"):
        make_bz2(tmp_path, f"Seed.{ext}.bz2", f"{ext}\n".encode())
    steps = [{"bunzip2": f"*.{ext}.bz2"} for ext in ("mmn", "uHu", "amn", "spn", "sHu", "sIu")]
    produced = run_prepare(steps, tmp_path)
    assert len(produced) == 6
    assert (tmp_path / "Seed.uHu").read_bytes() == b"uHu\n"


def test_a_pattern_matching_nothing_is_not_an_error(tmp_path):
    """Tests share a prepare list shape; absent optional inputs are normal."""
    assert run_prepare([{"bunzip2": "*.spn.bz2"}], tmp_path) == []


def test_an_unknown_step_is_rejected(tmp_path):
    with pytest.raises(PrepareError, match="unknown prepare step 'tar'"):
        run_prepare([{"tar": "*.tar"}], tmp_path)


def test_a_malformed_step_is_rejected(tmp_path):
    with pytest.raises(PrepareError, match="single-key mapping"):
        run_prepare([{"bunzip2": "a", "chk_from_bz2": "b"}], tmp_path)


def test_chk_from_bz2_without_the_converter_explains_itself(tmp_path):
    make_bz2(tmp_path, "Fe.chk.fmt.bz2")
    with pytest.raises(PrepareError, match="needs w90chk2chk"):
        run_prepare([{"chk_from_bz2": "*.chk.fmt.bz2"}], tmp_path, w90chk2chk=None)


# --- the one inter-test dependency ------------------------------------------------------

def test_copy_from_dependency_copies_the_named_artefacts(tmp_path):
    produced_dir = tmp_path / "checkpoint0_write"
    produced_dir.mkdir()
    (produced_dir / "copper.chk").write_bytes(b"binary chk\n")
    work = tmp_path / "checkpoint1_read"
    work.mkdir()

    run_prepare(
        [{"copy_from_dependency": {"test": "checkpoint0_write", "files": ["copper.chk"]}}],
        work,
        resolve_dependency=lambda name: tmp_path / name,
    )
    assert (work / "copper.chk").read_bytes() == b"binary chk\n"


def test_copy_from_dependency_fails_clearly_when_the_artefact_is_absent(tmp_path):
    (tmp_path / "dep").mkdir()
    work = tmp_path / "work"
    work.mkdir()
    with pytest.raises(PrepareError, match="did not produce"):
        run_prepare(
            [{"copy_from_dependency": {"test": "dep", "files": ["copper.chk"]}}],
            work,
            resolve_dependency=lambda name: tmp_path / name,
        )


def test_copy_from_dependency_without_a_resolver_explains_itself(tmp_path):
    with pytest.raises(PrepareError, match="cannot run standalone"):
        run_prepare(
            [{"copy_from_dependency": {"test": "dep", "files": ["a"]}}],
            tmp_path,
            resolve_dependency=None,
        )


# --- against the real w90chk2chk.x -------------------------------------------------------

REPO_ROOT = __import__("pathlib").Path(__file__).resolve().parents[3]
W90CHK2CHK = REPO_ROOT / "w90chk2chk.x"
SAMPLE_CHK = REPO_ROOT / "test-suite/tests/testpostw90_fe_ahc/Fe.chk.fmt.bz2"


@pytest.mark.skipif(not W90CHK2CHK.is_file(), reason="w90chk2chk.x is not built")
@pytest.mark.skipif(not SAMPLE_CHK.is_file(), reason="sample checkpoint not present")
def test_chk_from_bz2_produces_a_binary_checkpoint(tmp_path):
    """The full decompress -> convert -> discard-intermediate cycle, on a real checkpoint."""
    import shutil
    shutil.copy2(SAMPLE_CHK, tmp_path / "Fe.chk.fmt.bz2")

    produced = run_prepare([{"chk_from_bz2": "*.chk.fmt.bz2"}], tmp_path, w90chk2chk=W90CHK2CHK)

    assert produced == [tmp_path / "Fe.chk"]
    assert (tmp_path / "Fe.chk").stat().st_size > 0
    # The formatted intermediate is deleted; the compressed original is kept.
    assert not (tmp_path / "Fe.chk.fmt").exists()
    assert (tmp_path / "Fe.chk.fmt.bz2").is_file()
