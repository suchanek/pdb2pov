"""The command line, including the 1993 grammar it has to keep honouring."""

from __future__ import annotations

import os

import pytest

from pypdb2pov.cli import EXIT_NO_ATOMS, EXIT_NO_BONDS, EXIT_PARSE_ARGS, main


def run(capsys, *argv) -> tuple[int, str, str]:
    status = main(list(argv))
    captured = capsys.readouterr()
    return status, captured.out, captured.err


def test_the_historical_invocation_still_works(capsys, crambin_pdb, tmp_path):
    stem = crambin_pdb[: -len(".pdb")]
    out = str(tmp_path / "crambin")

    status, stdout, _ = run(capsys, stem, out, "-s", "-h", "-b", "-d", "1.5", "-x", "90")

    assert status == 0
    assert os.path.exists(out + ".pov")
    assert "found" in stdout and "bonds" in stdout


def test_dash_h_is_the_checkered_ground_not_help(capsys, crambin_pdb, tmp_path):
    """
    -h has meant the checkerboard since 1993.  If argparse ever reclaims it
    for help, every existing script silently stops producing a ground plane.
    """
    stem = crambin_pdb[: -len(".pdb")]
    out = str(tmp_path / "crambin")

    status, _, _ = run(capsys, stem, out, "-h")

    assert status == 0
    assert "checker color Black color White" in open(out + ".pov").read()


def test_object_only_writes_an_inc(capsys, crambin_pdb, tmp_path):
    stem = crambin_pdb[: -len(".pdb")]
    out = str(tmp_path / "crambin")

    status, _, _ = run(capsys, stem, out, "-o", "-b", "-d", "1.9")

    assert status == 0
    assert os.path.exists(out + ".inc")
    assert not os.path.exists(out + ".pov")


def test_an_output_path_with_an_extension_is_taken_as_given(capsys, crambin_pdb, tmp_path):
    stem = crambin_pdb[: -len(".pdb")]
    out = str(tmp_path / "scene.pov")

    status, _, _ = run(capsys, stem, out)

    assert status == 0
    assert os.path.exists(out)
    assert not os.path.exists(out + ".pov")
    assert "#declare scene = object { scene_obj }" in open(out).read()


def test_writing_to_stdout(capsys, data_dir):
    status, stdout, stderr = run(capsys, os.path.join(data_dir, "altloc.pdb"), "-", "-o")

    assert status == 0
    assert stdout.startswith("//\n")
    assert "#declare" in stdout
    # Progress must not contaminate the scene.
    assert "Scanning atom file" in stderr


def test_quiet_says_nothing_on_success(capsys, data_dir, tmp_path):
    status, stdout, stderr = run(
        capsys, os.path.join(data_dir, "altloc.pdb"), str(tmp_path / "x"), "--quiet"
    )
    assert status == 0
    assert stdout == ""
    assert stderr == ""


def test_a_missing_file_is_reported(capsys, tmp_path):
    status, _, stderr = run(capsys, str(tmp_path / "nope"), str(tmp_path / "out"))

    assert status != 0
    assert "nope.pdb" in stderr


def test_no_atoms_after_a_chain_filter_says_so(capsys, data_dir, tmp_path):
    status, _, stderr = run(
        capsys,
        os.path.join(data_dir, "models.pdb"),
        str(tmp_path / "out"),
        "--chain",
        "Z",
    )

    assert status == EXIT_NO_ATOMS
    assert "chain" in stderr


def test_no_bonds_suggests_a_larger_threshold(capsys, data_dir, tmp_path):
    status, _, stderr = run(
        capsys,
        os.path.join(data_dir, "models.pdb"),
        str(tmp_path / "out"),
        "-b",
        "-d",
        "0.1",
    )

    assert status == EXIT_NO_BONDS
    assert "-d" in stderr


def test_an_output_file_is_required(capsys, data_dir):
    status, _, stderr = run(capsys, os.path.join(data_dir, "altloc.pdb"))
    assert status == EXIT_PARSE_ARGS
    assert "output" in stderr


def test_info_reports_without_writing(capsys, data_dir, tmp_path):
    status, stdout, _ = run(capsys, os.path.join(data_dir, "mini.cif"), "--info")

    assert status == 0
    assert "format:   cif" in stdout
    assert "atoms:    4" in stdout
    assert "models:   [1, 2]" in stdout
    assert not list(tmp_path.iterdir())


def test_include_dir_points_at_the_bundled_files(capsys):
    status, stdout, _ = run(capsys, "--include-dir")

    assert status == 0
    directory = stdout.strip()
    assert os.path.exists(os.path.join(directory, "atoms2.inc"))
    assert os.path.exists(os.path.join(directory, "atoms_vdw.inc"))


def test_list_elements(capsys):
    status, stdout, _ = run(capsys, "--list-elements")

    assert status == 0
    assert "ZN  Atom_Zn" in stdout
    assert "33 elements" in stdout


def test_the_scene_identifier_comes_from_the_output_name(capsys, data_dir, tmp_path):
    out = str(tmp_path / "1crn")
    run(capsys, os.path.join(data_dir, "altloc.pdb"), out)

    text = open(out + ".pov").read()
    assert "#declare _1crn_obj = union {" in text  # a digit cannot lead


def test_an_explicit_name_overrides_it(capsys, data_dir, tmp_path):
    out = str(tmp_path / "1crn")
    run(capsys, os.path.join(data_dir, "altloc.pdb"), out, "--name", "widget")

    assert "#declare widget_obj = union {" in open(out + ".pov").read()


def test_no_timestamp_makes_runs_reproducible(capsys, data_dir, tmp_path):
    first, second = str(tmp_path / "a"), str(tmp_path / "b")
    run(capsys, os.path.join(data_dir, "altloc.pdb"), first, "--no-timestamp", "--name", "m")
    run(capsys, os.path.join(data_dir, "altloc.pdb"), second, "--no-timestamp", "--name", "m")

    assert open(first + ".pov").read() == open(second + ".pov").read()


def test_covalent_bonding_is_selectable(capsys, crambin_pdb, tmp_path):
    stem = crambin_pdb[: -len(".pdb")]
    distance, covalent = str(tmp_path / "d"), str(tmp_path / "c")

    run(capsys, stem, distance, "-b", "-d", "2.2")
    run(capsys, stem, covalent, "-b", "--bonds", "covalent")

    # 2.2 A bonds every carbon to its neighbour's neighbour; radii do not.
    assert open(covalent + ".pov").read().count("cylinder {") < open(
        distance + ".pov"
    ).read().count("cylinder {")


@pytest.mark.parametrize("flag", ["--keep-altlocs", "--altloc=all"])
def test_altloc_flags_agree(capsys, data_dir, tmp_path, flag):
    out = str(tmp_path / "out")
    run(capsys, os.path.join(data_dir, "altloc.pdb"), out, flag, "--name", "m")

    assert open(out + ".pov").read().count("object { Atom_") == 7
