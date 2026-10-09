"""Command line behaviour shared by all LayNii programs."""
import pytest

import laynii_testing as lt


def _param(program):
    marks = [pytest.mark.wip] if program in lt.WIP_PROGRAMS else []
    return pytest.param(program, marks=marks, id=program)


ALL = [_param(p) for p in lt.ALL_PROGRAMS]

# LN_PHYSIO_PARS uses positional arguments instead of -flags
FLAG_PROGRAMS = [p for p in ALL if p.values[0] != "LN_PHYSIO_PARS"]

# These programs print an error for unknown options but still exit with 0.
# The test is strict: once a program is fixed, remove it from this list.
EXIT_ZERO_ON_UNKNOWN_OPTION = {
    "LN2_DEVEIN", "LN2_LAYER_SMOOTH", "LN2_RIMIFY", "LN2_RIM_POLISH",
    "LN_3DCOLUMNS", "LN_COLUMNAR_DIST", "LN_DIRECT_SMOOTH", "LN_FLOAT_ME",
    "LN_GRADSMOOTH", "LN_IMAGIRO", "LN_INT_ME", "LN_LAYER_SMOOTH",
    "LN_LOITUMA", "LN_MP2RAGE_DNOISE", "LN_RAGRUG", "LN_SHORT_ME",
    "LN_TEMPSMOOTH",
}

# Name of the option that takes the main input file
MAIN_INPUT_OPTION = {
    "LN2_LAYERS": "-rim", "LN2_COLUMNS": "-rim", "LN2_MULTILATERATE": "-rim",
    "LN2_RIM_BORDERIZE": "-rim", "LN2_RIM_POLISH": "-rim", "LN3_LAYERS": "-rim",
    "LN_GROW_LAYERS": "-rim", "LN_LEAKY_LAYERS": "-rim",
    "LN2_PATCH_FLATTEN": "-values", "LN2_PATCH_FLATTEN_2D": "-values",
    "LN2_PATCH_UNFLATTEN": "-values", "LN2_UVD_FILTER": "-values",
    "LN2_UVD_LSTSQR": "-values", "LN2_LAYERDIMENSION": "-values",
    "LN2_PEAK_DETECT": "-values", "LN2_ZERO_CROSSING": "-values",
    "LN2_GEODISTANCE": "-domain", "LN2_IFPOINTS": "-domain",
    "LN2_VORONOI": "-domain", "LN2_HEXBIN": "-coord_uv",
    "LN2_CHOLMO": "-layers", "LN_3DCOLUMNS": "-layers",
    "LN_COLUMNAR_DIST": "-layers", "LN_IMAGIRO": "-layers",
    "LN_CONLAY": "-layers", "LN2_MASK": "-scores", "LN_INTPRO": "-image",
    "LN_BOCO": "-Nulled", "LN_MP2RAGE_DNOISE": "-INV1",
    "LN_CORREL2FILES": "-file1", "LN_LOITUMA": "-equidist",
    "LN2_REGRESS_OUT": "-input1",
}


def test_makefile_lists_programs():
    assert len(lt.RELEASED_PROGRAMS) > 50
    for program in lt.ALL_PROGRAMS:
        assert (lt.REPO_ROOT / "src" / (program + ".cpp")).exists(), program


@pytest.mark.parametrize("program", ALL)
def test_binary_was_built(program):
    assert lt.exe_path(program).exists(), (
        "%s is listed in the Makefile but was not built" % program)


@pytest.mark.parametrize("program", ALL)
def test_no_arguments_prints_help(program, tmp_path):
    result = lt.run_program(program, [], tmp_path, check=False, timeout=30)
    assert "Usage" in result.output, result.describe()
    # LN_PHYSIO_PARS treats missing arguments as an error
    expected = 1 if program == "LN_PHYSIO_PARS" else 0
    assert result.returncode == expected, result.describe()
    assert result.new_files == []


@pytest.mark.parametrize("program", FLAG_PROGRAMS)
def test_help_flag(program, tmp_path):
    result = lt.run_program(program, ["-help"], tmp_path, check=False,
                            timeout=30)
    assert result.returncode == 0, result.describe()
    assert "Usage" in result.output, result.describe()
    assert result.new_files == []


@pytest.mark.parametrize("program", FLAG_PROGRAMS)
def test_missing_input_file_fails(program, tmp_path):
    option = MAIN_INPUT_OPTION.get(program, "-input")
    result = lt.run_program(program, [option, "does_not_exist.nii.gz"],
                            tmp_path, check=False, timeout=30)
    assert result.returncode != 0, result.describe()
    # A negative exit code means the program was killed by a signal (crash)
    assert result.returncode > 0, "crashed\n" + result.describe()
    assert result.new_files == []


@pytest.mark.parametrize("program", FLAG_PROGRAMS)
def test_unknown_option_fails(program, tmp_path, request):
    if program in EXIT_ZERO_ON_UNKNOWN_OPTION:
        request.applymarker(pytest.mark.xfail(
            reason="exits with 0 on unknown options", strict=True))
    result = lt.run_program(program, ["-this_option_does_not_exist"],
                            tmp_path, check=False, timeout=30)
    assert result.returncode > 0, result.describe()
    assert result.new_files == []


@pytest.mark.xfail(reason="ifstream::bad() does not detect a missing file and "
                   "the read loop `while (input != 5003 || eof())` never ends",
                   strict=True)
def test_physio_pars_missing_file(tmp_path):
    result = lt.run_program("LN_PHYSIO_PARS", ["nope.puls", "out.txt"],
                            tmp_path, check=False, timeout=10)
    assert result.returncode != 0, result.describe()
