from io import StringIO
from os import path

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose

from open_petro_elastic.__main__ import (
    CalibrationError,
    InputError,
    calculate_results,
    calibration_error,
    make_input,
    parse_arguments,
    run,
    run_calibration,
)

EXAMPLE_DIR = path.join(path.dirname(__file__), "..", "examples")


def run_example(example_yaml, example_data_file, example_calibration, shape=100):
    with open(example_yaml, encoding="utf-8") as input_file:
        inp = make_input(input_file, example_data_file)
    results = calculate_results(inp)
    path_to_calibration_data = path.join(EXAMPLE_DIR, example_calibration)
    calibration_data = pd.read_csv(path_to_calibration_data)
    for key in ["kmin", "kdry", "ksat", "kmin_fls"]:
        result = results[key]
        calibration = calibration_data[key]
        error = calibration_error(result, calibration)
        assert_allclose(np.zeros(shape), error, atol=2e-4)


def example_yaml(f):
    return path.join(EXAMPLE_DIR, f + ".yaml")


def example_csv(f):
    return path.join(EXAMPLE_DIR, f + ".csv")


run_example_parameters = (
    "yaml_file, csv_file, calibration_file, shape",
    [
        [
            example_yaml("example1"),
            example_csv("example1"),
            example_csv("example1_calibration"),
            (),
        ],
        [
            example_yaml("example2"),
            example_csv("example2"),
            example_csv("example2_calibration"),
            (),
        ],
        [
            example_yaml("example3"),
            example_csv("example3"),
            example_csv("example3_calibration"),
            (),
        ],
        [
            example_yaml("example4"),
            example_csv("example4"),
            example_csv("example4_calibration"),
            (),
        ],
        [
            example_yaml("friable_sand_example"),
            None,
            example_csv("friable_sand_calibration"),
            1,
        ],
        [
            example_yaml("patchy_cement_example"),
            None,
            "patchy_cement_calibration.csv",
            1,
        ],
        [
            example_yaml("carbon_dioxide_example"),
            example_csv("carbon_dioxide_example"),
            example_csv("carbon_dioxide_calibration"),
            (),
        ],
    ],
)


@pytest.mark.parametrize(*run_example_parameters)
def test_run_example(yaml_file, csv_file, calibration_file, shape):
    run_example(yaml_file, csv_file, calibration_file, shape)


@pytest.mark.parametrize(*run_example_parameters)
def test_success_threshold_equivalent(
    yaml_file, csv_file, calibration_file, shape, tmpdir
):
    if csv_file is None:
        return
    out1 = path.join(tmpdir, "out1.csv")
    out2 = path.join(tmpdir, "out2.csv")
    code = run(MockArgs(config_file=yaml_file, data_file=csv_file, output_file=out1))
    assert code == 0
    code = run(
        MockArgs(
            config_file=yaml_file,
            data_file=csv_file,
            output_file=out2,
            success_threshold=99,
        )
    )
    assert code == 0
    assert_allclose(pd.read_csv(out1), pd.read_csv(out2))


@pytest.mark.parametrize(*run_example_parameters)
def test_success_no_data_error(yaml_file, csv_file, calibration_file, shape, capsys):
    if csv_file is None:
        return
        code = run(
            MockArgs(
                config_file=yaml_file,
                data_file=None,
                success_threshold=99,
            )
        )
        assert code == -5
        captured = capsys.readouterr()
        assert "data must be given in csv" in captured.out


@pytest.fixture
def failing_data():
    "This data has 50% failing saturations"
    return StringIO(
        """
,brine_saturation,gas_saturation,oil_saturation
0,0.9,0.1,0.0
1,0.5,0.4,0.2
2,0.5,0.4,0.1
3,0.7,0.3,0.1
"""
    )


@pytest.fixture
def failing_config():
    return StringIO(
        """
fluids:
    constituents:
        - material:
            type: "oil"
            gas_gravity: 0.5 # ratio to air (air gas gravity is 1)
            reference_density: 800.0 # kg/m3
            gas_oil_ratio: 153.0 # ratio of gas to oil
          fraction: # Ratio of brine to remaining fluid
              column: "oil_saturation"
        - material:
            type: "brine"
            salinity: 45000 # ppm
          fraction: # Ratio of brine to remaining fluid
              column: "brine_saturation"
        - material:
            type: "gas"
            gas_gravity: 0.5 # ratio to air (air gas gravity is 1)
          fraction: # Ratio of gas to remaining fluid
              column: "gas_saturation"
"""
    )


def test_success_threshold_code(failing_config, failing_data):
    assert (
        run(
            MockArgs(
                config_file=failing_config,
                data_file=failing_data,
                success_threshold=49,
            )
        )
        == 0
    )


def test_success_threshold_fail_code(failing_config, failing_data):
    assert (
        run(
            MockArgs(
                config_file=failing_config,
                data_file=failing_data,
                success_threshold=51,
            )
        )
        == -5
    )


def test_make_input_malformed_yaml():
    malformed_yaml = StringIO("'\n")
    good_data = StringIO("depth,value\n1,1")
    with pytest.raises(InputError, match="parsing the config file"):
        make_input(malformed_yaml, good_data)


def test_make_input_bad_format_yaml():
    bad_format_yaml = StringIO("[1,2,3]\n")
    good_data = StringIO("depth,value\n1,1")
    with pytest.raises(InputError, match="format of config"):
        make_input(bad_format_yaml, good_data)


def test_make_input_empty_data():
    good_yaml = StringIO("dry_rock:\nminerals:\nfluids:\n\n")
    empty_data = StringIO("\n\n\n")
    with pytest.raises(InputError, match="empty"):
        make_input(good_yaml, empty_data)


def test_make_input_bad_data():
    good_yaml = StringIO("dry_rock:\nminerals:\nfluids:\n\n")
    bad_data = StringIO('"')
    with pytest.raises(InputError, match="read"):
        make_input(good_yaml, bad_data)


def test_run_calibration_empty_csv():
    empty_data = StringIO("\n\n\n")
    with pytest.raises(CalibrationError, match="empty"):
        run_calibration(empty_data, None)


def test_run_calibration_bad_data():
    bad_data = StringIO('"')
    with pytest.raises(CalibrationError, match="read"):
        run_calibration(bad_data, None)


def test_run_calibration_missing_columns():
    missing_data = StringIO("a,b,c\n1,2,3\n")
    results = pd.DataFrame(
        {
            "ksat": [],
            "kmin": [],
            "kdry": [],
            "mysat": [],
            "rsat": [],
            "kmin_fls": [],
        },
    )
    with pytest.raises(CalibrationError, match="getting column"):
        run_calibration(missing_data, results)


class MockArgs:
    def __init__(
        self,
        config_file=path.join(EXAMPLE_DIR, "example1.yaml"),
        data_file=path.join(EXAMPLE_DIR, "example1.csv"),
        calibration_data=path.join(EXAMPLE_DIR, "example1_calibration.csv"),
        output_file=None,
        success_threshold=None,
    ):
        self.data_file = data_file
        if isinstance(self.data_file, str):
            self.data_file = open(self.data_file, encoding="utf-8")  # noqa: SIM115
        self.config_file = config_file
        if isinstance(self.config_file, str):
            self.config_file = open(self.config_file, encoding="utf-8")  # noqa: SIM115
        self.calibration_data = calibration_data
        if isinstance(self.calibration_data, str):
            self.calibration_data = open(self.calibration_data, encoding="utf-8")  # noqa: SIM115
        self.output_file = output_file
        self.success_threshold = success_threshold


def test_run_bad_config():
    args = MockArgs(config_file=StringIO("[1,2,3]\n"))
    assert run(args) == -2


def test_run_bad_calibration():
    args = MockArgs(calibration_data=StringIO('"\n'))
    assert run(args) == -4


def test_commandline_parsing_config_is_file():
    config_file = path.join(EXAMPLE_DIR, "example1.yaml")
    parsed_args = parse_arguments(["open_petro_elastic", config_file])
    assert path.abspath(parsed_args.config_file.name) == path.abspath(config_file)


def test_commandline_parsing_data_is_file():
    config_file = path.join(EXAMPLE_DIR, "example1.yaml")
    data_file = path.join(EXAMPLE_DIR, "example1.csv")
    parsed_args = parse_arguments(
        ["open_petro_elastic", "--data-file", data_file, config_file]
    )
    assert path.abspath(parsed_args.data_file.name) == path.abspath(data_file)


def test_commandline_parsing_calibration_is_file():
    config_file = path.join(EXAMPLE_DIR, "example1.yaml")
    calibration_file = path.join(EXAMPLE_DIR, "example1_calibration.csv")
    parsed_args = parse_arguments(
        ["open_petro_elastic", "--calibration", calibration_file, config_file]
    )
    assert path.abspath(parsed_args.calibration_data.name) == path.abspath(
        calibration_file
    )


def test_help_string(capsys):
    with pytest.raises(SystemExit) as sysexit:
        parse_arguments(["open_petro_elastic", "-h"])
    assert sysexit.value.code == 0
    captured = capsys.readouterr()
    assert "petro-elastic modelling tool" in captured.out
