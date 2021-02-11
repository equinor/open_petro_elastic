import argparse
import sys
import warnings
from functools import wraps
from os import path
from shutil import copytree
from traceback import print_tb

import numpy as np
import pandas as pd
import yaml
from pydantic import parse_obj_as

import open_petro_elastic

from .config.input import ColumnInsertionError, Input
from .material import fluid_substitution


def warn_without_traceback(message, category, filename, lineno, file=None, line=None):

    log = file if hasattr(file, "write") else sys.stderr
    log.write(f"Warning: {message}\n")


def calibration_error(calculated, calibrated):
    return (calculated - calibrated) / calibrated


def print_calibration_error(name, error):
    print(
        f"""{name}:
min: {np.min(error)}
mean: {np.mean(error)}
max: {np.max(error)}
"""
    )


def calculate_results(inp):
    mixed_mineral = inp.minerals.as_mixture
    mixed_fluid = inp.fluids.as_mixture(inp.pressure)
    dry_rock_material = inp.dry_rock.material(mixed_mineral, inp.pressure)
    dense_packing = inp.dry_rock.nonporous(mixed_mineral)
    saturated_rock = fluid_substitution(
        dry_rock_material,
        dense_packing,
        mixed_fluid,
        inp.dry_rock.porosity,
    )
    index = [0]
    if inp.data is not None:
        index = inp.data.index
    return pd.DataFrame(
        {
            "ksat": saturated_rock.bulk_modulus,
            "kmin": mixed_mineral.bulk_modulus,
            "kdry": dry_rock_material.bulk_modulus,
            "mysat": saturated_rock.shear_modulus,
            "rsat": saturated_rock.density,
            "kmin_fls": dense_packing.bulk_modulus,
        },
        index=index,
    )


def takes_stream(i, mode):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            if len(args) > i and args[i] is not None and isinstance(args[i], str):
                with open(args[i], mode) as f:
                    return func(*args[:i], f, *args[i + 1 :], **kwargs)
            else:
                return func(*args, **kwargs)

        return wrapper

    return decorator


class InputError(ValueError):
    pass


def load_config_from_yaml(config_file):
    try:
        config = yaml.safe_load(config_file)
    except yaml.YAMLError as exc:
        if hasattr(exc, "problem_mark"):
            m = "Error while parsing the config file:\n"
            if exc.context is not None:
                m += f"{exc.problem_mark}\n{exc.problem}\n{exc.context}"
            else:
                m += f"\n{exc.problem_mark}\n{exc.problem}"
        else:
            m = "Something went wrong while parsing the config file."
        raise InputError(m) from exc
    if not isinstance(config, dict):
        raise InputError("Unexpected format of config file")
    return config


def load_data_from_csv(data_file):
    try:
        if data_file is not None:
            data = pd.read_csv(data_file)
        else:
            data = None
    except pd.errors.ParserError as e:
        raise InputError("Error while reading data file: {e}") from e
    except pd.errors.EmptyDataError as e:
        raise InputError("Data file seems to be empty") from e
    return data


@takes_stream(0, "r")
def make_input(config_file, data_file):
    """
    :param config_file: Either path to or opened file to yaml config.
    :param data_file: Either path to or opened file to csv data file.
    :return: Input read from config and data
    """
    config = load_config_from_yaml(config_file)
    data = load_data_from_csv(data_file)
    try:
        config["data"] = data
        return parse_obj_as(Input, config)
    except ColumnInsertionError as e:
        raise InputError(f"Error while inserting column data from csv: {e}") from e
    except ValueError as e:
        raise InputError(f"Inconsistent values in input file: {e}") from e


def make_per_row_input(config_file, data_file):
    config = load_config_from_yaml(config_file)
    data = load_data_from_csv(data_file)
    inputs = []
    for i, row in data.iterrows():
        try:
            config["data"] = row.to_frame().transpose()
            inputs.append(parse_obj_as(Input, config))
        except ColumnInsertionError as e:
            print(f"Error while inserting column data from csv: {e}")
            inputs.append(None)
        except ValueError as e:
            print(f"Inconsistent values in input file for line {i}: {e}")
            inputs.append(None)
        except Exception as e:
            print(f"Encountered error while creating input for line {i}: {e}")
            inputs.append(None)
    return inputs


class CalibrationError(ValueError):
    pass


def run_calibration(calibration_file, results):
    try:
        calibration_data = pd.read_csv(calibration_file)
        for key in ["kmin", "kdry", "ksat", "kmin_fls"]:
            result = results[key]
            calibration = calibration_data[key]
            error = calibration_error(result, calibration)
            print_calibration_error(key, error)
    except pd.errors.ParserError as e:
        raise CalibrationError(f"Could not read calibration file: {e}") from e
    except pd.errors.EmptyDataError as e:
        raise CalibrationError("Calibration file seems to be empty") from e
    except (KeyError, IndexError) as e:
        raise CalibrationError(
            f"Error while getting column from calibration data: {e}"
        ) from e


def run_with_threshold(args):
    if args.data_file is None:
        print("In order to use success threshold, data must be given in csv file.")
        return -5
    inputs = make_per_row_input(args.config_file, args.data_file)
    num_failures = 0
    results = pd.DataFrame()
    for idx, inp in enumerate(inputs):
        if inp is None:
            num_failures += 1
            continue
        try:
            results = results.append(calculate_results(inp))
        except Exception as e:
            print(f"Encountered error while calculating rows {id}: {e}")
            num_failures += 1
    if len(inputs) == 0:
        success_percentage = 100
    else:
        success_percentage = (1 - num_failures / len(inputs)) * 100
    if success_percentage < args.success_threshold:
        print(
            f"Percentage of successes {success_percentage} did not reach threshold {args.success_threshold}"
        )
        return -5
    try:
        results.to_csv(args.output_file)
    except Exception as e:
        print(f"Encountered unexpected error while writing output: {e}")
        print_tb(e.__traceback__)
        return -1

    if args.calibration_data is not None:
        try:
            run_calibration(args.calibration_data, results)
        except CalibrationError as e:
            print(f"Error encountered in calibration: {e}")
            return -4
        except Exception as e:
            print(f"Encountered unexpected error during calibration: {e}")
            print_tb(e.__traceback__)
            return -1
    return 0


def run(args):
    if args.success_threshold is not None:
        return run_with_threshold(args)
    try:
        inp = make_input(args.config_file, args.data_file)
    except InputError as e:
        print(f"Could not read input: {e}")
        return -2
    except Exception as e:
        print(f"Encountered unexpected error while reading input: {e}")
        print_tb(e.__traceback__)
        return -1

    try:
        results = calculate_results(inp)
    except ValueError as e:
        print(f"Encountered invalid value in petro elastic calculation: {e}")
        return -3
    except Exception as e:
        print(f"Encountered unexpected error in petro elastic calculation: {e}")
        print_tb(e.__traceback__)
        return -1

    try:
        results.to_csv(args.output_file)
    except Exception as e:
        print(f"Encountered unexpected error while writing output: {e}")
        print_tb(e.__traceback__)
        return -1

    if args.calibration_data is not None:
        try:
            run_calibration(args.calibration_data, results)
        except CalibrationError as e:
            print(f"Error encountered in calibration: {e}")
            return -4
        except Exception as e:
            print(f"Encountered unexpected error during calibration: {e}")
            print_tb(e.__traceback__)
            return -1
    return 0


def generate_tutorial(out_dir):
    here = path.dirname(__file__)
    tutorial_path = path.join(here, "tutorial_config")
    copytree(tutorial_path, out_dir)


class GenerateTutorialAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        try:
            generate_tutorial(values)
        except Exception as err:
            print(f"Could not write tutorial config: {err}")
            parser.exit(-4)
        parser.exit()


def percentage(value):
    try:
        value = float(value)
    except (ValueError, TypeError):
        raise argparse.ArgumentTypeError(f"Value must be float {value}")
    if 0.0 <= value <= 100.0:
        return value
    raise argparse.ArgumentTypeError(f"Value must be in range [0, 100] {value}")


TOOL_DESCRIPTION = """
The Equinor petro-elastic modelling tool. Models elastic properties
of sandstone saturated with oil, gas, and brine. See github.com/equinor/open_petro_elastic
for documentation.
"""


def parse_arguments(argv):
    ap = argparse.ArgumentParser(prog=argv[0], description=TOOL_DESCRIPTION)
    ap.add_argument(
        "config_file", type=argparse.FileType("r"), help="petro elastic yaml input file"
    )
    ap.add_argument(
        "--data-file",
        type=argparse.FileType("r"),
        dest="data_file",
        help="Optional additional csv file for data to be inserted into config.",
    )
    ap.add_argument(
        "--calibration",
        type=argparse.FileType("r"),
        dest="calibration_data",
        help="Optional calibration data, if given, an error statistic is calculated and printed.",
    )
    ap.add_argument(
        "--output_file",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="Optional output file, defaults to stdout.",
    )
    ap.add_argument(
        "--version", action="version", version=open_petro_elastic.__version__
    )
    ap.add_argument(
        "--generate-tutorial",
        action=GenerateTutorialAction,
        help="Write tutorial config. Takes directory to write to as argument.",
    )
    ap.add_argument(
        "--success_threshold",
        type=percentage,
        default=None,
        help="Optional percentage required for success."
        + " Default is 100 percent with no partial results given."
        + " Unsuccessful rows will be given NaN. Note that "
        + "this option can deteriorate performance.",
    )
    return ap.parse_args(argv[1:])


def main():
    warnings.showwarning = warn_without_traceback
    sys.exit(run(parse_arguments(sys.argv)))


if __name__ == "__main__":
    main()
