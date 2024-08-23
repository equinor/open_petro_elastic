import argparse
import sys

import numpy as np

from open_petro_elastic.material.span_wagner.carbon_dioxide import (
    carbon_dioxide_density,
)
from open_petro_elastic.material.span_wagner.tables.lookup_table import (
    generate_lookup_table,
)

t_density = 273.15 + np.hstack(
    (np.arange(0, 50, step=0.05), np.arange(50, 200.1, step=0.5))
)
p_density = np.hstack((np.arange(0.5, 10, step=0.01), np.arange(10, 100.1, step=0.1)))


def run(args):
    generate_lookup_table(
        lambda _t, _p: carbon_dioxide_density(
            _t, _p, interpolate=False, force_vapor="auto"
        ),
        t_density,
        p_density,
        args.output_file,
    )
    return 0


TOOL_DESCRIPTION = """
Used to generate the tables included in open_petro_elastic.material.span_wagner.tables.
"""


def make_parser(prog="generate_sp_tables"):
    ap = argparse.ArgumentParser(prog=prog, description=TOOL_DESCRIPTION)
    ap.add_argument(
        "output_file",
        type=argparse.FileType("wb"),
        help="The output file to write the density data to. Data is stored on the uncompressed numpy npz format and it"
        " is recommended to use file suffix .npz",
    )
    return ap


def parse_arguments(argv):
    ap = make_parser(argv[0])
    return ap.parse_args(argv[1:])


def main():
    sys.exit(run(parse_arguments(sys.argv)))


if __name__ == "__main__":
    main()
