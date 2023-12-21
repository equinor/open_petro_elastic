# open_petro_elastic

A Python library for petro-elastic modelling. It contains a `Material` class for representing rocks and fluids, as well as various rock physics models and algorithms such as Hashin-Shtrikman bounds and Gassmann fluid substitution.


## Installation

```
pip install open_petro_elastic
```

Developers and contributors can download the repository and do `pip install ".[dev,test,docs]"` to install the package with all its dependencies for development, testing, and building the docs.


## Usage

The tool can be used as a library, or from the command line using [YAML](https://en.wikipedia.org/wiki/YAML) configuration files. For example, see the [tutorial config file](tutorial_config/tutorial_config.yaml).

Get help on the command line interface:

```
open_petro_elastic --help
```

See [the docs](docs/index.rst) for more usage instructions.


## Run tests

Developers and contributors should install everything `test_requirements.txt`. Then tests can be run with:

```
pytest
```

Developers should also intall everything in `doc_requirements.txt` and read the [Code of Conduct](CODE_OF_CONDUCT.md).
