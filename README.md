# open_petro_elastic

A Python library for petro-elastic modelling. It contains a `Material` class for representing rocks and fluids, as well as various rock physics models and algorithms such as Hashin-Shtrikman bounds and Gassmann fluid substitution.

[![Build and test](https://github.com/equinor/open_petro_elastic/actions/workflows/python-build-test.yml/badge.svg)](https://github.com/equinor/open_petro_elastic/actions/workflows/python-build-test.yml)
[![Build documentation](https://github.com/equinor/open_petro_elastic/actions/workflows/python-sphinx-doc.yml/badge.svg)](https://github.com/equinor/open_petro_elastic/actions/workflows/python-sphinx-doc.yml)
[![PyPI version](https://img.shields.io/pypi/v/open_petro_elastic.svg)](https://pypi.org/project/open_petro_elastic//)
[![PyPI versions](https://img.shields.io/pypi/pyversions/open_petro_elastic.svg)](https://pypi.org/project/open_petro_elastic//)
[![PyPI license](https://img.shields.io/pypi/l/open_petro_elastic.svg)](https://pypi.org/project/open_petro_elastic/)


## Installation

```shell
pip install open_petro_elastic
```

Developers and contributors can download the repository and do `pip install ".[dev,test,docs]"` to install the package with all its dependencies for development, testing, and building the docs.


## Usage

The tool can be used as a library, or from the command line using [YAML](https://en.wikipedia.org/wiki/YAML) configuration files. For example, see the [tutorial config file](tutorial_config/tutorial_config.yaml).

Generally, the `open_petro_elastic.material` module is used for interfacing with `open_petro_elastic` as a Python library:

```python
>>> from open_petro_elastic.material import Material
>>> from open_petro_elastic.material.sandstone import hertz_mindlin
>>> mineral = Material(bulk_modulus=1e9, shear_modulus=1e9, density=1000)
>>> sand = hertz_mindlin(mineral, porosity=0.4, pressure=1e6)
>>> print(sand.density)
600.0
```

To get help on the command line interface:

```shell
open_petro_elastic --help
```

See [the docs](https://equinor.github.io/open_petro_elastic/) for more usage instructions.


## Run tests

Developers and contributors should install everything `test_requirements.txt`. Then tests can be run with:

```
pytest
```

Developers should also intall everything in `doc_requirements.txt` and read the [Code of Conduct](CODE_OF_CONDUCT.md).
