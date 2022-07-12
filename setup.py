from pathlib import Path

from setuptools import find_packages, setup


def get_long_description() -> str:
    return Path("README.md").read_text(encoding="utf8")


setup(
    name="open_petro_elastic",
    author="Equinor",
    author_email="fg_sib-scout@equinor.com",
    license="LGPL",
    url="https://github.com/equinor/open_petro_elastic",
    description="Utility for calculating elastic properties of petroleum fields.",
    long_description=get_long_description(),
    long_description_content_type="text/markdown",
    use_scm_version=True,
    setup_requires=["setuptools_scm"],
    entry_points={
        "console_scripts": ["open_petro_elastic = open_petro_elastic.__main__:main"],
        "open_petro_elastic.fluid_model_providers": [
            "batzle_wang = open_petro_elastic.config.fluid_model_providers:BatzleWangFluidModelProvider",
            "span_wagner = open_petro_elastic.config.fluid_model_providers:SpanWagnerFluidModelProvider",
        ],
    },
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    install_requires=[
        "numpy",
        "scipy",
        "pyyaml",
        "pandas",
        "pydantic",
        "sympy",
        "typing_extensions",
    ],
    platforms="any",
    classifiers=[
        "Development Status :: 1 - Planning",
        "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    python_requires=">=3.6",
    include_package_data=True,
    package_data={
        "open_petro_elastic": ["tutorial_config/*"],
        "open_petro_elastic.material.span_wagner.tables": [
            "material/span_wagner/tables/carbon_dioxide_density.npz"
        ],
    },
)
