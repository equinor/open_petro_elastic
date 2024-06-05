project = "open_petro_elastic"
copyright = "2021, Equinor"
author = "Equinor"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinxarg.ext",
    "sphinx.ext.viewcode",
    "sphinx.ext.doctest",
    "sphinx.ext.githubpages",
    "matplotlib.sphinxext.plot_directive",
    "recommonmark",
]

autosummary_generate = True
master_doc = "index"
html_theme = "sphinx_rtd_theme"
# templates_path = ["_templates"]
