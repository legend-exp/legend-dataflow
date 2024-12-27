# Configuration file for the Sphinx documentation builder.
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, Path(__file__).parents[2].resolve().as_posix())

project = "legend-dataflow"
copyright = "2024, the LEGEND Collaboration"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx_copybutton",
    "sphinx_inline_tabs",
    "myst_parser",
    "IPython.sphinxext.ipython_console_highlighting",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}
master_doc = "index"

# Furo theme
html_theme = "furo"
html_theme_options = {
    "source_repository": "https://github.com/legend-exp/legend-dataflow",
    "source_branch": "main",
    "source_directory": "docs/source",
}
html_title = f"{project}"

# sphinx-napoleon
# enforce consistent usage of NumPy-style docstrings
napoleon_numpy_docstring = True
napoleon_google_docstring = False
napoleon_use_ivar = True
napoleon_use_rtype = False

# intersphinx
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "awkward": ("https://awkward-array.org/doc/stable", None),
    "numba": ("https://numba.readthedocs.io/en/stable", None),
    "pandas": ("https://pandas.pydata.org/docs", None),
    "h5py": ("https://docs.h5py.org/en/stable", None),
    "pint": ("https://pint.readthedocs.io/en/stable", None),
    "hist": ("https://hist.readthedocs.io/en/latest", None),
    "dspeed": ("https://dspeed.readthedocs.io/en/stable", None),
    "daq2lh5": ("https://legend-daq2lh5.readthedocs.io/en/stable", None),
    "lgdo": ("https://legend-pydataobj.readthedocs.io/en/stable", None),
    "dbetto": ("https://dbetto.readthedocs.io/en/stable", None),
    "pylegendmeta": ("https://pylegendmeta.readthedocs.io/en/stable", None),
}  # add new intersphinx mappings here

# sphinx-autodoc
autodoc_default_options = {"ignore-module-all": True}
# Include __init__() docstring in class docstring
autoclass_content = "both"
autodoc_typehints = "description"
autodoc_typehints_description_target = "documented_params"
autodoc_typehints_format = "short"
