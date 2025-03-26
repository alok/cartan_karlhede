"""Sphinx configuration file for the Cartan-Karlhede documentation."""

import os
import sys

# Add the project root directory to the Python path
sys.path.insert(0, os.path.abspath(".."))

# Project information
project = "Cartan-Karlhede Algorithm"
copyright = "2023, Alok Singh"
author = "Alok Singh"
version = "0.1.0"
release = "0.1.0"

# General Sphinx configuration
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# HTML output configuration
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

# LaTeX configuration
latex_elements = {
    "preamble": r"""
    \usepackage{amsmath}
    \usepackage{amssymb}
    \usepackage{bm}
    """,
}

# Extension configuration
autoclass_content = "both"
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_docstring_signature = True

# Enable math notation in docstrings
mathjax_path = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"
mathjax3_config = {
    "tex": {
        "inlineMath": [["$", "$"], ["\\(", "\\)"]],
        "displayMath": [["$$", "$$"], ["\\[", "\\]"]],
    }
}
