import time

import sphinx_rtd_theme


def _get_version():
    import bgen_reader

    return bgen_reader.__version__


def _get_name():
    import bgen_reader

    return bgen_reader.__name__


project = _get_name()
copyright = "2018, Danilo Horta"
author = "Danilo Horta"

version = _get_version()
release = version
today = time.strftime("%B %d, %Y")

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
]

autodoc_default_flags = ["members"]
autodoc_mock_imports = ["_tkinter"]
autosummary_generate = True
napoleon_numpy_docstring = True
templates_path = ["_templates"]

source_suffix = ".rst"

master_doc = "index"
man_pages = [(master_doc, _get_name(), "{} documentation".format(project), [author], 1)]
language = None

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "conf.py"]

pygments_style = "default"

html_theme = "sphinx_rtd_theme"
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
html_sidebars = {"**": ["relations.html", "searchbox.html"]}
htmlhelp_basename = "{}doc".format(project)

intersphinx_mapping = {
    "https://docs.python.org/": None,
    "numpy": ("https://docs.scipy.org/doc/numpy/", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable/", None),
}
