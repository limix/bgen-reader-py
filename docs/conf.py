import time


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
    "sphinx.ext.autosectionlabel",
]

autodoc_default_flags = ["members"]
autodoc_mock_imports = ["_tkinter"]
autosummary_generate = True
napoleon_numpy_docstring = True
templates_path = ["_templates"]
autosectionlabel_prefix_document = False

source_suffix = ".rst"

master_doc = "index"
man_pages = [(master_doc, _get_name(), "{} documentation".format(project), [author], 1)]
language = None

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "conf.py"]

pygments_style = "default"

html_theme = "bootstrap-limix"
html_theme_options = {
    "logo_only": False,
    "display_version": True,
    "style_external_links": True,
}
htmlhelp_basename = "{}doc".format(project)

intersphinx_mapping = {
    "https://docs.python.org/": None,
    "numpy": ("https://docs.scipy.org/doc/numpy/", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable/", None),
}
