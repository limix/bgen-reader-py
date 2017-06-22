from __future__ import unicode_literals

import os

import sphinx_rtd_theme

try:
    import bgen_reader
    version = bgen_reader.__version__
except ImportError:
    version = 'unknown'

extensions = [
    'sphinx.ext.autodoc', 'sphinx.ext.doctest', 'sphinx.ext.intersphinx',
    'sphinx.ext.coverage', 'sphinx.ext.mathjax', 'sphinx.ext.viewcode'
]
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
project = 'bgen-reader-py'
copyright = '2017, Danilo Horta'
author = 'Danilo Horta'
release = version
language = None
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'conf.py']
pygments_style = 'sphinx'
todo_include_todos = False
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
htmlhelp_basename = 'bgen-reader-pydoc'
latex_elements = {}
latex_documents = [
    (master_doc, 'bgen-reader-py.tex', 'bgen-reader-py Documentation',
     'Danilo Horta', 'manual'),
]
man_pages = [(master_doc, 'bgen-reader-py', 'bgen-reader-py Documentation',
              [author], 1)]
texinfo_documents = [
    (master_doc, 'bgen-reader-py', 'bgen-reader-py Documentation', author,
     'bgen-reader-py', 'One line description of project.', 'Miscellaneous'),
]
intersphinx_mapping = {'https://docs.python.org/': None}
