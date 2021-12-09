# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('../'))
#TOP_DIR = os.path.abspath(os.path.join('..', '..'))
#sys.path.insert(0, TOP_DIR)
# -- Project information -----------------------------------------------------

project = 'multiflap'
copyright = '2020, Gianmarco Ducci'
author = 'Gianmarco Ducci'

# The full version, including alpha/beta/rc tags
release = '1.1'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
#extensions = [ 'sphinx.ext.napoleon','sphinx.ext.coverage', 'sphinxcontrib.bibtex', 'sphinxcontrib.programoutput', 'sphinx.ext.todo', 'sphinx.ext.viewcode','sphinx.ext.autodoc']
extensions = ["sphinx.ext.napoleon", "sphinx.ext.viewcode", 'sphinx.ext.autodoc']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

master_doc = 'index'
# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

pygments_style = 'sphinx'
# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

html_theme = 'sphinx_rtd_theme'

#extensions = ['sphinxcontrib.contentui']
#if html_theme != 'sphinx_rtd_theme':
 #   def setup(app):
  #      app.add_stylesheet('custom.css')
# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
# The paper size ('letterpaper' or 'a4paper').
#'papersize': 'letterpaper',

'extraclassoptions': 'openany,oneside'
# The font size ('10pt', '11pt' or '12pt').
#'pointsize': '10pt',

# Additional stuff for the LaTeX preamble.
#'preamble': '',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).

latex_documents = [
  ('index', 'ReadtheDocsTemplate.tex', u'Read the Docs Template Documentation',
   u'Read the Docs', 'manual'),
]

def setup(app):
    app.add_css_file('my_width.css')

# The name of an image file (relative to this directory) to place at the top of
# the title page.
#latex_logo = None

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
#latex_use_parts = False

# If true, show page references after internal links.
#latex_show_pagerefs = False

# If true, show URL addresses after external links.
#latex_show_urls = False

# Documents to append as an appendix to all manuals.
#latex_appendices = []

# If false, no module index is generated.
#latex_domain_indices = True
