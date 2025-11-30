"""Sphinx configuration for English documentation."""

import os
import sys
from pathlib import Path

# Add the project root to the path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

# Project information
project = "ORCA Descriptors"
copyright = "2024, MassonNN"
author = "MassonNN"
release = "0.1.0"
version = "0.1.0"

# Extensions
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
]

# Templates
templates_path = ["../_templates"]

# Language
language = "en"

# HTML output
html_theme = "furo"
html_static_path = ["../_static"]
html_logo = "../_static/logo.jpg"
html_theme_options = {
    "navigation_with_keys": True,
    "sidebar_hide_name": False,
    "light_css_variables": {
        "color-brand-primary": "#2980b9",
        "color-brand-content": "#3498db",
        "color-background-primary": "#ffffff",
        "color-sidebar-background": "#ffffff",
        "color-sidebar-item-background--hover": "#e9ecef",
    },
    "dark_css_variables": {
        "color-brand-primary": "#3498db",
        "color-brand-content": "#5dade2",
        "color-background-primary": "#ffffff",
        "color-sidebar-background": "#ffffff",
        "color-sidebar-item-background--hover": "#2d2d2d",
    },
}

# Add language switcher script and CSS
html_js_files = [
    'language_switcher.js',
]
html_css_files = [
    'language_switcher.css',
]

# Autodoc settings
autodoc_member_order = "bysource"
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True

