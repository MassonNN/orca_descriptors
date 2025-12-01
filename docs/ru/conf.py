"""Sphinx configuration for Russian documentation."""

import os
import sys
from pathlib import Path

# Add the project root to the path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

# Read version from pyproject.toml
def get_version():
    """Read version from pyproject.toml."""
    project_root = Path(__file__).parent.parent.parent
    pyproject_path = project_root / "pyproject.toml"
    
    if pyproject_path.exists():
        import re
        content = pyproject_path.read_text()
        match = re.search(r'version\s*=\s*["\']([^"\']+)["\']', content)
        if match:
            return match.group(1)
    
    return "0.1.0"  # Fallback version

# Project information
project = "ORCA Descriptors"
copyright = "2024, MassonNN"
author = "MassonNN"
version = get_version()
release = version

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
language = "ru"

# HTML output
html_theme = "furo"
html_static_path = ["../_static"]
html_logo = "../_static/logo.png"
html_css_files = ["language_switcher.css"]
html_js_files = ["language_switcher.js"]
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
        "color-background-primary": "#1e1e1e",
        "color-sidebar-background": "#1e1e1e",
        "color-sidebar-item-background--hover": "#2d2d2d",
    },
}

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

# RST substitutions for version
rst_epilog = f"""
.. |version| replace:: {version}
.. |release| replace:: {release}
"""

