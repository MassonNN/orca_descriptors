# Building Documentation

This directory contains Sphinx documentation for ORCA Descriptors in two languages: English and Russian.

## Structure

- `en/` - English documentation
- `ru/` - Russian documentation
- `_static/` - Static files (CSS, images, etc.)
- `_templates/` - Custom templates

## Building Documentation

### Using Make (recommended)

Build English documentation:
```bash
cd docs/en
make html
```

Build Russian documentation:
```bash
cd docs/ru
make html
```

Build both:
```bash
make -C docs/en html
make -C docs/ru html
```

### Using Poetry

Build English documentation:
```bash
cd docs/en
poetry run sphinx-build -b html . _build/html
```

Build Russian documentation:
```bash
cd docs/ru
poetry run sphinx-build -b html . _build/html
```

## Viewing Documentation

After building, open the HTML files in your browser:

- English: `docs/en/_build/html/index.html`
- Russian: `docs/ru/_build/html/index.html`

## Theme

The documentation uses the **Furo** theme - a modern, clean theme for Sphinx documentation. Furo provides:
- Clean, modern design
- Excellent mobile responsiveness
- Dark mode support
- Fast navigation
- Accessible design

## Language Switching

The documentation includes a language switcher that allows users to switch between English and Russian versions. The switcher appears as:
- A button in the sidebar (Furo theme)
- A link in the top navigation bar (if available)

The switcher automatically detects the current language and provides a button to switch to the other language.

## Requirements

Documentation requires:
- Sphinx >= 7.0.0
- furo (modern documentation theme)

Install via Poetry:
```bash
poetry install --with dev
```

The `furo` theme is automatically installed as a development dependency.

