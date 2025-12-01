.PHONY: help docs docs-en docs-ru docs-update clean

help:
	@echo "Available targets:"
	@echo "  docs-en     - Build English documentation"
	@echo "  docs-ru     - Build Russian documentation"
	@echo "  docs        - Build both English and Russian documentation"
	@echo "  docs-update - Update translation files (gettext)"
	@echo "  clean       - Clean build files"

docs-en:
	cd docs/source && poetry run sphinx-build -b html -D language=en . _build/html/en

docs-ru:
	cd docs/source && poetry run sphinx-build -b html -D language=ru . _build/html/ru

docs: docs-en docs-ru

docs-update:
	cd docs/source && poetry run sphinx-build -b gettext . _build/gettext
	cd docs/source && poetry run sphinx-intl update -p _build/gettext -d ../locale -l ru

clean:
	rm -rf docs/source/_build
	rm -rf docs/locale/*/LC_MESSAGES/*.mo

