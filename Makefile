.PHONY: help docs docs-en docs-ru clean

help:
	@echo "Available targets:"
	@echo "  docs-en    - Build English documentation"
	@echo "  docs-ru    - Build Russian documentation"
	@echo "  docs       - Build both English and Russian documentation"
	@echo "  clean      - Clean build files"

docs-en:
	cd docs/en && poetry run sphinx-build -b html . _build/html

docs-ru:
	cd docs/ru && poetry run sphinx-build -b html . _build/html

docs: docs-en docs-ru

clean:
	rm -rf docs/en/_build
	rm -rf docs/ru/_build

