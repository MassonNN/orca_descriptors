SPHINXOPTS    ?=
SPHINXBUILD  ?= sphinx-build
SOURCEDIR    = .
BUILDDIR     = _build

# Languages
LANGUAGES    = en ru

.PHONY: help Makefile gettext update update-po html html-en html-ru clean

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Generate gettext translation files
gettext:
	@$(SPHINXBUILD) -b gettext "$(SOURCEDIR)" "$(BUILDDIR)/gettext"

# Update translation files
update: gettext
	@echo "Updating translation files..."
	@sphinx-intl update -p "$(BUILDDIR)/gettext" -d ../locale -l ru

# Update PO files (alias for update)
update-po: update

# Build HTML for all languages
html: html-en html-ru

# Build HTML for English
html-en:
	@$(SPHINXBUILD) -b html -D language=en "$(SOURCEDIR)" "$(BUILDDIR)/html/en" $(SPHINXOPTS) $(O)

# Build HTML for Russian
html-ru:
	@$(SPHINXBUILD) -b html -D language=ru "$(SOURCEDIR)" "$(BUILDDIR)/html/ru" $(SPHINXOPTS) $(O)

# Default target
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
