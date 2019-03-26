################################### PREAMBLE ###################################

# set directories
CDIR := code
DDIR := docs
RDIR := results
PKGDIR := rpackage

# set R script options
ROPTS := --no-save --no-restore --verbose

# retrieve package name and version
PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" $(PKGDIR)/DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" $(PKGDIR)/DESCRIPTION)

#################################### RULES #####################################

# building the R package
build:
	R CMD build $(PKGDIR)
	
# checking the R package build
check: build
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz

# run knitr if a knitr file changed
$(DDIR)/%.tex: $(DDIR)/%.Rnw
	Rscript -e "knitr::knit('$<', '$(DDIR)/%.tex')"

# recompile pdf if tex or bib files changed
$(DDIR)/%.pdf: $(DDIR)/%.tex $(DDIR)/refs.bib $(DDIR)/author_short3.bst
	pdflatex -output-directory=$(DDIR) $(DDIR)/%.tex