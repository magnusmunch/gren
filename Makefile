################################### PREAMBLE ###################################

# set directories
CDIR := code
DDIR := docs
PKGDIR := rpackage

# set R package compile options
ROPTS := --no-save --no-restore --verbose

# retrieve package name and version
PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" $(PKGDIR)/DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" $(PKGDIR)/DESCRIPTION)

#################################### RULES #####################################


