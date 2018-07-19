First submission of package

## Test environments
* local macOS High sierra 10.13.5 install, R 3.4.4
* remote server Ubuntu xenial 16.04.4 install, R 3.4.4
* win-builder

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs on the local macOS

There was 1 WARNING on the remote Ubuntu:
* checking PDF version of manual ... WARNING
LaTeX errors when creating PDF version.
This typically indicates Rd problems.

There was no pdflatex installed on the remote

There was 1 ERROR on the remote Ubuntu:
* checking PDF version of manual without hyperrefs or index ... ERROR
Re-running with no redirection of stdout/stderr.
Hmm ... looks like a package
Error in texi2dvi(file = file, pdf = TRUE, clean = clean, quiet = quiet,  : 
  pdflatex is not available
Error in texi2dvi(file = file, pdf = TRUE, clean = clean, quiet = quiet,  : 
  pdflatex is not available
Error in running tools::texi2pdf()

There was no pdflatex installed on the remote

There were no ERRORs, WARNINGs or NOTEs on win-builder

## Downstream dependencies
There are currently no downstream dependencies for this package