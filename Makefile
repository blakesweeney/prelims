SUBDIR+=neocles
SUBDIR+=bullerjohn
SUBDIR+=larsen
SUBDIR+=hans
EXTERNAL_BIB=/Users/bsweene/Documents/bibliography/library.bib

include /Users/bsweene/.local/share/latex-mk/latex.subdir.gmk

all.bib: $(EXTERNAL_BIB)
	cp $^ $@
