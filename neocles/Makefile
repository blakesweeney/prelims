NAME=writeup
TEXSRCS=$(patsubst %.mkd,%.tex,$(wildcard questions/*.mkd))
BIBTEXSRCS=all.bib
PANDOC_OPTS=--from=markdown+tex_math_dollars+tex_math_single_backslash+implicit_figures --to=latex --biblatex --bibliography='$(BIBTEXSRCS)'

include /Users/bsweene/.local/share/latex-mk/latex.gmk

writeup.tex: figs

figs: figs/properties.latex figs/diffusion.ps

questions/%.tex: questions/%.mkd
	pandoc $< $(PANDOC_OPTS) -o $@

tmp/properties.pickle: bin/compute-table data/molecules.json libs/molecular_properties.py
	$^ > $@

figs/properties.latex: bin/table2latex tmp/properties.pickle
	$^ > $@

tmp/properties.csv: bin/table2csv tmp/properties.pickle
	$^ > $@

figs/diffusion.ps: bin/plot-diffusion tmp/properties.csv libs/diffusion.R
	$^ $@
