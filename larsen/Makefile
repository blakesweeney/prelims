NAME=writeup
TEXSRCS=$(patsubst %.mkd,%.tex,$(wildcard questions/*.mkd))
BIBTEXSRCS=all.bib
PANDOC_OPTS=--from=markdown+tex_math_dollars+tex_math_single_backslash+implicit_figures --to=latex --biblatex --bibliography='$(BIBTEXSRCS)'

include /Users/bsweene/.local/share/latex-mk/latex.gmk

questions/%.tex: questions/%.mkd
	pandoc $< $(PANDOC_OPTS) -o $@
