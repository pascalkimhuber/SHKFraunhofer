BASENAME = Tremolo_Manual
LATEX = pdflatex

pdf:
	-./clean.sh
	-$(LATEX) $(BASENAME)
	bibtex $(BASENAME)
	-$(LATEX) $(BASENAME)
	$(LATEX) -interaction=nonstopmode $(BASENAME)
