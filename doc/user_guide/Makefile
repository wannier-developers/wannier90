TEXFILES=$(wildcard *.tex)

user_guide.pdf: $(TEXFILES)
	pdflatex user_guide 
	bibtex user_guide
	pdflatex user_guide
	bibtex user_guide 
	pdflatex user_guide
#	dvips -Ppdf -o user_guide.ps user_guide.dvi
#	ps2pdf user_guide.ps user_guide.pdf

clean :
	rm -f *.aux *.log *.lot *.lof *.toc *.dvi *.bak *~ *.bbl *.blg user_guide.ps *.out

veryclean : clean
	rm -f user_guide.pdf

.PHONY: clean veryclean


