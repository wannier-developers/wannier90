tutorial.pdf: tutorial.tex
	pdflatex tutorial.tex 
	bibtex tutorial
	pdflatex tutorial.tex
	bibtex tutorial
	pdflatex tutorial.tex
#	mv tutorial.pdf ../
#	dvips -Ppdf -o tutorial.ps tutorial.dvi
#	ps2pdf tutorial.ps tutorial.pdf

clean :
	/bin/rm -f *.aux *.log *.lot *.lof *.toc *.dvi *.bak *~ *.ps

veryclean : clean
	rm -f tutorial.pdf

.PHONY:	clean veryclean
