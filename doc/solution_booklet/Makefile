tutorial.pdf: solution_booklet.tex
	pdflatex solution_booklet.tex 
	bibtex   solution_booklet
	pdflatex solution_booklet.tex
	bibtex   solution_booklet
	pdflatex solution_booklet.tex
#	mv solution_booklet.pdf ../
#	dvips -Ppdf -o solution_booklet.ps solution_booklet.dvi
#	ps2pdf solution_booklet.ps solution_booklet.pdf

clean :
	/bin/rm -f *.aux *.log *.lot *.lof *.toc *.dvi *.bak *~ *.ps *.out *.bbl *.blg

veryclean : clean
	rm -f solution_booklet.pdf

.PHONY:	clean veryclean
