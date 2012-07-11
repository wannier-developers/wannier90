tutorial :
	latex tutorial 
	latex tutorial 
	dvips -Ppdf -o tutorial.ps tutorial.dvi
	ps2pdf tutorial.ps tutorial.pdf

clean :
	/bin/rm -f *.aux *.log *.lot *.lof *.toc *.dvi *.bak *~ *.ps

veryclean : clean
        rm -f user_guide.pdf

.PHONY:	clean veryclean
