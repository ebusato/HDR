
main: 
	pdflatex hdr.tex
	mpost fgraph-4topsSM
	mpost fgraph-ContactGluFu
	mpost fgraph-2UEDRPP
	bibtex hdr
	pdflatex hdr.tex
	pdflatex hdr.tex

clean:
	rm -f *.aux *.log *.bbl *.blg fgraph*
