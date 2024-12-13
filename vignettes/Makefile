LATEXFILES = *.aux\
        *.bbl\
        *.blg\
        *.ilg\
        *.log\
        *.nlo\
        *.nls\
        *.toc\
        *.aux\
        mzR.tex\
        Rplots.pdf\
	mzR-plotspectrum.pdf\
        *.dvi\
        *.map\
        *.out\

all:
	make vignette
	make clean

setvars:
ifeq (${R_HOME},)
R_HOME= $(shell R RHOME)
endif

vignette:
	"$(R_HOME)/bin/R" CMD Sweave spCCA.Rnw
	pdflatex spCCA.tex
	bibtex spCCA
	pdflatex spCCA.tex
	pdflatex spCCA.tex

clean:
	rm -f $(LATEXFILES)
	rm -f *~

