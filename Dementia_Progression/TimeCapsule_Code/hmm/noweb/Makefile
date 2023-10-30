PARTS=  part1.Rnw \
	part2.Rnw\
	part3.Rnw  part3b.Rnw\
	part4.Rnw \
	response.Rnw \
	score.Rnw \
        tail

CFUN  = cdecomp.c 

SFUN =  hmm.R \
	hmmdata.R \
	hmmscore.R \
	derivtest.R \
	response.R  

.PHONY: fun clean all
RDIR= ../R
CCODE= ../src
RFUN = $(SFUN:%=$(RDIR)/%)
CFUN2= $(CFUN:%=../src/%)

all: fun hmmcode.pdf

.PHONY: fun clean 
fun: $(RFUN) $(CFUN2) noweb.sty

$(RFUN): hmmcode.nw
$(CFUN2): hmmcode.nw

hmmcode.pdf: hmmcode.tex noweb.sty
	pdflatex hmmcode.tex
	pdflatex hmmcode.tex

hmmcode.nw: $(PARTS)
	cat $(PARTS) > hmmcode.nw

hmmcode.tex: hmmcode.nw
	echo "require(noweb); noweave('hmmcode.nw', out='$@')" | R --slave

%.R:	
	echo "# Automatically generated from the hmmcode" > $@
	echo "require(noweb); notangle('hmmcode.nw', target='$(*F)', out='zz')" | R --slave
	cat zz >> $@
	rm zz

%.c: 
	echo "/* Automatically generated from the noweb directory */" > $@
	echo "require(noweb); notangle('hmmcode.nw', target='$(*F)', out='zz')" | R  --slave
	cat zz >> $@
	rm zz

noweb.sty:
	echo 'library(noweb); data(noweb); cat(noweb.sty, sep="\n", file="noweb.sty")' |  R --slave

clean: 
	-rm *.log *.aux *.tex *.out
	-rm noweb.sty
