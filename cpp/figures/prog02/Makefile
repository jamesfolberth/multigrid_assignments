
.PHONY: all force clean cleandata deepclean

all: 2_3_wj.tikz 2_3_gs.tikz 2_3_rbgs.tikz

#2_3_wj.tikz: 2_3_wj.plt
#	gnuplot 2_3_wj.plt

%.tikz: %.gpi
	gnuplot $<

force:

clean:
	@rm -f *.tikz* # gnuplot output
	@rm -f *.pdf *.log *.md5 *.dpth # tikz-externalize files

cleandata:
	@rm -f *.txt 

deepclean: clean cleandata
