
.PHONY: all force clean cleandata deepclean

all: 2_3_wj.tikz 2_3_gs.tikz 2_3_rbgs.tikz mix_1.tikz mix_2.tikz

#2_3_wj.tikz: force 2_3_wj.gpi
#	gnuplot 2_3_wj.gpi
#
#2_3_gs.tikz: force 2_3_gs.gpi
#	gnuplot 2_3_gs.gpi
#
#2_3_rbgs.tikz: force 2_3_rbgs.gpi
#	gnuplot 2_3_rbgs.gpi
#
#mix_1.tikz: force mix_1.gpi
#	gnuplot mix_1.gpi

%.tikz: %.gpi force
	gnuplot $<

force:

clean:
	@rm -f *.tikz* # gnuplot output
	@rm -f *.pdf *.log *.md5 *.dpth # tikz-externalize files

cleandata:
	@rm -f *.txt 

deepclean: clean cleandata
