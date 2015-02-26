
.PHONY: all force clean cleandata deepclean

all: 

%.tikz: %.gpi force
	gnuplot $<

force:

clean:
	@rm -f *.tikz* # gnuplot output
	@rm -f *.pdf *.log *.md5 *.dpth # tikz-externalize files

cleandata:
	@rm -f *.txt 

deepclean: clean cleandata
