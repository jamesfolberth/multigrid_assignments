
.PHONY: prog02 force clean

prog02:
	make -C prog02

prog03:
	make -C prog03

prog04:
	make -C prog04

force:

clean:
	@rm -f *.pdf *.log *.dpth *.md5 # tikz-externalize stuff

cleandata:
	make -C prog02 cleandata
	make -C prog03 cleandata
	make -C prog04 cleandata

deepclean: clean cleandata
	make -C prog02 deepclean
	make -C prog03 deepclean
	make -C prog04 deepclean

