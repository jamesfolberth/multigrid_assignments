
.PHONY: prog02 force clean

prog02:
	make -C prog02

prog03:
	make -C prog03

force:

clean:
	@rm -f *.pdf *.log *.dpth *.md5 # tikz-externalize stuff

cleandata:
	make -C prog02 cleandata
	make -C prog03 cleandata

deepclean: clean cleandata
	make -C prog02 deepclean
	make -C prog03 deepclean
