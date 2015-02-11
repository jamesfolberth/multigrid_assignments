
.PHONY: prog02 force clean

prog02:
	make -C prog02

force:

clean:
	@rm -f *.pdf *.log *.dpth *.md5 # tikz-externalize stuff

cleandata:
	make -C prog02 cleandata

deepclean: clean cleandata
	make -C prog02 deepclean
