
CPP=g++
SHELL=/usr/bin/bash

CPPFLAGS=-Wall -std=c++11 -O2 -march=native

LDFLAGS=
LDLIBS=-lc -lm

.PHONY: all clean

all: test

test: matrix_coo.o test.cpp force
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c test.cpp
	$(CPP) $(CPPFLAGS) -flto $(LDLIBS) $(LDFLAGS) matrix_coo.o test.o -o test
	./test

matrix_coo.o: matrix_coo.cpp
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c matrix_coo.cpp

force:

clean:
	rm -f *.o
	rm -f test
