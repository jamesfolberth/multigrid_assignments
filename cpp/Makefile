
CPP=g++
SHELL=/usr/bin/bash

CPPFLAGS=-Wall -Wpedantic -std=c++11 -O2 -march=native

LDFLAGS=
LDLIBS=-lc -lm

.PHONY: all clean

all: test

test: matrix_coo.o matrix_crs.o test.cpp force
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c test.cpp
	$(CPP) $(CPPFLAGS) -flto $(LDLIBS) $(LDFLAGS) matrix_coo.o matrix_crs.o test.o -o test
	./test

matrix_coo.o: matrix_coo.cpp matrix_coo.hpp matrix_crs.hpp utils.hpp
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c matrix_coo.cpp

matrix_crs.o: matrix_crs.cpp matrix_crs.hpp matrix_coo.hpp utils.hpp
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c matrix_crs.cpp

force:

clean:
	rm -f *.o
	rm -f test
