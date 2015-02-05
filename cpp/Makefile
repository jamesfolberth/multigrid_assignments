
CPP=g++
SHELL=/usr/bin/bash

CPPFLAGS=-Wall -Wpedantic -g -std=c++11 -O2 -march=native

LDFLAGS=
LDLIBS=-lc -lm

MATLIB_COMMON= matrix_coo.cpp\
					matrix_coo.hpp\
					matrix_crs.cpp\
					matrix_crs.hpp\
					utils.hpp

.PHONY: all clean

all: test

test: matrix_coo.o matrix_crs.o model_problems.o test.cpp force
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c test.cpp
	$(CPP) $(CPPFLAGS) -flto $(LDLIBS) $(LDFLAGS) matrix_coo.o\
		matrix_crs.o model_problems.o test.o -o test
	./test

matrix_coo.o: $(MATLIB_COMMON)
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c matrix_coo.cpp

matrix_crs.o: $(MATLIB_COMMON)
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c matrix_crs.cpp

model_problems.o: $(MATLIB_COMMON) model_problems.cpp model_problems.hpp
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c model_problems.cpp

force:

clean:
	rm -f *.o
	rm -f test
