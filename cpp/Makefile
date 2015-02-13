# GNU Makefile

# use normal g++ if you don't have GCC version >= 4.9
#CPP=g++ 
CPP=g++ -fdiagnostics-color=auto
SHELL=/usr/bin/bash

CPPFLAGS=-Wall -Wpedantic -std=c++11 -g -Og -march=native

LDFLAGS=
LDLIBS=-lc -lm

MATLIB_COMMON= matrix_coo.cpp\
					matrix_coo.hpp\
					matrix_crs.cpp\
					matrix_crs.hpp\
					utils.hpp\
					utils.cpp\

.PHONY: all clean run prog02_fig

all: test run

run: test
	./test

# Code
test: matrix_coo.o matrix_crs.o utils.o model_problems.o classical_solvers.o test.cpp
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c test.cpp
	$(CPP) $(CPPFLAGS) -flto $(LDLIBS) $(LDFLAGS) matrix_coo.o\
		matrix_crs.o utils.o model_problems.o classical_solvers.o test.o -o test

prog02: matrix_coo.o matrix_crs.o utils.o model_problems.o classical_solvers.o prog02.cpp
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c prog02.cpp
	$(CPP) $(CPPFLAGS) -flto $(LDLIBS) $(LDFLAGS) matrix_coo.o\
		matrix_crs.o utils.o model_problems.o classical_solvers.o \
		prog02.o -o prog02

matrix_coo.o: $(MATLIB_COMMON)
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c matrix_coo.cpp

matrix_crs.o: $(MATLIB_COMMON)
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c matrix_crs.cpp

utils.o: utils.cpp utils.hpp
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c utils.cpp

model_problems.o: $(MATLIB_COMMON) model_problems.cpp model_problems.hpp
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c model_problems.cpp

classical_solvers.o: $(MATLIB_COMMON) classical_solvers.cpp classical_solvers.hpp
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c classical_solvers.cpp


# figures
prog02_fig: prog02
	./prog02
	make -C figures/ prog02


force:

clean:
	@rm -f *.o
	@rm -f test prog02
	@rm -f *.out # gprof

cleandata:
	make -C figures/ cleandata

deepclean: clean cleandata
	make -C figures/ deepclean

