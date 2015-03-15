# GNU Makefile

# use normal g++ if you don't have GCC version >= 4.9
#CPP=g++ 
CPP=g++ -fdiagnostics-color=auto
SHELL=/usr/bin/bash

CPPFLAGS=-Wall -Wpedantic -std=c++11 -pg -Og -march=native
#CPPFLAGS=-Wall -Wpedantic -std=c++11 -O2 -march=native

LDFLAGS=
LDLIBS=-lc -lm

COMMON_CPP= matrix_coo.cpp\
				matrix_crs.cpp\
				utils.cpp\
				model_problems.o\
				classical_solvers.cpp\
				multigrid.cpp

COMMON_HPP= $(COMMON_CPP:%.cpp=%.hpp)
COMMON_O  = $(COMMON_CPP:%.cpp=%.o)

.PHONY: all clean run prog02_fig prog03_fig prog04_fig

all: test run

run: test
	./test

# Test
test: $(COMMON_O) test.cpp
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c test.cpp
	$(CPP) $(CPPFLAGS) -flto $(LDLIBS) $(LDFLAGS) $(COMMON_O) test.o -o test


# Common files
matrix_coo.o: matrix_coo.cpp matrix_crs.cpp utils.cpp
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c matrix_coo.cpp

matrix_crs.o: matrix_coo.cpp matrix_crs.cpp utils.cpp
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c matrix_crs.cpp

utils.o: utils.cpp utils.hpp
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c utils.cpp

model_problems.o: matrix_coo.cpp matrix_crs.cpp utils.cpp model_problems.cpp model_problems.hpp
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c model_problems.cpp

classical_solvers.o: matrix_coo.cpp matrix_crs.cpp utils.cpp classical_solvers.cpp classical_solvers.hpp
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c classical_solvers.cpp

multigrid.o: matrix_crs.cpp utils.cpp classical_solvers.cpp multigrid.cpp multigrid.hpp
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c multigrid.cpp


# programming assignments
prog02: $(COMMON_O) prog02.cpp
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c prog02.cpp
	$(CPP) $(CPPFLAGS) -flto $(LDLIBS) $(LDFLAGS) matrix_coo.o\
		matrix_crs.o utils.o model_problems.o classical_solvers.o \
		prog02.o -o prog02

prog02_fig: prog02
	./prog02
	make -C figures/ prog02


prog03: $(COMMON_O) prog03.cpp
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c prog03.cpp
	$(CPP) $(CPPFLAGS) -flto $(LDLIBS) $(LDFLAGS) $(COMMON_O) prog03.o -o prog03

prog03_fig: prog03
	./prog03
	make -C figures/ prog03


prog04: $(COMMON_O) prog04.cpp
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c prog04.cpp
	$(CPP) $(CPPFLAGS) -flto $(LDLIBS) $(LDFLAGS) $(COMMON_O) prog04.o -o prog04

prog04_fig: prog04
	./prog04
	make -C figures/ prog04


force:

clean:
	@rm -f *.o
	@rm -f test prog02 prog03
	@rm -f *.out # gprof

cleandata:
	make -C figures/ cleandata

deepclean: clean cleandata
	make -C figures/ deepclean

