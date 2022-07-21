FFLAGS=-O3 -fcheck=all -Wall -Werror -g -fPIC --std=f95
CC=gcc
CXX=g++
FC=gfortran
CPPFLAGS=-O3 -Wall -Werror -g -fPIC --std=c++17 
.PHONY: all clean

all: arar.so

clean:
	rm -rf *.o *.so

arar.so: arar_impl.o arar.o
	$(CC) $^ -lstdc++ -lgfortran -shared -o $@

%.o: %.cc
	$(CXX) $(CPPFLAGS) -c $< -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

