FFLAGS=-O3 -fcheck=all -Wall -Werror -g -fPIC --std=f95
CC=gcc
CXX=g++
FC=gfortran
CPPFLAGS=-O3 -Wall -Werror -g -fPIC --std=c++17 
.PHONY: all clean

all: hehe.so hexe.so xexe.so

clean:
	rm -rf *.o *.so

hehe.so: hehe_impl.o hehe.o
	$(CC) $^ -lstdc++ -lgfortran -shared -o $@

hexe.so: hexe_impl.o hexe.o
	$(CC) $^ -lstdc++ -lgfortran -shared -o $@

xexe.so: xexe_impl.o xexe.o
	$(CC) $^ -lstdc++ -lgfortran -shared -o $@

%.o: %.cc
	$(CXX) $(CPPFLAGS) -c $< -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

