FFLAGS=-O3 -fcheck=all -Wall -Werror -g -fPIC 
CC=gcc
CXX=g++
FC=gfortran
CPPFLAGS=-O3 -Wall -Werror -g -fPIC --std=c++11 -I../external/include
OBJECTS=omega.o pec.o 
OBJDIR ?= ../build
SRC := alpha.cc beta.cc global.cc main.cc \
			 omega.F \
			 alpha_impl.F90 beta_impl.F90 inversion.F90 pec.F90 \
			 userdefpec11.F90 userdefpec12.F90 userdefpec22.F90
OBJ := $(SRC)
OBJ := $(OBJ:%.cc=%.o)
OBJ := $(OBJ:%.F=%.o)
OBJ := $(OBJ:%.F90=%.o)

VPATH = ../src

.PHONY: clean all

all: main.exe omega11.so omega12.so omega22.so 

test:
	./main.exe > tmp
	diff l answer || echo ""

main.o omega.o alpha.o alpha_impl.o beta.o beta_impl.o: param.h
main.o global.o alpha.o beta.o: global.h

main.exe: main.o beta.o inversion.o global.o beta_impl.o alpha.o alpha_impl.o
	$(CC) -o $@ $^ -ldl -lgfortran -lstdc++ -lm -L../external/lib -lfmt 

omega11.so: $(OBJECTS) userdefpec11.o 
	$(FC) $(FFLAGS) -shared -o $@ $^

omega12.so: $(OBJECTS) userdefpec12.o 
	$(FC) $(FFLAGS) -shared -o $@ $^

omega22.so: $(OBJECTS) userdefpec22.o 
	$(FC) $(FFLAGS) -shared -o $@ $^


clean:
	rm -f *.exe *.o *.so fort.*

%.o: %.cc
	$(CXX) $(CPPFLAGS) -c $< -o $@

%.o: %.F
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.F90
	$(FC) $(FFLAGS) -c $< -o $@

