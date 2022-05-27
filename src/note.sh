gfortran pot11.f90 -c -fPIC -shared
g++ pot1d.cc -c -fPIC
g++ pot1d.o pot11.o -lgfortran -fPIC -shared -o pot11.so
g++ test.cc -I../external/include --std=c++17 -ldl -g

