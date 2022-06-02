# Build the potential as a C libaray
# gfortran pot11.f90 -c -fPIC -shared
# g++ pot1d.fortran.cc -c -fPIC
# g++ pot1d.fortran.o pot11.o -lgfortran -fPIC -shared -o pot11.so

# The test program
# g++ test.cc atompair.cc -I../external/include --std=c++17 -ldl -g

# Potential 1D Optimization
g++ -O3 -Wall pot1d.cc toms424.cc -I ../external/include -L ../external/lib -lnlopt -Wl,-rpath,"../external/lib"

