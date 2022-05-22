# Dilute: A software package for the thermophysics of dilute gases.

Thermophysics of gases is an important field in chemical physics, 
and computation is one of the most important technology in the research therein,
just like in other fields.
No need to say, software packages play roles in the computational research.

However, many codes used in this field is out-of-date not in the sense of science but in programming.
FORTRAN programs were written in around 1960s to 1980s, when most programs were designed specifically for a type of mechine.
We actually need a software which can work out-of-the-box.

In this repo I present a version of code provided by H. O'Hara and F. J. Smith, 
which was published with the paper "Transport collision integrals for a dilute gas" (*Comp. Phys. Comm.* **2**, 47-54 (1971), [doi:10.1016/0010-4655(71)90014-2](https://doi.org/10.1016/0010-4655(71)90014-2) ).

In this contribution, the code is cleaned and is made possible to used in a modern mechine.

The files are listed below
1. `README.md`: This file.
2. `acqn.F`: The source code for the ACQN program, which can be built with `gfortran 10.2.0` and `9.3.0`, which, in 2022, can be thought as "modern" compilers.
3. `fort.1`: The example input file.
4. `fort.2`: The corresponding output file, which agree with the original work.
5. `Makefile`: Makefile for a typical building.

