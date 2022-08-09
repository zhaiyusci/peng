# <span style="font-variant: small-caps;">Peng</span>: A software package for the thermophysics of dilute gases.

Authors: [Yu Zhai](https://www.zhaiyusci.net/), You Li, [Hui Li](https://huiligroup.org/), and Frederick R. W. McCourt

## Introduction

Thermophysics of gases is an important field in chemical physics, 
and computation is one of the most important technologies in the research therein,
just like in other fields.
No need to say, software packages play roles in the computational research.

However, many codes used in this field are out-of-date not in the sense of science but in programming.
FORTRAN programs were written in around 1960s to 1980s, when most programs were designed specifically for a type of machine.
We struggled to make them work in new machines, and often some modification is required simply because we used different compilers.
We actually need a software which is written following the standard and can work out-of-the-box.

In this repo we present the **Platform of ENergetic Gasses**, <span style="font-variant: small-caps; font-weight: bold;">Peng</span>, 
which compute the collision integrals and the transport properties of binary dilute gases.
By now, we focus on the case where both kinds of particles are atoms, 
which limits the code to of rare gases mixtures.
However, you **can** test your effective potential energy curve for polyatomic molecules with the help of this project.

## Build

To build the project, simply run `make` in the root directory.
On a machine with Internet connection, 
the building script will download and build the required libraries automatically, 
and then compile the dilute project.

If Internet is not available, you can follow `external/prerequest.sh` and build the prerequests manually.

Please find the executable and libraries in `build` directory.  
We do not provided an installation script because some functions of the project are still under construction.
You can call the executable with absolute path.

## Usage

### The potential energy curves

User should also provide the potential energy curves (PECs) of the three interatomic interaction.
Due to [name mangling](https://en.wikipedia.org/wiki/Name_mangling), 
the lite way to take user provided PECs must be using C interface.
Most programming languages has a C interface, 
and the user can load the PECs at runtime.

At least one function providing the value of PEC should be written, of which the signature is (in C)
```c
double value(double r);
```
for a real system, `r` should be in Angstrom and the return value in Kelvin.

It is often the case where the analytical derivative of the PEC can be provided.
The API is like
```c
double derivative(double r);
```
Do **not** provide a wrong derivative.
Perhaps you do not have the spare time to do the maths to give the analytical derivative. 
Just do not provide it.
If the program did not find the `derivative` function, it will compute the numerical derivative instead.

Also note that, 
we always evaluate of the derivative is always clung to the one of value, 
i.e., instead of
```cpp
for (double r = 3.0; r <= 10.0; r += 0.1) {
  double v = value(r);
  // ... something else 0
}
for (double r = 3.0; r <= 10.0; r += 0.1) {
  double dv = derivative(r);
  // ... something else 1
}
```
we prefer
```cpp
for (double r = 3.0; r <= 10.0; r += 0.1) {
  double v = value(r);
  double dv = derivative(r);
  // ... something else
}
```
It means you can do some optimization by storing intermediate variables
to accelerate the computation of derivatives right after the evaluation of value, 
or vice versa.

Historically, people provide potential energy as FORTRAN subroutines.
Users should note the FORTRAN subroutine parameters are actually pointers in C-like languages,
in our simple case.
Thus, a simple interface for "f2c" can be written following the provided examples.
In the example provided, 
we write the working part of PEC in FORTRAN, while write the C interface in C++.

Compile the PEC as dynamic shared libraries (on GNU/Linux, `.so` file).
In the simplest case, we have `pec00.c` file, which contains the `value` function.
We can build the dynamic shared library as
```sh
gcc -o pec00.so -shared -fPIC pec00.c
```
We will have `pec00.so` in the working directory if everything is good.
Put the path in the JSON file.

### The JSON input file

Input files of <span style="font-variant: small-caps;">Peng</span> executable obey the [JSON format](https://en.wikipedia.org/wiki/JSON).
The users can find examples of input file in `examples/hexe/hexe.json`.

Although the sample JSON file is self-explained, we list the keys below.
If unmentioned keys are written in the input file, 
they will **not** affect the computation, thus, can be used as comments.
- `atoms`: Array of the two types of atoms, name is not really used, and mass in atomic mass unit [amu, 1/12 m(C-12)].
- `potentials`: Array of the PECs, in the order of interactions between `atoms[0]`-`atoms[0]`, `atoms[0]`-`atoms[1]`, and `atoms[1]`-`atoms[1]`. 
- `accuracy`: Maximum of allowed integration relative error.
- `temperatures`: Array of the temperatures at which the properties are computed.
- `molefractions0`: Array of the mole fractions of `atoms[0]`, and the mole fraction of `atoms[1]` is computed accordingly.
- `propertyorder`: The maximum order when compute the transport properties.


This project is supported by
- 2020-JCJQ project (GFJQ2126-007)
- National Natural Science Foundation of China (22073035)
- The Program for JLU Computational Interdiscipline Innovative Platform
