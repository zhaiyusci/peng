
// Define the global variables
#include"global.h"

std::vector<double> temperatures;
double temps[60]={}; // TTT in GETOMEGA in omega.F
double x[2]={}; 
std::vector<std::string> elements;
std::vector< std::vector<double> > xs;
int ntemp=0; 
double acc=1.;
int maxpq=1;
// Arrays which take result from omega.F 
// Dimensions here should be same as those in fortran file
double om11[60][MAXORD][MAXORD];
double om12[60][MAXORD][MAXORD];
double om22[60][MAXORD][MAXORD];
// Vectors takes the results from alpha.F90 and beta.F90
// Every 1D vector in the 2D vectors 
// is the result for the same system (mole fractions) and temperature
// but with different order of accuracy
// The outer index, is ordered with the following snippet
//
// for(int i = 0; i!=ntemp; ++i){
//   for(auto&& x : xs){
//   ...
//   }
// }
//
std::vector<std::vector<double> > etas;
std::vector<std::vector<double> > D12s;
std::vector<std::vector<double> > DTs;
std::vector<std::vector<double> > lambdas;
