#ifndef __GLOBAL_H__
#define __GLOBAL_H__
#include "param.hh"
#include <string>
#include <vector>
extern std::vector<double> temperatures;
extern double temps[60];
extern double x[2]; 
extern std::vector<std::string> elements;
extern std::vector< std::vector<double> > xs;
extern int ntemp;
extern double acc;
extern int maxpq;
extern double om11[60][MAXORD][MAXORD];
extern double om12[60][MAXORD][MAXORD];
extern double om22[60][MAXORD][MAXORD];
extern std::vector<std::vector<double> > etas;
extern std::vector<std::vector<double> > D12s;   
extern std::vector<std::vector<double> > DTs;    
extern std::vector<std::vector<double> > lambdas;
#endif // __GLOBAL_H__
