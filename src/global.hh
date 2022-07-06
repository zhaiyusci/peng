#ifndef __GLOBAL_H__
#define __GLOBAL_H__
#include "param.hh"
#include <iostream>
#include <string>
#include <vector>
// extern std::vector<double> temperatures;
// extern double temps[60];
// extern double x[2];
// extern std::vector<std::string> elements;
// extern std::vector< std::vector<double> > xs;
// extern int ntemp;
// extern double acc;
// extern double om11[60][MAXORD][MAXORD];
// extern double om12[60][MAXORD][MAXORD];
// extern double om22[60][MAXORD][MAXORD];
// extern std::vector<std::vector<double> > etas;
// extern std::vector<std::vector<double> > D12s;
// extern std::vector<std::vector<double> > DTs;
// extern std::vector<std::vector<double> > lambdas;

inline void omega_cpp2fort(std::vector<double> Omega, double om[MAXORD][MAXORD],
                           size_t omegaorder) {
  size_t i = 0;
  for (size_t l = 0; l != omegaorder; ++l) {
    for (size_t s = l; s != omegaorder; ++s) {
      om[s][l] = Omega[i];
      ++i;
    }
  }
}
inline size_t omegaorder(size_t maxpq) { return 2 * maxpq + 2; }
#endif // __GLOBAL_H__
