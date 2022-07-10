#ifndef __GLOBAL_H__
#define __GLOBAL_H__
#include "param.hh"
#include <iostream>
#include <string>
#include <vector>

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
