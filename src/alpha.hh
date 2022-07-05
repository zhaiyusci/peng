#ifndef _DILUTE_ALPHA_HH_
#define _DILUTE_ALPHA_HH_
#include "global.hh"
#include "param.hh"
#include <iostream>
#include <tuple>
#include <vector>

extern std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
alpha(double t, double xs, std::vector<double> Omega00,
      std::vector<double> Omega01, std::vector<double> Omega11, int maxpq);

extern void calcalpha();

#endif
