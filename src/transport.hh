#ifndef _DILUTE_ALPHA_HH_
#define _DILUTE_ALPHA_HH_
#include "global.hh"
#include <iostream>
#include <tuple>
#include <vector>

namespace dlt {
extern std::tuple<std::vector<double> /*D12*/, std::vector<double> /*DT*/,
                  std::vector<double> /*lambda*/, std::vector<double> /*eta*/>
transport(double t, double x0, std::vector<double> Omega00,
          std::vector<double> Omega01, std::vector<double> Omega11,
          double mass0, double mass1, int propertyorder);
}

#endif
