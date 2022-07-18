#ifndef _DILUTE_TRANSPORT_HH_
#define _DILUTE_TRANSPORT_HH_
#include <iostream>
#include <tuple>
#include <vector>

/**
 * @brief The top level namespace for the Dilute project.
 */
namespace dlt {
/**
 * @brief Compute the transport properties based on computed Omega's.
 *
 * @param T: temperature in Kelvin.
 * @param x0: mole fraction of the first (zeroth, as in C++) specie.
 * @param Omega00, Omega01, Omega11: computed Omega, stored in (l,s) = (1,1),
 * (1,2), (1,3), ... , (max, max).
 * @param mass0, mass1: masses in amu.
 * @param propertyorder: highest order of the property computation in
 * Chapman-Enskog solution.
 */
extern std::tuple<std::vector<double> /*D12*/, std::vector<double> /*DT*/,
                  std::vector<double> /*lambda*/, std::vector<double> /*eta*/>
transport(double T, double x0, std::vector<double> Omega00,
          std::vector<double> Omega01, std::vector<double> Omega11,
          double mass0, double mass1, int propertyorder);
}

#endif
