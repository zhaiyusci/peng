#ifndef __DILUTE_TOMS424_HH__
#define __DILUTE_TOMS424_HH__
#include <functional>
#include <tuple>
#include <vector>

extern std::tuple<double, double, size_t, std::vector<double>>
ccquad(const std::function<double(double)> &f, double a, double b,
       double tolerr, size_t limit);

#endif
