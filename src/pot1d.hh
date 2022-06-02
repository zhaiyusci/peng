#ifndef __DILUTE_POT1D_HH__
#define __DILUTE_POT1D_HH__
#include <DataFrame/DataFrame.h>
#include <functional>
/**
 * The 1-dimension function with its 1st order derivative.
 */
typedef std::function<std::tuple<double, double>(double)> FuncDeriv1D;

class CachedFuncDeriv1D {
  public:
    const std::function<std::tuple<double, double>(double)> *const func_;
    mutable std::vector<std::tuple<double, double, double, bool>> cache_;
    const double ftol_;

    std::tuple<double, double, bool, size_t> cubic_spline_(double x) const;

  public:
    CachedFuncDeriv1D(
        const std::function<std::tuple<double, double>(double)> &func,
        double ftol = 1.0e-8)
        : func_(&func), ftol_(ftol) {}
    std::tuple<double, double> operator()(double x) const;
    void add_to_cache(const std::vector<double> &xs) const;
};

#endif
