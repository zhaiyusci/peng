#ifndef _DILUTE_CACHEDFUNC_HH_
#define _DILUTE_CACHEDFUNC_HH_
#include "mathtools.hh"
#include <functional>

namespace dlt {

class CachedFuncDeriv1D : public FuncDeriv1D {
public:
  FuncDeriv1D *const pfunc_;
  mutable std::vector<std::tuple<double, double, double, bool>> cache_;
  const double ftol_;

  std::tuple<double, double, bool, size_t> cubic_spline_(double x) const;

public:
  CachedFuncDeriv1D(FuncDeriv1D &func, double ftol = 1.0e-8)
      : pfunc_(&func), ftol_(ftol) {}
  double value(double x) const override;
  double derivative(double x) const override;
  // We always set this true because we can get it from interpolation.
  bool provide_derivative() const override { return true; }
  void add_to_cache(const std::vector<double> &xs) const;
};
} // namespace dlt

#endif
