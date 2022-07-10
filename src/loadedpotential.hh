#ifndef _DILUTE_LOADEDPOTENTIAL_HH_
#define _DILUTE_LOADEDPOTENTIAL_HH_
#include "mathtools.hh"
#include <dlfcn.h>

namespace dlt {
class LoadedPotential : public dlt::FuncDeriv1D {
protected:
  std::string libpath_;
  void *plib_;
  double (*pvalue_)(double);
  double (*pderivative_)(double);

public:
  LoadedPotential(std::string libpath);
  ~LoadedPotential() { dlclose(plib_); }
  double value(double r) const override { return pvalue_(r); }
  bool provide_derivative() const override { return pderivative_ != nullptr; }
  double derivative(double r) const override;
};
} // namespace dlt
#endif
