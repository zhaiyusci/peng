#ifndef _PENG_LOADEDPOTENTIAL_HH_
#define _PENG_LOADEDPOTENTIAL_HH_
#include "mathtools.hh"
#include <dlfcn.h>

namespace peng {
/**
 * @brief Load external potential as FuncDeriv1D.
 */
class LoadedPotential : public peng::FuncDeriv1D {
protected:
  std::string libpath_;
  void *plib_;
  double (*pvalue_)(double);
  double (*pderivative_)(double);

public:
  /**
   * @brief Constructor.
   *
   * @param libpath: The path of the library in which the required functions are
   * located.
   */
  LoadedPotential(std::string libpath);
  ~LoadedPotential() { dlclose(plib_); }
  double value(double r) const override { return pvalue_(r); }
  bool provide_derivative() const override { return pderivative_ != nullptr; }
  double derivative(double r) const override;
};
} // namespace peng
#endif
