#ifndef _DILUTE_DILUTE_HH_
#define _DILUTE_DILUTE_HH_
#include "mathtools.hh"
#include "pot1dquad.hh"
#include <dlfcn.h>

class LoadedPotential : public dlt::FuncDeriv1D {
protected:
  std::string libpath_;
  void *plib_;
  double (*pvalue_)(double);
  double (*pderivative_)(double);

public:
  LoadedPotential(std::string libpath)
      : libpath_(libpath), plib_(nullptr), pvalue_(nullptr),
        pderivative_(nullptr) {
    if ((plib_ = dlopen(libpath.c_str(), RTLD_LAZY))) {
      std::cout << "Open user-defined potential function as "
                << static_cast<void *>(plib_) << std::endl;
      if ((pvalue_ =
               reinterpret_cast<double (*)(double)>(dlsym(plib_, "value")))) {
        pderivative_ =
            reinterpret_cast<double (*)(double)>(dlsym(plib_, "derivative"));
      }
    }
    if (pvalue_ == nullptr) {
      std::cerr
          << "Fail to open user-defined potential function. Library path is "
          << libpath << ".\n"
          << "Use abstract path to solve this problem." << std::endl;
    }
  }
  ~LoadedPotential() { dlclose(plib_); }
  double value(double r) const override { return pvalue_(r); }
  double derivative(double r) const override {
    if (provide_derivative()) {
      return pderivative_(r);
    } else {
      return this->dlt::FuncDeriv1D::derivative(r);
    }
  }
  bool provide_derivative() const override { return pderivative_ != nullptr; }
};
#endif
