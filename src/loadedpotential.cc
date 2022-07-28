#include "loadedpotential.hh"
#include <dlfcn.h>
#include <string>
namespace peng{
LoadedPotential::LoadedPotential(std::string libpath)
    : libpath_(libpath), plib_(nullptr), pvalue_(nullptr),
      pderivative_(nullptr) {
  if ((plib_ = dlopen(libpath.c_str(), RTLD_LAZY))) {
    std::cerr << "Open user-defined potential function as "
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
    std::cout
        << "Fail to open user-defined potential function. Library path is "
        << libpath << ".\n"
        << "Use abstract path to solve this problem." << std::endl;
  }
}
double LoadedPotential::derivative(double r) const {
  if (provide_derivative()) {
    return pderivative_(r);
  } else {
    return this->peng::FuncDeriv1D::derivative(r);
  }
}
}
