#include "atompair.hh"
#include <dlfcn.h>
#include <iostream>

auto getpot(std::string soname) {
  void *so = dlopen(soname.c_str(), RTLD_NOW);
  std::cerr << "Open user-defined dimer as " << static_cast<void *>(so)
            << std::endl;
  auto pot = reinterpret_cast<void (*)(double, double *, double *, double *)>(
      dlsym(so, "pot"));
  return [=](double r) -> std::tuple<double, double> {
    double v, dv, d2v;
    pot(r, &v, &dv, &d2v);
    return std::make_tuple(v, dv);
  };
  // for (double r = 1.0; r <= 10.0; ++r) {
  // double v, dv, d2v;
  // std::tie(v, dv, d2v) = pot(r);
  // std::cout << r << "\t" << v << std::endl;
  // }
  // dlclose(so);
  // return pot;
}

int main() {
  Particle Ar("Ar", 40.0);
  Particle Xe("Xe", 128.0);
  PotentialLibrary::getInstance()[std::make_pair(Ar, Xe)] =
      getpot("./pot11.so");
  for (double r = 1.0; r <= 10.0; ++r) {
    double v, dv;
    std::tie(v, dv) =
        PotentialLibrary::getInstance()[std::make_pair(Ar, Xe)](r);
    std::cout << r << "\t" << v << std::endl;
  }
  return 0;
}
