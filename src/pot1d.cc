#include "atompair.hh"
#include <functional>
#include <iostream>
#include <nlopt.hpp>
#include <tuple>

auto lj(double r) {
  double r_6(pow(r, -6));
  double v(4 * (r_6 * r_6 - r_6));
  double dv(4 * ((-12) * r_6 * r_6 / r - (-6) * r_6 / r));
  double d2v(4 * ((12 * 13) * r_6 * r_6 / r / r - (6 * 7) * r_6 / r / r));

  return std::make_tuple(v, dv, d2v);
}


double mlopt_wrapper(const std::vector<double> &x, std::vector<double> &grad,
                     void *f_data) {
  Pot1d pot = *reinterpret_cast<Pot1d *>(f_data);
  double r(x[0]);
  double v, dv, d2v;
  std::tie(v, dv, d2v) = pot(r);
  // grad.resize(1);
  // grad[0] = dv;
  return v;
}

int main() {
  nlopt::opt opt(nlopt::LN_COBYLA, 1);
  // opt.nlopt_set_min_objective(mlopt_wrapper, &lj);
  Pot1d pot(lj);
  opt.set_min_objective(mlopt_wrapper, &pot);
  opt.set_xtol_rel(1.0e-8);
  std::vector<double> r0{0.5};
  double f;
  // std::cerr << __LINE__ << std::endl;
  opt.optimize(r0, f);
  std::cout << "r0[0]"
            << "\t"
            << "f" << std::endl;
  std::cout << r0[0] << "\t" << f << std::endl;
  return 0;
}
