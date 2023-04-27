#include <fmt/core.h>
#include <peng.hh>

class LJ : public peng::FuncDeriv1D {
public:
  /* Some local cache */
  mutable double rm6 = 0.0;
  mutable double oldr = -1.0;

  void update_r(double r) const {
    if (oldr != r) {
      rm6 = pow(r, -6);
      oldr = r;
    }
  }

  LJ() : FuncDeriv1D() {}
  /* V(r)=4(r^-12-r^-6) */
  double value(double r) const {
    update_r(r);
    return 4 * (rm6 * rm6 - rm6);
  }

  /* V'(r)=4(-12*r^-13+6r^-7) */
  double derivative(double r) const {
    update_r(r);
    return 4 * (-12 * rm6 * rm6 / r + 6 * rm6 / r);
  }
};

int main(){
  LJ lj;
  peng::Pot1DFeatures pf(lj);
  peng::ReducedPotentialQuadrature rpq(pf.reduced_potential());
  rpq.set_algorithm<peng::ChiCG, peng::QCG, peng::OmegaGL>();
  for (double T = 1.0; T <= 100.0; T += 1.0) {
    std::cout << fmt::format("{:6.1f}", T) << ' ';
    for (int l = 1; l <= 8; ++l) {
      for (int s = l; s <= 8; ++s) {
        std::cout << fmt::format("{:10.5f}", rpq.Omega(l, s, T, 1.0e-5)) << ' ';
      }
    }
    std::cout << '\n';
  }

  return 0;
}
