#include "pot1dquad.hh"

namespace dlt {

ReducedPotentialQuadrature::PhiEff::PhiEff(FuncDeriv1D &pot, double b, double E)
    : ppot_(&pot), b_(b), E_(E) {}
double ReducedPotentialQuadrature::PhiEff::value(double r) const {
  double res;
  res = ppot_->value(r);
  res += b_ * b_ * E_ / r / r;
  return res;
}
double ReducedPotentialQuadrature::PhiEff::derivative(double r) const {
  double res;
  res = ppot_->derivative(r);
  res += -2 * b_ * b_ * E_ / r / r / r;
  return res;
}

double ReducedPotentialQuadrature::Y::value(double r) const {
  double res;
  res = ppot_->value(r);
  res += 0.5 * r * ppot_->derivative(r);
  return res;
}

double ReducedPotentialQuadrature::chi(double E, double r_m, double rtol) {
  return p_chi_impl_->chi(E, r_m, rtol);
}

ReducedPotentialQuadrature::ReducedPotentialQuadrature(FuncDeriv1D &reduced_pot)
    : p_reduced_pot_(&reduced_pot), p_chi_impl_(nullptr), p_q_impl_(nullptr), p_omega_impl_(nullptr) {

  y_.reset(new Y(*p_reduced_pot_));

  std::tie(r_C_, E_C_) = find_local_maximum(*y_, 0.5, 5);
  std::cerr << "r_C = " << r_C_ << '\n';
  std::cerr << "E_C = " << E_C_ << '\n';

  y_root_.reset(new LocalRoot(*y_, r_C_, 15.0, 21, 1.0e-12));
  v_root_.reset(new LocalRoot(*p_reduced_pot_, 0.1, 1.0, 21, 1.0e-12));
  // Maybe we want to try (0.0,1.0)... but consider the really high energy part
  // is rarely used...
}



std::tuple<double, double> ReducedPotentialQuadrature::r_range(double E) const {
  double r_O, r_Op, b_O;
  if (E >= E_C_) {
    r_O = r_C_;
    r_Op = r_C_;
    b_O = r2b(r_O, E);
  } else {
    r_O = (*y_root_)(E);
    b_O = r2b(r_O, E);
    PhiEff phieff(*p_reduced_pot_, b_O, E);
    r_Op = find_local_root(phieff, E, 1.0, r_C_, 21, 1.0e-12);
    if (fabs(phieff(r_Op) - E) >= 1.0e-6) {
      double _;
      std::tie(r_Op, _) = find_local_maximum(phieff, 1.0, r_C_);
      r_Op = find_local_root(phieff, E, 0.1, r_Op, 21, 1.0e-12);
    }
  }
  return std::make_tuple(r_O, r_Op);
}

double ReducedPotentialQuadrature::r2b(double r, double E) const {
  double v = p_reduced_pot_->value(r);
  if (r > 1.0e8) {
    v = 0.0;
  }
  double b = r * sqrt((E - v) / E);
  return b;
}

/** Compute Q.
 *
 * It is faster to keep the r_E unchanged and scan the l.
 */
double ReducedPotentialQuadrature::Q(size_t l, double r_E, double E,
                                     double rtol) {
  return p_q_impl_->Q(l, r_E, E, rtol);
}

///
/// This class inheriated from the Chebyshev-Gauss Quarduture class, compute the
/// integration required by the computation of Omega.
///

double ReducedPotentialQuadrature::Omega(size_t l, size_t s, double T,
                                         double rtol) {
  return p_omega_impl_->Omega(l,s,T,rtol);
}

} // namespace dlt
