#include "pot1dquad.hh"
#define CG_MAXORDER 10

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

ReducedPotentialQuadrature::ChiCG::ChiCG(ReducedPotentialQuadrature *rpq)
    : CGIntegrator(true), rpq_(rpq) {
  clean_cache();
  E_ = -2.0;
  r_m_ = -1.0;
}

void ReducedPotentialQuadrature::ChiCG::set_param(double r_m, double E,
                                                  double b) {
  if (E_ != E) {
    clean_workspace();
  }
  if (r_m_ != r_m) {
    clean_cache();
  }
  if (E_ != E || r_m_ != r_m) {
    if (b < 0) {
      b_ = rpq_->r2b(r_m, E);
    } else {
      b_ = b;
    }
  }
  E_ = E;
  r_m_ = r_m;
}

/** Compute the integrands with computed values cached.
 */
void ReducedPotentialQuadrature::ChiCG::calculate_integrands(size_t ordersize,
                                                             double rtol) {
  // See if we need an update based on the "flag".
  if (cache_ordersize_ < ordersize) {
    // We only need the positive half of the integration.
    CubicIter ci(cache_ordersize_, ordersize, true, true);
    // Reserve the memory... which help std::vector work faster
    vs_.reserve(ci.size_from_0());
    for (auto &&i : ci) {
      double y = map_pm1(CGIntegratorBackend::instance()->coss(i));
      double r = r_m_ / y;
      if (!(r >= 1.0e-5)) {
        r = 1.0e-5;
      }
      double v = (*rpq_->p_reduced_pot_)(r);
      vs_.push_back(v);
    }
    cache_ordersize_ = ordersize;
  }
  if (ordersize_ < ordersize) {
    // We only need the positive half of the integration.
    CubicIter ci(ordersize_, ordersize, true, true);
    // Reserve the memory... which help std::vector work faster
    integrands_.reserve(ci.size_from_0());
    for (auto &&i : ci) {
      double y = map_pm1(CGIntegratorBackend::instance()->coss(i));
      double r = r_m_ / y;
      // double v = (*rpq_->ppot_)(r);
      double F = (1.0 - vs_[i] / E_ - b_ * b_ / r / r);
      if (y <= 1.0e-8) { // r --> inf
        F = 1.0;
      }
      if (!(F >= 1.0e-18)) {
        // should never run here if Chebyshev-Gauss Quarduture is used
        // because y does not equal 1 in CG Quad
        // std::cerr << " ~~F~~ " << F << " ~~r~~ " << r << std::endl;
        // throw std::runtime_error(__FILE__ "  F < 0");
        F = 1.0e-18;
      }
      double res = 1.0 / sqrt(F);
      integrands_.push_back(res * sqrt(1.0 - y * y));
    }
    // Update the "flag".
    ordersize_ = ordersize;
  }

  return;
}

double ReducedPotentialQuadrature::chi(double E, double r_m, double rtol) {
  // double E = ppot_->value(r_E);
  double b = r2b(r_m, E);

  chicg.set_param(r_m, E, b);
  double quadrature, err;
  bool converged;
  std::tie(quadrature, err, converged) = chicg.integrate(rtol, CG_MAXORDER);
  if (!converged) {
    std::cerr << "Line " << __LINE__ << " ChiCG not converged with E = " << E
              << ", r_m = " << r_m << " and b = " << b << "." << '\n';
    std::cerr << "quadrature   " << quadrature << ' ' << M_PI - b * quadrature
              << '\n';
    // chicg.show_integrands();
  }
  return M_PI - b / r_m * quadrature;
}

ReducedPotentialQuadrature::ReducedPotentialQuadrature(FuncDeriv1D &reduced_pot)
    : p_reduced_pot_(&reduced_pot), chicg(this), qcg1_(this), qcg2_(this),
      omegacg_(this), omegagl_(this) {

  y_.reset(new Y(*p_reduced_pot_));

  std::tie(r_C_, E_C_) = find_local_maximum(*y_, 0.5, 5);
  std::cerr << "r_C = " << r_C_ << '\n';
  std::cerr << "E_C = " << E_C_ << '\n';

  y_root_.reset(new LocalRoot(*y_, r_C_, 15.0, 21, 1.0e-12));
  v_root_.reset(new LocalRoot(*p_reduced_pot_, 0.1, 1.0, 21, 1.0e-12));
  // Maybe we want to try (0.0,1.0)... but consider the really high energy part
  // is really used...
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
  // std::cerr << "E = " << E << ", r_O = " << r_O << ", r_Op = " << r_Op <<
  // '\n';
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

ReducedPotentialQuadrature::QCG1::QCG1(ReducedPotentialQuadrature *rpq)
    : CGIntegrator(false), rpq_(rpq) {
  clean_cache();
  l_ = 0;
  r_E_ = -1.0;
  r_Op_ = -1.0;
}

void ReducedPotentialQuadrature::QCG1::set_param(size_t l, double r_E,
                                                 double r_Op, double E) {
  if (l_ != l) {
    clean_workspace();
    l_ = l;
  }
  if (r_E_ != r_E || r_Op_ != r_Op) {
    clean_cache();
    set_a_b(r_E, r_Op);
    r_E_ = r_E;
    r_Op_ = r_Op;
    if (E < 0.0) {
      E_ = rpq_->p_reduced_pot_->value(r_E);
    } else {
      E_ = E;
    }
  }
  return;
}

/** Compute the integrands with computed values cached.
 */
void ReducedPotentialQuadrature::QCG1::calculate_integrands(size_t ordersize,
                                                            double rtol) {
  // See if we need an update based on the "flag".
  if (cache_ordersize_ < ordersize) {
    CubicIter ci(cache_ordersize_, ordersize, false, true);
    size_t num = ci.size_from_0();
    integrands_.reserve(num);
    coschis_.reserve(num);
    fct2_.reserve(num);
    for (auto &&i : ci) {
      double y = CGIntegratorBackend::instance()->coss(i);
      double r_m = map_pm1(y);
      double v, dv;
      v = rpq_->p_reduced_pot_->value(r_m);
      dv = rpq_->p_reduced_pot_->derivative(r_m);
      double coschi = cos(rpq_->chi(E_, r_m, rtol));
      // Here I put everything else in fct2, include the weight
      double fct2 = (2.0 * (E_ - v) - r_m * dv) * r_m * sqrt(1.0 - y * y);
      coschis_.push_back(coschi);
      fct2_.push_back(fct2);
    }
    cache_ordersize_ = ordersize;
  }
  if (ordersize_ < ordersize) {
    CubicIter ci(ordersize_, ordersize, false, false);
    size_t num = ci.size_from_0();
    integrands_.reserve(num);
    for (auto &&i : ci) {
      double res = (1.0 - pow(coschis_[i], l_)) * fct2_[i];
      integrands_.push_back(res);
    }
    ordersize_ = ordersize;
  }
  return;
}

ReducedPotentialQuadrature::QCG2::QCG2(ReducedPotentialQuadrature *rpq)
    : CGIntegrator(true), rpq_(rpq) {
  clean_cache();
  l_ = 0;
  r_E_ = -1.0;
  r_O_ = -1.0;
}

/** Set the parameters, clear the inner storage if it is needed.
 */
void ReducedPotentialQuadrature::QCG2::set_param(size_t l, double r_E,
                                                 double r_O, double E) {
  if (l_ != l) {
    clean_workspace();
    l_ = l;
  }
  if (r_E_ != r_E || r_O_ != r_O) {
    clean_cache();
    r_E_ = r_E;
    r_O_ = r_O;
    if (E < 0.0) {
      E_ = rpq_->p_reduced_pot_->value(r_E);
    } else {
      E_ = E;
    }
  }
  return;
}

/** Compute the integrands with computed values cached.
 */
void ReducedPotentialQuadrature::QCG2::calculate_integrands(size_t ordersize,
                                                            double rtol) {
  // See if we need an update based on the "flag".
  if (cache_ordersize_ < ordersize) {
    CubicIter ci(cache_ordersize_, ordersize, true, true);
    size_t num = ci.size_from_0();
    integrands_.reserve(num);
    coschis_.reserve(num);
    fct2_.reserve(num);
    for (auto &&i : ci) {
      double y = CGIntegratorBackend::instance()->coss(i);
      if (y < 1.0e-8) {
        y = 1.0e-8; // prevent div by 0
      }
      double r_m = r_O_ / y;
      double v, dv;
      v = rpq_->p_reduced_pot_->value(r_m);
      dv = rpq_->p_reduced_pot_->derivative(r_m);
      double coschi = cos(rpq_->chi(E_, r_m, rtol));
      // Here I put everything else in fct2, include the weight
      double fct2 =
          (2.0 * (E_ - v) - r_m * dv) * r_m * r_m / y * sqrt(1.0 - y * y);
      coschis_.push_back(coschi);
      fct2_.push_back(fct2);
    }
    cache_ordersize_ = ordersize;
  }
  if (ordersize_ < ordersize) {
    CubicIter ci(ordersize_, ordersize, true, true);
    size_t num = ci.size_from_0();
    integrands_.reserve(num);
    for (auto &&i : ci) {
      double res = (1.0 - pow(coschis_[i], l_)) * fct2_[i];
      integrands_.push_back(res);
    }
    ordersize_ = ordersize;
  }
  return;
}

/** Compute Q.
 *
 * It is faster to keep the r_E unchanged and scan the l.
 */
double ReducedPotentialQuadrature::Q(size_t l, double r_E, double E,
                                     double rtol) {
  // double old_r_E = 0.0;
  double r_O, r_Op;

  if (E < 0.0) {
    E = p_reduced_pot_->value(r_E); // This is common
  } else if (r_E < 0.0) {
    r_E = (*v_root_)(E);
  }
  std::tie(r_O, r_Op) = r_range(E);
  // std::cerr << "(r_Op , r_O) = " << r_Op << ' ' << r_O << std::endl;

  double coeff = 1.0 / (1.0 - (1.0 + pow(-1, l)) / 2.0 / (1.0 + l)) / E;
  double esterr;
  bool converged;
  // std::cerr << "l = " << l << ", E = " << E << std::endl;
  if (E <= 2 * E_C()) {
    double quadrature1;
    qcg1_.set_param(l, r_E, r_Op, E);
    std::tie(quadrature1, esterr, converged) =
        qcg1_.integrate(rtol, CG_MAXORDER);
    if (!converged) {
      std::cerr << "Line " << __LINE__ << " QCG1 not converged with E = " << E
                << "." << std::endl;
    }
    double quadrature2;
    qcg2_.set_param(l, r_E, r_O, E);
    std::tie(quadrature2, esterr, converged) =
        qcg2_.integrate(rtol, CG_MAXORDER);
    if (!converged) {
      std::cerr << "Line " << __LINE__ << " QCG2 not converged with E = " << E
                << "." << std::endl;
    }
    quadrature2 /= 2;
    return coeff * (quadrature1 + quadrature2);
  } else {
    double quadrature2;
    qcg2_.set_param(l, r_E, r_E, E);
    std::tie(quadrature2, esterr, converged) =
        qcg2_.integrate(rtol, CG_MAXORDER);
    if (!converged) {
      std::cerr << "Line " << __LINE__ << " QCG2 not converged with E = " << E
                << "." << std::endl;
    }
    quadrature2 /= 2;
    return coeff * quadrature2;
  }
}

// Omega stuff

///
/// This class inheriated from the Chebyshev-Gauss Quarduture class, compute the
/// integration required by the computation of Omega.
///
ReducedPotentialQuadrature::OmegaCG::OmegaCG(ReducedPotentialQuadrature *rpq)
    : CGIntegrator(true, -0.5, 1.0), rpq_(rpq) {
  // The integration range should be [-1,1],
  // however, -0.5 is good enough to avoid some numerical error.
  l_ = 0;
  s_ = 0;
  T_ = 0.0;
  clean_cache_Q();
  clean_cache_v();
}
/** Set the parameters, clear the inner storage if it is needed.
 */
void ReducedPotentialQuadrature::OmegaCG::set_param(size_t l, size_t s,
                                                    double T) {
  if (s_ != s || T_ != T) {
    clean_workspace();
    s_ = s;
    T_ = T;
  }
  if (l_ != l) {
    clean_cache_Q();
    l_ = l;
  }
  return;
}

/** Compute the integrands with computed values cached.
 */
void ReducedPotentialQuadrature::OmegaCG::calculate_integrands(size_t ordersize,
                                                               double rtol) {
  // See if we need an update based on the "flag".
  if (ordersize_Q_ < ordersize) {
    CubicIter ci(ordersize_Q_, ordersize, true, true);
    size_t num = ci.size_from_0();
    Qs_.reserve(num);
    for (auto &&i : ci) {
      double y = CGIntegratorBackend::instance()->coss(i);
      if (y < 1.0e-8) {
        y = 1.0e-8; // prevent div by 0
      }
      double r_E = map_pm1(y);
      Qs_.push_back(rpq_->Q(l_, r_E, -1.0, rtol));
    }
    ordersize_Q_ = ordersize;
  }
  if (ordersize_v_ < ordersize) {
    CubicIter ci(ordersize_v_, ordersize, true, true);
    size_t num = ci.size_from_0();
    Qs_.reserve(num);
    for (auto &&i : ci) {
      double y = CGIntegratorBackend::instance()->coss(i);
      double r_E = map_pm1(y);
      if (r_E < 1.0e-8) {
        r_E = 1.0e-8; // prevent div by 0
      }

      vs_.push_back(rpq_->p_reduced_pot_->value(r_E));
      // weight is included in dvs
      dvs_.push_back(rpq_->p_reduced_pot_->derivative(r_E) * sqrt(1.0 - y * y));
    }
    ordersize_v_ = ordersize;
  }
  if (ordersize_ < ordersize) {
    CubicIter ci(ordersize_, ordersize, true, true);
    size_t num = ci.size_from_0();
    vs_.reserve(num);
    dvs_.reserve(num);

    for (auto &&i : ci) {
      double x = vs_[i] / T_;
      double res = exp(-x) * pow(x, s_ + 1) * Qs_[i] * dvs_[i];
      integrands_.push_back(res);
    }
    ordersize_ = ordersize;
  }
  return;
}

ReducedPotentialQuadrature::OmegaGL::OmegaGL(ReducedPotentialQuadrature *rpq)
    : GLIntegrator(), rpq_(rpq) {
  // The integration range should be [-1,1],
  // however, -0.5 is good enough to avoid some numerical error.
  l_ = 0;
  s_ = 0;
  T_ = 0.0;
  set_alpha(s_ + 1);
  clean_workspace();
}

/** Set the parameters, clear the inner storage if it is needed.
 */
void ReducedPotentialQuadrature::OmegaGL::set_param(size_t l, size_t s,
                                                    double T) {
  s_ = s;
  set_alpha(s + 1.0);
  T_ = T;
  l_ = l;
  return;
}

/** Compute the integrands with computed values cached.
 */
void ReducedPotentialQuadrature::OmegaGL::calculate_integrands(double rtol) {
  for (size_t i = 0; i != xs_.size(); ++i) {
    double Qrtol = ws_[i] > 2 ? rtol / 2 : rtol / ws_[i];
    double value = rpq_->Q(l_, -1.0, xs_[i] * T_, Qrtol);
    // std::cerr << "Q(" << l_ << ", " << x * T_ << ") = " << value <<
    // std::endl;
    integrands_.push_back(value);
  }
  return;
}

/*
// This is the Chebyshev-Gauss version.
double ReducedPotentialQuadrature::Omega(size_t l, size_t s, double T) {
  double coeff = -1.0 / T / std::tgamma(s + 2);
  double esterr;
  double quadrature;
  double converged;
  omegacg_.set_param(l, s, T);
  std::tie(quadrature, esterr, converged) = omegacg_.integrate(tol_,
CG_MAXORDER); if (!converged) { std::cerr << "Line " << __LINE__ << " OmegaCG
not converged with l = " << l
              << " s = " << s << " and T = " << T << "." << std::endl;
  }
  quadrature /= 2;
  return coeff * quadrature;
}
*/

// This is the Gauss-Laguerre version.
double ReducedPotentialQuadrature::Omega(size_t l, size_t s, double T,
                                         double rtol) {
  double coeff = 1.0 / std::tgamma(s + 2);
  double esterr;
  double quadrature;
  double converged;
  omegagl_.set_param(l, s, T);
  std::tie(quadrature, esterr, converged) = omegagl_.integrate(rtol, 128);
  if (!converged) {
    std::cerr << "Line " << __LINE__ << " OmegaGL not converged with l = " << l
              << " s = " << s << " and T = " << T << "." << std::endl;
    // omegacg_.show_integrands();
  }
  // std::cerr << "l = " << l << ", s= " << s << ", T = " << T << std::endl;
  return coeff * quadrature;
}

} // namespace dlt
