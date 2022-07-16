#include "omega.hh"
#include "pot1dquad.hh"

namespace dlt {

// Omega stuff

///
/// This class inheriated from the Chebyshev-Gauss Quarduture class, compute the
/// integration required by the computation of Omega.
///
OmegaCG::OmegaCGInt::OmegaCGInt(ReducedPotentialQuadrature &rpq)
    : CGIntegrator(true, -0.5, 1.0), rpq_(&rpq) {
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
void OmegaCG::OmegaCGInt::set_param(size_t l, size_t s, double T) {
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
void OmegaCG::OmegaCGInt::calculate_integrands(size_t ordersize) {
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
      Qs_.push_back(rpq_->Q(l_, r_E, -1.0, q_rtol));
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

      vs_.push_back(rpq_->potential_value(r_E));
      // weight is included in dvs
      dvs_.push_back(rpq_->potential_derivative(r_E) * sqrt(1.0 - y * y));
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

OmegaGL::OmegaGLInt::OmegaGLInt(ReducedPotentialQuadrature &rpq)
    : GLIntegrator(), rpq_(&rpq) {
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
void OmegaGL::OmegaGLInt::set_param(size_t l, size_t s,
                                                    double T) {
  s_ = s;
  set_alpha(s + 1.0);
  T_ = T;
  l_ = l;
  return;
}

/** Compute the integrands with computed values cached.
 */
void OmegaGL::OmegaGLInt::calculate_integrands(size_t ngridsize) {
  if (ngridsize !=xs_.size()){
    throw std::runtime_error("Numbers of grids in GL quad mismatch.");
  }
  for (size_t i = 0; i != xs_.size(); ++i) {
    double grid_q_rtol = ws_[i] > 2 ? q_rtol / 2 : q_rtol / ws_[i];
    double value = rpq_->Q(l_, -1.0, xs_[i] * T_, grid_q_rtol);
    // std::cerr << "Q(" << l_ << ", " << x * T_ << ") = " << value <<
    // std::endl;
    integrands_.push_back(value);
  }
  return;
}

// This is the Chebyshev-Gauss version.
double OmegaCG::Omega(size_t l, size_t s, double T, double rtol) {
  double coeff = -1.0 / T / std::tgamma(s + 2);
  double esterr;
  double quadrature;
  double converged;
  omegacgint_.set_param(l, s, T);
  std::tie(quadrature, esterr, converged) =
      omegacgint_.integrate(rtol,CG_INT_ORDER_MAX);
  if (!converged) {
    std::cerr << "Line " << __LINE__ << " OmegaCG not converged with l = " << l
              << " s = " << s << " and T = " << T << "." << std::endl;
  }
  quadrature /= 2;
  return coeff * quadrature;
}

// This is the Gauss-Laguerre version.
double OmegaGL::Omega(size_t l, size_t s, double T, double rtol) {
  double coeff = 1.0 / std::tgamma(s + 2);
  double esterr;
  double quadrature;
  double converged;
  omegaglint_.set_param(l, s, T);
  omegaglint_.q_rtol = rtol;
  std::tie(quadrature, esterr, converged) = omegaglint_.integrate(rtol, GL_INT_ORDER_MAX);
  if (!converged) {
    std::cerr << "Line " << __LINE__ << " OmegaGL not converged with l = " << l
              << " s = " << s << " and T = " << T << "." << std::endl;
    // omegacg_.show_integrands();
  }
  // std::cerr << "l = " << l << ", s= " << s << ", T = " << T << std::endl;
  return coeff * quadrature;
}

} // namespace dlt
