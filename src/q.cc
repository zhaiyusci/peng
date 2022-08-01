#include "q.hh"
#include "pot1dquad.hh"
#include <cmath>

namespace peng {
QCG::QCGInt1::QCGInt1(ReducedPotentialQuadrature &rpq)
    : CGIntegrator(false), rpq_(&rpq) {
  clean_cache();
  l_ = 0;
  r_E_ = -1.0;
  r_Op_ = -1.0;
}

void QCG::QCGInt1::set_param(size_t l, double r_E, double r_Op, double E) {
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
      E_ = rpq_->potential_value(r_E);
    } else {
      E_ = E;
    }
  }
  return;
}

/** Compute the integrands with computed values cached.
 */
void QCG::QCGInt1::calculate_integrands(size_t ordersize) {
  // See if we need an update based on the "flag".
  if (cache_ordersize_ < ordersize) {
    CubicIter ci(cache_ordersize_, ordersize, false, true);
    size_t num = ci.size_from_0();
    integrands_.reserve(num);
    chis_.reserve(num);
    fct2_.reserve(num);
    for (auto &&i : ci) {
      double y = CGIntegratorBackend::instance()->coss(i);
      double r_m = map_pm1(y);
      double v, dv;
      v = rpq_->potential_value(r_m);
      dv = rpq_->potential_derivative(r_m);
      double chi = rpq_->chi(E_, r_m, chi_rtol);
      // Here I put everything else in fct2, include the weight
      double fct2 = (2.0 * (E_ - v) - r_m * dv) * r_m * sqrt(1.0 - y * y);
      chis_.push_back(chi);
      fct2_.push_back(fct2);
    }
    cache_ordersize_ = ordersize;
  }
  if (ordersize_ < ordersize) {
    CubicIter ci(ordersize_, ordersize, false, false);
    size_t num = ci.size_from_0();
    integrands_.reserve(num);
    for (auto &&i : ci) {
      double res;
      if (std::fabs(chis_[i]) <= 20 * M_PI) {
        res = (1.0 - pow(cos(chis_[i]), l_)) * fct2_[i];
      } else {
        if (l_ % 2 == 1) {
          res = fct2_[i];
        } else {
          res = (1.0 - std::tgamma(l_ / 2.0 + 0.5) /
                           std::tgamma(l_ / 2.0 + 1.0) / std::sqrt(M_PI)) *
                fct2_[i];
        }
      }
      integrands_.push_back(res);
    }
    ordersize_ = ordersize;
  }
  return;
}

QCG::QCGInt2::QCGInt2(ReducedPotentialQuadrature &rpq)
    : CGIntegrator(true), rpq_(&rpq) {
  clean_cache();
  l_ = 0;
  r_E_ = -1.0;
  r_O_ = -1.0;
}

/** Set the parameters, clear the inner storage if it is needed.
 */
void QCG::QCGInt2::set_param(size_t l, double r_E, double r_O, double E) {
  if (l_ != l) {
    clean_workspace();
    l_ = l;
  }
  if (r_E_ != r_E || r_O_ != r_O) {
    clean_cache();
    r_E_ = r_E;
    r_O_ = r_O;
    if (E < 0.0) {
      E_ = rpq_->potential_value(r_E);
    } else {
      E_ = E;
    }
  }
  return;
}

/** Compute the integrands with computed values cached.
 */
void QCG::QCGInt2::calculate_integrands(size_t ordersize) {
  // See if we need an update based on the "flag".
  if (cache_ordersize_ < ordersize) {
    CubicIter ci(cache_ordersize_, ordersize, true, true);
    size_t num = ci.size_from_0();
    integrands_.reserve(num);
    chis_.reserve(num);
    fct2_.reserve(num);
    for (auto &&i : ci) {
      double y = CGIntegratorBackend::instance()->coss(i);
      if (y < 1.0e-8) {
        y = 1.0e-8; // prevent div by 0
      }
      double r_m = r_O_ / y;
      double v, dv;
      v = rpq_->potential_value(r_m);
      dv = rpq_->potential_derivative(r_m);
      double chi = rpq_->chi(E_, r_m, chi_rtol);
      // Here I put everything else in fct2, include the weight
      double fct2 =
          (2.0 * (E_ - v) - r_m * dv) * r_m * r_m / y * sqrt(1.0 - y * y);
      chis_.push_back(chi);
      fct2_.push_back(fct2);
    }
    cache_ordersize_ = ordersize;
  }
  if (ordersize_ < ordersize) {
    CubicIter ci(ordersize_, ordersize, true, true);
    size_t num = ci.size_from_0();
    integrands_.reserve(num);
    for (auto &&i : ci) {
      double res;
      if (std::fabs(chis_[i]) <= 20 * M_PI) {
        res = (1.0 - pow(cos(chis_[i]), l_)) * fct2_[i];
      } else {
        if (l_ % 2 == 1) {
          res = fct2_[i];
        } else {
          res = (1.0 - std::tgamma(l_ / 2.0 + 0.5) /
                           std::tgamma(l_ / 2.0 + 1.0) / std::sqrt(M_PI)) *
                fct2_[i];
        }
      }
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
double QCG::Q(size_t l, double r_E, double E,
                                     double rtol) {
  // double old_r_E = 0.0;
  double r_O, r_Op;

  if (E < 0.0) {
    E = rpq_->potential_value(r_E); // This is common
  } else if (r_E < 0.0) {
    r_E = rpq_->v_root(E);
  }
  std::tie(r_O, r_Op) = rpq_->r_range(E);


  double coeff = 1.0 / (1.0 - (1.0 + pow(-1, l)) / 2.0 / (1.0 + l)) / E;
  double esterr;
  bool converged;
  if (E <= 2 * rpq_->E_C()) {
    double quadrature1;
    qcgint1_.set_param(l, r_E, r_Op, E);
    qcgint1_.chi_rtol=rtol;

    std::tie(quadrature1, esterr, converged) =
        qcgint1_.integrate(rtol,CG_INT_ORDER_MAX);
    if (!converged) {
      std::cerr << "Line " << __LINE__ << " QCG1 not converged with E = " << E
                << "." << std::endl;
    }
    double quadrature2;
    qcgint2_.set_param(l, r_E, r_O, E);
    qcgint2_.chi_rtol=rtol;
    std::tie(quadrature2, esterr, converged) =
        qcgint2_.integrate(rtol,CG_INT_ORDER_MAX);
    if (!converged) {
      std::cerr << "Line " << __LINE__ << " QCG2 not converged with E = " << E
                << "." << std::endl;
    }
    quadrature2 /= 2;
    return coeff * (quadrature1 + quadrature2);
  } else {
    double quadrature2;
    qcgint2_.set_param(l, r_E, r_E, E);
    qcgint2_.chi_rtol=rtol;
    std::tie(quadrature2, esterr, converged) =
        qcgint2_.integrate(rtol,CG_INT_ORDER_MAX);
    if (!converged) {
      std::cerr << "Line " << __LINE__ << " QCG2 not converged with E = " << E
                << "." << std::endl;
    }
    quadrature2 /= 2;
    return coeff * quadrature2;
  }
}


} // namespace peng
