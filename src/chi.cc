#include "chi.hh"
#include "pot1dquad.hh"

namespace peng {
ChiCG::ChiCGInt::ChiCGInt(ReducedPotentialQuadrature &rpq)
    : CGIntegrator(true), rpq_(&rpq) {
  clean_cache();
  E_ = -2.0;
  r_m_ = -1.0;
}

void ChiCG::ChiCGInt::set_param(double r_m, double E, double b) {
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
void ChiCG::ChiCGInt::calculate_integrands(size_t ordersize) {
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
      double v = rpq_->potential_value(r);
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

double ChiCG::chi(double E, double r_m, double rtol) {
  // double E = ppot_->value(r_E);
  double b = rpq_->r2b(r_m, E);

  chicgint.set_param(r_m, E, b);
  double quadrature, err;
  bool converged;
  std::tie(quadrature, err, converged) = chicgint.integrate(1.1, 5);
  double res = M_PI - b / r_m * quadrature;

  if (std::fabs(res) <= 20 * M_PI) {
    // Do this only when chi is not so large,
    // i.e., we do not want waste or computational resources on
    // orbiting ...
    // Following Barker et al., Phys Fluids, 7, 897 (1964)

    std::tie(quadrature, err, converged) =
        chicgint.integrate(rtol, CG_INT_ORDER_MAX);
    if (!converged) {
      std::cerr << "Line " << __LINE__ << " ChiCG not converged with E = " << E
                << ", r_m = " << r_m << " and b = " << b << "." << '\n';
      std::cerr << "quadrature   " << quadrature << ' ' << M_PI - b * quadrature
                << '\n';
      // chicg.show_integrands();
    }

    res = M_PI - b / r_m * quadrature;
  }
  return res;
}

} // namespace peng
