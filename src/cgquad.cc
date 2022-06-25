#include <chrono>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <vector>

#include "cgquad.hh"

namespace dlt {
CGIntegratorBackend *CGIntegratorBackend::instance_ = nullptr;
CGIntegrator::CGIntegrator(size_t order) {
  auto backend = CGIntegratorBackend::instance();
  backend->angles_.push_back(0.5); // pi/2
  backend->coss_.push_back(0.0);
  backend->maxorder_ = 0;
  backend->size_ = 1;
  backend->gap_ = 1.0;
  backend->allocate(order);
}
void CGIntegratorBackend::allocate(size_t order) {
  if (order <= maxorder_) {
    return;
  }
  size_t num_ = maxidx(order);
  angles_.reserve(num_);
  coss_.reserve(num_);
  for (; maxorder_ != order; ++maxorder_) {
    gap_ /= 3.0;
    for (size_t i = 1; i != maxidx(maxorder_); ++i) {
      angles_.push_back(angles_[i] - gap_);
      angles_.push_back(angles_[i] + gap_);
      coss_.push_back(cos(M_PI * (angles_[i] - gap_)));
      coss_.push_back(cos(M_PI * (angles_[i] + gap_)));
    }
    angles_.push_back(angles_[0] - gap_);
    coss_.push_back(cos(M_PI * (angles_[0] - gap_)));
  }
  return;
}
std::tuple<double, double>
CGIntegrator::integrate(std::function<double(double)> integrand, double a,
                        double b, double tol, size_t maxorder, bool negative) {
  auto backend = CGIntegratorBackend::instance();
  // init
  double k0 = 0.5 * (b + a);
  double k1 = 0.5 * (b - a);
  double res = (negative ? 1 : 0.5) *
               integrand(backend->map_pm1(backend->coss_[0], k0, k1));
  double oldint = 0.0;
  double newint = 0.0;
  double err = 0.0;
  for (size_t order = 1; order <= maxorder; ++order) {
    // std::cerr << "Order =    " << order << std::endl;
    backend->allocate(order);
    // std::cout << coss_.size() << std::endl;
    // std::cout << maxidx(order) << std::endl;
    for (size_t i = CGIntegratorBackend::maxidx(order - 1);
         i != CGIntegratorBackend::maxidx(order); ++i) {
      // std::cout << __FILE__ << __LINE__ << (map_pm1(coss_[i], k0, k1))
      // << std::endl;
      res += integrand(CGIntegratorBackend::map_pm1(backend->coss_[i], k0, k1));
    }
    if (negative) {
      for (size_t i = CGIntegratorBackend::maxidx(order - 1);
           i != CGIntegratorBackend::maxidx(order); ++i) {
        res +=
            integrand(CGIntegratorBackend::map_pm1(-backend->coss_[i], k0, k1));
      }
    }
    integrand(CGIntegratorBackend::map_pm1(backend->coss_[0], k0, k1));
    oldint = newint;
    newint = (negative ? 1 : 2) * res * M_PI / pow(3, order) * k1;
    err = fabs(oldint - newint);
    // std::cerr << "ERR for order" << order << " : " << err << "   ";
    err /= fabs(newint); // relative error...
    std::cerr << err << std::endl;
    // std::cout << newint << std::endl;
    if (err < tol) {
      std::cerr << "Meet the errtol requirement   " << order << std::endl;
      break;
    }
  }
  return std::make_tuple(newint, err);
}
double CGIntegratorBackend::map_pm1(double x, double k0, double k1) {
  return k0 + k1 * x;
}
size_t CGIntegratorBackend::maxidx(size_t order) {
  return (pow(3, order) + 1) / 2;
}
} // namespace dlt
