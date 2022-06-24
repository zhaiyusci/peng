#include <chrono>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <vector>

#include "cgquad.hh"

namespace dlt {
CGIntegrator::CGIntegrator(size_t order) {
  angles_.push_back(0.5); // pi/2
  coss_.push_back(0.0);
  maxorder_ = 0;
  order_ = 0;
  size_ = 1;
  gap_ = 1.0;
  allocate(order);
}
void CGIntegrator::allocate(size_t order) {
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
  // init
  double k0 = 0.5 * (b + a);
  double k1 = 0.5 * (b - a);
  double res = (negative ? 1 : 0.5) * integrand(map_pm1(coss_[0], k0, k1));
  double oldint = 0.0;
  double newint = 0.0;
  double err = 0.0;
  for (size_t order = 1; order <= maxorder; ++order) {
    // std::cerr << "Order =    " << order << std::endl;
    allocate(order);
    // std::cout << coss_.size() << std::endl;
    // std::cout << maxidx(order) << std::endl;
    for (size_t i = maxidx(order - 1); i != maxidx(order); ++i) {
      // std::cout << __FILE__ << __LINE__ << (map_pm1(coss_[i], k0, k1))
      // << std::endl;
      res += integrand(map_pm1(coss_[i], k0, k1));
    }
    if (negative) {
      for (size_t i = maxidx(order - 1); i != maxidx(order); ++i) {
        res += integrand(map_pm1(-coss_[i], k0, k1));
      }
    }
    integrand(map_pm1(coss_[0], k0, k1));
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
double CGIntegrator::map_pm1(double x, double k0, double k1) {
  return k0 + k1 * x;
}
size_t CGIntegrator::maxidx(size_t order) { return (pow(3, order) + 1) / 2; }
} // namespace dlt
