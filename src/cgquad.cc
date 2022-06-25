#include <chrono>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <vector>

#include "cgquad.hh"

namespace dlt {

size_t maxidx_symm(size_t order) { return (pow(3, order) + 1) / 2; }
// This class used to iterate over the n^3 data...
class CubicIter {
private:
  size_t border_;
  size_t eorder_;
  size_t symm_; // 0: only positive half
                // 1: also negative half

public:
  CubicIter() = delete;
  CubicIter(size_t eorder) : border_(0), eorder_(eorder), symm_(1) {}
  CubicIter(size_t border, size_t eorder)
      : border_(border), eorder_(eorder), symm_(1) {}
  CubicIter(size_t border, size_t eorder, size_t symm)
      : border_(border), eorder_(eorder), symm_(symm) {}
  class iterator {
  public:
    CubicIter *ci_;
    size_t order_;
    size_t negative_;
    size_t idx_;
    size_t maxidx_;
    size_t totidx_;

  public:
    iterator &operator++() {
      ++totidx_;
      ++idx_;
      if (order_ == 0) {
        idx_ = 0;
        negative_ = 0;
        order_ = 1;
        maxidx_ = maxidx_symm(order_ - 1);
        return *this;
      }
      if (idx_ != maxidx_) {
        return *this;
      }
      idx_ = 0;
      ++negative_;
      if (ci_->symm_ == 1 && negative_ == 1)
        totidx_ -= maxidx_;
      if (negative_ != ci_->symm_ + 1) {
        return *this;
      }
      idx_ = 0;
      negative_ = 0;
      ++order_;
      maxidx_ = maxidx_symm(order_ - 1);
      return *this;
    }
    int operator*() const { return totidx_ * (negative_ ? -1 : 1); }
    bool operator==(iterator &rhs) const {
      return order_ == rhs.order_ && negative_ == rhs.negative_ &&
             idx_ == rhs.idx_;
    }
    bool operator!=(iterator &rhs) const { return !(operator==(rhs)); }
    iterator(CubicIter *ci, size_t order, size_t negative, size_t idx)
        : ci_(ci), order_(order), negative_(negative), idx_(idx) {
      if (order == 0) {
        totidx_ = 0;
        maxidx_ = 1;
      } else {
        totidx_ = maxidx_symm(order - 1) + idx;
        maxidx_ = maxidx_symm(order - 1);
      }
    }
  };
  iterator begin() { return iterator(this, border_, 0, 0); }
  iterator end() { return iterator(this, eorder_ + 1, 0, 0); }
};

CGIntegratorBackend *CGIntegratorBackend::instance_ = nullptr;
CGIntegrator::CGIntegrator(double a, double b) : a_(a), b_(b) {
  k0_ = 0.5 * (b + a);
  k1_ = 0.5 * (b - a);
}

void CGIntegratorBackend::allocate(size_t order) {
  if (order <= maxorder_) {
    return;
  }
  size_t num_ = maxidx_symm(order);
  angles_.reserve(num_);
  coss_.reserve(num_);
  for (; maxorder_ != order; ++maxorder_) {
    gap_ /= 3.0;
    for (size_t i = 1; i != maxidx_symm(maxorder_); ++i) {
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

std::tuple<double, double> CGIntegrator::integrate(double tol, size_t maxorder,
                                                   bool negative) {
  auto backend = CGIntegratorBackend::instance();
  // init
  double res = (negative ? 1 : 0.5) * integrand(map_pm1(backend->coss_[0]));
  double oldint = 0.0;
  double newint = 0.0;
  double err = 0.0;
  for (size_t order = 1; order <= maxorder; ++order) {
    // std::cerr << "Order =    " << order << std::endl;
    backend->allocate(order);
    // std::cout << coss_.size() << std::endl;
    // std::cout << maxidx(order) << std::endl;
    for (size_t i = CGIntegratorBackend::maxidx_symm(order - 1);
         i != CGIntegratorBackend::maxidx_symm(order); ++i) {
      // std::cout << __FILE__ << __LINE__ << (map_pm1(coss_[i], k0, k1))
      // << std::endl;
      res += integrand(map_pm1(backend->coss_[i]));
    }
    if (negative) {
      for (size_t i = CGIntegratorBackend::maxidx_symm(order - 1);
           i != CGIntegratorBackend::maxidx_symm(order); ++i) {
        res += integrand(map_pm1(-backend->coss_[i]));
      }
    }
    integrand(map_pm1(backend->coss_[0]));
    oldint = newint;
    newint = (negative ? 1 : 2) * res * M_PI / pow(3, order) * k1_;
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
double CGIntegrator::map_pm1(double x) { return k0_ + k1_ * x; }
CGIntegratorBackend *CGIntegratorBackend::instance() {
  if (instance_ == nullptr) {
    instance_ = new CGIntegratorBackend;
  }
  return instance_;
}
CGIntegratorBackend::CGIntegratorBackend() {
  angles_.push_back(0.5); // pi/2
  coss_.push_back(0.0);
  maxorder_ = 0;
  size_ = 1;
  gap_ = 1.0;
};
} // namespace dlt
