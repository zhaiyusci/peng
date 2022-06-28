#include <chrono>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <vector>

#include "cgquad.hh"

namespace dlt {

CGIntegratorBackend *CGIntegratorBackend::instance_ = nullptr;
CGIntegrator::CGIntegrator(double a, double b, bool symm)
    : a_(a), b_(b), symm_(symm) {
  k0_ = 0.5 * (b + a);
  k1_ = 0.5 * (b - a);
}

void CGIntegratorBackend::allocate(size_t ordersize) {
  if (ordersize <= ordersize_) {
    return;
  }
  size_t num = (pow(3, ordersize - 1) + 1) / 2;
  angles_.reserve(num);
  coss_.reserve(num);
  for (; ordersize_ != ordersize; ++ordersize_) {
    CubicIter ci(1, ordersize_, true);
    gap_ /= 3.0;
    for (auto &&i : ci) {
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

std::tuple<double, double> CGIntegrator::integrate(double tol,
                                                   size_t maxordersize) {
  auto backend = CGIntegratorBackend::instance();
  backend->allocate(1);
  calculate_integrands(1);
  // init
  double res = (symm_ ? 0.5 : 1) * integrands_[0];
  double oldint = 0.0;
  double newint = 0.0;
  double err = 0.0;
  for (size_t ordersize = 2; ordersize != maxordersize; ++ordersize) {
    backend->allocate(ordersize);
    calculate_integrands(ordersize);
    // std::cout << integrands_.size() << std::endl;
    CubicIter ci(ordersize - 1, ordersize, symm_, symm_);
    for (auto &&i : ci) {
      // std::cout << __LINE__ << ' ' << i << ' ' << integrands_[i] <<
      // std::endl;
      res += integrands_[i];
    }
    oldint = newint;
    newint = (symm_ ? 2 : 1) * res * M_PI / pow(3, ordersize - 1) * k1_;
    err = fabs(oldint - newint);
    // std::cerr << "ERR for order" << order << " : " << err << "   ";
    err /= fabs(newint); // relative error...
    // std::cerr << err << std::endl;
    // std::cout << newint << std::endl;
    if (err < tol) {
      return std::make_tuple(newint, err);
    }
  }
  std::cerr << "WARNING\n"
            << __FILE__ << ' ' << __LINE__ << ": "
            << "The Chebyshev-Gauss Quadrature did NOT meet the errtol "
               "requirement with order = "
            << maxordersize << "." << std::endl;
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
  ordersize_ = 1;
  gap_ = 1.0;
}

double CGIntegratorBackend::coss(int i) {
  return (i < 0 ? -1 : 1) * coss_[abs(i)];
}

CubicIter::iterator::iterator(CubicIter *ci, size_t order, size_t negative,
                              size_t idx)
    : ci_(ci), order_(order), negative_(negative), idx_(idx) {
  if (order == 0) {
    totidx_ = 0;
    idxsize_ = 1;
  } else {
    idxsize_ = pow(3, order - 1);
    if (ci_->half_storage_) {
      totidx_ = (pow(3, order - 1) + 1) / 2 + idx;
    } else {
      totidx_ = pow(3, order - 1) + negative * idxsize_ + idx;
    }
  }
}
CubicIter::iterator &CubicIter::iterator::operator++() {
  ++totidx_;
  ++idx_;
  if (order_ == 0) {
    idx_ = 0;
    negative_ = 0;
    order_ = 1;
    idxsize_ = 1;
    return *this;
  }
  if (idx_ != idxsize_) {
    return *this;
  }
  idx_ = 0;
  ++negative_;
  if ((!ci_->half_visit_) && negative_ == 1 && ci_->half_storage_)
    totidx_ -= idxsize_;
  if (negative_ != (ci_->half_visit_ ? 0 : 1) + 1) {
    return *this;
  }
  idx_ = 0;
  negative_ = 0;
  ++order_;
  idxsize_ = pow(3, order_ - 1);
  return *this;
}

} // namespace dlt

/*
int main() {
  {
    dlt::CubicIter ic(2, 3);
    for (auto &&i : ic) {
      std::cout << i << '\n';
    }
    std::cout << "=======" << std::endl;
  }
  {
    dlt::CubicIter ic(2, 3, true);
    for (auto &&i : ic) {
      std::cout << i << '\n';
    }
    std::cout << "=======" << std::endl;
  }
  if (false) {
    dlt::CubicIter ic(2, 3, true, false);
    for (auto &&i : ic) {
      std::cout << i << '\n';
    }
    std::cout << "=======" << std::endl;
  }
  {
    dlt::CubicIter ic(2, 3, false, false);
    for (auto &&i : ic) {
      std::cout << i << '\n';
    }
    std::cout << "=======" << std::endl;
  }
  {
    dlt::CubicIter ic(2, 3, false, false);
    std::cout << *dlt::CubicIter::iterator(&ic, 2, 1, 1) << std::endl;
  }
  return 0;
}
*/
