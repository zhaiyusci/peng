#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

class CCIntegrator {
public:
  std::vector<double> coss_;
  std::vector<double> sins_;
  // Storage
  size_t maxorder_;
  size_t maxnum_;
  // Status
  size_t order_;
  size_t num_;

public:
  CCIntegrator() : CCIntegrator(5) {}
  CCIntegrator(size_t order) {
    sins_.push_back(0.0);
    sins_.push_back(1.0);
    coss_.push_back(1.0);
    coss_.push_back(0.0);
    maxorder_ = 0;
    allocate(order);
    order_ = 0;
  }
  void set_order(size_t order) {
    if (order > maxorder_) {
      allocate(order);
    }
    order_ = order;
    num_ = 1 << order >> 1;
    return;
  }
  size_t idx(size_t r) {
    if (r == 0)
      return 0;
    size_t b = (num_ >> 1);
    while (!(r & 1)) {
      r >>= 1;
      b >>= 1;
    }
    r >>= 1;
    r += b + 1;
    return r;
  }
  void allocate(size_t order) {
    std::cerr << __LINE__ << "    " << order << std::endl;
    if (order <= maxorder_) {
      return;
    }
    size_t tmpord = order_;
    set_order(order);
    sins_.reserve(maxorder_);
    coss_.reserve(maxorder_);
    for (; maxorder_ != order; ++maxorder_) {
      set_order(maxorder_);
      maxnum_ = num_;
      double lc = coss_[0];
      double ls = sins_[0];
      for (size_t i = 0; i != maxnum_; ++i) {
        size_t r = idx(i + 1);
        std::cerr << r << std::endl;
        double rc = coss_[r];
        double rs = sins_[r];
        double c = lc * rc - ls * rs;
        double s = sqrt((1 - c) / 2);
        c = sqrt((1 + c) / 2);
        sins_.push_back(s);
        coss_.push_back(c);
        lc = rc;
        ls = rs;
      }
    }
    set_order(tmpord);
    return;
  }
};

int main() {
  std::cout << std::setprecision(8);
  // Only consider [0,pi/2]
  auto t0 = std::chrono::high_resolution_clock::now();
  auto dt = t0 - t0;
  size_t maxorder = 3;
  CCIntegrator cc;
  cc.set_order(maxorder);
  // init
  auto t1 = std::chrono::high_resolution_clock::now();
  dt += t1 - t0;
  std::cout << "Algo4: " << dt.count() << std::endl;

  // std::vector<int> ord{0, 8, 4, 2, 6, 1, 3, 5, 7}; // For test 0-8
  std::vector<int> ord{0, 4, 2, 1, 3}; // For test 0-8
  for (int i = 0; i <= pow(2, maxorder - 1); ++i) {
    std::cout << cc.coss_[i] << "  " << cos(M_PI / pow(2, maxorder) * ord[i])
              << "  " << cc.sins_[i] << "  "
              << sin(M_PI / pow(2, maxorder) * ord[i]) << std::endl;
  }

  return 0;
}