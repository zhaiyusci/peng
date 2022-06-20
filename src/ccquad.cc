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
    order_ = 0;
    allocate(order);
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
  size_t idxrev(size_t r) {
    if (r == 0)
      return 0;
    size_t b = 0;
    size_t rt = r - 1;
    while (rt) {
      rt >>= 1;
      ++b;
    }
    // std::cout << "b " << b << std::endl;
    size_t r0 = (1 << b >> 1) + 1;
    rt = r - r0;
    return ((rt + 1) << (order_ - b)) - (1 << (order_ - b) >> 1);
  }
  double sins(size_t r) {
    return r > num_ ? sins_[idx(2 * num_ - r)] : sins_[idx(r)];
  }
  double coss(size_t r) {
    return r > num_ ? -coss_[idx(2 * num_ - r)] : coss_[idx(r)];
  }
  void allocate(size_t order) {
    if (order <= maxorder_) {
      return;
    }
    size_t tmpord = order_;
    maxnum_ = 1 << order >> 1;
    sins_.reserve(2 * maxnum_ + 1);
    coss_.reserve(2 * maxnum_ + 1);
    for (; maxorder_ != order; ++maxorder_) {
      set_order(maxorder_);
      maxnum_ = num_;
      double lc = coss_[0];
      double ls = sins_[0];
      for (size_t i = 0; i != maxnum_; ++i) {
        size_t r = idx(i + 1);
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
    maxnum_ *= 2;
    set_order(tmpord);
    return;
  }
};

int main() {
  std::cout << std::setprecision(8);
  // Only consider [0,pi/2]
  auto t0 = std::chrono::high_resolution_clock::now();
  auto dt = t0 - t0;
  CCIntegrator cc(12);
  // cc.set_order(maxorder);
  // init
  auto t1 = std::chrono::high_resolution_clock::now();
  dt += t1 - t0;
  std::cout << "Algo4: " << dt.count() << std::endl;

  for (size_t maxorder = 3; maxorder != 5; ++maxorder) {
    std::cout << "= = = = = = = = = = = = = = = = = = = = " << std::endl;
    cc.set_order(maxorder);
    // std::vector<int> ord{0, 8, 4, 2, 6, 1, 3, 5, 7}; // For test 0-8
    std::vector<int> ord{0, 4, 2, 1, 3}; // For test 0-8
    for (int i = 0; i != pow(2, maxorder) + 1; ++i) {
      std::cout << cc.coss(i) << "  " << cos(M_PI / pow(2, maxorder) * i)
                << "  " << cc.sins(i) << "  "
                << sin(M_PI / pow(2, maxorder) * i) << std::endl;
    }
  }

  cc.set_order(4);
  // cc.idxrev(0);
  std::cout << cc.idxrev(0);
  std::cout << cc.idxrev(1);
  std::cout << cc.idxrev(2);
  std::cout << cc.idxrev(3);
  std::cout << cc.idxrev(4);
  std::cout << cc.idxrev(5);
  std::cout << cc.idxrev(6);
  std::cout << cc.idxrev(7);
  std::cout << cc.idxrev(8);
  std::cout << std::endl;

  return 0;
}
