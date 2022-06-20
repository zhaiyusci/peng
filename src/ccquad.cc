#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

class CCIntegrator {
private:
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
    order = 0;
  }
  void set_order(size_t order) {
    order_ = order;
    num_ = 1 << order >> 1;
    return;
  }
  void allocate(size_t order) {
    if (order <= maxorder_) {
      return;
    }
    size_t n = pow(2, maxorder - 1) + 1;
    sins_.reserve(n);
    coss_.reserve(n);
    for (size_t order = 1; order != maxorder; ++order) {
      size_t maxnum = pow(2, order - 1);
      double lc = coss_[0];
      double ls = sins_[0];
      for (size_t i = 0; i != maxnum; ++i) {
        size_t r = i + 1;
        size_t b = (maxnum >> 1);
        while (!(r & 1)) {
          r >>= 1;
          b >>= 1;
        }
        r >>= 1;
        r += b + 1;
        // std::cerr << r << std::endl;
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
  }
};

int main() {
  std::cout << std::setprecision(8);
  // 0 - pi
  {
    // Only consider [0,pi/2]
    auto t0 = std::chrono::high_resolution_clock::now();
    auto dt = t0 - t0;
    // init
    auto t1 = std::chrono::high_resolution_clock::now();
    dt += t1 - t0;
    std::cout << "Algo4: " << dt.count() << std::endl;
  }

  std::cout << coss_.size() << std::endl;
  // std::vector<int> ord{0, 8, 4, 2, 6, 1, 3, 5, 7}; // For test 0-8
  std::vector<int> ord{0, 4, 2, 1, 3}; // For test 0-8
  for (int i = 0; i <= pow(2, maxorder - 1); ++i) {
    std::cout << coss_[i] << "  " << cos(M_PI / pow(2, maxorder) * ord[i])
              << "  " << sins_[i] << "  "
              << sin(M_PI / pow(2, maxorder) * ord[i]) << std::endl;
  }

  return 0;
}
