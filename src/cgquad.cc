#include <chrono>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <vector>

class CGIntegrator {
public:
  std::vector<double> angles_;
  std::vector<double> coss_;
  // Storage
  size_t maxorder_;
  double gap_;
  size_t size_ = 1;
  // Status
  size_t order_;

public:
  CGIntegrator() : CGIntegrator(1) {}
  CGIntegrator(size_t order) {
    angles_.push_back(0.5); // pi/2
    coss_.push_back(0.0);
    maxorder_ = 0;
    order_ = 0;
    size_ = 1;
    gap_ = 1.0;
    allocate(order);
  }
  void allocate(size_t order) {
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
  std::tuple<double, double> integrate(std::function<double(double)> integrand,
                                       double a, double b, double tol,
                                       size_t maxorder, bool negative = true) {
    // init
    double k0 = 0.5 * (b + a);
    double k1 = 0.5 * (b - a);
    double res = (negative ? 1 : 0.5) * integrand(map_pm1(coss_[0], k0, k1));
    double oldint = 0.0;
    double newint = 0.0;
    double err;
    for (size_t order = 1; order <= maxorder; ++order) {
      allocate(order);
      std::cout << coss_.size() << std::endl;
      std::cout << maxidx(order) << std::endl;
      for (size_t i = maxidx(order - 1); i != maxidx(order); ++i) {
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
      if (err < tol) {
        std::cout << "Meet the errtol requirement   " << order << std::endl;
        break;
      }
    }
    return std::make_tuple(newint, err);
  }
  static double map_pm1(double x, double k0, double k1) { return k0 + k1 * x; }
  static size_t maxidx(size_t order) { return (pow(3, order) + 1) / 2; }
};

int main() {
  std::cout << std::setprecision(16);
  // Only consider [0,pi/2]
  auto t0 = std::chrono::high_resolution_clock::now();
  auto dt = t0 - t0;
  CGIntegrator cc(12);
  // cc.set_order(maxorder);
  // init
  auto t1 = std::chrono::high_resolution_clock::now();
  dt += t1 - t0;
  std::cout << "Algo4: " << dt.count() << std::endl;

  size_t maxorder = 3;
  // std::vector<int> ord{0, 8, 4, 2, 6, 1, 3, 5, 7}; // For test 0-8
  for (int i = 0; i != (pow(3, maxorder) + 1) / 2; ++i) {
    std::cout << cc.coss_[i] << "  " << cc.angles_[i] << std::endl;
  }

  double res, err;
  std::tie(res, err) =
      cc.integrate([](double x) { return (1.0 - x * x); }, -1, 1, 1e-8, 12);
  std::cout << res << "      " << err << std::endl;

  return 0;
}
