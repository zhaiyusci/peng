#include <chrono>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <vector>

namespace dlt {

class CGIntegrator {
public:
  // Status

public:
  CGIntegrator() : CGIntegrator(1) {}
  CGIntegrator(size_t order);
  std::tuple<double, double> integrate(std::function<double(double)> integrand,
                                       double a, double b, double tol,
                                       size_t maxorder, bool negative = true);
};

class CGIntegratorBackend {

protected:
  std::vector<double> angles_;
  std::vector<double> coss_;
  // Storage
  size_t maxorder_;
  double gap_;
  size_t size_ = 1;
  static CGIntegratorBackend *instance_;
  CGIntegratorBackend(){};
  friend CGIntegrator;

public:
  void allocate(size_t order);
  static CGIntegratorBackend *instance() {
    if (instance_ == nullptr) {
      instance_ = new CGIntegratorBackend;
    }
    return instance_;
  }
  static double map_pm1(double x, double k0, double k1);
  static size_t maxidx(size_t order);
};
} // namespace dlt
