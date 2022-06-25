#include <chrono>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <vector>

namespace dlt {

class CGIntegrator {
protected:
  double a_;
  double b_;
  double k0_;
  double k1_;
  bool is_symm_;

public:
  CGIntegrator(double a, double b);
  std::tuple<double, double> integrate(double tol, size_t maxorder,
                                       bool negative = true);
  double map_pm1(double x);
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
  CGIntegratorBackend();
  friend CGIntegrator;

public:
  void allocate(size_t order);
  static CGIntegratorBackend *instance();
};

class N3Allocator {
protected:
  std::vector<double> vals_;
  size_t maxorder_;

public:
  void allocate(size_t order);
};
} // namespace dlt
