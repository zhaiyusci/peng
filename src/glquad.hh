#ifndef _GLQUAD_HH_
#define _GLQUAD_HH_
#include <chrono>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <vector>

// Credit: The working part is taken from
// https://people.sc.fsu.edu/~jburkardt/cpp_src/gen_laguerre_rule/gen_laguerre_rule.html
// and Yu Zhai add a decent C++ interface

namespace dlt {

class GLIntegrator {
  // public:
protected:
  double alpha_;
  // Working space
  size_t ngrids_;
  std::vector<double> xs_;
  std::vector<double> ws_;
  std::vector<double> integrands_;

public:
  GLIntegrator(double alpha = 0.0);
  void set_alpha(double alpha) {
    alpha_ = alpha;
    return;
  }
  std::tuple<double, double, bool> integrate(double tol, size_t ngridsize);
  void clean_workspace();
  void reset_workspace(size_t ngrids);

  ///
  /// User implemented method. Basically, user should
  ///
  /// 1. Fill in integrands_ using xs_.
  ///
  virtual void calculate_integrands() = 0;
};

} // namespace dlt

#endif
