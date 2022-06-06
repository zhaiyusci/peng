#ifndef __DILUTE_MATHTOOLS_HH__
#define __DILUTE_MATHTOOLS_HH__
#include <iostream>
#include <nlopt.hpp>
#include <valarray>
namespace dlt {

class NegFuncDeriv1D;

class FuncDeriv1D {
private:
  NegFuncDeriv1D *pneg_ = nullptr;

public:
  double operator()(double x) const { return value(x); }
  virtual double value(double x) const = 0;
  virtual double derivative(double x) const;
  virtual bool provide_derivative() const { return false; };
  virtual FuncDeriv1D &operator-();
  virtual ~FuncDeriv1D();
};

class NegFuncDeriv1D : public FuncDeriv1D {
private:
  FuncDeriv1D *const pneg_;

public:
  NegFuncDeriv1D(FuncDeriv1D &func) : pneg_(&func) {}
  NegFuncDeriv1D(FuncDeriv1D *pfunc) : pneg_(pfunc) {}

  inline double value(double x) const { return -pneg_->value(x); }
  inline double derivative(double x) const { return -pneg_->derivative(x); }
  inline bool provide_derivative() const { return pneg_->provide_derivative(); }
  FuncDeriv1D &operator-() const { return *pneg_; };
  ~NegFuncDeriv1D() = default;
};

std::tuple<double, double> find_local_minimum(FuncDeriv1D &func, double lower,
                                              double upper, size_t nsamp = 11,
                                              double xtol_err = 1.0e-8);
std::tuple<double, double> find_local_maximum(FuncDeriv1D &func, double lower,
                                              double upper, size_t nsamp = 11,
                                              double xtol_err = 1.0e-8);
class LocalRoot {
private:
  FuncDeriv1D *const pfunc_;
  const double xtol_err_;
  std::valarray<double> xs_;
  std::valarray<double> fs_;

public:
  LocalRoot(FuncDeriv1D &func, double lower, double upper, size_t nsamp = 11,
            double xtol_err = 1.0e-8);
  double operator()(double target) const;
};

double find_local_root(FuncDeriv1D &func, double target, double lower,
                       double upper, size_t nsamp = 11,
                       double xtol_err = 1.0e-8);
} // namespace dlt
#endif
