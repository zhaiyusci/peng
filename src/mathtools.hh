#ifndef __DILUTE_MATHTOOLS_HH__
#define __DILUTE_MATHTOOLS_HH__
#include <iostream>
#include <memory>
#include <nlopt.hpp>
#include <valarray>
namespace dlt {

class NegFuncDeriv1D;

/**
 * @brief Abstract class for one-dimension fucntion, with its derivative provide
 * as requred.
 */
class FuncDeriv1D {
private:
  std::unique_ptr<NegFuncDeriv1D> pneg_;

public:
  /**
   * @brief The class is "callable". Returns value(x).
   */
  double operator()(double x) const { return value(x); }
  /**
   * @brief Get the function value when it takes x.
   */
  virtual double value(double x) const = 0;
  /**
   * @brief Get the function's derivative when it takes x.
   */
  virtual double derivative(double x) const;
  /**
   * @brief True if analytical derivative is provided.
   */
  virtual bool provide_derivative() const { return false; };
  /**
   * @brief Returns negative function.
   */
  virtual FuncDeriv1D &operator-();
  virtual ~FuncDeriv1D() {}
};

/**
 * @brief Abstract class representing the negative of a function.
 */
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
};

/**
 * @brief Get the local minimum, based on a grid scanning and optimization
 * algorithm.
 *
 * @param func: the function to be optimized.
 * @param lower, upper: the range of searching function to be optimized.
 * @param nsamp: scanned grid size.
 * @param xtol_err: allowed error.
 *
 * @return tuple of x and f(x) at local minimum.
 */
std::tuple<double, double> find_local_minimum(FuncDeriv1D &func, double lower,
                                              double upper, size_t nsamp = 11,
                                              double xtol_err = 1.0e-8);

/**
 * @brief Get the local maximum, based on a grid scanning and optimization
 * algorithm.
 *
 * @param func: the function to be optimized.
 * @param lower, upper: the range of searching function to be optimized.
 * @param nsamp: scanned grid size.
 * @param xtol_err: allowed error.
 *
 * @return tuple of x and f(x) at local maximum.
 */
std::tuple<double, double> find_local_maximum(FuncDeriv1D &func, double lower,
                                              double upper, size_t nsamp = 11,
                                              double xtol_err = 1.0e-8);

/**
 * @brief Get the local root of f(x) = target, a long-live version of
 * find_local_root.
 */
class LocalRoot {
private:
  FuncDeriv1D *const pfunc_;
  const double xtol_err_;
  std::valarray<double> xs_;
  std::valarray<double> fs_;

public:
  /**
   * @brief Constructor.
   *
   * @param func: the function to be optimized.
   * @param lower, upper: the range of searching function to be optimized.
   * @param nsamp: scanned grid size.
   * @param xtol_err: allowed error.
   */
  LocalRoot(FuncDeriv1D &func, double lower, double upper, size_t nsamp = 11,
            double xtol_err = 1.0e-8);
  /**
   * @brief Find the root for f(x) = target.
   *
   * @param target: allowed error.
   *
   * @return x.
   */
  double operator()(double target) const;
};

/**
 * @brief Get the local root of f(x) = target.
 *
 * @param func: the function to be optimized.
 * @param lower, upper: the range of searching function to be optimized.
 * @param nsamp: scanned grid size.
 * @param xtol_err: allowed error.
 *
 * @return tuple of x and f(x) at local maximum.
 */
double find_local_root(FuncDeriv1D &func, double target, double lower,
                       double upper, size_t nsamp = 11,
                       double xtol_err = 1.0e-8);

/**
 * @brief factorial.
 */
inline double fact(int n) { return std::tgamma(n + 1); }

/**
 * @brief Pochhammer symbol..
 */
inline double poch(int z, int n) { // eq.93
  return tgamma(z + n) / tgamma(z);
}

/**
 * @brief The delta symbol.
 */
inline int delta(int i, int j) { return i == j ? 1 : 0; }
} // namespace dlt
#endif
