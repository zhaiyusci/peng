#ifndef __PENG_MATHTOOLS_HH__
#define __PENG_MATHTOOLS_HH__
#include <iostream>
#include <memory>
#include <nlopt.hpp>
#include <valarray>
namespace peng {

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

/**
 * @brief Data structure to hold data that is grow by n^2.
 */
template <typename T> class Ext2D {
protected:
  std::size_t size_;
  std::vector<T> data_;

public:
  Ext2D(std::size_t r = 0, std::size_t c = 0, const T &init = T()) {
    resize(r, c, init);
  }

  std::size_t size() const { return size_ * size_; }
  // As if we take the size in dimension dim.
  std::size_t size(size_t dim) const { return size_; }

  void resize(std::size_t r, std::size_t c, const T &init = T()) {
    // We take two arguments to remind user that this is a 2D data struct.
    if (r == c) {
      size_ = r;
      data_.resize(r * c, init);
    } else {
      throw std::runtime_error("Ext2d only deals with square maxtrix.");
    }
    return;
  }

  T &operator()(std::size_t idx) { return data_[idx]; }
  T &operator[](std::size_t idx) { return data_[idx]; }
  T &operator()(std::size_t r, std::size_t c) { return data_[index(r, c)]; }

  static std::tuple<std::size_t, std::size_t> index(std::size_t idx) {
    std::size_t maxrc = std::floor(std::sqrt(idx));
    std::size_t r = idx - maxrc * maxrc;
    if (r > maxrc) {
      return std::make_tuple(maxrc, 2 * maxrc - r);
    } else {
      return std::make_tuple(r, maxrc);
    }
  }

  static std::size_t index(std::size_t r, std::size_t c) {
    std::size_t idx = std::max(r, c);
    idx *= idx;
    if (r <= c) {
      idx += r;
    } else {
      idx += 2 * r - c;
    }
    return idx;
  }

  class iterator {
  protected:
    Ext2D<T> *p_container_;
    size_t idx_;
    size_t r_;
    size_t c_;

  public:
    iterator(Ext2D<T> *p_container, size_t idx)
        : p_container_(p_container), idx_(idx), r_(0), c_(0) {
      std::tie(r_, c_) = p_container_->index(idx);
    }
    iterator(Ext2D<T> *p_container, size_t r, size_t c)
        : p_container_(p_container), idx_(0), r_(r), c_(c) {
      idx_ = p_container_->index(r, c);
    }
    T &operator*() { return (*p_container_)[idx_]; }

    iterator &operator++() {
      ++idx_;
      if (r_ < c_) {
        ++r_;
      } else if (c_ != 0) {
        --c_;
      } else {
        c_ = r_ + 1;
        r_ = 0;
      }
      return *this;
    }

    bool operator!=(const iterator &rhs) {
      return idx_ != rhs.idx_ || p_container_ != rhs.p_container_;
    }

    bool operator<(const iterator &rhs) { return idx_ < rhs.idx_; }

    std::size_t idx() { return idx_; }
    std::size_t row() { return r_; }
    std::size_t col() { return c_; }
  };

  iterator begin() { return iterator(this, 0); }
  iterator end() { return iterator(this, size()); }
};
} // namespace peng
#endif
