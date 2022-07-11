#ifndef _DILUTE_CGQUAD_HH_
#define _DILUTE_CGQUAD_HH_
#include <chrono>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <vector>

// Following
// https://mathworld.wolfram.com/Chebyshev-GaussQuadrature.html
// Note that we only need ONE copy of the cosines...
// so a CGIntegratorBackend class is prepared
namespace dlt {

///
/// This class is used to iterate over the n^3 data...
///
/// There are two fashions to visit the data.
///
/// 1. In case half_visit == true, only the positive half the data is visited,
/// This can often happen if the function we want to treat is an even function.
/// When this happens, the half_storage must be true because if the data is
/// stored fully, i.e., the negative part of the function is also recorded as
/// grids, It is unreasonable to only visit the positive half.
///
/// 2. When half_visit == false, both positive half and negative half is
/// visited. This case coresponds to two fashions of storage. The full storage
/// will cache all the results in a list, while the half storage only record the
/// positive half. Chances are that we generate the full grids from a
/// half-stored grids, like in the CGIntegrator case.
///
class CubicIter {
private:
  size_t border_;
  size_t eorder_;     // eorder is not included
  bool half_visit_;   // 0: only positive half
                      // 1: also negative half
  bool half_storage_; // in storage, only save the positive half

public:
  CubicIter() = delete;
  CubicIter(size_t eorder)
      : border_(0), eorder_(eorder), half_visit_(false), half_storage_(true) {}
  CubicIter(size_t border, size_t eorder, size_t half_visit = false,
            bool half_storage = true)
      : border_(border), eorder_(eorder), half_visit_(half_visit),
        half_storage_(half_storage) {
    if (half_visit_ && !half_storage_) {
      throw std::runtime_error("You cannot only visit the positive half of the "
                               "elements in the full storage.");
    }
  }
  class iterator {
  protected:
    CubicIter *ci_;
    size_t order_;
    size_t negative_;
    size_t idx_;
    size_t idxsize_;
    size_t totidx_;

  public:
    ///
    /// Used with ++iterator.
    ///
    iterator &operator++();

    ///
    /// Returns the order.
    ///
    size_t order() { return order_; }
    bool negative() { return negative_ == 1; }

    ///
    /// Returns the index, 0-based, of the current order and sign
    ///(if fully visit the space, i.e., half_visit_ == false).
    ///
    size_t idx() { return idx_; }

    ///
    /// Returns the total index, but the sign is also returned.
    /// You can seperate the sign and index easily.
    ///
    int operator*() const {
      if (ci_->half_storage_) {
        return totidx_ * (negative_ ? -1 : 1);
      } else {
        return totidx_;
      }
    }
    bool operator==(iterator &rhs) const {
      return order_ == rhs.order_ && negative_ == rhs.negative_ &&
             idx_ == rhs.idx_;
    }
    bool operator!=(iterator &rhs) const { return !(operator==(rhs)); }

    iterator(CubicIter *ci, size_t order, size_t negative, size_t idx);
  };

  ///
  /// Returns the begin iterator, following the STL fashion.
  ///
  iterator begin() { return iterator(this, border_, 0, 0); }

  ///
  /// Returns the end iterator, not really visited, following the STL fashion.
  ///
  iterator end() { return iterator(this, eorder_, 0, 0); }

private:
  size_t size_from_0_(size_t order) {
    if (order == 0) {
      return 0;
    }
    if (half_visit_) {
      return (pow(3, order - 1) + 1) / 2;
    } else {
      return pow(3, order - 1);
    }
  }

public:
  ///
  /// The size of the Iterable object if counting from order 0.
  ///
  size_t size_from_0() { return size_from_0_(eorder_); }

  ///
  /// The size of the Iterable object.
  /// In my context, I do not really need this method.
  /// TODO: check it.
  ///
  size_t size() { return size_from_0_(eorder_) - size_from_0_(border_); }
};

class CGIntegrator {
public:
  // protected:
  const bool symm_;
  double a_;
  double b_;
  double k0_;
  double k1_;

  // Working space
  size_t ordersize_;
  std::vector<double> integrands_;

public:
  CGIntegrator(bool symm, double a = -1.0, double b = 1.0);
  void set_a_b(double a, double b);
  std::tuple<double, double, bool> integrate(double rtol, size_t ordersize);
  void clean_workspace() {
    ordersize_ = 0;
    integrands_.clear();
  }

  ///
  /// User implemented method. Basically, user should
  ///
  /// 1. Fill in integrands_ according the rule.  The internal iterating
  /// logic is implemented with CubicIter class, the document of which is
  /// provided;
  ///
  /// 2. Update ordersize_ together with integrands_;
  ///
  /// 3. Design an algorithm with cache mechanism if needed to save the
  /// computational resource.
  ///
  virtual void calculate_integrands(size_t ordersize) = 0;

  ///
  /// Map the variable in [-1,1] to [a,b].
  ///
  double map_pm1(double x);
  void show_integrands() {
    std::cerr << "Integrands : " << ' ';
    for (auto &&v : integrands_) {
      std::cerr << v << ' ';
    }
    std::cerr << std::endl;
    return;
  }
};

class CGIntegratorBackend {

protected:
  std::vector<double> angles_;
  std::vector<double> coss_;
  // Storage
  size_t ordersize_; // ordersize is not included following C++ convention
  double gap_;       // keep track of the gap between neighbouring two pivots
  static CGIntegratorBackend *instance_;
  CGIntegratorBackend();
  friend CGIntegrator;

public:
  static CGIntegratorBackend *instance();
  void allocate(size_t order);
  double coss(int i);
};

} // namespace dlt

#endif
