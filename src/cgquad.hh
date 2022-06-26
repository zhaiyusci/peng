#include <chrono>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <vector>

namespace dlt {

// size_t half_stosize(size_t order) { return (pow(3, order) + 1) / 2; }
// This class used to iterate over the n^3 data...
class CubicIter {
private:
  size_t border_;
  size_t eorder_; // eorder is not included
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
    iterator &operator++() {
      ++totidx_;
      ++idx_;
      if (order_ == 0) {
        idx_ = 0;
        negative_ = 0;
        order_ = 1;
        idxsize_ = 1;
        return *this;
      }
      if (idx_ != idxsize_) {
        return *this;
      }
      idx_ = 0;
      ++negative_;
      if ((!ci_->half_visit_) && negative_ == 1 && ci_->half_storage_)
        totidx_ -= idxsize_;
      if (negative_ != (ci_->half_visit_ ? 0 : 1) + 1) {
        return *this;
      }
      idx_ = 0;
      negative_ = 0;
      ++order_;
      idxsize_ = pow(3, order_ - 1);
      return *this;
    }
    size_t order() { return order_; }
    bool negative() { return negative_ == 1; }
    size_t idx() { return idx_; }
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

    iterator(CubicIter *ci, size_t order, size_t negative, size_t idx)
        : ci_(ci), order_(order), negative_(negative), idx_(idx) {
      if (order == 0) {
        totidx_ = 0;
        idxsize_ = 1;
      } else {
        idxsize_ = pow(3, order - 1);
        if (ci_->half_storage_) {
          totidx_ = (pow(3, order - 1) + 1) / 2 + idx;
        } else {
          totidx_ = pow(3, order - 1) + negative * idxsize_ + idx;
        }
      }
    }
  };
  iterator begin() { return iterator(this, border_, 0, 0); }
  iterator end() { return iterator(this, eorder_, 0, 0); }
};

class CGIntegrator {
protected:
  double a_;
  double b_;
  double k0_;
  double k1_;
  bool symm_;
  std::vector<double> integrands_;

public:
  CGIntegrator(double a, double b, bool symm);
  std::tuple<double, double> integrate(double tol, size_t ordersize);
  virtual void calculate_integrands(size_t ordersize) = 0;
  double map_pm1(double x);
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
