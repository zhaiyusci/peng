#ifndef _DILUTE_POT1DFEAT_HH_
#define _DILUTE_POT1DFEAT_HH_
#include "mathtools.hh"

namespace dlt{
///
/// This class defined the reduced potential energy surface,
/// i.e., the "star" one.
///
class ReducedPotential : public FuncDeriv1D {
public:
private:
  double sigma_;
  double epsilon_;
  double r_min_;
  FuncDeriv1D *const p_pri_pot_;

public:
  ReducedPotential(double sigma, double epsilon, double r_min,
                   FuncDeriv1D *const p_pri_pot)
      : sigma_(sigma), epsilon_(epsilon), r_min_(r_min), p_pri_pot_(p_pri_pot) {
  }

  ///
  /// The collision radius.
  ///
  const double &origin_sigma() const { return sigma_; }

  ///
  /// The well depth.
  ///
  const double &origin_epsilon() const { return epsilon_; }

  ///
  /// The equlibrium nuclear seperation.
  ///
  const double &origin_r_min() const { return r_min_; }

  double value(double r) const override {
    return p_pri_pot_->value(r * sigma_) / epsilon_;
  }

  double derivative(double r) const override {
    return p_pri_pot_->derivative(r * sigma_) / epsilon_ * sigma_;
  }

  bool provide_derivative() const override {
    return p_pri_pot_->provide_derivative();
  }
};

///
/// This class find the features of a potential function.
///
class Pot1DFeatures {
public:
private:
  double sigma_;
  double epsilon_;
  double r_min_;
  FuncDeriv1D *const ppot_;
  std::unique_ptr<ReducedPotential> p_reduced_;

public:
  Pot1DFeatures(FuncDeriv1D &pot);

  ///
  /// The collision radius.
  ///
  const double &sigma() const { return sigma_; }

  ///
  /// The well depth.
  ///
  const double &epsilon() const { return epsilon_; }

  ///
  /// The equlibrium nuclear seperation.
  ///
  const double &r_min() const { return r_min_; }

  ///
  /// The potential.
  ///
  FuncDeriv1D &pot() const { return *ppot_; }

  ///
  /// Reduced potential energy curve.
  ///
  ReducedPotential &reduced_potential() {
    if (p_reduced_ == nullptr) {
      p_reduced_.reset(new ReducedPotential(sigma_, epsilon_, r_min_, ppot_));
    }
    return *p_reduced_;
  }
};

}

#endif
