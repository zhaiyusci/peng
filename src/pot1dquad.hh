#ifndef _DILUTE_POT1DQUAD_HH_
#define _DILUTE_POT1DQUAD_HH_
// #include "atompair.hh"
#include "cgquad.hh"
#include "glquad.hh"
#include "mathtools.hh"
#include <iostream>
#include <memory>
#include <tuple>
#include <vector>
#include "chi.hh"
#include "q.hh"
#include "omega.hh"

namespace dlt {

/// The quadrature algorithm for reduced potential.
class ReducedPotentialQuadrature {
  // 1. The functions need to find the numerical roots / minima / maxima

  ///
  /// Phi functions
  ///
  class PhiEff : public FuncDeriv1D { // {{{
  protected:
    FuncDeriv1D *ppot_;
    const double b_;
    const double E_;

  public:
    PhiEff(FuncDeriv1D &pot, double b, double E);
    double value(double r) const;
    double derivative(double r) const;
    bool provide_derivative() const { return ppot_->provide_derivative(); }
  };
  // }}}

  class Y : public FuncDeriv1D { // {{{
  protected:
    FuncDeriv1D *const ppot_;

  public:
    Y(FuncDeriv1D &pot) : ppot_(&pot) {}
    double value(double r) const;
  };
  // }}}

  // 2. Worker class that actually do the inegration, inheriated from
  // CGIntegrator

  // 3. The members of THE class
protected:
  FuncDeriv1D *const p_reduced_pot_;
  double r_C_;
  double E_C_;
  std::unique_ptr<LocalRoot> y_root_;
  std::unique_ptr<LocalRoot> v_root_;
  std::unique_ptr<Y> y_;
  std::unique_ptr<ChiImpl> p_chi_impl_;
  std::unique_ptr<QImpl> p_q_impl_;
  std::unique_ptr<OmegaImpl> p_omega_impl_;

public:
  ReducedPotentialQuadrature(FuncDeriv1D &reduced_pot);

  inline double potential_value(double r){
    return p_reduced_pot_->value(r);
  }
  inline double potential_derivative(double r){
    return p_reduced_pot_->derivative(r);
  }
  inline double v_root(double E) {return (*v_root_)(E);}
  inline double r_C() const { return r_C_; }
  inline double E_C() const { return E_C_; }

  std::tuple<double, double> r_range(double E) const;
  double r2b(double r, double E) const;

  ///
  /// Compute chi.
  ///
  double chi(double E, double r_m, double rtol);

  ///
  /// Compute Q.
  ///
  /// Tip: It is faster to keep the r_E unchanged and scan the l.
  ///
  double Q(size_t l, double r_E, double E, double rtol);

  ///
  /// Compute Omega.
  ///
  /// Tip: It is faster to keep the l and T unchanged and scan the s.
  ///
  double Omega(size_t l, size_t s, double T, double rtol);

  ///
  /// Set up algorithms.
  ///
  template <typename TConcreteChi, typename TConcreteQ, typename TConcreteOmega>
  void set_algorithm() {
    p_chi_impl_.reset(new TConcreteChi(*this));
    p_q_impl_.reset(new TConcreteQ(*this));
    p_omega_impl_.reset(new TConcreteOmega(*this));
    return;
  }
};
} // namespace dlt

#endif
