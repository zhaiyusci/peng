// #include "atompair.hh"
#include "cgquad.hh"
#include "mathtools.hh"
#include <iostream>
#include <memory>
#include <tuple>
#include <vector>

namespace dlt {

///
/// This class find the features of a potential function.
///
class Pot1DFeatures {
private:
  double sigma_;
  double epsilon_;
  double r_min_;
  FuncDeriv1D *const ppot_;
  void init_() {}

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
};

/// The quadrature algorithm for reduced potential.
class ReducedPotentialQuadrature {
  // 1. The functions need to find the numerical roots / minima / maxima

  ///
  /// Phi functions
  ///
  class PhiEff : public FuncDeriv1D {
  private:
    FuncDeriv1D *ppot_;
    const double b_;
    const double E_;

  public:
    PhiEff(FuncDeriv1D &pot, double b, double E);
    double value(double r) const;
    double derivative(double r) const;
    bool provide_derivative() const { return ppot_->provide_derivative(); }
  };

  class Y : public FuncDeriv1D {
  private:
    FuncDeriv1D *const ppot_;

  public:
    Y(FuncDeriv1D &pot) : ppot_(&pot) {}
    double value(double r) const;
  };

  // 2. Worker class that actually do the inegration, inheriated from
  // CGIntegrator

  ///
  /// Compute the integration required by the computation of chi.
  ///
  class ChiCG : public CGIntegrator { // {{{
  private:
    ReducedPotentialQuadrature *rpq_;
    // the following is for recording running status
    const double r_m_;
    const double E_;
    const double b_;

  public:
    ChiCG(ReducedPotentialQuadrature *rpq, double r_m, double E);
    /** Compute the integrands with computed values cached.
     */
    void calculate_integrands(size_t ordersize) override;

  }; // }}}

  ///
  /// Compute the integration required by the computation of Q, range 1.
  ///
  class QCG1 : public CGIntegrator { //{{{
  private:
    ReducedPotentialQuadrature *rpq_;
    size_t cache_ordersize_;
    size_t l_;
    std::vector<double> coschis_;
    std::vector<double> fct2_;
    // the following is for recording running status
    double r_E_;
    double r_Op_;
    double E_;

  public:
    QCG1(ReducedPotentialQuadrature *rpq);

    ///
    /// Set the parameters, clear the inner storage if it is needed.
    ///
    /// Here E should always greater than zero.  If not ptovided, the default
    /// value of E is -1.0, thus we can compute E automatically.
    ///
    void set_param(size_t l, double r_E, double r_Op, double E = -1.0);

    /** Compute the integrands with computed values cached.
     */
    void calculate_integrands(size_t ordersize) override;

    void clean_cache() {
      coschis_.clear();
      fct2_.clear();
      cache_ordersize_ = 0;
      // Always clean workspace of the base class.
      clean_workspace();
    }
  }; // }}}

  ///
  /// Compute the integration required by the computation of Q, range 2.
  ///
  class QCG2 : public CGIntegrator { //{{{
  private:
    ReducedPotentialQuadrature *rpq_;
    size_t cache_ordersize_;
    size_t l_;
    std::vector<double> coschis_;
    std::vector<double> fct2_;
    // the following is for recording running status
    double r_E_;
    double E_;
    double r_O_;

  public:
    QCG2(ReducedPotentialQuadrature *rpq);

    /** Set the parameters, clear the inner storage if it is needed.
     */
    void set_param(size_t l, double r_E, double r_O, double E = -1.0);

    /** Compute the integrands with computed values cached.
     */
    void calculate_integrands(size_t ordersize) override;

    void clean_cache() {
      coschis_.clear();
      fct2_.clear();
      cache_ordersize_ = 0;
      // Always clean workspace of the base class.
      clean_workspace();
    }
  }; // }}}

  ///
  /// Compute the integration required by the computation of Omega.
  ///
  class OmegaCG : public CGIntegrator { //{{{
  private:
    ReducedPotentialQuadrature *rpq_;
    size_t v_ordersize_;
    size_t Q_ordersize_;
    std::vector<double> Qs_;
    std::vector<double> vs_;
    std::vector<double> dvs_;
    // the following is for recording running status
    // const double r_E_;
    size_t l_;
    size_t s_;
    double T_;

  public:
    OmegaCG(ReducedPotentialQuadrature *rpq);
    /** Set the parameters, clear the inner storage if it is needed.
     */
    void set_l_s_T(size_t l, size_t s, double T);

    /** Compute the integrands with computed values cached.
     */
    void calculate_integrands(size_t ordersize) override;
  }; // }}}

private:
  const Pot1DFeatures *const pf_;
  FuncDeriv1D *const ppot_;
  double r_C_;
  double E_C_;
  std::unique_ptr<LocalRoot> y_root_;
  std::unique_ptr<Y> y_;
  OmegaCG omegacg_;
  QCG1 qcg1;
  QCG2 qcg2;

public:
  ReducedPotentialQuadrature(Pot1DFeatures &pf);
  const double &r_C() const { return r_C_; }
  const double &E_C() const { return E_C_; }

  std::tuple<double, double> r_range(double E) const;
  double r2b(double r, double E) const;

  ///
  /// Compute chi.
  ///
  double chi(double E, double r_m);

  ///
  /// Compute Q.
  ///
  /// Tip: It is faster to keep the r_E unchanged and scan the l.
  ///
  double Q(size_t l, double r_E);

  ///
  /// Compute Omega.
  ///
  /// Tip: It is faster to keep the l and T unchanged and scan the s.
  ///
  double Omega(size_t l, size_t s, double T);
};
} // namespace dlt
