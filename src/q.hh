#ifndef _DILUTE_Q_HH_
#define _DILUTE_Q_HH_
#include "cgquad.hh"

namespace dlt {
class ReducedPotentialQuadrature;
class QImpl {
protected:
  ReducedPotentialQuadrature *rpq_;

public:
  QImpl(ReducedPotentialQuadrature &rpq) : rpq_(&rpq) {}

  ///
  /// Compute chi. User can add some of their thought here.
  ///
  virtual double Q(size_t l, double r_E, double E, double rtol, size_t maxorder) =0;
};

class QCG : public QImpl{
  ///
  /// Compute the integration required by the computation of Q, range 1.
  ///
  class QCGInt1 : public CGIntegrator { 
    //{{{
  protected:
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
    QCGInt1(ReducedPotentialQuadrature &rpq);

    ///
    /// Set the parameters, clear the inner storage if it is needed.
    ///
    /// Here E should always greater than zero.  If not ptovided, the default
    /// value of E is -1.0, thus we can compute E automatically.
    ///
    void set_param(size_t l, double r_E, double r_Op, double E = -1.0);

    /** Compute the integrands with computed values cached.
     */
    void calculate_integrands(size_t ordersize, double rtol) override;

    void clean_cache() {
      coschis_.clear();
      fct2_.clear();
      cache_ordersize_ = 0;
      // Always clean workspace of the base class.
      clean_workspace();
    }
  // }}}
  }; 

  ///
  /// Compute the integration required by the computation of Q, range 2.
  ///
  class QCGInt2 : public CGIntegrator { 
    //{{{
  protected:
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
    QCGInt2(ReducedPotentialQuadrature &rpq);

    /** Set the parameters, clear the inner storage if it is needed.
     */
    void set_param(size_t l, double r_E, double r_O, double E = -1.0);

    /** Compute the integrands with computed values cached.
     */
    void calculate_integrands(size_t ordersize, double rtol) override;

    void clean_cache() {
      coschis_.clear();
      fct2_.clear();
      cache_ordersize_ = 0;
      // Always clean workspace of the base class.
      clean_workspace();
    }
  // }}}
  };

  QCGInt1 qcgint1_;
  QCGInt2 qcgint2_;

public:
  QCG(ReducedPotentialQuadrature & rpq): QImpl(rpq), qcgint1_(rpq), qcgint2_(rpq){}
  double Q(size_t l, double r_E, double E, double rtol, size_t maxorder) override;
};


} // namespace dlt

#endif
