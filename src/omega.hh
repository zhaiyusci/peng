#ifndef _DILUTE_OMEGA_HH_
#define _DILUTE_OMEGA_HH_
#include "cgquad.hh"
#include "glquad.hh"

namespace dlt {
class ReducedPotentialQuadrature;
class OmegaImpl {
protected:
  ReducedPotentialQuadrature *rpq_;

public:
  OmegaImpl(ReducedPotentialQuadrature &rpq) : rpq_(&rpq) {}

  ///
  /// Compute Omega. User can add some of their thought here.
  ///
  virtual double Omega(size_t l, size_t s, double T, double rtol, size_t maxorder) = 0;
};

class OmegaGL: public OmegaImpl{

  ///
  /// Compute the integration required by the computation of Omega.
  ///
  class OmegaGLInt : public GLIntegrator { 
    // {{{
  protected:
    ReducedPotentialQuadrature *rpq_;
    size_t l_;
    size_t s_;
    double T_;

  public:
    OmegaGLInt(ReducedPotentialQuadrature &rpq);
    /** Set the parameters, clear the inner storage if it is needed.
     */
    void set_param(size_t l, size_t s, double T);

    /** Compute the integrands with computed values cached.
     */
    void calculate_integrands(double rtol) override;

  // }}}
  }; 

  OmegaGLInt omegaglint_;
  

public:
  OmegaGL(ReducedPotentialQuadrature & rpq): OmegaImpl(rpq), omegaglint_(rpq){}
  double Omega(size_t l, size_t s, double T, double rtol, size_t maxorder) override;
};


class OmegaCG: public OmegaImpl{

  ///
  /// Compute the integration required by the computation of Omega.
  ///
  class OmegaCGInt : public CGIntegrator { //{{{
  protected:
    ReducedPotentialQuadrature *rpq_;
    size_t ordersize_v_;
    size_t ordersize_Q_;
    std::vector<double> Qs_;
    std::vector<double> vs_;
    std::vector<double> dvs_;
    // the following is for recording running status
    // const double r_E_;
    size_t l_;
    size_t s_;
    double T_;

  public:
    OmegaCGInt(ReducedPotentialQuadrature &rpq);
    /** Set the parameters, clear the inner storage if it is needed.
     */
    void set_param(size_t l, size_t s, double T);

    /** Compute the integrands with computed values cached.
     */
    void calculate_integrands(size_t ordersize, double rtol) override;

    void clean_cache_v() {
      vs_.clear();
      dvs_.clear();
      ordersize_v_ = 0;
      // Always clean workspace of the base class.
      clean_workspace();
    }

    void clean_cache_Q() {
      Qs_.clear();
      ordersize_Q_ = 0;
      // Always clean workspace of the base class.
      clean_workspace();
    }
  }; // }}}

  OmegaCGInt omegacgint_;
  

public:
  OmegaCG(ReducedPotentialQuadrature & rpq): OmegaImpl(rpq), omegacgint_(rpq){}
  double Omega(size_t l, size_t s, double T, double rtol, size_t maxorder) override;
};
} // namespace dlt

#endif
