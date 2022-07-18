#ifndef _DILUTE_OMEGA_HH_
#define _DILUTE_OMEGA_HH_
#include "cgquad.hh"
#include "glquad.hh"

namespace dlt {
class ReducedPotentialQuadrature;

/**
 * @brief Abstract class for computing Omega.
 */
class OmegaImpl {
protected:
  ReducedPotentialQuadrature *rpq_;

public:
  OmegaImpl(ReducedPotentialQuadrature &rpq) : rpq_(&rpq) {}

  /**
   * @brief Compute omega. User can add some of their thought here.
   *
   * @param l: Omega's order 1.
   * @param s: Omega's order 2.
   * @param T: temperature.
   * @param rtol: allowed relative error.
   */
  virtual double Omega(size_t l, size_t s, double T, double rtol) = 0;
};

/**
 * @brief Compute Omega using Gauss-Laguerre quadratures.
 */
class OmegaGL: public OmegaImpl{

  /**
   * @brief Compute the integration required by the computation of Omega.
   */
  class OmegaGLInt : public GLIntegrator { 
    // {{{
  protected:
    ReducedPotentialQuadrature *rpq_;
    size_t l_;
    size_t s_;
    double T_;

  public:
    double q_rtol;
    OmegaGLInt(ReducedPotentialQuadrature &rpq);
    /** Set the parameters, clear the inner storage if it is needed.
     */
    void set_param(size_t l, size_t s, double T);

    /** Compute the integrands with computed values cached.
     */
    void calculate_integrands(size_t ngridsize) override;

  // }}}
  }; 

  OmegaGLInt omegaglint_;
  

public:
  OmegaGL(ReducedPotentialQuadrature & rpq): OmegaImpl(rpq), omegaglint_(rpq){}
  double Omega(size_t l, size_t s, double T, double rtol) override;
};

/**
 * @brief Compute Omega using Chebyshev-Gauss quadratures.
 */
class OmegaCG: public OmegaImpl{

  /**
   * @brief Compute the integration required by the computation of Omega.
   */
  class OmegaCGInt : public CGIntegrator { 
    // {{{
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
    double q_rtol;
    OmegaCGInt(ReducedPotentialQuadrature &rpq);
    /** Set the parameters, clear the inner storage if it is needed.
     */
    void set_param(size_t l, size_t s, double T);

    /** Compute the integrands with computed values cached.
     */
    void calculate_integrands(size_t ordersize) override;

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
  // }}}
  }; 

  OmegaCGInt omegacgint_;
  

public:
  OmegaCG(ReducedPotentialQuadrature & rpq): OmegaImpl(rpq), omegacgint_(rpq){}
  double Omega(size_t l, size_t s, double T, double rtol) override;
};
} // namespace dlt

#endif
