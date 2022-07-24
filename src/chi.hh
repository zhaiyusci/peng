#ifndef _DILUTE_CHI_HH_
#define _DILUTE_CHI_HH_
#include "cgquad.hh"

namespace dlt {
class ReducedPotentialQuadrature;

/**
 * @brief Abstract class for computing chi.
 */
class ChiImpl {
protected:
  ReducedPotentialQuadrature *rpq_;

public:
  ChiImpl(ReducedPotentialQuadrature &rpq) : rpq_(&rpq) {}

  /**
   * @brief Compute chi.
   *
   * Users should implement this method, and add some of their thought here.
   *
   * @param E: Initial kinect energy for collision.
   * @param r_m: Closest distance of the collision.
   * @param rtol: allowed relative error.
   */
  virtual double chi(double E, double r_m, double rtol) = 0;
  virtual ~ChiImpl() {}
};

/**
 * @brief Compute chi using Chebyshev-Gauss quadrature.
 */
class ChiCG : public ChiImpl{
  /**
   * @brief Compute the integration required by the computation of chi.
   */
  class ChiCGInt : public CGIntegrator {
    // {{{
  public:
  private:
    ReducedPotentialQuadrature *rpq_;
    size_t cache_ordersize_;
    std::vector<double> vs_;
    // the following is for recording running status
    double r_m_;
    double E_;
    double b_;

  public:
    ChiCGInt(ReducedPotentialQuadrature &rpq);

    void set_param(double r_m, double E, double b = -1.0);

    /** Compute the integrands with computed values cached.
     */
  void calculate_integrands(size_t ordersize) override;

  void clean_cache() {
    cache_ordersize_ = 0;
    vs_.clear();
    clean_workspace();
    }
   // }}}
  };
  ChiCGInt chicgint;
public:
  ChiCG(ReducedPotentialQuadrature & rpq): ChiImpl(rpq), chicgint(rpq){}
  double chi(double E, double r_m, double rtol) override;
};


} // namespace dlt

#endif
