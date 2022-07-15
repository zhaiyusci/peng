#ifndef _DILUTE_CHI_HH_
#define _DILUTE_CHI_HH_
#include "cgquad.hh"

namespace dlt {
class ReducedPotentialQuadrature;
class ChiImpl {
protected:
  ReducedPotentialQuadrature *rpq_;

public:
  ChiImpl(ReducedPotentialQuadrature &rpq) : rpq_(&rpq) {}

  ///
  /// Compute chi. User can add some of their thought here.
  ///
  virtual double chi(double E, double r_m, double rtol, size_t maxorder) = 0;
};

class ChiCG : public ChiImpl{
  ///
  /// Compute the integration required by the computation of chi.
  ///
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
    void calculate_integrands(size_t ordersize, double rtol) override;

    void clean_cache() {
      cache_ordersize_ = 0;
      vs_.clear();
      clean_workspace();
    }
   // }}}
  };
  ChiCGInt chicgint;
  ChiCG(ReducedPotentialQuadrature & rpq): ChiImpl(rpq), chicgint(rpq){}

public:
  double chi(double E, double r_m, double rtol, size_t maxorder) override;
};


} // namespace dlt

#endif
