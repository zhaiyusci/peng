#ifndef _DILUTE_ATOMPAIR_HH_
#define _DILUTE_ATOMPAIR_HH_
// #include "cachedfunc.hh"
#include "mathtools.hh"
#include "pot1dfeat.hh"
#include "pot1dquad.hh"
#include <functional>
#include <map>
#include <string>
#include <valarray>
#include <vector>

namespace dlt {

/**
 * @brief Base class for all types of particles.
 */
class Particle {
public:
private:
  const std::string symbol_;
  const double mass_;

public:
  /**
   * @brief Constructor.
   * Constructor.
   * @param symbol: The unique identifier for the particle.
   * @param mass: Mass of the particle in atomic mass unit (amu).
   */
  Particle(const std::string &symbol, double mass)
      : symbol_(symbol), mass_(mass) {}
  /**
   * Get the symbol.
   */
  inline std::string symbol() const { return symbol_; }
  /**
   * Get the mass in amu.
   */
  inline double mass() const { return mass_; }
  /**
   * "Less then" used by the sorting algorithm of STL.
   */
  bool operator<(Particle rhs) const { return mass_ < rhs.mass(); }
};

/**
 * Monoatomic Molecule.
 */
class Atom : public Particle {
public:
  /**
   * Constructor.
   * @see Particle
   */
  Atom(const std::string &symbol, double mass) : Particle(symbol, mass) {}
};

/**
 * @brief Atomic pair, which hold the info for collision integrals.
 */
class AtomPair {
public:
private:
  const Atom atom0_;
  const Atom atom1_;
  FuncDeriv1D *const ppot_;
  double tol_;

  const double reduced_mass_;
  std::unique_ptr<Pot1DFeatures> pf_;
  std::unique_ptr<ReducedPotentialQuadrature> rpq_;

public:
  AtomPair(const Atom &atom0, const Atom &atom1, FuncDeriv1D &pot)
      : atom0_(atom0), atom1_(atom1), ppot_(&pot),
        reduced_mass_((atom0_.mass() * atom1.mass()) /
                      (atom0_.mass() + atom1_.mass())) {
    pf_.reset(new Pot1DFeatures(*ppot_));
    rpq_.reset(new ReducedPotentialQuadrature(pf_->reduced_potential()));
  }
  /**
   * @brief Set up algorithms.
   *
   * Pass all type info to ReducedPotentialQuadrature.
   * Note that this is a template method and the template parameters must be
   * definite in compile time.
   *
   * @param TConcreteChi: name for chi integrator concrete class.
   * @param TConcreteQ: name for Q integrator concrete class.
   * @param TConcreteOmega: name for Omega integrator concrete class.
   */
  template <typename TConcreteChi, typename TConcreteQ, typename TConcreteOmega>
  void set_algorithm() {
    rpq_->set_algorithm<TConcreteChi, TConcreteQ, TConcreteOmega>();
    return;
  }

  /**
   * @brief Returns the collision integral Omega in SI units.
   *
   * @param l: first order of Omega integral.
   * @param s: second order of Omega integral.
   * @param T: temperature in Kelvin.
   * @param rtol: allowed relative error.
   */
  double Omega(size_t l, size_t s, double T, double rtol) const;
};
} // namespace dlt

#endif
