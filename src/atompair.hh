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

///
/// Base class for all types of particles.
///
class Particle {
public:
private:
  const std::string symbol_;
  const double mass_;

public:
  Particle(const std::string &symbol, double mass)
      : symbol_(symbol), mass_(mass) {}
  inline std::string symbol() const { return symbol_; }
  inline double mass() const { return mass_; }
  bool operator<(Particle rhs) const { return mass_ < rhs.mass(); }
};

///
/// Monoatomic Molecule.
///
class Atom : public Particle {
public:
  Atom(const std::string &symbol, double mass) : Particle(symbol, mass) {}
};


///
/// Class for atomic pair, which is hold the info for collision integrals.
/// 
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
  ///
  /// Set up algorithms.
  ///
  template <typename TConcreteChi, typename TConcreteQ, typename TConcreteOmega>
  void set_algorithm() {
    rpq_->set_algorithm<TConcreteChi, TConcreteQ, TConcreteOmega>();
    return;
  }

  ///
  /// Returns the collision integral SI units.
  /// 
  double Omega(size_t l, size_t s, double T, double rtol) const;
};
} // namespace dlt

#endif
