#ifndef _DILUTE_ATOMPAIR_HH_
#define _DILUTE_ATOMPAIR_HH_
// #include "cachedfunc.hh"
#include "mathtools.hh"
#include "pot1dquad.hh"
#include <cppitertools/itertools.hpp>
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

/*
///
/// The Gaseous mixture
///
class DiluteGas {
  static std::valarray<double> normalize(std::valarray<double> vec) {
    vec /= vec.sum();
    return vec;
  }

public:
private:
  const std::vector<Particle> species_;
  const std::valarray<double> mole_fraction_;

public:
  DiluteGas(std::vector<Particle> species, std::valarray<double> mole_fraction)
      : species_(species), mole_fraction_(normalize(mole_fraction)) {}
  inline const std::vector<Particle> species() const { return species_; }
  inline const std::valarray<double> mole_fraction() const {
    return mole_fraction_;
  }
  auto species_info() { return iter::zip(species_, mole_fraction_); }
  auto pairs() { return iter::combinations(species_info(), 2); }
};

// extern std::map<std::pair<Particle, Particle>, Pot1d> POTLIB;
//
class PotentialLibrary {
public:
private:
  std::map<std::pair<Particle, Particle>, FuncDeriv1D> potlib_;
  static PotentialLibrary *instance_;
  PotentialLibrary() {}

public:
  static std::map<std::pair<Particle, Particle>, FuncDeriv1D> &potlib() {
    if (instance_ == nullptr) {
      instance_ = new PotentialLibrary();
    }
    return instance_->potlib_;
  }
};

*/

class AtomPair {
public:
private:
  const Atom atom0_;
  const Atom atom1_;

  FuncDeriv1D *const ppot_;

  const double reduced_mass_;
  std::unique_ptr<Pot1DFeatures> pf_;
  std::unique_ptr<ReducedPotentialQuadrature> rpq_;

public:
  AtomPair(Atom &atom0, Atom &atom1, FuncDeriv1D &pot)
      : atom0_(atom0), atom1_(atom1), ppot_(&pot),
        reduced_mass_((atom0_.mass() * atom1.mass()) /
                      (atom0_.mass() + atom1_.mass())) {
    pf_.reset(new Pot1DFeatures(*ppot_));
    rpq_.reset(new ReducedPotentialQuadrature(pf_->reduced_potential()));
  }
  double Omega(size_t l, size_t s, double T) const {

    const double kB = 1.380649e-23;      // BY DEFINITION
    const double amu = 1.6605390666e-27; // CODATA2018
    const double AA = 1.e-10;            // BY DEFINITION
    double Omegastar = rpq_->Omega(l, s, T / pf_->epsilon());
    std::cout << "Omega* = " << Omegastar << std::endl;
    return Omegastar * 0.5 * std::tgamma(s + 2) *
           (1.0 - 0.5 * (1.0 + pow(-1, l)) / (1. + l)) * M_PI *
           pow((pf_->sigma() * AA), 2) /
           sqrt(2 * M_PI * reduced_mass_ * amu / (kB * T));
  }
};
} // namespace dlt

#endif
