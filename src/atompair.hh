#ifndef __DILUTE_ATOMPAIR_HH__
#define __DILUTE_ATOMPAIR_HH__
#include "pot1d.hh"
#include <cppitertools/itertools.hpp>
#include <functional>
#include <map>
#include <string>
#include <valarray>
#include <vector>
/**
 * class Particle
 * Base class for all types of particles.
 */
class Particle {
  private:
    const std::string _symbol;
    const double _mass;

  public:
    Particle(const std::string &symbol, double mass)
        : _symbol(symbol), _mass(mass) {}
    inline std::string symbol() const { return _symbol; }
    inline double mass() const { return _mass; }
    bool operator<(Particle rhs) const { return _mass < rhs.mass(); }
};

/**
 * class Atom
 * Monoatomic Molecule.
 */
class Atom : public Particle {
  public:
    Atom(const std::string &symbol, double mass) : Particle(symbol, mass) {}
};

class DiluteGas {
    static std::valarray<double> normalize(std::valarray<double> vec) {
      vec /= vec.sum();
      return vec;
    }

  private:
    const std::vector<Particle> species_;
    const std::valarray<double> mole_fraction_;

  public:
    DiluteGas(std::vector<Particle> species,
              std::valarray<double> mole_fraction)
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
  private:
    std::map<std::pair<Particle, Particle>, FuncDeriv1D> potlib_;
    static PotentialLibrary *instance_;
    PotentialLibrary() {}

  public:
    static std::map<std::pair<Particle, Particle>, FuncDeriv1D> &getInstance() {
      if (instance_ == nullptr) {
        instance_ = new PotentialLibrary();
      }
      return instance_->potlib_;
    }
};

#endif
