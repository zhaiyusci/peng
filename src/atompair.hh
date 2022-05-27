#include <cppitertools/itertools.hpp>
#include <functional>
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
