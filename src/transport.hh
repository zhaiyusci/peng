#ifndef _PENG_TRANSPORT_HH_
#define _PENG_TRANSPORT_HH_
#include "atompair.hh"
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <tuple>
#include <vector>

/**
 * @brief The top level namespace for the Peng project.
 */
namespace peng {
/**
 * @brief Abstract class defining the interface of a Omega provider.
 *
 * User can provide their version of OmegaComp by inherite from this class.
 * E.g., you can compute the Omegas in other program and store them.
 * @see OmegaCache as an example, which compute Omega on-the-fly.
 */
class OmegaComp {
public:
  /**
   * @brief Get the value of Omega
   *
   * @param l, s : the index.
   * @temperature: temperature in Kelvin.
   */
  virtual double operator()(size_t l, size_t s, double temperature) const = 0;
};

/**
 * @brief A concrete class of OmegaComp.
 */
class OmegaCache : public OmegaComp {
protected:
  AtomPair *const p_atompair_;
  double rtol_;
  mutable std::map<std::tuple<size_t, size_t, double>, double> cache_;

public:
  /**
   * @brief Constructor.
   *
   * @param atompair: The AtomPair object compute Omega directly.
   * @param rtol: rtol for AtomPair.Omega() method.
   */
  OmegaCache(AtomPair &atompair, double rtol);
  double operator()(size_t l, size_t s, double temperature) const override;
};

class TransportProperties {

protected:
  double mass0_;
  double mass1_;

  // Workspace
  OmegaComp *p_Omega00_;
  OmegaComp *p_Omega01_;
  OmegaComp *p_Omega11_;
  double temperature_;
  size_t propertyorder_;

  struct Sexdec {
    double data[16];
    double &operator()(size_t i0, size_t i1, size_t i2, size_t i3) {
      return data[i0 + i1 * 2 + i2 * 4 + i3 * 8];
    }
  };

  Ext2D<Sexdec> Hint;
  Ext2D<Sexdec> Lint;

  double molefraction0_;
  double molefraction1_;

  double mtot_;
  double m0_;
  double m1_;

  inline double mass_(size_t i) const {
    switch (i) {
    case 0:
      return mass0_;
    case 1:
      return mass1_;
    default:
      return 0;
    }
  }

  const OmegaComp &Omega_(size_t i0, size_t i1) const {
    if (i0 != i1) {
      return *p_Omega01_;
    } else if (i0 == 0) {
      return *p_Omega00_;
    } else {
      return *p_Omega11_;
    }
  }

  void update_temperature_(double temperature);

  void update_Hint_Lint_(size_t propertyorder);

  Eigen::MatrixXd compute_raw_D_mat_();

  Eigen::MatrixXd compute_raw_B_mat();

  double bracket_int_H_(int p, int q, size_t i0, size_t i1, size_t oi0,
                        size_t oi1);

  double H12_(int p, int q, size_t i0, size_t i1, const OmegaComp &Omega);

  double H1_(int p, int q, size_t i0, [[maybe_unused]] size_t i1,
             const OmegaComp &Omega);

  double HSG_(int p, int q, size_t i0, size_t i1, const OmegaComp &Omega);

  double A_(int p, int q);

  double bracket_int_L_(int p, int q, size_t i0, size_t i1, size_t oi0,
                        size_t oi1);

  double L12_(int p, int q, size_t i0, size_t i1, const OmegaComp &Omega);

  double L1_(int p, int q, size_t i0, [[maybe_unused]] size_t i1,
             const OmegaComp &Omega);

  double LSG_(int p, int q, size_t i0, size_t i1, const OmegaComp &Omega);

  double B_(int p, int q);

public:
  /**
   * @brief Constructor.
   *
   * @param mass0, mass1: mass of the two species, in amu.
   * @param Omega00, Omega01, Omega11: objects returns Omegas, if l, s, and
   * temperature are provided.
   */
  TransportProperties(double mass0, double mass1, OmegaComp &Omega00,
                      OmegaComp &Omega01, OmegaComp &Omega11);

  /**
   * @brief Compute the transport properties.
   *
   * @param temperature: temperature in Kelvin.
   * @param molefraction0: mole fraction of the first (as in C++, 0th) specie.
   * @param po: properties will be computed in po-th order.
   *
   * @return in tuple, H12, HT, lambda, and eta, all in SI units.
   */
  std::tuple<double, double, double, double>
  evaluate(double temperature, double molefraction0, size_t po);
};

} // namespace peng

#endif
