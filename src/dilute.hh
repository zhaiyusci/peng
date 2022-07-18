#ifndef _DILUTE_DILUTE_HH_
#define _DILUTE_DILUTE_HH_
#include "atompair.hh"
#include "pot1dquad.hh"
#include "transport.hh"
#include <fmt/core.h>
// #include <json.hpp>

namespace dlt {
const char DILUTE_VERSION[] = "0.3.0";

namespace {
/**
 * @brief Class to save the transport properties, in an elegent way.
 */
class ComputedData {
protected:
  const std::string name_;
  const double coeff_;
  const std::string unit_;
  // const size_t molefractions_size_;
  const size_t temperatures_size_;
  const size_t propertyorder_;
  double *const data_;

public:
  ComputedData(const std::string &name, double coeff, const std::string &unit,
               size_t molefractions_size, size_t temperatures_size,
               size_t propertyorder)
      : name_(name), coeff_(coeff), unit_(unit),
        temperatures_size_(temperatures_size), propertyorder_(propertyorder),
        data_(new double[molefractions_size * temperatures_size *
                         propertyorder]) {}

  double &at(size_t molefractions_size, size_t temperatures_size,
             size_t propertyorder) {
    return data_[molefractions_size * temperatures_size_ * propertyorder_ +
                 temperatures_size * propertyorder_ + propertyorder];
  }
  ~ComputedData() { delete[] data_; }
  friend class Task;
};

} // namespace

/**
 * @brief Class that hold the data of a Task.
 */
class Task {
  const dlt::Atom atom0_;
  const dlt::Atom atom1_;
  const std::vector<double> temperatures_;
  const std::vector<double> molefractions0_;
  const size_t propertyorder_;
  const double accuracy_;
  FuncDeriv1D *const pot00_;
  FuncDeriv1D *const pot01_;
  FuncDeriv1D *const pot11_;

  // Data...
  dlt::AtomPair pair00_;
  dlt::AtomPair pair11_;
  dlt::AtomPair pair01_;

  ComputedData D12s_;
  ComputedData DTs_;
  ComputedData lambdas_;
  ComputedData etas_;

  ComputedData &get_property_(size_t i);

  void chant_(size_t ix, size_t m);

public:
  /**
   * @brief Constructor.
   *
   * @param atom0, atom1: atoms.
   * @param temperatures: in Kelvin.
   * @param molefractions0: mole fractions for the first specie.
   * @param propertyorder: order of the preperties need computing in
   * Chapman-Enskog solution.
   * @param rtol: allowed relative error.
   * @param pot00, pot01, pot11: Interatomic interaction.
   */
  Task(const Atom &atom0, const Atom &atom1,
       const std::vector<double> &temperatures,
       const std::vector<double> &molefractions0, size_t propertyorder,
       double rtol, FuncDeriv1D &pot00, FuncDeriv1D &pot01, FuncDeriv1D &pot11);

  /**
   * @brief Do the computation.
   */
  void execute();

  /**
   * @brief Print out the final results, i.e., the properties.
   *
   * We name it "chant" to show our respect to Barker, Fock and Smith,
   * in whose work [Phys. Fluids, 7, 897 (1964)]
   * the output subroutine was named "SHOUT".
   */
  void chant();

};

} // namespace dlt

#endif
