#ifndef _DILUTE_DILUTE_HH_
#define _DILUTE_DILUTE_HH_
#include "atompair.hh"
#include "pot1dquad.hh"
#include "transport.hh"
#include <fmt/core.h>
#include <json.hpp>

namespace dlt {
const char DILUTE_VERSION[] = "0.1.1";
// Class to save the transport properties, in an elegent way.
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

// Class that hold the data of a Task.
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

  ComputedData &get_property(size_t i) {
    switch (i) {
    case (0):
      return this->D12s_;
    case (1):
      return this->DTs_;
    case (2):
      return this->lambdas_;
    case (3):
      return this->etas_;
    }
    throw std::runtime_error("No such property.");
  }

public:
  Task(const Atom &atom0, const Atom &atom1,
       const std::vector<double> &temperatures,
       const std::vector<double> &molefractions0, size_t propertyorder,
       double accuracy, FuncDeriv1D &pot00, FuncDeriv1D &pot01,
       FuncDeriv1D &pot11);

  void execute() {
    size_t maxls = omegaorder(propertyorder_);
    size_t omegasize = (1 + maxls) * maxls / 2;
    for (size_t ti = 0; ti != temperatures_.size(); ++ti) {
      std::vector<double> Omega00;
      std::vector<double> Omega01;
      std::vector<double> Omega11;
      Omega00.reserve(omegasize);
      Omega01.reserve(omegasize);
      Omega11.reserve(omegasize);
      for (size_t l = 0; l != maxls; ++l) {
        for (size_t s = l; s != maxls; ++s) {
          Omega00.push_back(pair00_.Omega(l + 1, s + 1, temperatures_[ti]));
          Omega01.push_back(pair01_.Omega(l + 1, s + 1, temperatures_[ti]));
          Omega11.push_back(pair11_.Omega(l + 1, s + 1, temperatures_[ti]));
        }
      }
      std::vector<double> D12;
      std::vector<double> DT;
      std::vector<double> lambda;
      std::vector<double> eta;
      for (size_t xi = 0; xi != molefractions0_.size(); ++xi) {
        std::tie(D12, DT, lambda, eta) =
            transport(temperatures_[ti], molefractions0_[xi], Omega00, Omega01,
                      Omega11, atom0_.mass(), atom1_.mass(), propertyorder_);
        for (size_t mi = 0; mi != propertyorder_; ++mi) {
          // clang-format off
          D12s_   .at(xi, ti, mi) = D12   [mi];
          DTs_    .at(xi, ti, mi) = DT    [mi];
          lambdas_.at(xi, ti, mi) = lambda[mi];
          etas_   .at(xi, ti, mi) = eta   [mi];
          // clang-format on
        }
      }
    }
    return;
  }

  void chant() {
    for (size_t ix = 0; ix != molefractions0_.size(); ++ix) {
      for (size_t m = 0; m != propertyorder_; ++m) {
        chant(ix, m);
      }
    }
    return;
  }

  void chant(size_t ix, size_t m) {
    std::cout << std::endl;
    const size_t number_of_properties = 4;

    size_t linelen = number_of_properties * 26 + 15;
    std::cout << format(fmt::runtime(fmt::format("{{:^{}s}}", linelen - 15) +
                                     std::string(" (of order {:>3d})")),
                        "Results", m + 1)
              << '\n';
    std::cout << std::string(linelen, '-') << '\n';
    std::cout << fmt::format(
                     "For binary mixture formed by {} and {}, with mole "
                     "fraction {:10.5f} and {:10.5f}.",
                     atom0_.symbol(), atom1_.symbol(), molefractions0_[ix],
                     1.0 - molefractions0_[ix])
              << "\n";
    std::cout << std::string(linelen, '=') << '\n';

    std::cout << fmt::format("{:>15s}", "Temperature [K]");
    for (size_t ip = 0; ip != number_of_properties; ++ip) {
      std::cout << fmt::format("{:>26s}", get_property(ip).name_ + " [" +
                                              get_property(ip).unit_ + "]");
    }
    std::cout << '\n';
    std::cout << std::string(linelen, '-') << std::endl;
    for (size_t it = 0; it != temperatures_.size(); ++it) {
      std::cout << fmt::format("{:>15.2f}", temperatures_[it]);
      for (size_t ip = 0; ip != number_of_properties; ++ip) {
        std::cout << fmt::format("{:>26.8g}", get_property(ip).at(ix, it, m) /
                                                  get_property(ip).coeff_);
      }
      std::cout << '\n';
    }
    std::cout << std::string(linelen, '-') << std::endl;
    std::cout << std::endl;

    return;
  }
};

} // namespace dlt

#endif
