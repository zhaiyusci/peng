#include "dilute.hh"
#include "atompair.hh"
#include "mathtools.hh"
#include "transport.hh"
#include <fmt/core.h>
#include <fstream>
#include <iomanip>
#include <json.hpp>
#include <string>
#include <vector>

namespace dlt {
Task::Task(const Atom &atom0, const Atom &atom1,
           const std::vector<double> &temperatures,
           const std::vector<double> &molefractions0, size_t maxpq,
           double accuracy, FuncDeriv1D &pot00, FuncDeriv1D &pot01,
           FuncDeriv1D &pot11)
    : atom0_(atom0), atom1_(atom1), temperatures_(temperatures),
      molefractions0_(molefractions0), propertyorder_(maxpq),
      accuracy_(accuracy), pot00_(&pot00), pot01_(&pot01), pot11_(&pot11),
      pair00_(atom0, atom0, pot00), pair11_(atom1, atom1, pot11),
      pair01_(atom0, atom1, pot01),
      //
      // clang-format off
        D12s_   ("Diffusion",       1e-4, "10⁻⁴m²/s", molefractions0.size(), temperatures.size(), maxpq),
        DTs_    ("Therm. Diff.",    1e-4, "10⁻⁴m²/s", molefractions0.size(), temperatures.size(), maxpq),
        lambdas_("Therm. Conduct.", 1e-3, "mW/m⋅K",   molefractions0.size(), temperatures.size(), maxpq),
        etas_   ("Viscosity",       1e-6, "μPa⋅s",    molefractions0.size(), temperatures.size(), maxpq)
// clang-format on
{}

ComputedData &Task::get_property_(size_t i) {
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

void Task::execute() {
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
        Omega00.push_back(
            pair00_.Omega(l + 1, s + 1, temperatures_[ti], accuracy_));
        Omega01.push_back(
            pair01_.Omega(l + 1, s + 1, temperatures_[ti], accuracy_));
        Omega11.push_back(
            pair11_.Omega(l + 1, s + 1, temperatures_[ti], accuracy_));
      }
    }
    std::vector<double> D12;
    std::vector<double> DT;
    std::vector<double> lambda;
    std::vector<double> eta;
    for (size_t xi = 0; xi != molefractions0_.size(); ++xi) {
      std::cout << __FILE__ << ' ' << __LINE__ << std::endl;
      std::tie(D12, DT, lambda, eta) = dlt::transport(
          temperatures_[ti], molefractions0_[xi], Omega00, Omega01, Omega11,
          atom0_.mass(), atom1_.mass(), propertyorder_);
      std::cout << __FILE__ << ' ' << __LINE__ << std::endl;
      for (size_t mi = 0; mi != propertyorder_; ++mi) {
        // clang-format off
          D12s_   .at(xi, ti, mi) = D12   [mi];
          DTs_    .at(xi, ti, mi) = DT    [mi];
          lambdas_.at(xi, ti, mi) = lambda[mi];
          etas_   .at(xi, ti, mi) = eta   [mi];
        // clang-format on
      }
      std::cout << __FILE__ << ' ' << __LINE__ << std::endl;
    }
  }
  return;
}
void Task::chant() {
  for (size_t ix = 0; ix != molefractions0_.size(); ++ix) {
    for (size_t m = 0; m != propertyorder_; ++m) {
      chant_(ix, m);
    }
  }
  return;
}

void Task::chant_(size_t ix, size_t m) {
  std::cout << std::endl;
  const size_t number_of_properties = 4;

  size_t linelen = number_of_properties * 26 + 15;
  std::cout << format(fmt::runtime(fmt::format("{{:^{}s}}", linelen - 15) +
                                   std::string(" (of order {:>3d})")),
                      "Results", m + 1)
            << '\n';
  std::cout << std::string(linelen, '-') << '\n';
  std::cout << fmt::format("For binary mixture formed by {} and {}, with mole "
                           "fraction {:10.5f} and {:10.5f}.",
                           atom0_.symbol(), atom1_.symbol(),
                           molefractions0_[ix], 1.0 - molefractions0_[ix])
            << "\n";
  std::cout << std::string(linelen, '=') << '\n';

  std::cout << fmt::format("{:>15s}", "Temperature [K]");
  for (size_t ip = 0; ip != number_of_properties; ++ip) {
    std::cout << fmt::format("{:>26s}", get_property_(ip).name_ + " [" +
                                            get_property_(ip).unit_ + "]");
  }
  std::cout << '\n';
  std::cout << std::string(linelen, '-') << std::endl;
  for (size_t it = 0; it != temperatures_.size(); ++it) {
    std::cout << fmt::format("{:>15.2f}", temperatures_[it]);
    for (size_t ip = 0; ip != number_of_properties; ++ip) {
      std::cout << fmt::format("{:>26.8g}", get_property_(ip).at(ix, it, m) /
                                                get_property_(ip).coeff_);
    }
    std::cout << '\n';
  }
  std::cout << std::string(linelen, '-') << std::endl;
  std::cout << std::endl;

  return;
}
} // namespace dlt
