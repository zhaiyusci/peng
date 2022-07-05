#include "dilute.hh"
#include "alpha.hh"
#include "atompair.hh"
#include "mathtools.hh"
#include "param.hh"
#include <fmt/core.h>
#include <fstream>
#include <json.hpp>
#include <string>
#include <vector>

int main() {
  nlohmann::json config;
  std::cin >> config;
  dlt::Atom atom0(config["atoms"][0]["name"].get<std::string>(),
                  config["atoms"][0]["mass"].get<double>());
  dlt::Atom atom1(config["atoms"][1]["name"].get<std::string>(),
                  config["atoms"][1]["mass"].get<double>());
  std::cout << atom0.symbol() << ' ' << atom0.mass() << '\n';
  std::cout << atom1.symbol() << ' ' << atom1.mass() << '\n';
  std::cout << std::flush;

  std::vector<double> temperatures(
      config["temperatures"].get<std::vector<double>>());
  size_t maxpq(config["maxpq"].get<size_t>());
  double accuracy(config["accuracy"].get<size_t>());
  std::cout << "maxpq << ' '<< accuracy" << std::endl;
  std::cout << maxpq << ' ' << accuracy << std::endl;
  std::vector<double> molefractions0(
      config["molefractions0"].get<std::vector<double>>());
  std::vector<double> molefractions1(molefractions0.size());
  for (size_t i = 0; i != molefractions0.size(); ++i) {
    molefractions1[i] = 1.0 - molefractions0[i];
  }

  LoadedPotential pot00(config["potentials"][0]["path"].get<std::string>());
  LoadedPotential pot01(config["potentials"][1]["path"].get<std::string>());
  LoadedPotential pot11(config["potentials"][2]["path"].get<std::string>());
  std::cout << "pot00(3.0) = " << pot00(3.0) << std::endl;
  std::cout << "pot01(3.0) = " << pot01(3.0) << std::endl;
  std::cout << "pot11(3.0) = " << pot11(3.0) << std::endl;

  dlt::AtomPair pair00(atom0, atom0, pot00);
  dlt::AtomPair pair11(atom1, atom1, pot11);
  dlt::AtomPair pair01(atom0, atom1, pot01);

  size_t oo = omegaorder(maxpq);
  size_t omegasize = (1 + oo) * oo / 2;
  for (auto &&T : temperatures) {
    std::vector<double> Omega00;
    std::vector<double> Omega01;
    std::vector<double> Omega11;
    Omega00.reserve(omegasize);
    Omega01.reserve(omegasize);
    Omega11.reserve(omegasize);
    for (size_t l = 0; l != oo; ++l) {
      for (size_t s = l; s != oo; ++s) {
        Omega00.push_back(pair00.Omega(l + 1, s + 1, T));
        Omega01.push_back(pair01.Omega(l + 1, s + 1, T));
        Omega11.push_back(pair11.Omega(l + 1, s + 1, T));
      }
    }
    fmt::print("maxpq = {} \n", maxpq);
    auto alpha_res =
        alpha(T, molefractions0[0], Omega00, Omega01, Omega11, maxpq);
    break;
  }

  // Atom atom0;
  return 0;
}
