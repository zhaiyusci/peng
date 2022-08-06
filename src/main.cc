#include "peng.hh"
#include "loadedpotential.hh"
#include <fmt/core.h>
#include <iostream>
#include <json.hpp>
// using namespace peng;

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

/**
 * @brief Class that hold the data of a Task.
 */
class Task {
  const peng::Atom atom0_;
  const peng::Atom atom1_;
  const std::vector<double> temperatures_;
  const std::vector<double> molefractions0_;
  const size_t propertyorder_;
  const double accuracy_;

  // Data...
  peng::AtomPair pair00_;
  peng::AtomPair pair11_;
  peng::AtomPair pair01_;

  ComputedData D12s_;
  ComputedData DTs_;
  ComputedData lambdas_;
  ComputedData etas_;

  ComputedData &get_property_(size_t i) {
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

  void chant_(size_t ix, size_t m) {
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
  Task(const peng::Atom &atom0, const peng::Atom &atom1,
       const std::vector<double> &temperatures,
       const std::vector<double> &molefractions0, size_t maxpq, double accuracy,
       peng::FuncDeriv1D &pot00, peng::FuncDeriv1D &pot01,
       peng::FuncDeriv1D &pot11)
      : atom0_(atom0), atom1_(atom1), temperatures_(temperatures),
        molefractions0_(molefractions0), propertyorder_(maxpq),
        accuracy_(accuracy), pair00_(atom0, atom0, pot00),
        pair11_(atom1, atom1, pot11), pair01_(atom0, atom1, pot01),
        //
        // clang-format off
        D12s_   ("Diffusion",       1e-4, "10⁻⁴m²/s", molefractions0.size(), temperatures.size(), maxpq),
        DTs_    ("Therm. Diff.",    1e-4, "10⁻⁴m²/s", molefractions0.size(), temperatures.size(), maxpq),
        lambdas_("Therm. Conduct.", 1e-3, "mW/m⋅K",   molefractions0.size(), temperatures.size(), maxpq),
        etas_   ("Viscosity",       1e-6, "μPa⋅s",    molefractions0.size(), temperatures.size(), maxpq)
  // clang-format on
  {
    pair00_.set_algorithm<peng::ChiCG, peng::QCG, peng::OmegaGL>();
    pair01_.set_algorithm<peng::ChiCG, peng::QCG, peng::OmegaGL>();
    pair11_.set_algorithm<peng::ChiCG, peng::QCG, peng::OmegaGL>();
  }

  /**
   * @brief Do the computation.
   */
  void execute() {
    peng::OmegaCache om00(pair00_, accuracy_);
    peng::OmegaCache om01(pair01_, accuracy_);
    peng::OmegaCache om11(pair11_, accuracy_);
    peng::TransportProperties tp(atom0_.mass(), atom1_.mass(), om00, om01,
                                 om11);

    for (size_t ti = 0; ti != temperatures_.size(); ++ti) {
      std::cerr << "T = " << temperatures_[ti] << std::endl;
      for (size_t xi = 0; xi != molefractions0_.size(); ++xi) {
        std::cerr << "x_0 = " << molefractions0_[xi] << std::endl;
        for (size_t mi = 1; mi <= propertyorder_; ++mi) {
          double D12, DT, lambda, eta;
          std::tie(D12, DT, lambda, eta) =
              tp.evaluate(temperatures_[ti], molefractions0_[xi], mi);
          // clang-format off
          D12s_   .at(xi, ti, mi - 1) = D12   ;
          DTs_    .at(xi, ti, mi - 1) = DT    ;
          lambdas_.at(xi, ti, mi - 1) = lambda;
          etas_   .at(xi, ti, mi - 1) = eta   ;
          // clang-format on
        }
      }
    }
    return;
  }

  /**
   * @brief Print out the final results, i.e., the properties.
   *
   * We name it "chant" to show our respect to Barker, Fock and Smith,
   * in whose work [Phys. Fluids, 7, 897 (1964)]
   * the output subroutine was named "SHOUT".
   */
  void chant() {
    for (size_t ix = 0; ix != molefractions0_.size(); ++ix) {
      for (size_t m = 0; m != propertyorder_; ++m) {
        chant_(ix, m);
      }
    }
    return;
  }
};
int main() {
  // clang-format off
  std::cout << "*=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=*\n" 
            << "{      ,ggggggggggg,                                          }\n"
            << "{     dP^^^88^^^^^^Y8,                               88       }\n"
            << "{     Yb,  88      `8b                              I88I      }\n"
            << "{      `^  88      ,8P                              i88i      }\n"
            << "{          88aaaad8P^                                88       }\n"
            << "{          88^^^^^,ggg,    ,ggg,,ggg,     ,gggg,gg   88       }\n"
            << "{          88    i8^ ^8i  ,8^ ^8P^ ^8,   dP^  ^Y8I   88       }\n"
            << "{          88    I8, ,8I  I8   8I   8I  i8'    ,8I   ii       }\n"
            << "{          88    `YbadP' ,dP   8I   Yb,,d8,   ,d8I            }\n"
            << "{          88   888P^Y8888P'   8I   `Y8P^Y8888P^888  db       }\n"
            << "{                                             ,d8I'  qp       }\n"
            << "{                                           ,dP'8I            }\n"
            << "{                                          ,8^  8I            }\n"
            << "{                                          I8   8I            }\n"
            << "{                                          `8, ,8I            }\n"
            << "{                                           `Y8P^             }\n"
            << "*=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=*\n"
            << "\n"
            << "Version: " << peng::PENG_VERSION << '\n'
            << "\n"
            << "Authors: \n"
            << "  Yu Zhai, You Li, Hui Li, and Frederick R. W. McCourt\n"
            << "\n"
            << "GitHub Repo: \n"
            << "  https://github.com/zhaiyusci/peng\n"
            << "\n"
            << "The authors appreciate the financial support from \n"
            << "  * 2020-JCJQ project (GFJQ2126-007)\n"
            << "  * National Natural Science Foundation of China (22073035) \n"
            << "  * The Program for JLU Computational Interdiscipline Innovative Platform\n"
            ;

  // clang-format on
  std::cout.flush();

  std::cout << std::setprecision(16);
  std::cerr << std::setprecision(16);
  nlohmann::json config;
  std::cin >> config;

  peng::Atom atom0(config["atoms"][0]["name"].get<std::string>(),
                   config["atoms"][0]["mass"].get<double>());
  peng::Atom atom1(config["atoms"][1]["name"].get<std::string>(),
                   config["atoms"][1]["mass"].get<double>());
  std::vector<double> temperatures(
      config["temperatures"].get<std::vector<double>>());
  std::vector<double> molefractions0(
      config["molefractions0"].get<std::vector<double>>());
  size_t propertyorder(config["propertyorder"].get<size_t>());
  double accuracy(config["accuracy"].get<double>());
  peng::LoadedPotential pot00(
      config["potentials"][0]["path"].get<std::string>());
  peng::LoadedPotential pot01(
      config["potentials"][1]["path"].get<std::string>());
  peng::LoadedPotential pot11(
      config["potentials"][2]["path"].get<std::string>());

  Task task(atom0, atom1, temperatures, molefractions0, propertyorder, accuracy,
            pot00, pot01, pot11);
  task.execute();
  task.chant();

  return 0;
}
