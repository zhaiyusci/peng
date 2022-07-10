#include "dilute.hh"
#include "loadedpotential.hh"
#include <iostream>

using namespace dlt;
int main() {
  // clang-format off
  std::cout << "   /) ,  /)          " << '\n'
            << " _(/    //    _/_  _ " << '\n'
            << "(_(__(_(/_(_(_(___(/_" << "Ver. " << DILUTE_VERSION << '\n';
  // clang-format on
  std::cout.flush();

  std::cout << std::setprecision(16);
  std::cerr << std::setprecision(16);
  nlohmann::json config;
  std::cin >> config;

  Atom atom0(config["atoms"][0]["name"].get<std::string>(),
             config["atoms"][0]["mass"].get<double>());
  Atom atom1(config["atoms"][1]["name"].get<std::string>(),
             config["atoms"][1]["mass"].get<double>());
  std::vector<double> temperatures(
      config["temperatures"].get<std::vector<double>>());
  std::vector<double> molefractions0(
      config["molefractions0"].get<std::vector<double>>());
  size_t maxpq(config["maxpq"].get<size_t>());
  double accuracy(config["accuracy"].get<double>());
  LoadedPotential pot00(config["potentials"][0]["path"].get<std::string>());
  LoadedPotential pot01(config["potentials"][1]["path"].get<std::string>());
  LoadedPotential pot11(config["potentials"][2]["path"].get<std::string>());

  dlt::Task task(atom0, atom1, temperatures, molefractions0, maxpq, accuracy,
                 pot00, pot01, pot11);
  task.execute();
  task.chant();

  return 0;
}
