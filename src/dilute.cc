#include "dilute.hh"
#include "atompair.hh"
#include "mathtools.hh"
#include <fstream>
#include <json.hpp>
#include <string>
#include <vector>

using namespace dlt;
using namespace std;
using namespace nlohmann;

int main() {
  json config;
  cin >> config;
  Atom atom0(config["atoms"][0]["name"].get<string>(),
             config["atoms"][0]["mass"].get<double>());
  Atom atom1(config["atoms"][1]["name"].get<string>(),
             config["atoms"][1]["mass"].get<double>());
  cout << atom0.symbol() << ' ' << atom0.mass() << '\n';
  cout << atom1.symbol() << ' ' << atom1.mass() << '\n';
  cout << flush;

  vector<double> temperatures(config["temperatures"].get<vector<double>>());
  size_t maxpq(config["maxpq"].get<size_t>());
  double accuracy(config["accuracy"].get<size_t>());
  cout << "maxpq << ' '<< accuracy" << endl;
  cout << maxpq << ' ' << accuracy << endl;
  vector<double> molefractions0(config["molefractions0"].get<vector<double>>());
  vector<double> molefractions1(molefractions0.size());
  for (size_t i = 0; i != molefractions0.size(); ++i) {
    molefractions1[i] = 1.0 - molefractions0[i];
  }

  LoadedPotential pot00(config["potentials"][0]["path"].get<string>());
  LoadedPotential pot01(config["potentials"][1]["path"].get<string>());
  LoadedPotential pot11(config["potentials"][2]["path"].get<string>());
  std::cout << "pot00(3.0) = " << pot00(3.0) << std::endl;
  std::cout << "pot01(3.0) = " << pot01(3.0) << std::endl;
  std::cout << "pot11(3.0) = " << pot11(3.0) << std::endl;

  AtomPair pair00(atom0, atom0, pot00);
  AtomPair pair11(atom1, atom1, pot11);
  AtomPair pair01(atom0, atom1, pot01);

  // fstream fs;
  // fs.open("out.txt", std::ios::out);
  // for (double r = 0.5; r <= 3; r += 0.05) {
  // // clang-format off
  // cout << r
  // << '\t' << pair00.rpq_->y_->value(r)
  // << '\t' << pair00.rpq_->p_reduced_pot_->value(r)
  // << '\t' << pair00.rpq_->p_reduced_pot_->derivative(r)
  // << '\t' << pair00.rpq_->p_reduced_pot_->dlt::FuncDeriv1D::derivative(r)
  // << '\t' << pair00.ppot_->derivative(r)
  // << '\t' << pair00.ppot_->dlt::FuncDeriv1D::derivative(r)
  // << '\n';
  // // clang-format on
  // }
  // cout << flush;
  // fs.close();

  cout << pair00.Omega(1, 1, 300) << endl;

  // Atom atom0;
  return 0;
}
