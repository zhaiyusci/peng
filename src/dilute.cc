#include "dilute.hh"
#include "atompair.hh"
#include "mathtools.hh"
#include "param.hh"
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

  double Omega00[MAXORD * MAXORD];
  double Omega01[MAXORD * MAXORD];
  double Omega11[MAXORD * MAXORD];
  for (auto &&T : temperatures) {
    for (size_t l = 0; l != maxpq; ++l) {
      for (size_t s = 0; s != maxpq; ++s) {
        Omega00[l * MAXORD + s] = (l + 1) * 10 + s + 1;
        Omega01[l * MAXORD + s] = (l + 1) * 10 + s + 1;
        Omega11[l * MAXORD + s] = (l + 1) * 10 + s + 1;
        // << pair00.Omega(l + 1, s + 1, T) << endl;
      }
    }
  }

  // Atom atom0;
  return 0;
}
