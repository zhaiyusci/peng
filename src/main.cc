#include <dlfcn.h>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <fmt/core.h>

#include "global.h"
#include "param.h"

using fmt::format;

using std::cout;
using std::endl;
using std::fstream;
using std::stringstream;

using std::map;
using std::tuple;
using std::tie;
using std::get;
using std::make_tuple;

using std::vector;
using std::string;

vector<double> comp2chant(const vector<vector<double>> &compdata, size_t ix, size_t m) {
  vector<double> res;
  for (size_t it = 0; it != static_cast<size_t>(ntemp); ++it)
    res.push_back(compdata[it * xs.size() + ix][m - 1]);
  cout << res[0] << endl;
  return res;
}

void chant(vector<string> properties, int ix, int m = maxpq) {
  //           Property             0. data         1. scaling factor 2. unit
  map<string, tuple<vector<double>, double, string>> formatdb = {
      {"Viscosity", make_tuple(comp2chant(etas, ix, m), 1e-6, "μPa⋅s")},
      {"Diffusion", make_tuple(comp2chant(D12s, ix, m), 1e-4, "10⁻⁴m²/s")},
      {"Therm. Diff.", make_tuple(comp2chant(DTs, ix, m), 1e-4, "10⁻⁴m²/s")},
      {"Therm. Conduct.",
       make_tuple(comp2chant(lambdas, ix, m), 1e-3, "mW/m⋅K")},
  };

  auto linelen = properties.size() * 26 + 15;
  cout << format(fmt::runtime(format("{{:^{}s}}", linelen - 15) +
                     string(" (of order {:>3d})")),
                 "Results", m)
       << endl;
  cout << string(linelen, '-') << endl;
  cout << format("For binary mixture formed by {} and {}, with mole "
                 "fraction {:10.5f} and {:10.5f}.",
                 elements[0], elements[1], xs[ix][0], xs[ix][1])
       << endl;
  cout << string(linelen, '=') << endl;

  cout << format("{:>15s}", "Temperature [K]");
  for (auto &&property : properties) {
    cout << format("{:>26s}",
                   property + " [" + get<2>(formatdb[property]) + "]");
  }
  cout << endl;
  cout << string(linelen, '-') << endl;
  for (size_t i = 0; i != temperatures.size(); ++i) {
    cout << format("{:>15.2f}", temperatures[i]);
    for (auto &&property : properties) {
      cout << format("{:>26.8g}", get<0>(formatdb[property])[i] /
                                      get<1>(formatdb[property]));
      //                                              data scaling factor
    }
    cout << endl;
  }
  cout << string(linelen, '-') << endl;

  return;
}

extern double gasinfo_[2]; // mass[2] in beta.F90 and  alpha.F90

double *const mass(gasinfo_);

extern void calcbeta();
extern void calcalpha();

void getom(double om[60][MAXORD][MAXORD], string soname, double mass1,
           double mass2) {
  void *so = dlopen(soname.c_str(), RTLD_NOW);
  cout << "Open user-defined dimer as " << static_cast<void *>(so) << endl;
  auto getomega =
      reinterpret_cast<void (*)(double[60][MAXORD][MAXORD], double[60], int *,
                                double *, double *)>(dlsym(so, "getomega_"));
  getomega(om, temps, &ntemp, &mass1, &mass2);
  dlclose(so);
  return;
}

tuple<string, vector<string>> understand(const string &s, char delim = ' ') {
  vector<string> results;
  vector<string> result;
  string title;
  stringstream ss(s);
  string item;

  while (getline(ss, item, delim)) {
    if (item.length() != 0)
      results.push_back(item);
  }

  if (results.size() > 0) {
    title = results[0];
    result = vector<string>(results.begin() + 1, results.end());
  }

  return tie(title, result);
}

int main() {
  fstream inputfile;
  inputfile.open("in.txt", std::ios::in);
  ntemp = 0;
  string line;

  string title;
  vector<string> words;
  while (getline(inputfile, line)) {
    tie(title, words) = understand(line);
    // for ( auto&& i : words) cout << i << endl;
    //  cout << "== This line ends. ==" <<  endl;
    if (title == "temperatures") {
      for (size_t i = 0; i != words.size(); ++i) {
        temps[i] = stod(words[i]);
        temperatures.push_back(stod(words[i]));
      }
      ntemp = words.size();
    } else if (title == "masses") {
      for (size_t i = 0; i != words.size(); ++i) {
        mass[i] = stod(words[i]);
      }
    } else if (title == "x") {
      vector<double> x = {};
      for (size_t i = 0; i != words.size(); ++i) {
        x.push_back(stod(words[i]));
      }
      xs.push_back(x);
    } else if (title == "accuracy") {
      acc = stod(words[0]);
    } else if (title == "maxpq") {
      maxpq = stoi(words[0]);
    } else if (title == "elements") {
      elements = words;
    }
  }

  // for(auto&&i : x) cout << i << endl;
  // for(auto&&i : mass) cout << i << endl;

  getom(om11, "./omega11.so", mass[1 - 1], mass[1 - 1]);
  getom(om12, "./omega12.so", mass[1 - 1], mass[2 - 1]);
  getom(om22, "./omega22.so", mass[2 - 1], mass[2 - 1]);

  calcbeta();
  calcalpha();

  for (size_t ix = 0; ix != xs.size(); ++ix) {
    chant(vector<string>(
              {"Viscosity", "Diffusion", "Therm. Diff.", "Therm. Conduct."}),
          ix, maxpq);
  }

  return 0;
}
