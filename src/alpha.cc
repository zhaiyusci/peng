#include "alpha.hh"
#include "global.hh"
#include "param.hh"
#include <iostream>
#include <tuple>
#include <vector>

extern "C" {
void alpha_(double *T, int *maxpq, double x[2], double omega11[MAXORD][MAXORD],
            double omega12[MAXORD][MAXORD], double omega22[MAXORD][MAXORD],
            double *mass1, double *mass2, double D12[], double DT[], double lambda[]);
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
alpha(double t, double x0, std::vector<double> Omega00,
      std::vector<double> Omega01, std::vector<double> Omega11, double mass0, double mass1, int maxpq) {
  // Some dirty work: turn C++ omega to fortran 2D array
  // The following arrays should be "FORTRAN-ready"

  double om11[MAXORD][MAXORD];
  double om12[MAXORD][MAXORD];
  double om22[MAXORD][MAXORD];
  for(size_t i = 0; i!= MAXORD; ++i)
    for(size_t j = 0 ; j != MAXORD; ++j){
      om11[i][j] = 0.99999999999;
      om12[i][j] = 0.99999999999;
      om22[i][j] = 0.99999999999;
    }
  omega_cpp2fort(Omega00, om11, omegaorder(maxpq));
  omega_cpp2fort(Omega01, om12, omegaorder(maxpq));
  omega_cpp2fort(Omega11, om22, omegaorder(maxpq));
  double xs[2];
  xs[0] = x0;
  xs[1] = 1.0 - x0;

  std::vector<double> D12(maxpq);
  std::vector<double> DT(maxpq);
  std::vector<double> lambda(maxpq);
  std::cout << "maxpq = " << maxpq << std::endl;
  alpha_(&t, &maxpq, xs, om11, om12, om22, &mass0, &mass1, D12.data(), DT.data(), lambda.data());
  return make_tuple(D12, DT, lambda);
}

/*
void calcalpha() {

  for (int i = 0; i != ntemp; ++i) {
    for (auto &&x : xs) {
      vector<double> D12;
      vector<double> DT;
      vector<double> lambda;
      // Results at different orders
      tie(D12, DT, lambda) =
          alpha(temperatures[i], x.data(), om11[i], om12[i], om22[i]);
      // Add results to global variables
      D12s.push_back(D12);
      DTs.push_back(DT);
      lambdas.push_back(lambda);
    }
  }
}
*/
