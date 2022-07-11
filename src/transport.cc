#include "transport.hh"
#include "global.hh"
#include "param.hh"
#include <iostream>
#include <tuple>
#include <vector>

// This file call FORTRAN files alpha_impl.F90 and beta_impl.F90,
// in which the algorithms are located.
extern "C" {
void alpha_(double *T, int *maxpq, double x[2], double omega11[MAXORD][MAXORD],
            double omega12[MAXORD][MAXORD], double omega22[MAXORD][MAXORD],
            double *mass1, double *mass2, double D12[], double DT[], double lambda[]);
}

extern "C" {
void beta_(double *T, int *maxpq, double x[2], double omega11[MAXORD][MAXORD],
           double omega12[MAXORD][MAXORD], double omega22[MAXORD][MAXORD],
           double *mass1, double *mass2, double eta[]);
}

std::tuple<std::vector<double> /*D12*/, std::vector<double> /*DT*/,
           std::vector<double> /*lambda*/, std::vector<double> /*eta*/>
transport(double t, double x0, std::vector<double> Omega00,
          std::vector<double> Omega01, std::vector<double> Omega11,
          double mass0, double mass1, int propertyorder) {
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
  omega_cpp2fort(Omega00, om11, omegaorder(propertyorder));
  omega_cpp2fort(Omega01, om12, omegaorder(propertyorder));
  omega_cpp2fort(Omega11, om22, omegaorder(propertyorder));
  double xs[2];
  xs[0] = x0;
  xs[1] = 1.0 - x0;

  std::vector<double> D12(propertyorder);
  std::vector<double> DT(propertyorder);
  std::vector<double> lambda(propertyorder);
  std::vector<double> eta(propertyorder);
  // std::cout << "maxpq = " << maxpq << std::endl;
  alpha_(&t, &propertyorder, xs, om11, om12, om22, &mass0, &mass1, D12.data(),
         DT.data(), lambda.data());
  beta_(&t, &propertyorder, xs, om11, om12, om22, &mass0, &mass1, eta.data());
  return make_tuple(D12, DT, lambda, eta);
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
