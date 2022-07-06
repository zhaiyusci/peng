#include "global.hh"
#include "param.hh"
#include <iostream>
#include <vector>

std::vector<double> beta(double t, double x[2], double om11[MAXORD][MAXORD],
                         double om12[MAXORD][MAXORD],
                         double om22[MAXORD][MAXORD]) {
  std::vector<double> eta(maxpq);
  beta_(&t, &maxpq, x, om11, om12, om22, eta.data());
  return eta;
}

void calcbeta(){

  for(int i = 0; i!=ntemp; ++i){
    for(auto&& x : xs){
      auto eta=beta(temperatures[i], x.data(), om11[i], om12[i], om22[i]); // eta with different order
      etas.push_back(eta);
    }
  }

}

