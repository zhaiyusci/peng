#include<vector>
#include"param.h"
#include"global.h"
#include<iostream>

using std::vector;
using std::string;

extern "C"{
  void beta_(double*, int*, double [2], double [MAXORD][MAXORD], double [MAXORD][MAXORD], double [MAXORD][MAXORD], double []);
  //            T     MAXPQ      X          Omega11                  omega12                    omega22              eta
}

vector<double> beta(double t, double x[2], double om11[MAXORD][MAXORD], double om12[MAXORD][MAXORD], double om22[MAXORD][MAXORD]){
  vector<double> eta(maxpq);
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

