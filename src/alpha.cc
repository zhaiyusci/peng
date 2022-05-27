#include "global.hh"
#include "param.hh"
#include <iostream>
#include <tuple>
#include <vector>

using std::tuple;
using std::tie;
using std::make_tuple;

using std::vector;
using std::string;

extern "C"{
  void alpha_(double*, int*, double [2], double [MAXORD][MAXORD], double [MAXORD][MAXORD], double [MAXORD][MAXORD], double [], double [], double []);
  //            T     MAXPQ      X          Omega11                  omega12                    omega22              D12        DT        lambda   
}

tuple<vector<double>, vector<double>, vector<double>> alpha(double t, double x[2], double om11[MAXORD][MAXORD], double om12[MAXORD][MAXORD], double om22[MAXORD][MAXORD]){
  vector<double> D12(maxpq);
  vector<double> DT(maxpq);
  vector<double> lambda(maxpq);
  alpha_(&t, &maxpq, x, om11, om12, om22, D12.data(), DT.data(), lambda.data());
  return make_tuple(D12, DT, lambda);
}

void calcalpha(){

  for(int i = 0; i!=ntemp; ++i){
    for(auto&& x : xs){
      vector<double> D12;
      vector<double> DT;
      vector<double> lambda;
      // Results at different orders
      tie(D12, DT, lambda) = alpha(temperatures[i], x.data(), om11[i], om12[i], om22[i]); 
      // Add results to global variables
      D12s.push_back(D12);
      DTs.push_back(DT);
      lambdas.push_back(lambda);
    }
  }
}

