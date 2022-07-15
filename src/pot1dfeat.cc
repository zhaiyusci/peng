#include "pot1dfeat.hh"

namespace dlt {

Pot1DFeatures::Pot1DFeatures(FuncDeriv1D &pot)
    : ppot_(&pot), p_reduced_(nullptr) {
  std::tie(r_min_, epsilon_) = find_local_minimum(*ppot_, 1.0, 5.0);
  // TODO: I believe this range is safe...
  epsilon_ *= -1;
  std::cerr << "r_min  = " << r_min_ << '\n';
  std::cerr << "epsilon  = " << epsilon_ << '\n';
  sigma_ = find_local_root(*ppot_, 0.0, 0.5 * r_min_, r_min_);
  std::cerr << "sigma  = " << sigma_ << '\n';
  // TODO: I believe this range is safe...
  return;
}

}
