#include "atompair.hh"

namespace dlt {
double AtomPair::Omega(size_t l, size_t s, double T, double rtol) const {

  const double kB = 1.380649e-23;      // BY DEFINITION
  const double amu = 1.6605390666e-27; // CODATA2018
  const double AA = 1.e-10;            // BY DEFINITION
  double Omegastar = rpq_->Omega(l, s, T / pf_->epsilon(), rtol);
  // std::cout << "Omega* = " << Omegastar << std::endl;
  return Omegastar * 0.5 * std::tgamma(s + 2) *
         (1.0 - 0.5 * (1.0 + pow(-1, l)) / (1. + l)) * M_PI *
         pow((pf_->sigma() * AA), 2) /
         sqrt(2 * M_PI * reduced_mass_ * amu / (kB * T));
}
} // namespace dlt
// PotentialLibrary *PotentialLibrary::instance_ = nullptr;
