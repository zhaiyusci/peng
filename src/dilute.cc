#include "dilute.hh"
#include "atompair.hh"
#include "mathtools.hh"
#include "param.hh"
#include "transport.hh"
#include <fmt/core.h>
#include <fstream>
#include <iomanip>
#include <json.hpp>
#include <string>
#include <vector>

namespace dlt {
Task::Task(const Atom &atom0, const Atom &atom1,
           const std::vector<double> &temperatures,
           const std::vector<double> &molefractions0, size_t maxpq,
           double accuracy, FuncDeriv1D &pot00, FuncDeriv1D &pot01,
           FuncDeriv1D &pot11)
    : atom0_(atom0), atom1_(atom1), temperatures_(temperatures),
      molefractions0_(molefractions0), maxpq_(maxpq), accuracy_(accuracy),
      pot00_(&pot00), pot01_(&pot01), pot11_(&pot11),
      pair00_(atom0, atom0, pot00, accuracy),
      pair11_(atom1, atom1, pot11, accuracy),
      pair01_(atom0, atom1, pot01, accuracy),
      //
      // clang-format off
        D12s_   ("Diffusion",       1e-4, "10⁻⁴m²/s", molefractions0.size(), temperatures.size(), maxpq),
        DTs_    ("Therm. Diff.",    1e-4, "10⁻⁴m²/s", molefractions0.size(), temperatures.size(), maxpq),
        lambdas_("Therm. Conduct.", 1e-3, "mW/m⋅K",   molefractions0.size(), temperatures.size(), maxpq),
        etas_   ("Viscosity",       1e-6, "μPa⋅s",    molefractions0.size(), temperatures.size(), maxpq)
// clang-format on
{}
} // namespace dlt
