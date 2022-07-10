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

// Class to save the transport properties, in an elegent way.
class ComputedData {
protected:
  const std::string name_;
  const double coeff_;
  const std::string unit_;
  // const size_t molefractions_size_;
  const size_t temperatures_size_;
  const size_t maxpq_;
  double *const data_;

public:
  ComputedData(const std::string &name, double coeff, const std::string &unit,
               size_t molefractions_size, size_t temperatures_size,
               size_t maxpq)
      : name_(name), coeff_(coeff), unit_(unit),
        temperatures_size_(temperatures_size), maxpq_(maxpq),
        data_(new double[molefractions_size * temperatures_size * maxpq]) {}

  double &at(size_t molefractions_size, size_t temperatures_size,
             size_t maxpq) {
    return data_[molefractions_size * temperatures_size_ * maxpq_ +
                 temperatures_size * maxpq_ + maxpq];
  }
  ~ComputedData() { delete[] data_; }
  friend class Task;
};

// Class that hold the data of a Task.
class Task {
  const dlt::Atom atom0;
  const dlt::Atom atom1;
  const std::vector<double> temperatures;
  const std::vector<double> molefractions0;
  const size_t maxpq;
  const double accuracy;
  LoadedPotential pot00;
  LoadedPotential pot01;
  LoadedPotential pot11;

  // Data...
  dlt::AtomPair pair00;
  dlt::AtomPair pair11;
  dlt::AtomPair pair01;

  ComputedData D12s;
  ComputedData DTs;
  ComputedData lambdas;
  ComputedData etas;

  ComputedData &get_property(size_t i) {
    switch (i) {
    case (0):
      return this->D12s;
    case (1):
      return this->DTs;
    case (2):
      return this->lambdas;
    case (3):
      return this->etas;
    }
    throw std::runtime_error("No such property.");
  }

public:
  Task(nlohmann::json config)
      : atom0(config["atoms"][0]["name"].get<std::string>(),
              config["atoms"][0]["mass"].get<double>()),
        atom1(config["atoms"][1]["name"].get<std::string>(),
              config["atoms"][1]["mass"].get<double>()),
        temperatures(config["temperatures"].get<std::vector<double>>()),
        molefractions0(config["molefractions0"].get<std::vector<double>>()),
        maxpq(config["maxpq"].get<size_t>()),
        accuracy(config["accuracy"].get<size_t>()),
        pot00(config["potentials"][0]["path"].get<std::string>()),
        pot01(config["potentials"][1]["path"].get<std::string>()),
        pot11(config["potentials"][2]["path"].get<std::string>()),
        pair00(atom0, atom0, pot00, config["accuracy"].get<double>()),
        pair11(atom1, atom1, pot11, config["accuracy"].get<double>()),
        pair01(atom0, atom1, pot01, config["accuracy"].get<double>()),
        //
        // clang-format off
        D12s   ("Diffusion",       1e-4, "10⁻⁴m²/s", molefractions0.size(), temperatures.size(), maxpq),
        DTs    ("Therm. Diff.",    1e-4, "10⁻⁴m²/s", molefractions0.size(), temperatures.size(), maxpq),
        lambdas("Therm. Conduct.", 1e-3, "mW/m⋅K",   molefractions0.size(), temperatures.size(), maxpq),
        etas   ("Viscosity",       1e-6, "μPa⋅s",    molefractions0.size(), temperatures.size(), maxpq)
  // clang-format on
  {}

  void execute() {
    size_t maxls = omegaorder(maxpq);
    size_t omegasize = (1 + maxls) * maxls / 2;
    for (size_t ti = 0; ti != temperatures.size(); ++ti) {
      std::vector<double> Omega00;
      std::vector<double> Omega01;
      std::vector<double> Omega11;
      Omega00.reserve(omegasize);
      Omega01.reserve(omegasize);
      Omega11.reserve(omegasize);
      for (size_t l = 0; l != maxls; ++l) {
        for (size_t s = l; s != maxls; ++s) {
          Omega00.push_back(pair00.Omega(l + 1, s + 1, temperatures[ti]));
          Omega01.push_back(pair01.Omega(l + 1, s + 1, temperatures[ti]));
          Omega11.push_back(pair11.Omega(l + 1, s + 1, temperatures[ti]));
        }
      }
      // fmt::print("maxpq = {} \n", maxpq);
      std::vector<double> D12;
      std::vector<double> DT;
      std::vector<double> lambda;
      std::vector<double> eta;
      for (size_t xi = 0; xi != molefractions0.size(); ++xi) {
        std::tie(D12, DT, lambda, eta) =
            transport(temperatures[ti], molefractions0[xi], Omega00, Omega01,
                      Omega11, atom0.mass(), atom1.mass(), maxpq);
        for (size_t mi = 0; mi != maxpq; ++mi) {
          // clang-format off
          D12s   .at(xi, ti, mi) = D12   [mi];
          DTs    .at(xi, ti, mi) = DT    [mi];
          lambdas.at(xi, ti, mi) = lambda[mi];
          etas   .at(xi, ti, mi) = eta   [mi];
          // clang-format on
        }
      }
    }
    return;
  }

  void chant() {
    for (size_t ix = 0; ix != molefractions0.size(); ++ix) {
      for (size_t m = 0; m != maxpq; ++m) {
        chant(ix, m);
      }
    }
    return;
  }

  void chant(size_t ix, size_t m) {
    std::cout << std::endl;
    const size_t number_of_properties = 4;

    size_t linelen = number_of_properties * 26 + 15;
    std::cout << format(fmt::runtime(fmt::format("{{:^{}s}}", linelen - 15) +
                                     std::string(" (of order {:>3d})")),
                        "Results", m + 1)
              << '\n';
    std::cout << std::string(linelen, '-') << '\n';
    std::cout << fmt::format(
                     "For binary mixture formed by {} and {}, with mole "
                     "fraction {:10.5f} and {:10.5f}.",
                     atom0.symbol(), atom1.symbol(), molefractions0[ix],
                     1.0 - molefractions0[ix])
              << "\n";
    std::cout << std::string(linelen, '=') << '\n';

    std::cout << fmt::format("{:>15s}", "Temperature [K]");
    for (size_t ip = 0; ip != number_of_properties; ++ip) {
      std::cout << fmt::format("{:>26s}", get_property(ip).name_ + " [" +
                                              get_property(ip).unit_ + "]");
    }
    std::cout << '\n';
    std::cout << std::string(linelen, '-') << std::endl;
    for (size_t it = 0; it != temperatures.size(); ++it) {
      std::cout << fmt::format("{:>15.2f}", temperatures[it]);
      for (size_t ip = 0; ip != number_of_properties; ++ip) {
        std::cout << fmt::format("{:>26.8g}", get_property(ip).at(ix, it, m) /
                                                  get_property(ip).coeff_);
      }
      std::cout << '\n';
    }
    std::cout << std::string(linelen, '-') << std::endl;
    std::cout << std::endl;

    return;
  }
};

int main() {
  // clang-format off
  std::cout << "   /) ,  /)          " << '\n'
            << " _(/    //    _/_  _ " << '\n'
            << "(_(__(_(/_(_(_(___(/_" << "Ver. " << _VERSION_ << '\n';
  // clang-format on
  std::cout.flush();

  std::cout << std::setprecision(16);
  std::cerr << std::setprecision(16);
  nlohmann::json config;
  std::cin >> config;

  if (1) {
    Task task(config);
    task.execute();
    task.chant();
  }

  if (0) {
    LoadedPotential hehe("./hehe.so");
    dlt::Pot1DFeatures hehepf(hehe);
    dlt::ReducedPotentialQuadrature heherpq(hehepf.reduced_potential());
    std::cout << heherpq.chi(5.286831024184222e-05, 6.543886890852044)
              << std::endl;

    LoadedPotential hexe("./hexe.so");
    dlt::Pot1DFeatures hexepf(hexe);
    dlt::ReducedPotentialQuadrature hexerpq(hexepf.reduced_potential());
    std::cout << hexerpq.chi(5.286831024184222e-05, 6.543886890852044)
              << std::endl;

    LoadedPotential xexe("./xexe.so");
    dlt::Pot1DFeatures xexepf(xexe);
    dlt::ReducedPotentialQuadrature xexerpq(xexepf.reduced_potential());
    std::cout << xexerpq.chi(5.286831024184222e-05, 6.543886890852044)
              << std::endl;

    std::cout << 1 << ' ' << 1 << ' ' << heherpq.Omega(1, 1, 4.5514316911066173)
              << std::endl;
    std::cout << 2 << ' ' << 2 << ' ' << heherpq.Omega(2, 2, 4.5514316911066173)
              << std::endl;
    std::cout << 1 << ' ' << 6 << ' ' << heherpq.Omega(1, 6, 4.5514316911066173)
              << std::endl;
    std::cout << 4 << ' ' << 4 << ' ' << heherpq.Omega(4, 4, 4.5514316911066173)
              << std::endl;

    std::cout << 1 << ' ' << 1 << ' ' << heherpq.Omega(1, 1, 31.860021837746324)
              << std::endl;
    std::cout << 2 << ' ' << 2 << ' ' << heherpq.Omega(2, 2, 31.860021837746324)
              << std::endl;
    std::cout << 1 << ' ' << 6 << ' ' << heherpq.Omega(1, 6, 31.860021837746324)
              << std::endl;
    std::cout << 4 << ' ' << 4 << ' ' << heherpq.Omega(4, 4, 31.860021837746324)
              << std::endl;
    return 0;

    std::cout << "===============0.54823988e-2" << std::endl;
    for (size_t l = 1; l != 7; ++l) {
      double Q = heherpq.Q(l, -1.0, 0.54823988e-2);
      Q *= M_PI * (1 - (1 + pow(-1, l)) / (2 * (1 + l)));
      std::cout << l << ' ' << Q << std::endl;
    }

    std::cout << "===============0.71356795e-1" << std::endl;
    for (size_t l = 1; l != 7; ++l) {
      double Q = heherpq.Q(l, -1.0, 0.71356795e-1);
      Q *= M_PI * (1 - (1 + pow(-1, l)) / (2 * (1 + l)));
      std::cout << l << ' ' << Q << std::endl;
    }

    std::cout << "===============0.43800684" << std::endl;
    for (size_t l = 1; l != 7; ++l) {
      double Q = heherpq.Q(l, -1.0, 0.43800684);
      Q *= M_PI * (1 - (1 + pow(-1, l)) / (2 * (1 + l)));
      std::cout << l << ' ' << Q << std::endl;
    }

    std::cout << "===============0.51775559D+02" << std::endl;
    for (size_t l = 1; l != 7; ++l) {
      double Q = heherpq.Q(l, -1.0, 0.51775559e+02);
      Q *= M_PI * (1 - (1 + pow(-1, l)) / (2 * (1 + l)));
      std::cout << l << ' ' << Q << std::endl;
    }

    std::cout << "==============================" << std::endl;
    std::cout << "===============0.12105661D+02" << std::endl;
    for (size_t l = 1; l != 7; ++l) {
      double Q = hexerpq.Q(l, -1.0, 0.12105661e+02);
      Q *= M_PI * (1 - (1 + pow(-1, l)) / (2 * (1 + l)));
      std::cout << l << ' ' << Q << std::endl;
    }
  }

  if (0) {
    /** The Lennerd-Jones Potential function.
     */
    class : public dlt::FuncDeriv1D {
    private:
      mutable double r_6_;
      mutable double old_r_ = -1.0;

      void prepare_(double r) const {
        r_6_ = std::pow(r, -6);
        // std::this_thread::sleep_for(std::chrono::milliseconds(1000));
        // fake time cost
        old_r_ = r;
        return;
      }

    public:
      double value(double r) const {
        if (r != old_r_) {
          prepare_(r);
        }
        return 4 * (r_6_ * r_6_ - r_6_);
      }
      double derivative(double r) const {
        if (r != old_r_) {
          prepare_(r);
        }
        return 4 * ((-12) * r_6_ * r_6_ / r + 6 * r_6_ / r);
      }
      bool provide_derivative() const { return true; };

    } lj;

    const double T = 10.0;
    std::cout << T << " <== T \n";
    std::cout << "Omega...\n";
    dlt::Pot1DFeatures pf(lj);
    dlt::ReducedPotentialQuadrature rpq(pf.reduced_potential());
    for (size_t l = 1; l != 7; ++l) {
      for (size_t s = l; s != 7; ++s) {
        double om = rpq.Omega(l, s, T);
        std::cout << l << ' ' << s << ' ' << om << std::endl;
      }
    }
    double om = rpq.Omega(1, 1, T);
    std::cout << 1 << ' ' << 1 << ' ' << om << std::endl;

    std::cout << "Q...\n";
    for (size_t l = 1; l != 7; ++l) {
      double Q = rpq.Q(l, -1.0, T);
      Q *= M_PI * (1 - (1 + pow(-1, l)) / (2 * (1 + l)));
      std::cout << l << ' ' << Q << std::endl;
    }
  }

  if (0) {
    // HFDHE2 Potential for helium
    class : public dlt::FuncDeriv1D {
    private:
      mutable double old_r_ = -1.0;
      mutable double v;
      mutable double dv;
      void prepare_(double r) const {
        const double a = 0.544850e6;
        const double alpha = 13.353384;
        const double c6 = 1.3732412;
        const double c8 = 0.4253785;
        const double c10 = 0.178100;
        const double d = 1.241314;
        double x = r; // / 2.9673;
        double vr = (x <= 7.6 ? a * exp(-alpha * x) : 0.0);
        double dvr = -alpha * vr;
        double x2 = x * x;
        double x4 = x2 * x2;
        double x6 = x4 * x2;
        double x8 = x4 * x4;
        double x10 = x8 * x2;
        double va = -(c6 / x6 + c8 / x8 + c10 / x10);
        double dva =
            6 * c6 / (x6 * x) + 8 * c8 / (x8 * x) + 10 * c10 / (x10 * x);
        if (x < d) {
          double fx = exp(-pow(d / x - 1, 2));
          double dfx = 2 * fx * (d * d / (x2 * x) - d / x2);
          v = vr + va * fx;
          dv = (dvr + dva * fx + va * dfx); // / 2.9673;
          old_r_ = r;
          return;
        } else {
          v = vr + va;
          dv = (dvr + dva); // / 2.9673;
          old_r_ = r;
          return;
        }
      }

    public:
      double value(double r) const {
        if (r != old_r_) {
          prepare_(r);
        }
        return v;
      }
      double derivative(double r) const {
        if (r != old_r_) {
          prepare_(r);
        }
        return dv;
      }
      bool provide_derivative() const { return true; };
    } hfdhe2;

    /*
    std::ofstream f;
    f.open("pec.txt", std::ios::out);
    for (double r = 0.5; r <= 10.0; r += 0.05) {
      f << r << ' ' << hfdhe2(r) << ' ' << hfdhe2.derivative(r) << ' '
        << hfdhe2.dlt::FuncDeriv1D::derivative(r) << '\n';
    }
    f.flush();
    f.close();
    return 0;
    */

    dlt::Pot1DFeatures pf(hfdhe2);
    dlt::ReducedPotentialQuadrature rpq(pf.reduced_potential());
    std::cout << rpq.chi(1.0, 1.0) << std::endl;
    for (size_t l = 1; l != 7; ++l) {
      for (size_t s = l; s != 7; ++s) {
        double om = rpq.Omega(l, s, 1.0);
        std::cout << l << ' ' << s << ' ' << om << std::endl;
      }
    }
  }

  return 0;
}
