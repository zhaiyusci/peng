// #include "atompair.hh"
#include "mathtools.hh"
#include "toms424.hh"
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <nlopt.hpp>
#include <tuple>
#include <vector>
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

namespace dlt {

class Pot1DFeatures {
private:
  double sigma_;
  double epsilon_;
  double r_min_;
  FuncDeriv1D *const ppot_;
  void init_() {}

public:
  Pot1DFeatures(FuncDeriv1D &pot) : ppot_(&pot) {
    std::tie(r_min_, epsilon_) = find_local_minimum(*ppot_, 1.0, 5.0);
    // TODO: I believe this range is safe...
    epsilon_ *= -1;
    std::cout << "epsilon\tr_min" << std::endl;
    std::cout << epsilon_ << "\t" << r_min_ << std::endl;

    sigma_ = find_local_root(*ppot_, 0.0, 0.5 * r_min_, r_min_);
    // TODO: I believe this range is safe...
    std::cout << "sigma_" << std::endl;
    std::cout << sigma_ << std::endl;
    return;
  }
  const double &sigma() const { return sigma_; }
  const double &epsilon() const { return epsilon_; }
  const double &r_min() const { return r_min_; }
  FuncDeriv1D &pot() const { return *ppot_; }
};

class PhiEff : public FuncDeriv1D {
private:
  FuncDeriv1D *ppot_;
  const double b_;
  const double E_;

public:
  PhiEff(FuncDeriv1D &pot, double b, double E) : ppot_(&pot), b_(b), E_(E) {}
  double value(double r) const {
    double res;
    res = ppot_->value(r);
    res += b_ * b_ * E_ / r / r;
    return res;
  }
  double derivative(double r) const {
    double res;
    res = ppot_->derivative(r);
    res += -2 * b_ * b_ * E_ / r / r / r;
    return res;
  }
  bool provide_derivative() const { return ppot_->provide_derivative(); }
};

class Y : public FuncDeriv1D {
private:
  FuncDeriv1D *const ppot_;

public:
  Y(FuncDeriv1D &pot) : ppot_(&pot) {}
  double value(double r) const {
    double res;
    res = ppot_->value(r);
    res += 0.5 * r * ppot_->derivative(r);
    return res;
  }
};

class IntegralRange {
private:
  const Pot1DFeatures *const pf_;
  FuncDeriv1D *const ppot_;
  double r_C_;
  double E_C_;
  std::unique_ptr<LocalRoot> y_root_;
  std::unique_ptr<Y> y;

public:
  IntegralRange(Pot1DFeatures &pf) : pf_(&pf), ppot_(&(pf.pot())) {

    y.reset(new Y(*ppot_));

    std::tie(r_C_, E_C_) = find_local_maximum(*y, 1.0, 2 * pf_->r_min());
    // TODO: fit reduced potential
    std::cout << "r_C\tE_C" << std::endl;
    std::cout << r_C_ << "\t" << E_C_ << std::endl;
    y_root_.reset(new LocalRoot(*y, r_C_, 10.0, 21));
  }
  const double &r_C() const { return r_C_; }
  const double &E_C() const { return E_C_; }

  std::tuple<double, double> r_range(double E) const {
    if (E >= E_C()) {
      return std::make_tuple(r_C(), r_C());
    } else {
      double r_O, r_Op;
      r_O = (*y_root_)(E);
      double b_O = r2b(r_O, E);
      PhiEff phieff(*ppot_, b_O, E);
      r_Op = find_local_root(phieff, E, 1.0, r_C_);
      return std::make_tuple(r_O, r_Op);
    }
  }

  double chi(double E, double r_m) {
    // double E = ppot_->value(r_E);
    double b = r2b(r_m, E);
    // std::cout << "b    " << b << std::endl;
    auto integrated = [&](double y) {
      double r = r_m / y;
      double v = (*ppot_)(r);
      if (y <= 1.0e-8) {
        v = 0.0;
      }
      if (y >= 1 - 1.0e-8) {
        return 0.0;
      }
      double F = (1.0 - v / E - b * b / r / r);
      if (F < 0.0) {
        // std::cout << " ~~F~~ " << F << " ~~r~~ " << r << std::endl;
      }
      double res = 1.0 / sqrt(F) / r_m;

      // std::cout << "integrated    " << res << " y    " << y << std::endl;
      return res;
    };
    double quadrature, esterr;
    size_t used;
    std::vector<double> _;
    std::tie(quadrature, esterr, used, _) =
        ccquad(integrated, 0.0, 1.0, 1.0e-8, 10000); // TODO
    // for (auto &&x : _)
    // std::cout << x << std::endl;

    return M_PI - 2 * b * quadrature;
  }
  double r2b(double r, double E) const {
    double v = ppot_->value(r);
    double b = r * sqrt((E - v) / E);
    // std::cout << "b = " << b << std::endl;
    return b;
  }

  double Q(int l, double r_E) {

    double E = ppot_->value(r_E);
    std::cout << "E" << std::endl;
    std::cout << E << std::endl;
    double r_O, r_Op;
    std::tie(r_O, r_Op) = r_range(E);
    std::cout << "r_O,   r_Op" << std::endl;
    std::cout << r_O << "    " << r_Op << std::endl;
    double coeff = 1.0 / (1.0 - (1.0 + pow(-1, l)) / 2.0 / (1.0 + l)) / E;
    auto integrated1 = [&](double r_m) {
      double v, dv;
      v = ppot_->value(r_m);
      dv = ppot_->derivative(r_m);
      double chival = chi(E, r_m);
      std::cout << "YO " << r_m << "    " << chival << std::endl;
      double res =
          (1.0 - pow(cos(chival), l)) * (2.0 * (E - v) - r_m * dv) * r_m;
      return res;
    };
    double quadrature1, quadrature2;
    double esterr;
    int used;
    std::vector<double> _;
    std::cout << "r_E << "
                 " << r_Op"
              << std::endl;
    std::cout << r_E << "    " << r_Op << std::endl;
    std::tie(quadrature1, esterr, used, _) =
        ccquad(integrated1, r_E, r_Op, 1.0e-3, 100000); // TODO
    auto integrated2 = [&](double y) {
      double x;
      x = r_O / y;
      if (y <= 1.0e-8) {
        return 0.0;
      }
      // if (y >= 1 - 1.0e-8) {
      // return 0.0;
      // }
      return x / y * integrated1(x);
    };
    // std::tie(quadrature2, esterr, used, _) =
    // ccquad(integrated2, 0.0, 1.0, 1.0e-3, 100000); // TODO

    for (double r = 90; r <= 100; r += 1.0) {
      std::cout << "dbg==  " << r << "           " << integrated1(r)
                << std::endl;
    }

    std::tie(quadrature2, esterr, used, _) =
        ccquad(integrated1, r_O, 1000.0, 1.0e-3, 100000); // TODO
    std::cout << "quadrature1 << "
                 " << quadrature2"
              << std::endl;
    std::cout << quadrature1 << "    " << quadrature2 << std::endl;
    return coeff * (quadrature1 + quadrature2);
  }
};
} // namespace dlt

int main() {
  std::cerr << __FILE__ << " : " << __LINE__ << " : " << __FUNCTION__
            << std::endl;
  std::cout << std::setprecision(8);

  dlt::Pot1DFeatures pf(lj);
  std::cout << " " << pf.r_min() << " " << pf.sigma() << " " << pf.epsilon()
            << " " << std::endl;
  dlt::IntegralRange ir(pf);
  for (double E = 0.0; E <= 0.8; E += 0.1) {
    std::cout << "r_O ==> " << E << " " << std::get<0>(ir.r_range(E)) << " "
              << std::get<1>(ir.r_range(E)) << std::endl;
  }

  std::cout << "ir.chi(lj(0.9), 0.9)" << std::endl;
  std::cout << ir.chi(lj(0.9), 0.9) << std::endl;
  std::cout << "ir.chi(lj(0.9), 99.9)" << std::endl;
  std::cout << ir.chi(lj(0.9), 99.9) << std::endl;
  std::cout << "E = " << lj(0.999) << std::endl;
  std::cout << "ir.Q(1, 0.999)" << std::endl;
  std::cout << ir.Q(1, 0.999) << std::endl;
  std::cout << "ir.Q(2, 0.999)" << std::endl;
  std::cout << ir.Q(2, 0.999) << std::endl;

  // double r_O, r_Op;
  // std::tie(r_O, r_Op) = ir.r_range(0.7);
  // std::cout << r_O << "  " << r_Op << std::endl;
  // std::cout << "ir.chi(0.7, 1.2)" << std::endl;
  // std::cout << ir.chi(0.7, 1.2) << std::endl;
  // for (double rm = 0.9; rm <= 1.3; rm += 0.10) {
  // std::cout << ir.chi(lj(0.9), rm) << std::endl;
  // }

  return 0;
}
