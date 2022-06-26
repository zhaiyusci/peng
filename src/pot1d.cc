// #include "atompair.hh"
#include "cgquad.hh"
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
  std::unique_ptr<Y> y_;

  class ChiCG : public CGIntegrator {
  private:
    IntegralRange *ir_;
    size_t ordersize_;

    // the following is for recording running status
    double r_m_;
    double E_;
    double b_;

  public:
    ChiCG(IntegralRange *ir)
        : CGIntegrator(-1, 1, true), ir_(ir), ordersize_(1) {
      ordersize_ = 0;
      integrands_.clear();
      integrands_.resize(0);
      // calculate_integrands(1);
    }
    void set_r_m_E_b(double r_m, double E, double b) {
      if (r_m_ != r_m || E_ != E) {
        r_m_ = r_m;
        E_ = E;
        b_ = b;
        ordersize_ = 0;
        integrands_.clear();
        integrands_.resize(0);
        // calculate_integrands(1);
        // std::cerr << " integrands_.size -- " << integrands_.size() <<
        // std::endl;
      }
      return;
    }
    void calculate_integrands(size_t ordersize) override {
      if (ordersize_ >= ordersize) {
        return;
      }
      CubicIter ci(ordersize_, ordersize, true, true);
      size_t num = (pow(3, ordersize - 1) + 1) / 2;
      integrands_.reserve(num);
      for (auto &&i : ci) {
        double y = map_pm1(CGIntegratorBackend::instance()->coss(i));
        double r = r_m_ / y;
        double v = (*ir_->ppot_)(r);
        double F = (1.0 - v / E_ - b_ * b_ / r / r);
        if (y <= 1.0e-8) { // r --> inf
          F = 1.0;
        }
        if (F < 0.0) {
          // should never run here if Chebyshev-Gauss Quarduture is used
          // because y does not equal 1 in CG Quad
          F = 0.0;
          std::cout << " ~~F~~ " << F << " ~~r~~ " << r << std::endl;
        }
        double res = 1.0 / sqrt(F) / r_m_;

        std::cout << __LINE__ << ' ' << i << ' ' << (res * sqrt(1.0 - y * y))
                  << std::endl;
        integrands_.push_back(res * sqrt(1.0 - y * y));
      }
      ordersize_ = ordersize;
      return;
    }
  };

  ChiCG chicg;

public:
  double chi(double E, double r_m) {
    // double E = ppot_->value(r_E);
    double b = r2b(r_m, E);

    chicg.set_r_m_E_b(r_m, E, b);
    // std::cout << "b    " << b << " " << r_m << " " << E << std::endl;
    double quadrature, err;
    std::tie(quadrature, err) = chicg.integrate(1.0e-4, 4);
    return M_PI - b * quadrature;
  }

  IntegralRange(Pot1DFeatures &pf) : pf_(&pf), ppot_(&(pf.pot())), chicg(this) {

    y_.reset(new Y(*ppot_));

    std::tie(r_C_, E_C_) = find_local_maximum(*y_, 1.0, 2 * pf_->r_min());
    // TODO: fit reduced potential
    std::cout << "r_C\tE_C" << std::endl;
    std::cout << r_C_ << "\t" << E_C_ << std::endl;
    y_root_.reset(new LocalRoot(*y_, r_C_, 10.0, 21));
  }

  const double &r_C() const { return r_C_; }
  const double &E_C() const { return E_C_; }

  std::tuple<double, double> r_range(double E) const {
    double r_O, r_Op;
    if (E >= E_C_) {
      r_O = r_C_;
      r_Op = r_C_;
    } else {
      r_O = (*y_root_)(E);
      double b_O = r2b(r_O, E);
      PhiEff phieff(*ppot_, b_O, E);
      r_Op = find_local_root(phieff, E, 1.0, r_C_);
    }
    return std::make_tuple(r_O, r_Op);
  }

  double r2b(double r, double E) const {
    double v = ppot_->value(r);
    if (r > 1.0e8) {
      v = 0.0;
    }
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
      double res =
          (1.0 - pow(cos(chival), l)) * (2.0 * (E - v) - r_m * dv) * r_m;
      // std::cout << __LINE__ << "   " << res << std::endl;
      std::cout << __LINE__ << "   " << (1.0 - pow(cos(chival), l)) << "   "
                << (2.0 * (E - v) - r_m * dv) << "   " << r_m << std::endl;
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
        ccquad(integrated1, r_E, r_Op, 1.0e-3 / coeff, 100000); // TODO
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

    // for (double r = 90; r <= 100; r += 1.0) {
    // std::cout << "dbg==  " << r << "           " << integrated1(r)
    // << std::endl;
    // }

    std::cout << "quadrature2" << std::endl;
    std::tie(quadrature2, esterr, used, _) =
        // ccquad(integrated1, r_O, 2000.0, 1.0e-3 / coeff, 100000); // TODO
        ccquad(integrated1, 10.0, 20.0, 1.0e-3 / coeff, 100000); // TODO
    std::cout << " USED " << used << std::endl;
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
  std::cout << std::setprecision(18);

  dlt::Pot1DFeatures pf(lj);
  std::cout << " " << pf.r_min() << " " << pf.sigma() << " " << pf.epsilon()
            << " " << std::endl;
  dlt::IntegralRange ir(pf);
  /*
  for (double E = 0.0; E <= 1.0; E += 0.1) {
    std::cout << "r_O ==> " << E << " " << std::get<0>(ir.r_range(E)) << " "
              << std::get<1>(ir.r_range(E)) << std::endl;
  }
  */

  std::cout << "ir.E_C()" << ir.E_C() << std::endl;
  std::cout << "lj(0.9)" << lj(0.9) << std::endl;
  double E_O, E_Op;
  std::tie(E_O, E_Op) = ir.r_range(lj(0.9));
  std::cout << "E_O" << E_O << std::endl;

  std::cout << ir.chi(10.0, 0.9) << std::endl;

  std::cout << " ========= ========== =====" << std::endl;

  std::cout << "ir.chi(lj(0.9), 0.9)" << std::endl;
  std::cout << ir.chi(lj(0.9), 0.9) << std::endl;
  std::cout << "ir.chi(lj(0.9), 99.9)" << std::endl;
  std::cout << ir.chi(lj(0.9), 99.9) << std::endl;
  std::cout << "ir.chi(lj(0.9), 9999.9)" << std::endl;
  std::cout << ir.chi(lj(0.9), 9999.9) << std::endl;
  std::cout << "E = " << lj(0.999) << std::endl;
  // for (double rmm = 0.9; rmm <= 50.0; rmm += 5.0) {
  // std::cout << rmm << "   " << ir.chi(lj(0.9), rmm) << std::endl;
  // }

  /*
  std::cout << "ir.Q(1, 0.999)" << std::endl;
  std::cout << ir.Q(1, 0.999) << std::endl;
  std::cout << "ir.Q(2, 0.999)" << std::endl;
  std::cout << ir.Q(2, 0.999) << std::endl;
  */

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
