// #include "atompair.hh"
#include "cgquad.hh"
#include "mathtools.hh"
#include "toms424.hh"
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <nlopt.hpp>
#include <tuple>
#include <vector>
/** The Lennerd-Jones Potential function.
 */
class : public dlt::FuncDeriv1D { // {{{
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

//}}}

namespace dlt {

/** This class find the features of a potential function.
 */
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
    // std::cout << "epsilon\tr_min" << std::endl;
    // std::cout << epsilon_ << "\t" << r_min_ << std::endl;

    sigma_ = find_local_root(*ppot_, 0.0, 0.5 * r_min_, r_min_);
    // TODO: I believe this range is safe...
    // std::cout << "sigma_" << std::endl;
    // std::cout << sigma_ << std::endl;
    return;
  }
  /** The collision radius.
   */
  const double &sigma() const { return sigma_; }

  /** The well depth.
   */
  const double &epsilon() const { return epsilon_; }

  /** The equlibrium nuclear seperation.
   */
  const double &r_min() const { return r_min_; }

  /** The potential.
   */
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

  /** This class inheriated from the Chebyshev-Gauss Quarduture class, compute
   * the integration required by the computation of chi.
   */
  class ChiCG : public CGIntegrator { // {{{
  private:
    IntegralRange *ir_;
    size_t ordersize_;
    // the following is for recording running status
    const double r_m_;
    const double E_;
    const double b_;

  public:
    ChiCG(IntegralRange *ir, double r_m, double E)
        : CGIntegrator(-1, 1, true), ir_(ir), r_m_(r_m), E_(E),
          b_(ir_->r2b(r_m, E)) {
      ordersize_ = 0;
      integrands_.clear();
    }

    /** Compute the integrands with computed values cached.
     */
    void calculate_integrands(size_t ordersize) override {
      // See if we need an update based on the "flag".
      if (ordersize_ >= ordersize) {
        return;
      }
      // We only need the positive half of the integration.
      CubicIter ci(ordersize_, ordersize, true, true);
      // Reserve the memory... which help std::vector work faster
      integrands_.reserve(ci.size_from_0());
      for (auto &&i : ci) {
        double y = map_pm1(CGIntegratorBackend::instance()->coss(i));
        double r = r_m_ / y;
        double v = (*ir_->ppot_)(r);
        double F = (1.0 - v / E_ - b_ * b_ / r / r);
        if (y <= 1.0e-8) { // r --> inf
          F = 1.0;
        }
        if (F <= 1.0e-12) {
          // should never run here if Chebyshev-Gauss Quarduture is used
          // because y does not equal 1 in CG Quad
          std::cerr << " ~~F~~ " << F << " ~~r~~ " << r << std::endl;
          throw std::runtime_error("178 F<0");
          F = 1.0e-12;
        }
        double res = 1.0 / sqrt(F) / r_m_;
        integrands_.push_back(res * sqrt(1.0 - y * y));
      }
      // Update the "flag".
      ordersize_ = ordersize;
      return;
    }
  }; // }}}

public:
  double chi(double E, double r_m) {
    // double E = ppot_->value(r_E);
    double b = r2b(r_m, E);

    ChiCG chicg(this, r_m, E);
    double quadrature, err;
    std::tie(quadrature, err) = chicg.integrate(1.0e-4, 15);
    return M_PI - b * quadrature;
  }

  IntegralRange(Pot1DFeatures &pf) : pf_(&pf), ppot_(&(pf.pot())) {

    y_.reset(new Y(*ppot_));

    std::tie(r_C_, E_C_) = find_local_maximum(*y_, 1.0, 2 * pf_->r_min());
    // TODO: fit reduced potential
    y_root_.reset(new LocalRoot(*y_, r_C_, 15.0, 21, 1.0e-12));
  }

  const double &r_C() const { return r_C_; }
  const double &E_C() const { return E_C_; }

  std::tuple<double, double> r_range(double E) const {
    double r_O, r_Op, b_O;
    if (E >= E_C_) {
      r_O = r_C_;
      r_Op = r_C_;
      b_O = r2b(r_O, E);
    } else {
      r_O = (*y_root_)(E);
      b_O = r2b(r_O, E);
      PhiEff phieff(*ppot_, b_O, E);
      r_Op = find_local_root(phieff, E, 1.0, r_C_, 21, 1.0e-12);
      // std::cout << "Zero check: " << phieff(r_Op) - E << std::endl;
      if (fabs(phieff(r_Op) - E) >= 1.0e-6) {
        double _;
        std::tie(r_Op, _) = find_local_maximum(phieff, 1.0, r_C_);
        r_Op = find_local_root(phieff, E, 0.0, r_Op, 21, 1.0e-12);
      }
    }
    // std::cout << "r_O, r_Op, r_C   " << r_O << ' ' << r_Op << ' ' << r_C()
    // << std::endl;

    // std::cout << "F(r_O) = "
    // << 1.0 - ppot_->value(r_O) / E - b_O * b_O / r_O / r_O << '\n';
    // std::cout << "F(r_Op) = "
    // << 1.0 - ppot_->value(r_Op) / E - b_O * b_O / r_Op / r_Op << '\n';
    return std::make_tuple(r_O, r_Op);
  }

  double r2b(double r, double E) const {
    double v = ppot_->value(r);
    if (r > 1.0e8) {
      v = 0.0;
    }
    double b = r * sqrt((E - v) / E);
    return b;
  }

  // Q stuff
  /** This class inheriated from the Chebyshev-Gauss Quarduture class, compute
   * the integration required by the computation of Q.
   */
  class QCG1 : public CGIntegrator { //{{{
  private:
    IntegralRange *ir_;
    size_t ordersize_;
    size_t cache_ordersize_;
    size_t l_;
    std::vector<double> coschis_;
    std::vector<double> fct2_;
    // the following is for recording running status
    // const double r_E_;
    const double E_;

  public:
    QCG1(IntegralRange *ir, double r_E, double r_Op)
        : CGIntegrator(r_E, r_Op, false), ir_(ir), E_(ir->ppot_->value(r_E)) {
      ordersize_ = 0;
      cache_ordersize_ = 0;
      l_ = 0;
      integrands_.clear();
      coschis_.clear();
      fct2_.clear();
    }
    /** Set the parameters, clear the inner storage if it is needed.
     */
    void set_l(size_t l) {
      if (l_ != l) {
        ordersize_ = 0;
        integrands_.clear();
        l_ = l;
      }
      return;
    }

    /** Compute the integrands with computed values cached.
     */
    void calculate_integrands(size_t ordersize) override {
      // See if we need an update based on the "flag".
      if (cache_ordersize_ < ordersize) {
        CubicIter ci(cache_ordersize_, ordersize, false, true);
        size_t num = ci.size_from_0();
        integrands_.reserve(num);
        coschis_.reserve(num);
        fct2_.reserve(num);
        for (auto &&i : ci) {
          double y = CGIntegratorBackend::instance()->coss(i);
          double r_m = map_pm1(y);
          double v, dv;
          v = ir_->ppot_->value(r_m);
          dv = ir_->ppot_->derivative(r_m);
          double coschi = cos(ir_->chi(E_, r_m));
          // Here I put everything else in fct2, include the weight
          double fct2 = (2.0 * (E_ - v) - r_m * dv) * r_m * sqrt(1.0 - y * y);
          coschis_.push_back(coschi);
          fct2_.push_back(fct2);
        }
        cache_ordersize_ = ordersize;
      }
      // std::cout << coschis_.size() << '\n' << fct2_.size() << std::endl;
      if (ordersize_ < ordersize) {
        CubicIter ci(ordersize_, ordersize, false, false);
        size_t num = ci.size_from_0();
        integrands_.reserve(num);
        for (auto &&i : ci) {
          double res = (1.0 - pow(coschis_[i], l_)) * fct2_[i];
          integrands_.push_back(res);
        }
        ordersize_ = ordersize;
      }
      return;
    }
  }; // }}}

  class QCG2 : public CGIntegrator { //{{{
  private:
    IntegralRange *ir_;
    size_t ordersize_;
    size_t cache_ordersize_;
    size_t l_;
    std::vector<double> coschis_;
    std::vector<double> fct2_;
    // the following is for recording running status
    // double r_E_;
    const double E_;
    double r_O_;

  public:
    QCG2(IntegralRange *ir, double r_E, double r_O)
        : CGIntegrator(-1, 1, true), ir_(ir),
          // r_E_(r_E),
          E_(ir->ppot_->value(r_E)), r_O_(r_O) {
      ordersize_ = 0;
      cache_ordersize_ = 0;
      l_ = 0;
      integrands_.clear();
      coschis_.clear();
      fct2_.clear();
    }
    /** Set the parameters, clear the inner storage if it is needed.
     */
    void set_l(size_t l) {
      if (l_ != l) {
        ordersize_ = 0;
        integrands_.clear();
        l_ = l;
      }
      return;
    }

    /** Compute the integrands with computed values cached.
     */
    void calculate_integrands(size_t ordersize) override {
      // See if we need an update based on the "flag".
      if (cache_ordersize_ < ordersize) {
        CubicIter ci(cache_ordersize_, ordersize, true, true);
        size_t num = ci.size_from_0();
        integrands_.reserve(num);
        coschis_.reserve(num);
        fct2_.reserve(num);
        for (auto &&i : ci) {
          double y = CGIntegratorBackend::instance()->coss(i);
          if (y < 1.0e-8) {
            y = 1.0e-8; // prevent div by 0
          }
          double r_m = r_O_ / y;
          double v, dv;
          v = ir_->ppot_->value(r_m);
          dv = ir_->ppot_->derivative(r_m);
          double coschi = cos(ir_->chi(E_, r_m));
          // Here I put everything else in fct2, include the weight
          double fct2 =
              (2.0 * (E_ - v) - r_m * dv) * r_m * r_m / y * sqrt(1.0 - y * y);
          coschis_.push_back(coschi);
          fct2_.push_back(fct2);
        }
        cache_ordersize_ = ordersize;
      }
      if (ordersize_ < ordersize) {
        CubicIter ci(ordersize_, ordersize, true, true);
        size_t num = ci.size_from_0();
        integrands_.reserve(num);
        for (auto &&i : ci) {
          double res = (1.0 - pow(coschis_[i], l_)) * fct2_[i];
          integrands_.push_back(res);
        }
        ordersize_ = ordersize;
      }
      return;
    }
  }; // }}}

  /** Compute Q.
   *
   * It is faster to keep the r_E unchanged and scan the l.
   */
  double Q(size_t l, double r_E) {
    static double old_r_E = 0.0;
    static double E, r_O, r_Op;
    static std::unique_ptr<QCG1> qcg1;
    static std::unique_ptr<QCG2> qcg2;

    if (old_r_E != r_E) {
      old_r_E = r_E;
      E = ppot_->value(r_E); // This is common
      std::tie(r_O, r_Op) = r_range(E);
      qcg1.reset(new QCG1(this, r_E, r_Op));
      qcg2.reset(new QCG2(this, r_E, r_O));
    }

    double coeff = 1.0 / (1.0 - (1.0 + pow(-1, l)) / 2.0 / (1.0 + l)) / E;
    double esterr;
    double quadrature1;
    qcg1->set_l(l);
    std::tie(quadrature1, esterr) = qcg1->integrate(1.e-4, 15);
    double quadrature2;
    qcg2->set_l(l);
    std::tie(quadrature2, esterr) = qcg2->integrate(1.e-4, 15);
    quadrature2 /= 2;
    return coeff * (quadrature1 + quadrature2);
  }

  // Omega stuff
  /** This class inheriated from the Chebyshev-Gauss Quarduture class, compute
   * the integration required by the computation of Omega.
   */
  class OmegaCG : public CGIntegrator { //{{{
  private:
    IntegralRange *ir_;
    size_t ordersize_;
    size_t v_ordersize_;
    size_t Q_ordersize_;
    std::vector<double> Qs_;
    std::vector<double> vs_;
    std::vector<double> dvs_;
    // the following is for recording running status
    // const double r_E_;
    size_t l_;
    size_t s_;
    double T_;

  public:
    OmegaCG(IntegralRange *ir) : CGIntegrator(-1.0, 1.0, true), ir_(ir) {
      ordersize_ = 0;
      v_ordersize_ = 0;
      Q_ordersize_ = 0;
      l_ = 0;
      s_ = 0;
      T_ = 0.0;
      integrands_.clear();
      Qs_.clear();
      vs_.clear();
      dvs_.clear();
    }
    /** Set the parameters, clear the inner storage if it is needed.
     */
    void set_l_s_T(size_t l, size_t s, double T) {
      if (s_ != s || l_ != l || T_ != T) {
        ordersize_ = 0;
        integrands_.clear();
      }
      if (l_ != l) {
        Q_ordersize_ = 0;
        Qs_.clear();
      }
      s_ = s;
      l_ = l;
      T_ = T;
      return;
    }

    /** Compute the integrands with computed values cached.
     */
    void calculate_integrands(size_t ordersize) override {
      // See if we need an update based on the "flag".
      if (Q_ordersize_ < ordersize) {
        CubicIter ci(Q_ordersize_, ordersize, true, true);
        size_t num = ci.size_from_0();
        Qs_.reserve(num);
        for (auto &&i : ci) {
          double y = CGIntegratorBackend::instance()->coss(i);
          if (y < 1.0e-8) {
            y = 1.0e-8; // prevent div by 0
          }
          double r_E = map_pm1(y);
          Qs_.push_back(ir_->Q(l_, r_E));
        }
        Q_ordersize_ = ordersize;
      }
      if (v_ordersize_ < ordersize) {
        CubicIter ci(v_ordersize_, ordersize, true, true);
        size_t num = ci.size_from_0();
        Qs_.reserve(num);
        std::ofstream fs;
        for (auto &&i : ci) {
          double y = CGIntegratorBackend::instance()->coss(i);
          if (y < 1.0e-8) {
            y = 1.0e-8; // prevent div by 0
          }
          double r_E = map_pm1(y);

          vs_.push_back(ir_->ppot_->value(r_E));
          // weight is included in dvs
          dvs_.push_back(ir_->ppot_->derivative(r_E) * sqrt(1.0 - y * y));
        }
        v_ordersize_ = ordersize;
      }
      // std::cout << "Line " << __LINE__ << std::endl;
      if (ordersize_ < ordersize) {
        CubicIter ci(ordersize_, ordersize, true, true);
        size_t num = ci.size_from_0();
        vs_.reserve(num);
        dvs_.reserve(num);
        // std::cout << "Line " << __LINE__ << std::endl;

        for (auto &&i : ci) {
          double x = vs_[i] / T_;
          double res = exp(-x) * pow(x, s_ + 1) * Qs_[i] * dvs_[i];
          integrands_.push_back(res);
        }
        ordersize_ = ordersize;
      }
      // std::cout << "Line " << __LINE__ << std::endl;
      return;
    }
  }; // }}}

  double Omega(size_t l, size_t s, double T) {
    static OmegaCG omegacg(this);
    double coeff = -1.0 / T / std::tgamma(s + 2);
    double esterr;
    double quadrature1;
    omegacg.set_l_s_T(l, s, T);
    std::tie(quadrature1, esterr) = omegacg.integrate(1.e-4, 7);
    // std::cout << "Line " << __LINE__ << std::endl;
    quadrature1 /= 2;
    return coeff * quadrature1;
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

  std::cout << "ir.Q(1, 0.999)" << std::endl;
  std::cout << ir.Q(1, 0.999) << std::endl;

  std::cout << "ir.Q(2, 0.999)" << std::endl;
  std::cout << ir.Q(2, 0.999) << std::endl;

  std::cout << "ir.Q(3, 0.999)" << std::endl;
  std::cout << ir.Q(3, 0.999) << std::endl;

  std::cout << "ir.Q(3, 0.99)" << std::endl;
  std::cout << ir.Q(3, 0.99) << std::endl;

  std::cout << "ir.Q(4, 0.999)" << std::endl;
  std::cout << ir.Q(3, 0.999) << std::endl;

  std::cout << "ir.Q(3, 0.99)" << std::endl;
  std::cout << ir.Q(3, 0.99) << std::endl;

  std::cout << "ir.Q(3, 0.999)" << std::endl;
  std::cout << ir.Q(3, 0.999) << std::endl;

  std::cout << "ir.Omega(1, 1, 30)" << std::endl;
  std::cout << ir.Omega(1, 1, 30) << std::endl;

  std::cout << "ir.Omega(1, 2, 30)" << std::endl;
  std::cout << ir.Omega(1, 2, 30) << std::endl;

  /*
  std::cout << "ir.Omega(2, 2, 30)" << std::endl;
  std::cout << ir.Omega(2, 2, 30) << std::endl;
  */

  std::cout << "ir.Omega(1, 1, 30)" << std::endl;
  std::cout << ir.Omega(1, 1, 30) << std::endl;

  return 0;
}
