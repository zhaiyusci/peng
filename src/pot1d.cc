#include "atompair.hh"
#include <functional>
#include <iomanip>
#include <iostream>
#include <nlopt.hpp>
#include <tuple>
#include <vector>

auto lj(double r) {
  double r_6(pow(r, -6));
  double v(4 * (r_6 * r_6 - r_6));
  double dv(4 * ((-12) * r_6 * r_6 / r - (-6) * r_6 / r));
  // double d2v(4 * ((12 * 13) * r_6 * r_6 / r / r - (6 * 7) * r_6 / r / r));

  return std::make_tuple(v, dv);
}

double wrapper(const std::vector<double> &x, std::vector<double> &grad,
               void *f_data) {
  Pot1d pot = *reinterpret_cast<Pot1d *>(f_data);
  double v, dv;
  std::tie(v, dv) = pot(x[0]);
  return v;
}

template <typename T>
std::tuple<T, size_t> min_and_index(const std::valarray<T> &fs) {
  size_t min_i(0);
  T min_f(fs[0]);
  for (size_t i = 1; i != fs.size(); ++i) {
    if (fs[i] <= min_f) {
      min_i = i;
      min_f = fs[i];
    }
  }
  return std::make_tuple(min_f, min_i);
}

auto find_local_minimum(const Pot1d &pot, double lower, double upper,
                        size_t nsamp = 11, double xtol_err = 1.0e-8) {
  std::valarray<double> xs(nsamp); // TODO
  std::valarray<double> fs(nsamp); // TODO

  std::vector<double> x{0.5 * (lower + upper)};
  // Compute the pointwise values to get the optimization range
  if (nsamp >= 3) {
    {
      double width = upper - lower;
      double step = width / (nsamp - 1);
      for (size_t i = 0; i != nsamp; ++i) {
        xs[i] = lower + i * step;
        double _;
        std::tie(fs[i], _) = pot(xs[i]);
        // std::cerr << xs[i] << "  " << fs[i] << std::endl;
      }
    }
    size_t min_i;
    double min_f;
    std::tie(min_f, min_i) = min_and_index(fs);
    lower = xs[min_i] - 1; // This can be dangerous because min_i is unsigned...
                           // relax, see below
    upper = xs[min_i] + 1;
    std::cerr << "lower << "
                 " << upper"
              << std::endl;
    std::cerr << lower << " " << upper << std::endl;
    if (min_i == 0 || min_i == nsamp - 1) {
      std::cerr << "WARNING:" << std::endl
                << __FILE__ << " " << __LINE__
                << ": Searching range may not be good, with index = " << min_i
                << " and val = " << min_f << std::endl;
      if (min_i == 0)
        lower = 0;
      if (min_i == nsamp - 1)
        upper = nsamp - 1;
    }
    x[0] = xs[min_i];
    std::cerr << "===================" << min_i << std::endl;
  }
  nlopt::opt opt(nlopt::LN_COBYLA, 1);
  opt.set_xtol_rel(xtol_err);
  opt.set_lower_bounds(lower);
  opt.set_upper_bounds(upper);
  Pot1d p(pot); // NLOpt not support a const data, so...
  opt.set_min_objective(wrapper, &p);
  double f;
  opt.optimize(x, f);
  return std::make_tuple(x[0], f);
}

class Pot1dFeatures {
  private:
    double sigma_;
    double epsilon_;
    double r_min_;
    Pot1d *const pot_;

  public:
    Pot1dFeatures(Pot1d &pot) : pot_(&pot) {
      std::tie(r_min_, epsilon_) = find_local_minimum(
          *pot_, 1.0, 5.0); // TODO: I believe this range is safe...
      epsilon_ *= -1;
      std::cout << "epsilon\tr_min" << std::endl;
      std::cout << epsilon_ << "\t" << r_min_ << std::endl;

      double _;
      Pot1d potabs = [=](double x) {
        double v, dv;
        std::tie(v, dv) = (*(this->pot_))(x);
        // std::cerr << x << "  " << v << std::endl;
        return std::make_tuple(fabs(v), (v > 0 ? 1 : -1) * dv);
      };
      std::tie(sigma_, _) =
          find_local_minimum(potabs, 0.5 * r_min_,
                             r_min_); // TODO: I believe this range is safe...
      std::cout << "sigma_" << std::endl;
      std::cout << sigma_ << std::endl;

    }
    const double &sigma() const { return sigma_; }
    const double &epsilon() const { return epsilon_; }
    const double &r_min() const { return r_min_; }
    Pot1d *pot() const { return pot_; }
};

class IntegralRange {
  private:
    const Pot1dFeatures *const pf_;
    Pot1d *const pot_;
    double r_C_;
    double E_C_;
    /**
     * This function is for finding r_C and E_C
     */
    static double wrapper_y_(const std::vector<double> &x,
                             std::vector<double> &grad, void *f_data) {
      Pot1d pot = *reinterpret_cast<Pot1d *>(f_data);
      double r(x[0]);
      double v, dv;
      std::tie(v, dv) = pot(r);
      // std::cout << " r = " << r << "   v = " << v << "    dv = " << dv
      // << std::endl;

      return -(2 * v + r * dv);
    }

  public:
    IntegralRange(Pot1dFeatures &pf) : pf_(&pf), pot_(pf.pot()) {

      nlopt::opt opt(nlopt::LN_COBYLA, 1);
      opt.set_xtol_rel(1.0e-8);  // TODO
      std::vector<double> r0(1); // TODO: Let user to give the init guess
      {
        opt.set_lower_bounds(pf_->r_min());
        opt.set_upper_bounds(20.0);
        r0[0] = 2 * pf_->r_min();
        opt.set_min_objective(wrapper_y_, this->pot_);
        // r0[0] = 10.0; // TODO
        opt.optimize(r0, E_C_);
        r_C_ = r0[0];
        E_C_ *= -1;
        std::cout << "r_C\tE_C" << std::endl;
        std::cout << r_C_ << "\t" << E_C_ << std::endl;
      }
    }
    const double &r_C() const { return r_C_; }
    const double &E_C() const { return E_C_; }
};
int main() {
  std::cout << std::setprecision(8);
  nlopt::opt opt(nlopt::LN_COBYLA, 1);
  Pot1d pot(lj);
  Pot1dFeatures pf(pot);
  std::cout << " " << pf.r_min() << " " << pf.sigma() << " " << pf.epsilon()
            << " " << std::endl;
  return 0;
}
