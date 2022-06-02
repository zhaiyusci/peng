#include "atompair.hh"
#include <functional>
#include <iomanip>
#include <iostream>
#include <nlopt.hpp>
#include <tuple>
#include <vector>


std::tuple<double, double, bool, size_t>
CachedFuncDeriv1D::cubic_spline_(double x) const {

  // find the right interval using bisec method...
  size_t li, ri, mi;
  size_t N_ = cache_.size();
  bool abinit(false);
  li = 0;
  ri = N_ - 1; // here we want the extropolate value be reasonable, as if
               // they stay in the first or last section
  if (N_ > 2 && x > std::get<0>(cache_[li]) && x <= std::get<0>(cache_[ri])) {
    mi = (li + ri) / 2;
    // mid point.  because they are all ints here, it is ok to do so
    while (ri - li != 1) {
      (x > std::get<0>(cache_[li]) && x <= std::get<0>(cache_[mi])) ? ri = mi
                                                                    : li = mi;
      mi = (li + ri) / 2;
      // std::cout << "liri  " << li << "   " << ri << std::endl;
    }
  } else {
    abinit = true;
    return std::make_tuple(std::nan("1"), std::nan("1"), abinit, N_);
  }
  abinit = !std::get<3>(cache_[li]); // if the interval is not qualified,
                                     // need an ab init computation
  // if (x - std::get<0>(cache_[li]) < ftol_ / fabs(std::get<2>(cache_[li]))
  // || std::get<0>(cache_[ri]) - x < ftol_ / fabs(std::get<2>(cache_[ri])))
  // { abinit = false;
  // }

  double w = std::get<0>(cache_[ri]) - std::get<0>(cache_[li]);

  double t = (x - std::get<0>(cache_[li])) / w;
  double t2 = t * t;
  double t3 = t2 * t;

  double h00 = 2 * t3 - 3 * t2 + 1;
  double h10 = t3 - 2 * t2 + t;
  double h01 = -2 * t3 + 3 * t2;
  double h11 = t3 - t2;

  double p00 = 6 * t2 - 6 * t;
  double p10 = 3 * t2 - 4 * t + 1;
  double p01 = -6 * t2 + 6 * t;
  double p11 = 3 * t2 - 2 * t;

  double e00 = 2 * t3;
  double e10 = t3;
  double e01 = -2 * t3;
  double e11 = t3;
  // cerr << "t   " << t<< endl;
  double v, dv;
  double esterr;
  v = h00 * std::get<1>(cache_[li]) + h10 * w * std::get<2>(cache_[li]) +
      h01 * std::get<1>(cache_[ri]) + h11 * w * std::get<2>(cache_[ri]);
  dv = p00 * std::get<1>(cache_[li]) + p10 * w * std::get<2>(cache_[li]) +
       p01 * std::get<1>(cache_[ri]) + p11 * w * std::get<2>(cache_[ri]);
  esterr = e00 * std::get<1>(cache_[li]) + e10 * w * std::get<2>(cache_[li]) +
           e01 * std::get<1>(cache_[ri]) + e11 * w * std::get<2>(cache_[ri]);
  abinit = !(fabs(esterr) < ftol_);
  return std::make_tuple(
      v, dv / (std::get<0>(cache_[ri]) - std::get<0>(cache_[li])), abinit, li);
}

std::tuple<double, double> CachedFuncDeriv1D::operator()(double x) const {
  double f, df;
  size_t li, ri;
  bool abinit, qq;
  std::tie(f, df, abinit, li) = cubic_spline_(x);
  qq = false;
  if (li != cache_.size()) {
    ri = li + 1;
    double lx = std::get<0>(cache_[li]);
    double rx = std::get<0>(cache_[ri]);
    double w = rx - lx;
    qq = (x >= lx + 0.3 * w && x <= rx - 0.3 * w);
  }
  // abinit = true;
  if (abinit) {
    double af, adf; // The ab init ones
    std::tie(af, adf) = (*func_)(x);
    if (fabs(f - af) < ftol_ && fabs(df - adf) < ftol_ && qq) {
      double x, f, df;
      bool _;
      std::tie(x, f, df, _) = cache_[li];
      cache_[li] = std::make_tuple(x, f, df, true);

      // std::cerr << "GOOD SPLINE!!!" << std::endl;
    } else {
      cache_.push_back(std::make_tuple(x, af, adf, false));
      std::sort(cache_.begin(), cache_.end());
    }
    f = af;
    df = adf;
  } else {
    // std::cerr << "GREAT!!!" << std::endl;
    std::fstream fs("test.txt", std::ios_base::out | std::fstream::app);
    fs << x << "   " << f << "  " << df << std::endl;
  }
  // cache_.push_back(std::make_tuple(x, f, df, false));
  // std::sort(cache_.begin(), cache_.end());
  // std::cout << f << "----" << df << std::endl;
  return std::make_tuple(f, df);
}
void CachedFuncDeriv1D::add_to_cache(const std::vector<double> &xs) const {
  for (auto &&x : xs) {
    double af, adf; // The ab init ones
    std::tie(af, adf) = (*func_)(x);
    cache_.push_back(std::make_tuple(x, af, adf, false));
  }
  std::sort(cache_.begin(), cache_.end());
  return;
}
auto lj(double r) {
  double r_6(pow(r, -6));
  double v(4 * (r_6 * r_6 - r_6));
  double dv(4 * ((-12) * r_6 * r_6 / r - (-6) * r_6 / r));

  return std::make_tuple(v, dv);
}

double wrapper(const std::vector<double> &x, std::vector<double> &grad,
               void *f_data) {
  FuncDeriv1D pot = *reinterpret_cast<FuncDeriv1D *>(f_data);
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

auto find_local_minimum(const FuncDeriv1D &pot, double lower, double upper,
                        size_t nsamp = 11, double xtol_err = 1.0e-8) {
  std::valarray<double> xs(nsamp); // TODO
  std::valarray<double> fs(nsamp); // TODO

  std::vector<double> x{0.5 * (lower + upper)};
  // Compute the pointwise values to get the optimization range
  if (nsamp >= 3) {
    {
      double width = upper - lower;
      double step = width / ((signed)nsamp - 1);
      for (size_t i = 0; i != nsamp; ++i) {
        xs[i] = lower + (signed)i * step;
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
    x[0] = xs[min_i];
    if (min_i == 0 || min_i == nsamp - 1) {
      std::cerr << "WARNING:" << std::endl
                << __FILE__ << " " << __LINE__
                << ": Searching range may not be good, with index = " << min_i
                << " and val = " << min_f << std::endl;
      if (min_i == 0)
        lower = xs[0];
      if (min_i == nsamp - 1)
        upper = xs[nsamp - 1];
      x[0] = 0.5 * (lower + upper);
    }
  }
  nlopt::opt opt(nlopt::LN_COBYLA, 1);
  opt.set_xtol_rel(xtol_err);
  opt.set_lower_bounds(lower);
  opt.set_upper_bounds(upper);
  FuncDeriv1D p(pot); // NLOpt not support a const data, so...
  opt.set_min_objective(wrapper, &p);
  double f;
  opt.optimize(x, f);
  return std::make_tuple(x[0], f);
}

auto find_local_maximum(const FuncDeriv1D &pot, double lower, double upper,
                        size_t nsamp = 11, double xtol_err = 1.0e-8) {
  FuncDeriv1D mpot = [&](double x) {
    double v, dv;
    std::tie(v, dv) = pot(x);
    return std::make_tuple(-v, -dv);
  };
  double x, f;
  std::tie(x, f) = find_local_minimum(mpot, lower, upper, nsamp, xtol_err);
  return std::make_tuple(x, -f);
}

/**
 * The 1-dimension function with its 1st order derivative.
 */
typedef FuncDeriv1D Pot1D;

class Pot1DFeatures {
  private:
    double sigma_;
    double epsilon_;
    double r_min_;
    Pot1D *const pot_;

  public:
    Pot1DFeatures(Pot1D &pot) : pot_(&pot) {
      std::tie(r_min_, epsilon_) = find_local_minimum(
          *pot_, 1.0, 5.0); // TODO: I believe this range is safe...
      epsilon_ *= -1;
      std::cout << "epsilon\tr_min" << std::endl;
      std::cout << epsilon_ << "\t" << r_min_ << std::endl;

      double _;
      Pot1D potabs = [=](double x) {
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
    Pot1D *pot() const { return pot_; }
};

class IntegralRange {
  private:
    const Pot1DFeatures *const pf_;
    Pot1D *const pot_;
    double r_C_;
    double E_C_;
    std::vector<double> xs_;
    std::vector<double> ys_;
    std::tuple<double, double> y_(double r) const {
      double v, dv;
      std::tie(v, dv) = (*pot_)(r);
      return std::make_tuple(v + 0.5 * r * dv,
                             0.0); // here we do not really need the derivative.
    }

  public:
    IntegralRange(Pot1DFeatures &pf) : pf_(&pf), pot_(pf.pot()) {

      auto y = [&](double r) { return this->y_(r); };
      std::tie(r_C_, E_C_) =
          find_local_maximum(y, pf_->r_min(), pf_->r_min() * 3);
      std::cout << "r_C\tE_C" << std::endl;
      std::cout << r_C_ << "\t" << E_C_ << std::endl;

      // Generate needed grid for y function
      const size_t gridsize(30);
      xs_.resize(gridsize);
      ys_.resize(gridsize);
      for (size_t i = 0; i != gridsize; ++i) {
        xs_[i] = r_C_ + 0.1 * i;
        double _;
        std::tie(ys_[i], _) = y_(xs_[i]);
      }
    }
    const double &r_C() const { return r_C_; }
    const double &E_C() const { return E_C_; }
    std::tuple<double, double> r_range(double E) const {
      double r_O, r_Op;
      double _;
      FuncDeriv1D potabs = [=](double x) {
        double v, dv;
        std::tie(v, dv) = this->y_(x);
        return std::make_tuple(fabs(v - E), (v - E > 0 ? 1 : -1) * dv);
      };
      size_t l, r, m;
      double lower, upper;
      l = 0;
      r = xs_.size() - 1;
      if (E > ys_[r]) {
        while (r - l != 1) {
          m = (l + r) / 2;
          (E > ys_[m] ? r : l) = m;
        }
        lower = xs_[l];
        upper = xs_[r];
      } else {
        lower = xs_[r];
        upper = 1.0 / 0.0; // inf
      }
      // std::cerr << lower << "  " << upper << std::endl;
      std::tie(r_O, _) = find_local_minimum(potabs, lower, upper, 1);
      r_Op = 0.0; // TODO
      return std::make_tuple(r_O, r_Op);
    }
};

int main() {
  std::cout << std::setprecision(8);
  nlopt::opt opt(nlopt::LN_COBYLA, 1);

  Pot1D pot(lj);
  Pot1DFeatures pf(pot);
  std::cout << " " << pf.r_min() << " " << pf.sigma() << " " << pf.epsilon()
            << " " << std::endl;
  IntegralRange ir(pf);
  for (double E = 0.0; E <= 0.8; E += 0.1)
    std::cout << "r_O " << E << " " << std::get<0>(ir.r_range(E)) << std::endl;
  return 0;
}
