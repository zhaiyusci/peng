#include "mathtools.hh"
#include <nlopt.hpp>
#include <tuple>
#include <valarray>
#include <vector>
namespace dlt {
double FuncDeriv1D::derivative(double x) const {
  std::cerr << "Fallback to the finite differential method." << std::endl;
  const double dx = 1.0e-6;
  return (value(x + dx) - value(x - dx)) / dx / 2;
}
FuncDeriv1D &FuncDeriv1D::operator-() {
  if (pneg_ == nullptr)
    pneg_.reset(new NegFuncDeriv1D(this));
  return *pneg_;
}

double wrapper(const std::vector<double> &x, std::vector<double> &grad,
               void *f_data) {
  FuncDeriv1D *pot = reinterpret_cast<FuncDeriv1D *>(f_data);
  double v;
  v = pot->value(x[0]);
  return v;
}

double wrapper_deriv(const std::vector<double> &x, std::vector<double> &grad,
                     void *f_data) {
  FuncDeriv1D *pot = reinterpret_cast<FuncDeriv1D *>(f_data);
  double v;
  v = pot->value(x[0]);
  grad.resize(1);
  grad[0] = pot->derivative(x[0]);
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

std::tuple<double, double> find_local_minimum(FuncDeriv1D &func, double lower,
                                              double upper, size_t nsamp,
                                              double xtol_err) {
  std::valarray<double> xs(nsamp);
  std::valarray<double> fs(nsamp);

  std::vector<double> x{0.5 * (lower + upper)};
  // Compute the pointwise values to get the optimization range
  if (nsamp >= 3) {
    {
      double width = upper - lower;
      double step = width / ((signed)nsamp - 1);
      for (size_t i = 0; i != nsamp; ++i) {
        xs[i] = lower + (signed)i * step;
        fs[i] = func(xs[i]);
      }
    }
    size_t min_i;
    double min_f;
    std::tie(min_f, min_i) = min_and_index(fs);
    size_t l, u;
    l = min_i - 1;
    u = min_i + 1;
    x[0] = xs[min_i];
    if (min_i == 0 || min_i == nsamp - 1) {
      std::cerr << "WARNING:\n"
                << __FILE__ << " " << __LINE__
                << ": Searching range may not be good, with index = " << min_i
                << " and val = " << min_f
                << ", which evaluate with x = " << xs[min_i] << ".\n";
      if (min_i == 0)
        l = 0;
      if (min_i == nsamp - 1)
        u = nsamp - 1;
    }
    lower = xs[l];
    upper = xs[u];
    x[0] = 0.5 * (lower + upper);
  }
  double f;
  if (func.provide_derivative()) {
    nlopt::opt opt(nlopt::LD_LBFGS, 1);
    opt.set_xtol_rel(xtol_err);
    opt.set_lower_bounds(lower);
    opt.set_upper_bounds(upper);
    opt.set_min_objective(wrapper_deriv, &func);
    opt.optimize(x, f);
  } else {
    nlopt::opt opt(nlopt::LN_COBYLA, 1);
    opt.set_xtol_rel(xtol_err);
    opt.set_lower_bounds(lower);
    opt.set_upper_bounds(upper);
    opt.set_min_objective(wrapper, &func);
    opt.optimize(x, f);
  }
  return std::make_tuple(x[0], f);
}

std::tuple<double, double> find_local_maximum(FuncDeriv1D &func, double lower,
                                              double upper, size_t nsamp,
                                              double xtol_err) {
  double x, f;
  auto &&mfunc = -func;
  std::tie(x, f) = find_local_minimum(mfunc, lower, upper, nsamp, xtol_err);
  return std::make_tuple(x, -f);
}

class AbsDiff : public FuncDeriv1D {
private:
  FuncDeriv1D *const func_;
  const double target_;

public:
  AbsDiff(FuncDeriv1D &func, double target) : func_(&func), target_(target) {

    func_->value(1.00);
  };
  double value(double x) const { return fabs(func_->value(x) - target_); }
  bool provide_derivative() const { return false; }
};

double find_local_root(FuncDeriv1D &func, double target, double lower,
                       double upper, size_t nsamp, double xtol_err) {
  double x, f;
  func(1.0);
  AbsDiff absfunc(func, target);
  std::tie(x, f) = find_local_minimum(absfunc, lower, upper, nsamp, xtol_err);
  return x;
}

LocalRoot::LocalRoot(FuncDeriv1D &func, double lower, double upper,
                     size_t nsamp, double xtol_err)
    : pfunc_(&func), xtol_err_(xtol_err) {
  if (nsamp >= 3) {
    double width = upper - lower;
    double step = width / ((signed)nsamp - 1);
    xs_.resize(nsamp);
    fs_.resize(nsamp);
    for (size_t i = 0; i != nsamp; ++i) {
      xs_[i] = lower + (signed)i * step;
      fs_[i] = func(xs_[i]);
    }
  }
}

inline bool inrange(double target, double lower, double upper) {
  return (target >= lower && target <= upper) ||
         (target <= lower && target >= upper);
}

double LocalRoot::operator()(double target) const {
  // We assume that the function is monotonic in the range of lower:upper
  size_t min_i;
  double min_f;
  std::tie(min_f, min_i) = min_and_index(fs_);
  double lower, upper;
  size_t l, r, m;
  l = 0;
  r = xs_.size() - 1;
  if (inrange(target, fs_[l], fs_[r])) {
    while (r - l != 1) {
      m = (l + r) / 2;
      (inrange(target, fs_[l], fs_[m]) ? r : l) = m;
    }
    lower = xs_[l];
    upper = xs_[r];
  } else {
    if ((fs_[l] < fs_[r]) == (target > fs_[r])) {
      lower = xs_[r] - 1.0e-6; //
      upper = xs_[r] + 1.0e-6; // +inf
    } else {
      lower = xs_[l] - 1.0e-6; // -inf
      upper = xs_[l] + 1.0e-6; //
    }
  }
  double res = find_local_root(*pfunc_, target, lower, upper, 1, xtol_err_);
  return res;
}

} // namespace dlt
