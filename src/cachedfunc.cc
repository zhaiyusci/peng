#include "cachedfunc.hh"
#include<cmath>
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
    // std::fstream fs("test.txt", std::ios_base::out | std::fstream::app);
    // fs << x << "   " << f << "  " << df << std::endl;
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
