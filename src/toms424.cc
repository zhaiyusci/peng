#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <tuple>
#include <utility>
#include <vector>

// using std::cout;
// using std::endl;
using std::function;
using std::make_tuple;
using std::tuple;
using std::vector;

double ccquad_fortran(double (*f)(double *), double *a, double *b,
                      double *tolerr, int *limit, double *esterr, int *used,
                      double *csxfrm) {
  auto f_fort = [&](double x) -> double { return f(&x); };
  tuple<double, double, size_t, vector<double>> ccquad(
      const function<double(double)> &f, double a, double b, double tolerr,
      size_t limit);
  double quadrature(0.0);
  vector<double> cs;
  std::tie(quadrature, *esterr, *used, cs) =
      ccquad(f_fort, *a, *b, *tolerr, *limit);
  for (size_t i = 0; i != (size_t)(*used); ++i) {
    csxfrm[i] = cs[i];
  }
  return quadrature;
}

extern "C" double ccquad_(double (*f)(double *), double *a, double *b,
                          double *tolerr, int *limit, double *esterr, int *used,
                          double *csxfrm) {
  auto res = ccquad_fortran(f, a, b, tolerr, limit, esterr, used, csxfrm);
  return res;
}

class Jterator {
  private:
    const std::vector<size_t> l;
    friend class iterator;

  public:
    Jterator(const std::vector<size_t> &l)
        : l(l), begin_(iterator(this, std::vector<size_t>(l.size(), 1))),
          end_(++iterator(this, l)){};
    class iterator {
      private:
        Jterator *jterator;
        std::vector<size_t> val;

      public:
        iterator(Jterator *jterator, const std::vector<size_t> &val)
            : jterator(jterator), val(val) {}
        iterator &operator++() {
          size_t s = val.size();

          for (size_t i = 0; i != s; ++i) { // i in [0, s)
            if (s - 1 - i != 0) {
              val[s - 1 - i] += jterator->l[s - 2 - i];
            } else {
              ++val[0];
            }
            if (val[s - 1 - i] <= jterator->l[s - 1 - i]) {
              for (size_t j = s - i; j != s; ++j) {
                val[j] = val[j - 1];
              }
              break;
            }
          }
          return *this;
        }
        const std::vector<size_t> &operator*() { return val; }
        bool operator==(iterator rhs) {
          // return val[0] == 2;
          return jterator == rhs.jterator && val == rhs.val;
        }
        bool operator!=(iterator rhs) { return !(*this == rhs); }
    };

    iterator begin_;
    iterator end_;

    const iterator &begin() const { return begin_; }
    const iterator &end() const { return end_; }
};

tuple<double, double, size_t, vector<double>>
ccquad(const function<double(double)> &f, double a, double b, double tolerr,
       size_t limit) {
  // return make_tuple(quadrature, used, csxfrm_container);
  size_t used;
  double esterr, oldint, newint;

  void r3pass(size_t n2, size_t m, size_t length, double *const x0,
              double *const x1, double *const x2);

  const double centre = (a + b) * 0.5;
  const double width = (b - a) * 0.5;
  // maxp = 2 * pow(3, (m_max + 1));
  // maxp = limit < maxp ? limit : maxp;
  size_t m_max = (int)ceil(log(limit / 2.0) / log(3.0)) - 1;
  vector<size_t> l(m_max, 1);
  // size_t maxp = limit;
  limit = 2 * pow(3, (m_max)) + 1;
  vector<double> csxfrm_container(limit);
  double *csxfrm =
      csxfrm_container.data(); // Use STL to do the memory management

  // Initialization, for 2*3^1+1=7 pivot quadrature
  {
    csxfrm[0] = f(b);
    csxfrm[6] = f(a);
    double shift = -width * (double)sqrt(3) * 0.5;
    csxfrm[1] = f(centre - shift);
    csxfrm[5] = f(centre + shift);
    shift = -width * 0.5;
    csxfrm[2] = f(centre - shift);
    csxfrm[4] = f(centre + shift);
    csxfrm[3] = f(centre);
    //  evaluate the factored n=6 cosine transform.;
    double t1 = csxfrm[0] + csxfrm[6];
    double t2 = csxfrm[0] - csxfrm[6];
    double t3 = 2.0 * csxfrm[3];
    double t4 = csxfrm[1] + csxfrm[5];
    double t5 = (csxfrm[1] - csxfrm[5]) * (double)sqrt(3);
    double t6 = csxfrm[2] + csxfrm[4];
    double t7 = csxfrm[2] - csxfrm[4];
    double t8 = t1 + 2.0 * t6;
    double t9 = 2.0 * t4 + t3;
    double t10 = t2 + t7;
    double t11 = t1 - t6;
    double t12 = t4 - t3;
    csxfrm[0] = t8 + t9;
    csxfrm[1] = t10 + t5;
    csxfrm[2] = t11 + t12;
    csxfrm[3] = t2 - 2.0 * t7;
    csxfrm[4] = t11 - t12;
    csxfrm[5] = t10 - t5;
    csxfrm[6] = t8 - t9;
    used = 7;
    oldint = (t1 + 2.0 * t3) / 3.0;
  }

  bool firstrun = true;
  size_t n = 6;
  do {
    double fund;
    if (!firstrun) {
      for (size_t j = 1; j != m_max; ++j) {
        l[j - 1] = l[j];
      }
      l[m_max - 1] = 3 * l[m_max - 2];
      size_t j = used;
      fund = (double)M_PI / (3.0 * n);
      const Jterator jter(l);
      for (auto &&jj : jter) {
        double angle = fund * (3.0 * jj[m_max - 1] - 2);
        double shift = -width * cos(angle);
        double t1 = f(centre - shift);
        double t3 = f(centre + shift);
        shift = -width * sin(angle);
        double t2 = f(centre + shift);
        double t4 = f(centre - shift);
        double t5 = t1 + t3;
        double t6 = t2 + t4;
        csxfrm[j++] = t5 + t6;
        csxfrm[j++] = t1 - t3;
        csxfrm[j++] = t5 - t6;
        csxfrm[j++] = t2 - t4;
      }
      //  do radix 3 passes of fast fourier transform.

      {
        size_t step = 4;
        do {
          size_t j1 = used + step;
          size_t j2 = used + 2 * step;
          r3pass(n * 2, step, n * 2 - 2 * step, csxfrm + used, csxfrm + j1,
                 csxfrm + j2);
          step = 3 * step;
        } while (step < n);
      }

      //
      //  combine results.
      //
      //  first do j=0 and j=n.
      //
      {
        double t1 = csxfrm[0];
        double t2 = csxfrm[used];
        csxfrm[0] = t1 + 2.e0 * t2;
        csxfrm[used] = t1 - t2;
        t1 = csxfrm[n];
        t2 = csxfrm[n * 2 + 1];
        csxfrm[n] = t1 + t2;
        csxfrm[n * 2 + 1] = t1 - 2.e0 * t2;
      }
      //  now do remaining values of j.
      for (size_t j = 1; j <= n - 1; ++j) {
        size_t j1 = n + j;
        size_t j2 = n * 3 - j;
        double angle = fund * double(j);
        double c = cos(angle);
        double s = sin(angle);
        double t1 = c * csxfrm[j1 + 1] - s * csxfrm[j2 + 1];
        double t2 = (s * csxfrm[j1 + 1] + c * csxfrm[j2 + 1]) * (double)sqrt(3);
        csxfrm[j1 + 1] = csxfrm[j] - t1 - t2;
        csxfrm[j2 + 1] = csxfrm[j] - t1 + t2;
        csxfrm[j] = csxfrm[j] + 2.e0 * t1;
      }
      //  now unscramble.
      {
        double t1 = csxfrm[n * 2];
        double t2 = csxfrm[n * 2 + 1];
        for (size_t j = 0; j < n - 1; ++j) {
          size_t j1 = used + j;
          size_t j2 = n * 2 + j;
          csxfrm[j2] = csxfrm[j1];
          csxfrm[j1] = csxfrm[j2 + 2];
        }
        n *= 3;
        csxfrm[n - 1] = t1;
        csxfrm[n] = t2;
        used = n + 1;
      }
    }

    {
      firstrun = false;
      newint = 0.5 * csxfrm[used - 1] / (1 - (signed)(n * n));
      for (size_t j = 1; j <= n - 3; j += 2) {
        size_t j_rev = n - j;
        // std::cerr << j_rev << std::endl;
        newint += csxfrm[j_rev - 1] / ((signed)j_rev * (2 - (signed)j_rev));
      }
      newint += 0.5 * csxfrm[0];
    }
    //
    //  test if done.
    //  test if estimated error adequate.
    esterr = fabs(oldint * 3.0 - newint);
    // std::cerr <<std::setprecision(15)<< "inloop =====" << std::endl;
    // std::cerr << oldint << std::endl;
    // std::cerr << newint << std::endl;
    // std::cerr << esterr << std::endl;
    //
    //  insert the following four statements to trace program flow.
    //     sclint=width*newint/double(n/2))
    //     sclerr=width*(oldint*3.e0-newint)/double(n/2)
    //     write(6,900) n,sclint,sclerr
    // 900  format ( 3h n=,i5,23h integral estimated as ,e15.8,
    //    *  7h error ,e15.8 )
    if (fabs(newint) * tolerr >= esterr)
      break;
    //  if estimated error too large, refine sampling if permitted.
    oldint = newint;
  } while (3 * n + 1 <= limit);
  //  if refinement not permitted, or if estimated error
  //  satisfactory, rescale answers and return.
  //  insert the following two statements to trace program flow.
  //     write (6,910)
  // 910 format ( 25h refinement not permitted )
  double quadrature = width * newint / (n / 2.0);
  // std::cerr << "esterr :>" << std::endl;
  // std::cerr << "raw esterr" << std::endl;
  // std::cerr << esterr << std::endl;
  esterr = width * esterr / (n / 2.0);
  // std::cerr << "ral esterr" << std::endl;
  // std::cerr << esterr << std::endl;
  // std::cerr << width << "   "<< n << std::endl;

  if (esterr > tolerr) {
    std::cerr
        << "=== WARNING ===" << std::endl
        << "The estimate relative error is GREATER than the one user assigned."
        << std::endl
        << "That is mainly because user choose a small number of "
        << "function evaluations (`limit`)." << std::endl
        << "Because of the structure of program, the actual maxium "
        << "number of function evaluations is " << limit << "." << std::endl
        << "To break the limitation, set limit >= " << (limit - 1) * 3 + 1
        << " in the next run." << std::endl;
  }

  return make_tuple(quadrature, esterr, used, std::move(csxfrm_container));
}

void r3pass(const size_t n2, const size_t m, const size_t length,
            double *const x0, double *const x1, double *const x2) {
  const double hafrt3(0.5 * sqrt(3.0));

  //  do all transforms for c hat = 0, i.e., twiddle factor unity.
  for (size_t k = 0; k < n2; k += m * 3) {
    double rsum = (x1[k] + x2[k]);
    double rdiff = (x1[k] - x2[k]) * hafrt3;
    x1[k] = x0[k] - rsum * 0.5;
    x2[k] = rdiff;
    x0[k] = x0[k] + rsum;
  }
  //  do all transforms for c hat = cap c/2, i.e., twiddle factor.
  //  e(b/6)
  for (size_t k = m / 2; k < n2; k += m * 3) {
    double rsum = (x1[k] + x2[k]) * hafrt3;
    double rdiff = (x1[k] - x2[k]);
    x1[k] = x0[k] - rdiff;
    x2[k] = rsum;
    x0[k] = x0[k] + rdiff * 0.5;
  }
  //  do all transforms for remaining values of c hat.  observe
  //  that c hat and cap c-c hat must be paired.
  //  choose a frequency index.
  for (size_t j = 1; j <= (m - 1) / 2; ++j) {
    size_t j0 = j + 1;
    size_t j1 = m - j + 1;
    //  compute the twiddle factor.
    double angle = 2 * M_PI * j / (m * 3);
    double c1 = cos(angle);
    double s1 = sin(angle);
    double c2 = c1 * c1 - s1 * s1;
    double s2 = 2.e0 * s1 * c1;
    //  choose the replication.
    for (size_t k0 = j0 - 1; k0 < n2; k0 += m * 3) {
      size_t k1 = k0 - j0 + j1;
      //  obtain twiddled values.
      double r0 = x0[k0];
      double i0 = x0[k1];
      double r1 = c1 * x1[k0] - s1 * x1[k1];
      double i1 = s1 * x1[k0] + c1 * x1[k1];
      double r2 = c2 * x2[k0] - s2 * x2[k1];
      double i2 = s2 * x2[k0] + c2 * x2[k1];
      //  compute transforms and return in place.;
      double rsum = r1 + r2;
      double rdiff = (r1 - r2) * hafrt3;
      double rsum2 = r0 - 0.5 * rsum;
      double isum = i1 + i2;
      double idiff = (i1 - i2) * hafrt3;
      double idiff2 = i0 - 0.5 * isum;
      x0[k0] = r0 + rsum;
      x0[k1] = rsum2 + idiff;
      x1[k0] = rsum2 - idiff;
      x1[k1] = rdiff + idiff2;
      x2[k0] = rdiff - idiff2;
      x2[k1] = i0 + isum;
    }
  }

  return;
}
