#define EIGEN_NO_DEBUG
#include "transport.hh"
#include "mathtools.hh"
#include <algorithm>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <tuple>
#include <vector>

static const double kB = 1.380649e-23;      // BY DEFINITION
static const double amu = 1.6605390666e-27; // CODATA2018
                                            //
namespace {
/**
 * @brief Store a 6-dimension array.
 */
class Hexa {
private:
  size_t size0_;
  size_t size1_;
  size_t size2_;
  size_t size3_;
  size_t size4_;
  size_t size5_;
  size_t sizet_;

  std::vector<double> data_;

public:
  size_t idx(size_t i0, size_t i1, size_t i2, size_t i3, size_t i4, size_t i5) {
    return
        // clang-format off
        i0 * size0_ 
      + i1 * size1_ 
      + i2 * size2_ 
      + i3 * size3_ 
      + i4 * size4_ 
      + i5 * size5_;
    // clang-format on
  }
  Hexa() { resize(0, 0, 0, 0, 0, 0); }
  Hexa(size_t i0, size_t i1, size_t i2, size_t i3, size_t i4, size_t i5) {
    resize(i0, i1, i2, i3, i4, i5);
  }
  void resize(size_t i0, size_t i1, size_t i2, size_t i3, size_t i4,
              size_t i5) {
    // clang-format off
    size0_ = i1 * i2 * i3 * i4 * i5 * 1;
    size1_ =      i2 * i3 * i4 * i5 * 1;
    size2_ =           i3 * i4 * i5 * 1;
    size3_ =                i4 * i5 * 1;
    size4_ =                     i5 * 1;
    size5_ =                          1;
    // clang-format on
    sizet_ = i0 * i1 * i2 * i3 * i4 * i5 * 1;
    data_.resize(sizet_);
  }
  double &operator()(size_t i0, size_t i1, size_t i2, size_t i3, size_t i4,
                     size_t i5) {
    return data_[idx(i0, i1, i2, i3, i4, i5)];
  }
  size_t size() const { return sizet_; }
};
} // namespace

namespace dlt {
class AlphaImpl {
private:
  size_t propertyorder_;
  double mass0_;
  double mass1_;

  // Workspace
  double temperature_;
  Eigen::MatrixXd Omega00_;
  Eigen::MatrixXd Omega01_;
  Eigen::MatrixXd Omega11_;
  Hexa Hint;

  double molefraction0_;
  double molefraction1_;

  double mtot_;
  double m0_;
  double m1_;

  std::vector<double> lambda_;
  std::vector<double> D12_;
  std::vector<double> DT_;

public:
  AlphaImpl(size_t propertyorder, double mass0, double mass1)
      : propertyorder_(propertyorder), mass0_(mass0), mass1_(mass1),
        temperature_(-1.0),
        Hint(propertyorder + 1, propertyorder + 1, 2, 2, 2, 2),
        lambda_(propertyorder), D12_(propertyorder), DT_(propertyorder) {
    mtot_ = mass0 + mass1;
    m0_ = mass0 / mtot_;
    m1_ = mass1 / mtot_;
  }

  const std::vector<double> &lambda() const { return lambda_; }
  const std::vector<double> &D12() const { return D12_; }
  const std::vector<double> &DT() const { return DT_; }

  double mass(size_t i) const {
    switch (i) {
    case 0:
      return mass0_;
    case 1:
      return mass1_;
    default:
      return 0;
    }
  }

  const Eigen::MatrixXd &Omega(size_t i0, size_t i1) const {
    if (i0 != i1) {
      return Omega01_;
    } else if (i0 == 0) {
      return Omega00_;
    } else {
      return Omega11_;
    }
  }

  void evaluate(double temperature, const Eigen::MatrixXd &Omega00,
                const Eigen::MatrixXd &Omega01, const Eigen::MatrixXd &Omega11,
                double molefraction0) {
    if (temperature_ != temperature) {
      Omega00_ = Omega00;
      Omega01_ = Omega01;
      Omega11_ = Omega11;
      temperature_ = temperature;
      for (size_t p = 0; p <= propertyorder_; ++p) {
        for (size_t q = 0; q <= propertyorder_; ++q) {
          for (size_t i0 = 0; i0 != 2; ++i0) {
            for (size_t i1 = 0; i1 != 2; ++i1) {
              for (size_t io0 = 0; io0 != 2; ++io0) {
                for (size_t io1 = 0; io1 != 2; ++io1) {
                  Hint(p, q, i0, i1, io0, io1) =
                      bracket_int_H(p, q, i0, i1, io0, io1);
                }
              }
            }
          }
        }
      }
    }

    molefraction0_ = molefraction0;
    molefraction1_ = 1 - molefraction0;

    Eigen::MatrixXd raw_D_mat(2 * propertyorder_ + 1, 2 * propertyorder_ + 1);
    for (int p = 1; p <= (int)propertyorder_; ++p) {
      for (int q = 1; q <= (int)propertyorder_; ++q) {
        // clang-format off
        raw_D_mat(propertyorder_ + p, propertyorder_ + q) = A( p,  q); // SE
        raw_D_mat(propertyorder_ + p, propertyorder_ - q) = A( p, -q); // SW
        raw_D_mat(propertyorder_ - p, propertyorder_ - q) = A(-p, -q); // NW
        raw_D_mat(propertyorder_ - p, propertyorder_ + q) = A(-p,  q); // NE
        // clang-format on
      }
    }
    for (int p = 1; p <= (int)propertyorder_; ++p) {
      // clang-format off
      raw_D_mat(propertyorder_ + 0, propertyorder_ + p) = A( 0,  p); // E
      raw_D_mat(propertyorder_ + p, propertyorder_ + 0) = A( p, -0); // S
      raw_D_mat(propertyorder_ + 0, propertyorder_ - p) = A(-0, -p); // W
      raw_D_mat(propertyorder_ - p, propertyorder_ + 0) = A(-p,  0); // N
      // clang-format on
    }
    raw_D_mat(propertyorder_, propertyorder_) = A(0, 0); // C
                                                         //
    std::cerr << "raw_D_mat:" << '\n';
    std::cerr << raw_D_mat << std::endl;

    for (size_t i = 1; i <= propertyorder_; ++i) {
      Eigen::MatrixXd A_mat(2 * i, 2 * i);
      // clang-format off
      A_mat.block(i, i, i, i) = raw_D_mat.block(propertyorder_ + 1, propertyorder_ + 1, i, i); // SE
      A_mat.block(i, 0, i, i) = raw_D_mat.block(propertyorder_ + 1, propertyorder_ - i, i, i); // SW
      A_mat.block(0, 0, i, i) = raw_D_mat.block(propertyorder_ - i, propertyorder_ - i, i, i); // NW
      A_mat.block(0, i, i, i) = raw_D_mat.block(propertyorder_ - i, propertyorder_ + 1, i, i); // NE
      // clang-format on
      Eigen::VectorXd alpha_vec = Eigen::VectorXd::Zero(2 * i);
      // clang-format off
      alpha_vec(i - 1) = -15.0 / 4.0 * molefraction1_ * sqrt(2 * kB * temperature / mass1_ / amu); // ALPHA_[-1]
      alpha_vec(i    ) = -15.0 / 4.0 * molefraction0_ * sqrt(2 * kB * temperature / mass0_ / amu); // ALPHA_[ 1]
      // clang-format on
      auto a_vec = A_mat.inverse() * alpha_vec;
      lambda_[i - 1] = -5.0 / 4.0 * kB *
                       sqrt(2 * kB * temperature / mtot_ / amu) *
                       (molefraction0_ / sqrt(m0_) * a_vec(i) +
                        molefraction1_ / sqrt(m1_) * a_vec(i - 1));
      Eigen::MatrixXd D_mat(2 * i + 1, 2 * i + 1);
      D_mat = raw_D_mat.block(propertyorder_ - i, propertyorder_ - i, 2 * i + 1,
                              2 * i + 1);
      Eigen::VectorXd delta_vec = Eigen::VectorXd::Zero(2 * i + 1);
      delta_vec(i) = 3.0 / 2.0 * sqrt(2 * kB * temperature / mtot_ / amu);
      Eigen::VectorXd d_vec = D_mat.inverse() * delta_vec;
      D12_[i - 1] = (kB * temperature) / 1.013e5 * 1.0 / 2.0 * molefraction0_ *
                    molefraction1_ * sqrt(2 * kB * temperature / mtot_ / amu) *
                    d_vec(i);
      DT_[i - 1] = (kB * temperature) / 1.013e5 * (-5.0) / 4.0 *
                   molefraction0_ * molefraction1_ *
                   sqrt(2 * kB * temperature / mtot_ / amu) *
                   (molefraction0_ / sqrt(m0_) * d_vec(i + 1) +
                    molefraction1_ / sqrt(m1_) * d_vec(i - 1));
    }
    return;
  }

  double bracket_int_H(int p, int q, size_t i0, size_t i1, size_t oi0,
                       size_t oi1) {

    const Eigen::MatrixXd &Omega = this->Omega(oi0, oi1);

    if (i0 != i1) {
      return H12(p, q, i0, i1, Omega);
    } else {
      if (oi0 != oi1) {
        return H1(p, q, i0, i1, Omega);
      } else {
        return HSG(p, q, i0, i1, Omega);
      }
    }
  }

  double H12(int p, int q, size_t i0, size_t i1,
             const Eigen::MatrixXd &Omega) { // eq 109
    double s_sum = 0.0;
    double m0 = mass(i0) / mtot_;
    double m1 = mass(i1) / mtot_;
    for (int l = 1; l <= std::min(p, q) + 1; ++l) {
      for (int r = l; r <= p + q + 2 - l; ++r) {
        double a_sum = 0.0;
        // remenber that 1-delta=0 cause zero value
        for (int i = l - 1; i <= std::min({p, q, r, p + q + 1 - r});
             ++i) { // eq 110
          a_sum += pow(8, i) * fact(p + q - 2 * i) /
                   (fact(p - i) * fact(q - i)) * pow(-1, l) /
                   (fact(l) * fact(i + 1 - l)) * pow(-1, r + i) /
                   (fact(r - i) * (fact(p + q + 1 - i - r))) * fact(r + 1) /
                   fact(2 * r + 2) * fact(2 * (p + q + 2 - i)) /
                   fact(p + q + 2 - i) * pow(2, 2 * r) / pow(4.0, p + q + 1) *
                   ((i + 1 - l) * (p + q + 1 - i - r) - (l * (r - i)));
        }

        double a = a_sum * pow(m1, p + 0.5) * pow(m0, q + 0.5);

        s_sum += a * Omega(l, r);
      }
    }
    return s_sum * 8;
  }

  double H1(int p, int q, size_t i0, [[maybe_unused]] size_t i1,
            const Eigen::MatrixXd &Omega) {
    double m0 = mass(i0) / mtot_;
    double m1 = mass(1 - i0) / mtot_; // if 0, get 1, and vice versa.
    double G = (m0 - m1) / m1;
    double F = (m0 * m0 + m1 * m1) / (2 * m0 * m1);
    double s_sum = 0.0;
    for (int l = 1; l <= std::min(p, q) + 1; ++l) {
      for (int r = l; r <= p + q + 2 - l; ++r) {

        double a = 0.0;
        double a_sum1 = 0.0;
        for (int i = l - 1; i <= std::min({p, q, r, (p + q + 1 - r)}); ++i) {
          a_sum1 = pow(8, i) * fact(p + q - 2 * i) /
                   (fact(p - i) * fact(q - i)) * 1.0 /
                   (fact(l) * fact(i + 1 - l)) * pow(-1, r + i) /
                   (fact(r - i) * fact(p + q + 1 - i - r)) * fact(r + 1) /
                   fact(2 * r + 2) * fact(2 * (p + q + 2 - i)) /
                   fact(p + q + 2 - i) * pow(2, 2 * r) / pow(4, p + q + 1);
          double a_sum2 = 0.0;
          for (int w = 0; w <= std::min({p, q, p + q + 1 - r}) - i; ++w) {
            a_sum2 += pow(F, i + 1 - l) * pow(G, w) / fact(w) *
                      poch(p + q + 2 - i - r - w, w) * poch(p + 1 - i - w, w) /
                      poch(2 * (p + q + 2 - i) - 2 * w + 1, w) *
                      poch(p + q + 3 - i - w, w) * poch(q + 1 - i - w, w) /
                      poch(2 * (p + q + 2 - i) - w + 1, w) * pow(2, 2 * w - 1) *
                      pow(m0, i) * pow(m1, i) * pow(m1, (p + q - 2 * i - w)) /
                      poch(p + q + 1 - 2 * i - w, w) *
                      (2.0 * m0 / F * ((i + 1 - l) * (p + q + 1 - i - r - w)) -
                       2.0 * m1 * l * (r - i));
          }
          a += a_sum1 * a_sum2;
        }

        s_sum += a * Omega(l, r);
      }
    }
    return s_sum * 8;
  }

  double HSG(int p, int q, size_t i0, size_t i1, const Eigen::MatrixXd &Omega) {
    double s_sum = 0.0;
    for (int l = 1; l <= std::min(p, q) + 1; ++l) {
      for (int r = l; r <= p + q + 2 - l; ++r) {
        double a_sum = 0.0;
        for (int i = l - 1; i <= std::min({p, q, r, (p + q + 1 - r)}); ++i) {
          a_sum += pow(8, i) * fact(p + q - 2 * i) /
                   (fact(p - i) * fact(q - i)) * (1 + pow(-1.0, l)) /
                   (fact(l) * fact(i + 1 - l)) * pow(-1, r + i) /
                   (fact(r - i) * fact(p + q + 1 - i - r)) * fact(r + 1) /
                   fact(2 * r + 2) * fact(2 * (p + q + 2 - i)) /
                   fact(p + q + 2 - i) * pow(2.0, 2 * r) / pow(4.0, p + q + 1) *
                   ((i + 1 - l) * (p + q + 1 - i - r) - l * (r - i));
        }

        double a = a_sum / pow(2, p + q + 1);

        s_sum += a * Omega(l, r);
      }
    }
    return s_sum * 8;
  }

  double A(int p, int q) {
    int ap = abs(p);
    int aq = abs(q);
    if (p > 0 && q > 0) {
      return molefraction0_ * molefraction0_ * Hint(ap, aq, 0, 0, 0, 0) +
             molefraction0_ * molefraction1_ * Hint(ap, aq, 0, 0, 0, 1);
    } else if (p > 0 && q < 0) {
      return molefraction0_ * molefraction1_ * Hint(ap, aq, 0, 1, 0, 1);
    } else if (p < 0 && q > 0) {
      return molefraction0_ * molefraction1_ * Hint(ap, aq, 1, 0, 1, 0);
    } else if (p < 0 && q < 0) {
      return molefraction1_ * molefraction1_ * Hint(ap, aq, 1, 1, 1, 1) +
             molefraction0_ * molefraction1_ * Hint(ap, aq, 1, 1, 1, 0);
    } else if (p > 0 && q == 0) {
      return molefraction0_ * molefraction1_ * sqrt(m0_) *
             Hint(ap, 0, 0, 0, 0, 1);
    } else if (p == 0 && q > 0) {
      return molefraction0_ * molefraction1_ * sqrt(m0_) *
             Hint(aq, 0, 0, 0, 0, 1);
    } else if (p < 0 && q == 0) {
      return -molefraction0_ * molefraction1_ * sqrt(m1_) *
             Hint(ap, 0, 1, 1, 1, 0);
    } else if (p == 0 && q < 0) {
      return -molefraction0_ * molefraction1_ * sqrt(m1_) *
             Hint(aq, 0, 1, 1, 1, 0);
    } else { // if (p == 0 &&q == 0) {
      return molefraction0_ * molefraction1_ * (8 * m0_ * m1_ * Omega01_(1, 1));
    }
  }
};

class BetaImpl {
private:
  size_t propertyorder_;
  double mass0_;
  double mass1_;

  // Workspace
  double temperature_;
  Eigen::MatrixXd Omega00_;
  Eigen::MatrixXd Omega01_;
  Eigen::MatrixXd Omega11_;
  Hexa Lint;

  double molefraction0_;
  double molefraction1_;

  double mtot_;
  double m0_;
  double m1_;

  std::vector<double> eta_;

public:
  BetaImpl(size_t propertyorder, double mass0, double mass1)
      : propertyorder_(propertyorder), mass0_(mass0), mass1_(mass1),
        temperature_(-1.0),
        Lint(propertyorder + 1, propertyorder + 1, 2, 2, 2, 2),
        eta_(propertyorder) {
    mtot_ = mass0 + mass1;
    m0_ = mass0 / mtot_;
    m1_ = mass1 / mtot_;
  }

  const std::vector<double> &eta() const { return eta_; }

  double mass(size_t i) const {
    switch (i) {
    case 0:
      return mass0_;
    case 1:
      return mass1_;
    default:
      return 0;
    }
  }

  const Eigen::MatrixXd &Omega(size_t i0, size_t i1) const {
    if (i0 != i1) {
      return Omega01_;
    } else if (i0 == 0) {
      return Omega00_;
    } else {
      return Omega11_;
    }
  }

  void evaluate(double temperature, const Eigen::MatrixXd &Omega00,
                const Eigen::MatrixXd &Omega01, const Eigen::MatrixXd &Omega11,
                double molefraction0) {

    if (temperature_ != temperature) {
      Omega00_ = Omega00;
      Omega01_ = Omega01;
      Omega11_ = Omega11;
      temperature_ = temperature;
      for (size_t p = 0; p <= propertyorder_; ++p) {
        for (size_t q = 0; q <= propertyorder_; ++q) {
          for (size_t i0 = 0; i0 != 2; ++i0) {
            for (size_t i1 = 0; i1 != 2; ++i1) {
              for (size_t io0 = 0; io0 != 2; ++io0) {
                for (size_t io1 = 0; io1 != 2; ++io1) {
                  Lint(p, q, i0, i1, io0, io1) =
                      bracket_int_L(p, q, i0, i1, io0, io1);
                }
              }
            }
          }
        }
      }
    }

    molefraction0_ = molefraction0;
    molefraction1_ = 1 - molefraction0;

    Eigen::MatrixXd raw_B_mat(2 * propertyorder_, 2 * propertyorder_);
    for (int p = 1; p <= (int)propertyorder_; ++p) {
      for (int q = 1; q <= (int)propertyorder_; ++q) {
        // clang-format off
        raw_B_mat(propertyorder_ - 1 + p, propertyorder_ - 1 + q) = B( p,  q); // SE
        raw_B_mat(propertyorder_ - 1 + p, propertyorder_     - q) = B( p, -q); // SW
        raw_B_mat(propertyorder_     - p, propertyorder_     - q) = B(-p, -q); // NW
        raw_B_mat(propertyorder_     - p, propertyorder_ - 1 + q) = B(-p,  q); // NE
        // clang-format on
      }
    }

    std::cerr << "raw_B_mat:" << '\n';
    std::cerr << raw_B_mat << std::endl;

    for (size_t i = 1; i <= propertyorder_; ++i) {
      Eigen::MatrixXd B_mat(2 * i, 2 * i);
      B_mat.block(0, 0, 2 * i, 2 * i) = raw_B_mat.block(
          propertyorder_ - i, propertyorder_ - i, 2 * i, 2 * i); // NW
      Eigen::VectorXd beta_vec = Eigen::VectorXd::Zero(2 * i);
      // clang-format off
      beta_vec(i - 1) = 5.0 / 2.0 * molefraction1_; // ALPHA_[-1]
      beta_vec(i    ) = 5.0 / 2.0 * molefraction0_; // ALPHA_[ 1]
      // clang-format on
      auto b_vec = B_mat.inverse() * beta_vec;
      eta_[i - 1] =
          (molefraction0_ * b_vec(i) + molefraction1_ * b_vec(i - 1)) *
          temperature * kB;
    }
    return;
  }

  double bracket_int_L(int p, int q, size_t i0, size_t i1, size_t oi0,
                       size_t oi1) {

    const Eigen::MatrixXd &Omega = this->Omega(oi0, oi1);

    if (i0 != i1) {
      return L12(p, q, i0, i1, Omega);
    } else {
      if (oi0 != oi1) {
        return L1(p, q, i0, i1, Omega);
      } else {
        return LSG(p, q, i0, i1, Omega);
      }
    }
  }

  double L12(int p, int q, size_t i0, size_t i1,
             const Eigen::MatrixXd &Omega) { // eq 109
    double s_sum = 0.0;
    double m0 = mass(i0) / mtot_;
    double m1 = mass(i1) / mtot_;
    for (int l = 1; l <= std::min(p, q) + 2; ++l) {
      for (int r = l; r <= p + q + 4 - l; ++r) {
        double b_sum = 0.0;
        // REMENBER THAT 1-DELTA=0 CAUSE ZERO VALUE
        if (r != p + q + 3) {
          for (int i = l - 2; i <= std::min({p, q, r, (p + q + 2 - r)}); ++i) {
            // REMENBER THAT 1-DELTA(I,-1)=0 CAUSE ZERO VALUE
            if (i != -1) {
              b_sum +=
                  pow(2, 2 * r) / pow(4, p + q + 2) * pow(8.0, i) *
                  fact(p + q - 2 * i) / (fact(p - i) * fact(q - i)) *
                  pow(-1, r + i) * (1.0 - delta(i, -1)) /
                  (fact(r - i) * (fact(p + q + 2 - i - r))) * fact(r + 1) /
                  fact(2 * r + 2) * fact(2 * (p + q + 3 - i)) /
                  fact(p + q + 3 - i) * pow(-1.0, l) /
                  (fact(l) * fact(i + 2 - l)) *
                  ((i + 1 - l) * (i + 2 - l) *
                       ((p + q + 1 - i - r) * (p + q + 2 - i - r) -
                        0.5 * ((r - i) * (r - i - 1))) +
                   1.5 * ((l - 1) * l * (r - i) * (r - i - 1)) -
                   2.0 * (l * (i + 2 - l) * (r - i) * (p + q + 2 - i - r)));
            }
          }
        }

        double b = b_sum * pow(m1, p + 1) * pow(m0, q + 1);

        s_sum += b * Omega(l, r);
      }
    }
    return s_sum * 16 / 3;
  }

  double L1(int p, int q, size_t i0, [[maybe_unused]] size_t i1,
            const Eigen::MatrixXd &Omega) {
    double m0 = mass(i0) / mtot_;
    double m1 = mass(1 - i0) / mtot_; // if 0, get 1, and vice versa.
    double G = (m0 - m1) / m1;
    double F = (m0 * m0 + m1 * m1) / (2 * m0 * m1);
    double s_sum = 0.0;
    for (int l = 1; l <= std::min(p, q) + 2; ++l) {
      for (int r = l; r <= p + q + 4 - l; ++r) {

        double b = 0.0;
        double b_sum1 = 0.0;
        // REMENBER THAT 1-DELTA=0 CAUSE ZERO VALUE
        if (r != p + q + 3) {
          for (int i = l - 2; i <= std::min({p, q, r, (p + q + 2 - r)}); ++i) {
            // REMENBER THAT 1-DELTA(I,-1)=0 CAUSE ZERO VALUE
            if (i != -1) {
              b_sum1 =
                  pow(2.0, 2 * r) / pow(4.0, p + q + 2) * pow(8.0, i) *
                  fact(p + q - 2 * i) / (fact(p - i) * fact(q - i)) *
                  pow(-1.0, r + i) * fact(r + 1) /
                  (fact(r - i) * fact(p + q + 2 - i - r) * fact(2 * r + 2)) *
                  fact(2 * (p + q + 3 - i)) / fact(p + q + 3 - i);
              double b_sum2 = 0.0;
              for (int w = 0; w <= std::min({p, q, p + q + 2 - r}) - i; ++w) {
                b_sum2 +=
                    poch(p + 1 - i - w, w) * poch(q + 1 - i - w, w) /
                    (fact(w) * poch(p + q + 1 - 2 * i - w, w)) *
                    poch(p + q + 3 - i - r - w, w) /
                    poch(2 * (p + q + 3 - i) - 2 * w + 1, w) *
                    pow(2.0, 2 * w - 2) * pow(G, w) *
                    poch(p + q + 4 - i - w, w) /
                    poch(2 * (p + q + 3 - i) - w + 1, w) * pow(m0, i) *
                    pow(m1, i) * pow(m1, p + q - 2 * i - w) *
                    pow(F, i + 2 - l) * m0 * m0 / (fact(l) * fact(i + 2 - l)) *
                    4.0 *
                    (1.5 * m1 * m1 / (m0 * m0) *
                         (l * (l - 1) * (r - i) * (r - i - 1)) -
                     2.0 / F * m1 / m0 *
                         (l * (i + 2 - l) * (r - i) * (p + q + 2 - i - r - w)) +
                     1.0 / F * F * ((i + 1 - l) * (i + 2 - l)) *
                         ((p + q + 1 - i - r - w) * (p + q + 2 - i - r - w) -
                          0.5 * m1 * m1 / (m0 * m0) * ((r - i) * (r - i - 1))));
              }
              b += b_sum1 * b_sum2;
            }
          }
        }

        s_sum += b * Omega(l, r);
      }
    }
    return s_sum * 16 / 3;
  }

  double LSG(int p, int q, size_t i0, size_t i1, const Eigen::MatrixXd &Omega) {
    double s_sum = 0.0;
    for (int l = 2; l <= std::min(p, q) + 2; ++l) {
      for (int r = l; r <= p + q + 4 - l; ++r) {
        double b_sum = 0.0;
        for (int i = l - 2; i <= std::min({p, q, r, (p + q + 2 - r)}); ++i) {
          if (i != -1) {
            b_sum += pow(-1.0, r + i) / fact(p + q + 2 - i - r) *
                     fact(p + q - 2 * i) /
                     (fact(p - i) * fact(q - i) * fact(r - i)) *
                     fact(2 * (p + q + 3 - i)) / fact(p + q + 3 - i) *
                     pow(8.0, i) / fact(i + 2 - l) *
                     (((i + 1 - l) * (i + 2 - l)) *
                          (((p + q + 1 - i - r) * (p + q + 2 - i - r)) -
                           0.5 * ((r - i) * (r - i - 1))) +
                      1.5 * ((l - 1) * l * (r - i) * (r - i - 1)) -
                      (2 * l * (i + 2 - l) * (r - i) * (p + q + 2 - i - r)));
          }
        }

        double b = b_sum * pow(0.5, p + q + 2) * pow(2.0, 2 * r) /
                   pow(4.0, p + q + 2) * fact(r + 1) / fact(2 * r + 2) *
                   (1.0 + pow(-1.0, l)) / fact(l);
        s_sum += b * Omega(l, r);
      }
    }
    return s_sum * 16 / 3;
  }

  double B(int p, int q) {
    int ap = abs(p);
    int aq = abs(q);
    if (p > 0 && q > 0) {
      return molefraction0_ * molefraction0_ *
                 Lint(ap - 1, aq - 1, 0, 0, 0, 0) +
             molefraction0_ * molefraction1_ * Lint(ap - 1, aq - 1, 0, 0, 0, 1);
    } else if (p > 0 && q < 0) {
      return molefraction0_ * molefraction1_ * Lint(ap - 1, aq - 1, 0, 1, 0, 1);
    } else if (p < 0 && q > 0) {
      return molefraction0_ * molefraction1_ * Lint(ap - 1, aq - 1, 1, 0, 1, 0);
    } else // if (p<0 && q<0) {
      return molefraction1_ * molefraction1_ *
                 Lint(ap - 1, aq - 1, 1, 1, 1, 1) +
             molefraction0_ * molefraction1_ * Lint(ap - 1, aq - 1, 1, 1, 1, 0);
  }
};

std::tuple<std::vector<double> /*D12*/, std::vector<double> /*DT*/,
           std::vector<double> /*lambda*/, std::vector<double> /*eta*/>
transport(double t, double x0, std::vector<double> Omega00,
          std::vector<double> Omega01, std::vector<double> Omega11,
          double mass0, double mass1, int propertyorder) {
  // Some dirty work: turn C++ omega to fortran 2D array
  // The following arrays should be "FORTRAN-ready"

  int maxord = 2 * propertyorder + 3;

  Eigen::MatrixXd om00(maxord + 1, maxord + 1);
  Eigen::MatrixXd om01(maxord + 1, maxord + 1);
  Eigen::MatrixXd om11(maxord + 1, maxord + 1);

  {
    size_t i = 0;
    for (int l = 0; l != maxord; ++l) {
      for (int s = l; s != maxord; ++s) {
        om00(l + 1, s + 1) = Omega00[i];
        om01(l + 1, s + 1) = Omega01[i];
        om11(l + 1, s + 1) = Omega11[i];
        ++i;
      }
    }
  }

  std::vector<double> D12;
  std::vector<double> DT;
  std::vector<double> lambda;
  std::vector<double> eta;

  AlphaImpl alpha(propertyorder, mass0, mass1);
  BetaImpl beta(propertyorder, mass0, mass1);
  alpha.evaluate(t, om00, om01, om11, x0);
  beta.evaluate(t, om00, om01, om11, x0);
  D12 = alpha.D12();
  DT = alpha.DT();
  lambda = alpha.lambda();
  eta = beta.eta();

  return std::make_tuple(D12, DT, lambda, eta);
}

} // namespace dlt
