#define EIGEN_NO_DEBUG
#include "transport.hh"
#include "atompair.hh"
#include "mathtools.hh"
#include <algorithm>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <map>
#include <tuple>
#include <vector>

static const double kB = 1.380649e-23;      // BY DEFINITION
static const double amu = 1.6605390666e-27; // CODATA2018

namespace peng {

OmegaCache::OmegaCache(AtomPair &atompair, double rtol)
    : p_atompair_(&atompair), rtol_(rtol) {}
double OmegaCache::operator()(size_t l, size_t s, double temperature) const {
  // Store the result after computation, and read the cache_ for a duplicated
  // call.
  if (cache_.count(std::make_tuple(l, s, temperature)) == 0) {
    double omega = p_atompair_->Omega(l, s, temperature, rtol_);
    cache_[std::make_tuple(l, s, temperature)] = omega;
    return omega;
  } else {
    return cache_[std::make_tuple(l, s, temperature)];
  }
}

TransportProperties::TransportProperties(double mass0, double mass1,
                                         OmegaComp &Omega00, OmegaComp &Omega01,
                                         OmegaComp &Omega11)
    : mass0_(mass0), mass1_(mass1), p_Omega00_(&Omega00), p_Omega01_(&Omega01),
      p_Omega11_(&Omega11), temperature_(-1.0), propertyorder_(0), Hint(),
      Lint() {
  mtot_ = mass0 + mass1;
  m0_ = mass0 / mtot_;
  m1_ = mass1 / mtot_;
}

void TransportProperties::update_temperature_(double temperature) {
  if (temperature != temperature_) {
    temperature_ = temperature;
    // Destroy all cached result.
    propertyorder_ = 0;
    // For order == 0, we need H and L (0,0) element prepared
    Hint.resize(1, 1);
    Lint.resize(1, 1);
    for (size_t i0 = 0; i0 != 2; ++i0) {
      for (size_t i1 = 0; i1 != 2; ++i1) {
        for (size_t io0 = 0; io0 != 2; ++io0) {
          for (size_t io1 = 0; io1 != 2; ++io1) {
            Hint(0, 0)(i0, i1, io0, io1) =
                bracket_int_H_(0, 0, i0, i1, io0, io1);
            Lint(0, 0)(i0, i1, io0, io1) =
                bracket_int_L_(0, 0, i0, i1, io0, io1);
          }
        }
      }
    }
  }
  return;
}

void TransportProperties::update_Hint_Lint_(size_t propertyorder) {
  Hint.resize(propertyorder + 1, propertyorder + 1);
  for (Ext2D<Sexdec>::iterator h(&Hint, 0, propertyorder_ + 1); h < Hint.end();
       ++h) {
    std::size_t p = h.row();
    std::size_t q = h.col();

    for (size_t i0 = 0; i0 != 2; ++i0) {
      for (size_t i1 = 0; i1 != 2; ++i1) {
        for (size_t io0 = 0; io0 != 2; ++io0) {
          for (size_t io1 = 0; io1 != 2; ++io1) {
            (*h)(i0, i1, io0, io1) = bracket_int_H_(p, q, i0, i1, io0, io1);
          }
        }
      }
    }
  }

  Lint.resize(propertyorder + 1, propertyorder + 1);
  for (Ext2D<Sexdec>::iterator l(&Lint, 0, propertyorder_ + 1); l < Lint.end();
       ++l) {
    std::size_t p = l.row();
    std::size_t q = l.col();

    for (size_t i0 = 0; i0 != 2; ++i0) {
      for (size_t i1 = 0; i1 != 2; ++i1) {
        for (size_t io0 = 0; io0 != 2; ++io0) {
          for (size_t io1 = 0; io1 != 2; ++io1) {
            (*l)(i0, i1, io0, io1) = bracket_int_L_(p, q, i0, i1, io0, io1);
          }
        }
      }
    }
  }

  propertyorder_ = propertyorder;
}

Eigen::MatrixXd TransportProperties::compute_raw_D_mat_() {
  Eigen::MatrixXd D_mat(2 * propertyorder_ + 1, 2 * propertyorder_ + 1);
  for (int p = 1; p <= (int)propertyorder_; ++p) {
    for (int q = 1; q <= (int)propertyorder_; ++q) {
      // clang-format off
        D_mat(propertyorder_ + p, propertyorder_ + q) = A_( p,  q); // SE
        D_mat(propertyorder_ + p, propertyorder_ - q) = A_( p, -q); // SW
        D_mat(propertyorder_ - p, propertyorder_ - q) = A_(-p, -q); // NW
        D_mat(propertyorder_ - p, propertyorder_ + q) = A_(-p,  q); // NE
      // clang-format on
    }
  }
  for (int p = 1; p <= (int)propertyorder_; ++p) {
    // clang-format off
      D_mat(propertyorder_ + 0, propertyorder_ + p) = A_( 0,  p); // E
      D_mat(propertyorder_ + p, propertyorder_ + 0) = A_( p, -0); // S
      D_mat(propertyorder_ + 0, propertyorder_ - p) = A_(-0, -p); // W
      D_mat(propertyorder_ - p, propertyorder_ + 0) = A_(-p,  0); // N
    // clang-format on
  }
  D_mat(propertyorder_, propertyorder_) = A_(0, 0); // C

  return D_mat;
}

Eigen::MatrixXd TransportProperties::compute_raw_B_mat() {
  Eigen::MatrixXd B_mat(2 * propertyorder_, 2 * propertyorder_);
  for (int p = 1; p <= (int)propertyorder_; ++p) {
    for (int q = 1; q <= (int)propertyorder_; ++q) {
      // clang-format off
        B_mat(propertyorder_ - 1 + p, propertyorder_ - 1 + q) = B_( p,  q); // SE
        B_mat(propertyorder_ - 1 + p, propertyorder_     - q) = B_( p, -q); // SW
        B_mat(propertyorder_     - p, propertyorder_     - q) = B_(-p, -q); // NW
        B_mat(propertyorder_     - p, propertyorder_ - 1 + q) = B_(-p,  q); // NE
      // clang-format on
    }
  }

  return B_mat;
}

std::tuple<double, double, double, double>
TransportProperties::evaluate(double temperature, double molefraction0,
                              size_t po) {
  update_temperature_(temperature);
  update_Hint_Lint_(po);

  molefraction0_ = molefraction0;
  molefraction1_ = 1 - molefraction0;

  // raw_D_mat
  auto &&raw_D_mat = compute_raw_D_mat_();
  std::cerr << "raw_D_mat:" << '\n';
  std::cerr << raw_D_mat << std::endl;

  auto &&raw_B_mat = compute_raw_B_mat();
  std::cerr << "raw_B_mat:" << '\n';
  std::cerr << raw_B_mat << std::endl;

  Eigen::MatrixXd A_mat(2 * po, 2 * po);
  // clang-format off
    A_mat.block(po, po, po, po) = raw_D_mat.block(propertyorder_ +  1, propertyorder_ +  1, po, po); // SE
    A_mat.block(po,  0, po, po) = raw_D_mat.block(propertyorder_ +  1, propertyorder_ - po, po, po); // SW
    A_mat.block( 0,  0, po, po) = raw_D_mat.block(propertyorder_ - po, propertyorder_ - po, po, po); // NW
    A_mat.block( 0, po, po, po) = raw_D_mat.block(propertyorder_ - po, propertyorder_ +  1, po, po); // NE
  // clang-format on
  Eigen::VectorXd alpha_vec = Eigen::VectorXd::Zero(2 * po);
  // clang-format off
    alpha_vec(po - 1) = -15.0 / 4.0 * molefraction1_ * sqrt(2 * kB * temperature / mass1_ / amu); // ALPHA_[-1]
    alpha_vec(po    ) = -15.0 / 4.0 * molefraction0_ * sqrt(2 * kB * temperature / mass0_ / amu); // ALPHA_[ 1]
  // clang-format on
  auto a_vec = A_mat.inverse() * alpha_vec;
  double lambda = -5.0 / 4.0 * kB * sqrt(2 * kB * temperature / mtot_ / amu) *
                  (molefraction0_ / sqrt(m0_) * a_vec(po) +
                   molefraction1_ / sqrt(m1_) * a_vec(po - 1));
  Eigen::MatrixXd D_mat(2 * po + 1, 2 * po + 1);
  D_mat = raw_D_mat.block(propertyorder_ - po, propertyorder_ - po, 2 * po + 1,
                          2 * po + 1);
  Eigen::VectorXd delta_vec = Eigen::VectorXd::Zero(2 * po + 1);
  delta_vec(po) = 3.0 / 2.0 * sqrt(2 * kB * temperature / mtot_ / amu);
  Eigen::VectorXd d_vec = D_mat.inverse() * delta_vec;
  double D12 = (kB * temperature) / 1.013e5 * 1.0 / 2.0 * molefraction0_ *
               molefraction1_ * sqrt(2 * kB * temperature / mtot_ / amu) *
               d_vec(po);
  double DT = (kB * temperature) / 1.013e5 * (-5.0) / 4.0 * molefraction0_ *
              molefraction1_ * sqrt(2 * kB * temperature / mtot_ / amu) *
              (molefraction0_ / sqrt(m0_) * d_vec(po + 1) +
               molefraction1_ / sqrt(m1_) * d_vec(po - 1));
  Eigen::MatrixXd B_mat(2 * po, 2 * po);
  B_mat.block(0, 0, 2 * po, 2 * po) = raw_B_mat.block(
      propertyorder_ - po, propertyorder_ - po, 2 * po, 2 * po); // NW
  Eigen::VectorXd beta_vec = Eigen::VectorXd::Zero(2 * po);
  // clang-format off
    beta_vec(po - 1) = 5.0 / 2.0 * molefraction1_; // ALPHA_[-1]
    beta_vec(po    ) = 5.0 / 2.0 * molefraction0_; // ALPHA_[ 1]
  // clang-format on
  auto b_vec = B_mat.inverse() * beta_vec;
  double eta = (molefraction0_ * b_vec(po) + molefraction1_ * b_vec(po - 1)) *
               temperature * kB;
  return std::make_tuple(D12, DT, lambda, eta);
}

double TransportProperties::bracket_int_H_(int p, int q, size_t i0, size_t i1,
                                           size_t oi0, size_t oi1) {

  const OmegaComp &Omega = this->Omega_(oi0, oi1);

  if (i0 != i1) {
    return H12_(p, q, i0, i1, Omega);
  } else {
    if (oi0 != oi1) {
      return H1_(p, q, i0, i1, Omega);
    } else {
      return HSG_(p, q, i0, i1, Omega);
    }
  }
}

double TransportProperties::H12_(int p, int q, size_t i0, size_t i1,
                                 const OmegaComp &Omega) { // eq 109
  double s_sum = 0.0;
  double m0 = mass_(i0) / mtot_;
  double m1 = mass_(i1) / mtot_;
  for (int l = 1; l <= std::min(p, q) + 1; ++l) {
    for (int r = l; r <= p + q + 2 - l; ++r) {
      double a_sum = 0.0;
      // remenber that 1-delta=0 cause zero value
      for (int i = l - 1; i <= std::min({p, q, r, p + q + 1 - r});
           ++i) { // eq 110
        a_sum += pow(8, i) * fact(p + q - 2 * i) / (fact(p - i) * fact(q - i)) *
                 pow(-1, l) / (fact(l) * fact(i + 1 - l)) * pow(-1, r + i) /
                 (fact(r - i) * (fact(p + q + 1 - i - r))) * fact(r + 1) /
                 fact(2 * r + 2) * fact(2 * (p + q + 2 - i)) /
                 fact(p + q + 2 - i) * pow(2, 2 * r) / pow(4.0, p + q + 1) *
                 ((i + 1 - l) * (p + q + 1 - i - r) - (l * (r - i)));
      }

      double a = a_sum * pow(m1, p + 0.5) * pow(m0, q + 0.5);

      s_sum += a * Omega(l, r, temperature_);
    }
  }
  return s_sum * 8;
}

double TransportProperties::H1_(int p, int q, size_t i0,
                                [[maybe_unused]] size_t i1,
                                const OmegaComp &Omega) {
  double m0 = mass_(i0) / mtot_;
  double m1 = mass_(1 - i0) / mtot_; // if 0, get 1, and vice versa.
  double G = (m0 - m1) / m1;
  double F = (m0 * m0 + m1 * m1) / (2 * m0 * m1);
  double s_sum = 0.0;
  for (int l = 1; l <= std::min(p, q) + 1; ++l) {
    for (int r = l; r <= p + q + 2 - l; ++r) {

      double a = 0.0;
      double a_sum1 = 0.0;
      for (int i = l - 1; i <= std::min({p, q, r, (p + q + 1 - r)}); ++i) {
        a_sum1 = pow(8, i) * fact(p + q - 2 * i) / (fact(p - i) * fact(q - i)) *
                 1.0 / (fact(l) * fact(i + 1 - l)) * pow(-1, r + i) /
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

      s_sum += a * Omega(l, r, temperature_);
    }
  }
  return s_sum * 8;
}

double TransportProperties::HSG_(int p, int q, size_t i0, size_t i1,
                                 const OmegaComp &Omega) {
  double s_sum = 0.0;
  for (int l = 1; l <= std::min(p, q) + 1; ++l) {
    for (int r = l; r <= p + q + 2 - l; ++r) {
      double a_sum = 0.0;
      for (int i = l - 1; i <= std::min({p, q, r, (p + q + 1 - r)}); ++i) {
        a_sum += pow(8, i) * fact(p + q - 2 * i) / (fact(p - i) * fact(q - i)) *
                 (1 + pow(-1.0, l)) / (fact(l) * fact(i + 1 - l)) *
                 pow(-1, r + i) / (fact(r - i) * fact(p + q + 1 - i - r)) *
                 fact(r + 1) / fact(2 * r + 2) * fact(2 * (p + q + 2 - i)) /
                 fact(p + q + 2 - i) * pow(2.0, 2 * r) / pow(4.0, p + q + 1) *
                 ((i + 1 - l) * (p + q + 1 - i - r) - l * (r - i));
      }

      double a = a_sum / pow(2, p + q + 1);

      s_sum += a * Omega(l, r, temperature_);
    }
  }
  return s_sum * 8;
}

double TransportProperties::A_(int p, int q) {
  int ap = abs(p);
  int aq = abs(q);
  if (p > 0 && q > 0) {
    return molefraction0_ * molefraction0_ * Hint(ap, aq)(0, 0, 0, 0) +
           molefraction0_ * molefraction1_ * Hint(ap, aq)(0, 0, 0, 1);
  } else if (p > 0 && q < 0) {
    return molefraction0_ * molefraction1_ * Hint(ap, aq)(0, 1, 0, 1);
  } else if (p < 0 && q > 0) {
    return molefraction0_ * molefraction1_ * Hint(ap, aq)(1, 0, 1, 0);
  } else if (p < 0 && q < 0) {
    return molefraction1_ * molefraction1_ * Hint(ap, aq)(1, 1, 1, 1) +
           molefraction0_ * molefraction1_ * Hint(ap, aq)(1, 1, 1, 0);
  } else if (p > 0 && q == 0) {
    return molefraction0_ * molefraction1_ * sqrt(m0_) *
           Hint(ap, 0)(0, 0, 0, 1);
  } else if (p == 0 && q > 0) {
    return molefraction0_ * molefraction1_ * sqrt(m0_) *
           Hint(aq, 0)(0, 0, 0, 1);
  } else if (p < 0 && q == 0) {
    return -molefraction0_ * molefraction1_ * sqrt(m1_) *
           Hint(ap, 0)(1, 1, 1, 0);
  } else if (p == 0 && q < 0) {
    return -molefraction0_ * molefraction1_ * sqrt(m1_) *
           Hint(aq, 0)(1, 1, 1, 0);
  } else { // if (p == 0 &&q == 0) {
    return molefraction0_ * molefraction1_ *
           (8 * m0_ * m1_ * (*p_Omega01_)(1, 1, temperature_));
  }
}

double TransportProperties::bracket_int_L_(int p, int q, size_t i0, size_t i1,
                                           size_t oi0, size_t oi1) {

  const OmegaComp &Omega = this->Omega_(oi0, oi1);

  if (i0 != i1) {
    return L12_(p, q, i0, i1, Omega);
  } else {
    if (oi0 != oi1) {
      return L1_(p, q, i0, i1, Omega);
    } else {
      return LSG_(p, q, i0, i1, Omega);
    }
  }
}

double TransportProperties::L12_(int p, int q, size_t i0, size_t i1,
                                 const OmegaComp &Omega) { // eq 109
  double s_sum = 0.0;
  double m0 = mass_(i0) / mtot_;
  double m1 = mass_(i1) / mtot_;
  for (int l = 1; l <= std::min(p, q) + 2; ++l) {
    for (int r = l; r <= p + q + 4 - l; ++r) {
      double b_sum = 0.0;
      // REMENBER THAT 1-DELTA=0 CAUSE ZERO VALUE
      if (r != p + q + 3) {
        for (int i = l - 2; i <= std::min({p, q, r, (p + q + 2 - r)}); ++i) {
          // REMENBER THAT 1-DELTA(I,-1)=0 CAUSE ZERO VALUE
          if (i != -1) {
            b_sum += pow(2, 2 * r) / pow(4, p + q + 2) * pow(8.0, i) *
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

      s_sum += b * Omega(l, r, temperature_);
    }
  }
  return s_sum * 16 / 3;
}

double TransportProperties::L1_(int p, int q, size_t i0,
                                [[maybe_unused]] size_t i1,
                                const OmegaComp &Omega) {
  double m0 = mass_(i0) / mtot_;
  double m1 = mass_(1 - i0) / mtot_; // if 0, get 1, and vice versa.
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
            b_sum1 = pow(2.0, 2 * r) / pow(4.0, p + q + 2) * pow(8.0, i) *
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
                  pow(2.0, 2 * w - 2) * pow(G, w) * poch(p + q + 4 - i - w, w) /
                  poch(2 * (p + q + 3 - i) - w + 1, w) * pow(m0, i) *
                  pow(m1, i) * pow(m1, p + q - 2 * i - w) * pow(F, i + 2 - l) *
                  m0 * m0 / (fact(l) * fact(i + 2 - l)) * 4.0 *
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

      s_sum += b * Omega(l, r, temperature_);
    }
  }
  return s_sum * 16 / 3;
}

double TransportProperties::LSG_(int p, int q, size_t i0, size_t i1,
                                 const OmegaComp &Omega) {
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
      s_sum += b * Omega(l, r, temperature_);
    }
  }
  return s_sum * 16 / 3;
}

double TransportProperties::B_(int p, int q) {
  int ap = abs(p);
  int aq = abs(q);
  if (p > 0 && q > 0) {
    return molefraction0_ * molefraction0_ * Lint(ap - 1, aq - 1)(0, 0, 0, 0) +
           molefraction0_ * molefraction1_ * Lint(ap - 1, aq - 1)(0, 0, 0, 1);
  } else if (p > 0 && q < 0) {
    return molefraction0_ * molefraction1_ * Lint(ap - 1, aq - 1)(0, 1, 0, 1);
  } else if (p < 0 && q > 0) {
    return molefraction0_ * molefraction1_ * Lint(ap - 1, aq - 1)(1, 0, 1, 0);
  } else // if (p<0 && q<0) {
    return molefraction1_ * molefraction1_ * Lint(ap - 1, aq - 1)(1, 1, 1, 1) +
           molefraction0_ * molefraction1_ * Lint(ap - 1, aq - 1)(1, 1, 1, 0);
}

} // namespace peng
