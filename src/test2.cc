#include "mathtools.hh"
#include <chrono>
#include <cmath>
#include <iostream>
#include <thread>

int main() {
  class : public dlt::FuncDeriv1D{
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

  std::cout << lj(1.0) << std::endl;
  std::cout << lj.derivative(1.0) << std::endl;
  std::cout << lj(2.0) << std::endl;
  std::cout << lj.derivative(2.0) << std::endl;

  double r, v;
  std::tie(r, v) = dlt::find_local_minimum(lj, 1.0, 3.0);
  std::cout << r << " " << v << std::endl;
  r = dlt::find_local_root(lj, 0.0, 0.5, 1.5);
  std::cout << r << std::endl;
  std::tie(r, v) = dlt::find_local_maximum(-lj, 1.0, 3.0);

  std::cout << r << std::endl;
  dlt::LocalRoot lr(lj, 0.5, 1.5);
  std::cout << lr(0.0) << std::endl;

  return 0;
}
