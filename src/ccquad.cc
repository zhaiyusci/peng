#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

int main() {
  std::cout << std::setprecision(8);
  // 0 - pi
  size_t maxorder = 9;
  // init
  if (0) {
    auto t0 = std::chrono::high_resolution_clock::now();
    std::vector<double> sins;
    std::vector<double> coss;
    sins.push_back(0.0);
    sins.push_back(0.0);
    coss.push_back(1.0);
    coss.push_back(-1.0);
    for (size_t order = 0; order != maxorder; ++order) {
      size_t maxnum = pow(2, order + 1);
      for (size_t i = 0; i != maxnum; i += 2) {
        double c = coss[i] * coss[i + 1] - sins[i] * sins[i + 1];
        double s = sqrt((1 - c) / 2);
        c = sqrt((1 + c) / 2) * (coss[i + 1] >= 0 ? 1 : -1);
        sins.insert(std::next(sins.begin(), i + 1), s);
        coss.insert(std::next(coss.begin(), i + 1), c);
        // sins.push_back(s);
        // coss.push_back(c);
      }
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    auto dt = t1 - t0;
    std::cout << "Algo1: " << dt.count() << std::endl;
  }

  {
    auto t0 = std::chrono::high_resolution_clock::now();
    size_t n = pow(2, maxorder) + 1;
    std::vector<double> sins(n);
    std::vector<double> coss(n);

    for (size_t i = 0; i != n; ++i) {
      sins[i] = (sin(M_PI / n * i));
      coss[i] = (cos(M_PI / n * i));
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    auto dt = t1 - t0;
    std::cout << "Algo2: " << dt.count() << std::endl;
  }

  {
    auto t0 = std::chrono::high_resolution_clock::now();
    auto dt = t0 - t0;
    // init
    std::vector<double> sins;
    std::vector<double> coss;
    size_t n = pow(2, maxorder) + 1;
    sins.reserve(n);
    coss.reserve(n);
    sins.push_back(0.0);
    sins.push_back(0.0);
    coss.push_back(1.0);
    coss.push_back(-1.0);
    for (size_t order = 0; order != maxorder; ++order) {
      size_t maxnum = pow(2, order);
      double lc = coss[0];
      double ls = sins[0];
      for (size_t i = 0; i != maxnum; ++i) {

        size_t r = i + 1;
        size_t b = (maxnum >> 1);
        while (!(r & 1)) {
          r >>= 1;
          b >>= 1;
        }
        r >>= 1;
        r += b + 1;
        double rc = coss[r];
        double rs = sins[r];
        double c = lc * rc - ls * rs;
        double s = sqrt((1 - c) / 2);
        c = sqrt((1 + c) / 2) * (rc >= 0 ? 1 : -1);
        sins.push_back(s);
        coss.push_back(c);
        lc = rc;
        ls = rs;
      }
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    dt += t1 - t0;
    std::cout << "Algo3: " << dt.count() << std::endl;
  }

  std::vector<double> sins;
  std::vector<double> coss;
  {
    // Only consider [0,pi/2]
    auto t0 = std::chrono::high_resolution_clock::now();
    auto dt = t0 - t0;
    // init
    size_t n = pow(2, maxorder - 1) + 1;
    sins.reserve(n);
    coss.reserve(n);
    sins.push_back(0.0);
    sins.push_back(1.0);
    coss.push_back(1.0);
    coss.push_back(0.0);
    for (size_t order = 1; order != maxorder; ++order) {
      size_t maxnum = pow(2, order - 1);
      double lc = coss[0];
      double ls = sins[0];
      for (size_t i = 0; i != maxnum; ++i) {
        size_t r = i + 1;
        size_t b = (maxnum >> 1);
        while (!(r & 1)) {
          r >>= 1;
          b >>= 1;
        }
        r >>= 1;
        r += b + 1;
        // std::cerr << r << std::endl;
        double rc = coss[r];
        double rs = sins[r];
        double c = lc * rc - ls * rs;
        double s = sqrt((1 - c) / 2);
        c = sqrt((1 + c) / 2);
        sins.push_back(s);
        coss.push_back(c);
        lc = rc;
        ls = rs;
      }
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    dt += t1 - t0;
    std::cout << "Algo4: " << dt.count() << std::endl;
  }

  std::cout << coss.size() << std::endl;
  // std::vector<int> ord{0, 8, 4, 2, 6, 1, 3, 5, 7}; // For test 0-8
  std::vector<int> ord{0, 4, 2, 1, 3}; // For test 0-8
  for (int i = 0; i <= pow(2, maxorder - 1); ++i) {
    std::cout << coss[i] << "  " << cos(M_PI / pow(2, maxorder) * ord[i])
              << "  " << sins[i] << "  "
              << sin(M_PI / pow(2, maxorder) * ord[i]) << std::endl;
  }

  return 0;
}
