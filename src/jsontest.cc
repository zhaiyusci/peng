#include <iostream>
#include <json.hpp>

int main() {
  nlohmann::json j;
  std::cin >> j;
  std::cout << j["temperatures"].get<std::vector<double>>()[0] << std::endl;
  return 0;
}
