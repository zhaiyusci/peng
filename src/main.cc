#include "dilute.hh"
#include "loadedpotential.hh"
#include <iostream>

using namespace dlt;
int main() {
  // clang-format off
  std::cout << "*=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=*\n" 
            << "{   ,gggggggggggg,                                            }\n"
            << "{  dP'''88''''''Y8b,       ,dPYb,              I8             }\n"
            << "{  Yb,  88       `8b,      IP'`Yb              I8             }\n"
            << "{   `'  88        `8b gg   I8  8I           88888888          }\n"
            << "{       88         Y8 ''   I8  8'              I8             }\n"
            << "{       88         d8 gg   I8 dP  gg      gg   I8    ,ggg,    }\n"
            << "{       88        ,8P 88   I8dP   I8      8I   I8   i8' '8i   }\n"
            << "{       88       ,8P' 88   I8P    I8,    ,8I  ,I8,  I8, ,8I   }\n"
            << "{       88______,dP'_,88,_,d8b,_ ,d8b,  ,d8b,,d88b, `YbadP'   }\n"
            << "{      888888888P'  8P''Y88P''Y888P''Y88P'`Y88P''Y8888P'Y888  }\n"
            << "*=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=~-=*\n"
            << "\n"
            << "Version: " << DILUTE_VERSION << '\n'
            << "\n"
            << "Authors: \n"
            << "  Yu Zhai, You Li, Hui Li, and Frederick R. W. McCourt\n"
            << "\n"
            << "GitHub Repo: \n"
            << "  https://github.com/zhaiyusci/dilute\n"
            << "\n"
            << "The authors appreciate the financial support from \n"
            << "  * National Natural Science Foundation of China\n"
            << "  * Jilin University, China\n" 
            << "    - The Program for JLU Computational Interdiscipline \n"
            << "      Innovative Platform\n"
            ;

  // clang-format on
  std::cout.flush();

  std::cout << std::setprecision(16);
  std::cerr << std::setprecision(16);
  nlohmann::json config;
  std::cin >> config;

  Atom atom0(config["atoms"][0]["name"].get<std::string>(),
             config["atoms"][0]["mass"].get<double>());
  Atom atom1(config["atoms"][1]["name"].get<std::string>(),
             config["atoms"][1]["mass"].get<double>());
  std::vector<double> temperatures(
      config["temperatures"].get<std::vector<double>>());
  std::vector<double> molefractions0(
      config["molefractions0"].get<std::vector<double>>());
  size_t propertyorder(config["propertyorder"].get<size_t>());
  double accuracy(config["accuracy"].get<double>());
  LoadedPotential pot00(config["potentials"][0]["path"].get<std::string>());
  LoadedPotential pot01(config["potentials"][1]["path"].get<std::string>());
  LoadedPotential pot11(config["potentials"][2]["path"].get<std::string>());

  dlt::Task task(atom0, atom1, temperatures, molefractions0, propertyorder,
                 accuracy, pot00, pot01, pot11);
  task.execute();
  task.chant();

  return 0;
}
