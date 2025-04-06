#ifndef INPUT_H
#define INPUT_H

#include <yaml-cpp/yaml.h>

using namespace std;

namespace IO {

extern YAML::Node config;

struct boundary {
  std::vector<double> velocity = std::vector<double>(2,0.0);
  double pressure;
  string type;
};

struct ConfigData {
  // domain settings
  double lx;
  double ly;
  int    nx;
  int    ny;
  // solver settings
  int    iter;
  double cfl;
  bool   fvflag;
  int    fvfreq;
  bool   resflag;
  int    resfreq;
  string resFile;
  string pMethod;
  double sorOmega; 
  // convergence criteria
  double toler;
  double cfli;
  double cflf;
  bool   dcfl;
  // IO parameters
  string foutDir;
  bool ghost;
  // BCs
  boundary bcLeft;
  boundary bcRight;
  boundary bcBottom;
  boundary bcTop;
  // initial conditions
  double uinit;
  double vinit;
  double pinit;

  static ConfigData fromYAMLConfig(const YAML::Node& config);
};

ConfigData parseInputDeck(const string& inFile);

}

#endif // INPUT_H
