#ifndef INPUT_H
#define INPUT_H

#include <yaml-cpp/yaml.h>
#include <argparse/argparse.hpp>

using namespace std;

namespace IO_input {

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
  int    pIter;
  double sorWeight; 
  // convergence criteria
  double toler;
  double cfli;
  double cflf;
  bool   dcfl;
  // IO parameters
  string foutDir;


  static ConfigData fromYAMLConfig(const YAML::Node& config);
};

extern ConfigData parseInputDeck(const string& inFile);

void getUserInput(int argc, char* argv[],
                  argparse::ArgumentParser& program);

}

#endif // INPUT_H
