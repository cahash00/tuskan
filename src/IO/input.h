#ifndef INPUT_H
#define INPUT_H

#include <yaml-cpp/yaml.h>

using namespace std;

namespace IO {

extern YAML::Node config;

struct levelset {
  bool reinit;
  int ireinit;
};
struct droplet {
  double r;
  double M;
  double x;
  double y;
};
struct fluid {
  double u;
  double v;
  double p;
  double rho;
  double mu;
};
struct boundary {
  std::vector<double> velocity = std::vector<double>(2,0.0);
  double pressure;
  string type;
};
struct Restart {
  bool save;
  bool load;
};
struct Solver {
  int iter;
  double toler;
  double cfli;
  double cflf;
  double omega;
};
struct Fv {
  bool enabled;
  int freq;
  string dir;
  bool ghost;
};
struct Mesh {
  double lx;
  double ly;
  int nx;
  int ny;
};
struct Res {
  bool enabled;
  int freq;
  string file;
};
struct ConfigData {
  // domain settings
  Mesh mesh;
  // solver settings
  Restart restart;
  Solver solver;
  // IO parameters
  Fv fv;
  Res res;
  // BCs
  boundary bcLeft;
  boundary bcRight;
  boundary bcBottom;
  boundary bcTop;
  // initial conditions
  fluid igas;
  fluid iliq;
  // droplet
  droplet drop;

  // level-set settings
  levelset levset;

  static ConfigData fromYAMLConfig(const YAML::Node& config);
};

ConfigData parseInputDeck(const string& inFile);

}

#endif // INPUT_H
