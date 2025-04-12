/**
 * @file inputDeck.cpp
 * @brief input deck reader library. Going to use YAML-cpp to parse the input 
 *        file.
 * @author Caleb Hash
 * @date 2025-02-22
 */

#include <input.h>
#include <iostream>
#include <yaml-cpp/yaml.h>

YAML::Node config;

namespace IO {

ConfigData ConfigData::fromYAMLConfig(const YAML::Node& iconfig) {
  ConfigData ideck;

  // ... domain settings
  ideck.mesh.lx = iconfig["domain"]["lengths"]["x"].as<double>();
  ideck.mesh.ly = iconfig["domain"]["lengths"]["y"].as<double>();
  ideck.mesh.nx = iconfig["domain"]["dimensions"]["x"].as<int>();
  ideck.mesh.ny = iconfig["domain"]["dimensions"]["y"].as<int>();

  // ... solver settings
  ideck.restart.load = iconfig["solver"]["restart"]["load"].as<bool>();
  ideck.restart.save = iconfig["solver"]["restart"]["save"].as<bool>();
  ideck.solver.iter    = iconfig["solver"]["iterations"].as<int>();
  ideck.solver.omega = iconfig["solver"]["pressure solver"]["omega"].as<double>(); 
  ideck.solver.toler = iconfig["convergence"]["residual"].as<double>();
  ideck.solver.cfli  = iconfig["dynamic CFL"]["cfli"].as<double>();
  ideck.solver.cflf  = iconfig["dynamic CFL"]["cflf"].as<double>();

  // ... IO parameters
  ideck.res.enabled = iconfig["output"]["residuals"]["enabled"].as<bool>();
  ideck.res.freq = iconfig["output"]["residuals"]["frequency"].as<int>();
  ideck.res.file = iconfig["output"]["residuals"]["file"].as<string>();
  ideck.fv.enabled  = iconfig["output"]["flowviz"]["enabled"].as<bool>();
  ideck.fv.freq  = iconfig["output"]["flowviz"]["frequency"].as<int>();
  ideck.fv.dir = iconfig["output"]["flowviz"]["directory"].as<string>();
  ideck.fv.ghost = iconfig["output"]["flowviz"]["ghost cells"].as<bool>();
  ideck.fv.mode = iconfig["output"]["flowviz"]["mode"].as<string>();

  // ... BC parameters
  ideck.bcLeft.type = iconfig["boundary conditions"]["left"]["type"].as<string>();
  if (ideck.bcLeft.type=="moving wall" || ideck.bcLeft.type=="inlet") {
    ideck.bcLeft.velocity[0] = iconfig["boundary conditions"]["left"]["velocity"]["u"].as<double>();
    ideck.bcLeft.velocity[1] = iconfig["boundary conditions"]["left"]["velocity"]["v"].as<double>();
  } else if (ideck.bcLeft.type=="periodic") {
    ideck.bcLeft.pressure = iconfig["boundary conditions"]["left"]["pressure"].as<double>();
  }
  ideck.bcRight.type = iconfig["boundary conditions"]["right"]["type"].as<string>();
  if (ideck.bcRight.type=="moving wall" || ideck.bcRight.type=="inlet") {
    ideck.bcRight.velocity[0] = iconfig["boundary conditions"]["right"]["velocity"]["u"].as<double>();
    ideck.bcRight.velocity[1] = iconfig["boundary conditions"]["right"]["velocity"]["v"].as<double>();
  } else if (ideck.bcRight.type=="periodic") {
    ideck.bcRight.pressure = iconfig["boundary conditions"]["right"]["pressure"].as<double>();
  }
  ideck.bcBottom.type = iconfig["boundary conditions"]["bottom"]["type"].as<string>();
  if (ideck.bcBottom.type=="moving wall" || ideck.bcBottom.type=="inlet") {
    ideck.bcBottom.velocity[0] = iconfig["boundary conditions"]["bottom"]["velocity"]["u"].as<double>();
    ideck.bcBottom.velocity[1] = iconfig["boundary conditions"]["bottom"]["velocity"]["v"].as<double>();
  } else if (ideck.bcBottom.type=="periodic") {
    ideck.bcBottom.pressure = iconfig["boundary conditions"]["bottom"]["pressure"].as<double>();
  }
  ideck.bcTop.type = iconfig["boundary conditions"]["top"]["type"].as<string>();
  if (ideck.bcTop.type=="moving wall" || ideck.bcTop.type=="inlet") {
    ideck.bcTop.velocity[0] = iconfig["boundary conditions"]["top"]["velocity"]["u"].as<double>();
    ideck.bcTop.velocity[1] = iconfig["boundary conditions"]["top"]["velocity"]["v"].as<double>();
  } else if (ideck.bcTop.type=="periodic") {
    ideck.bcTop.pressure = iconfig["boundary conditions"]["top"]["pressure"].as<double>();
  }

  // ... initial conditions
  ideck.igas.u = iconfig["initial conditions"]["gas"]["u"].as<double>();
  ideck.igas.v = iconfig["initial conditions"]["gas"]["v"].as<double>();
  ideck.igas.p = iconfig["initial conditions"]["gas"]["p"].as<double>();
  ideck.igas.rho = iconfig["initial conditions"]["gas"]["rho"].as<double>();
  ideck.igas.mu = iconfig["initial conditions"]["gas"]["mu"].as<double>();
  ideck.iliq.u = iconfig["initial conditions"]["liquid"]["u"].as<double>();
  ideck.iliq.v = iconfig["initial conditions"]["liquid"]["v"].as<double>();
  ideck.iliq.p = iconfig["initial conditions"]["liquid"]["p"].as<double>();
  ideck.iliq.rho = iconfig["initial conditions"]["liquid"]["rho"].as<double>();
  ideck.iliq.mu = iconfig["initial conditions"]["liquid"]["mu"].as<double>();

  // droplets
  ideck.drop.enabled = iconfig["droplet"]["enabled"].as<bool>();
  ideck.drop.r = iconfig["droplet"]["diameter"].as<double>() *0.5;
  ideck.drop.M = iconfig["droplet"]["M"].as<double>();
  ideck.drop.x = iconfig["droplet"]["center"]["x"].as<double>();
  ideck.drop.y = iconfig["droplet"]["center"]["y"].as<double>();

  // levelset options
  ideck.levset.reinit = iconfig["level set"]["reinitialization"]["enabled"].as<bool>();
  ideck.levset.ireinit = iconfig["level set"]["reinitialization"]["iterations"].as<int>();
  return ideck;
}

ConfigData parseInputDeck(const string& inFile){

  // ... create struct for the yaml node
  try {
    YAML::Node iconfig = YAML::LoadFile(inFile);
    return ConfigData::fromYAMLConfig(iconfig);
  } catch (const std::exception& e) {
    std::cerr << "Error loading YAML input deck: " << e.what() << std::endl;
    throw;
  }

} // end function parseInputDeck

} // end namespace IO
