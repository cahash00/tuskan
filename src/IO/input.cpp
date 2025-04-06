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

  // domain settings
  ideck.lx = iconfig["domain"]["lengths"]["x"].as<double>();
  ideck.ly = iconfig["domain"]["lengths"]["y"].as<double>();
  ideck.nx = iconfig["domain"]["dimensions"]["x"].as<int>();
  ideck.ny = iconfig["domain"]["dimensions"]["y"].as<int>();
  // solver settings
  ideck.iter    = iconfig["solver"]["iterations"].as<int>();
  ideck.fvflag  = iconfig["output"]["flowviz"]["enabled"].as<bool>();
  ideck.fvfreq  = iconfig["output"]["flowviz"]["frequency"].as<int>();
  ideck.resflag = iconfig["output"]["residuals"]["enabled"].as<bool>();
  ideck.resfreq = iconfig["output"]["residuals"]["frequency"].as<int>();
  ideck.resFile = iconfig["output"]["residuals"]["file"].as<string>();
  ideck.pMethod = iconfig["solver"]["pressure solver"]["method"].as<string>();
  if (ideck.pMethod == "SOR") {
    if (iconfig["solver"]["pressure solver"]["omega"]) {
      ideck.sorOmega = iconfig["solver"]["pressure solver"]["omega"].as<double>(); 
    } else {
      throw runtime_error("ERROR: INPUT DECK: SOR solver requires \"omega\".");
    }
  }
  // convergence criteria
  ideck.toler = iconfig["convergence"]["residual"].as<double>();
  ideck.cfli  = iconfig["dynamic CFL"]["cfli"].as<double>();
  ideck.cflf  = iconfig["dynamic CFL"]["cflf"].as<double>();
  // IO parameters
  ideck.foutDir = iconfig["output"]["flowviz"]["directory"].as<string>();
  ideck.ghost = iconfig["output"]["flowviz"]["ghost cells"].as<bool>();
  // BC parameters
  ideck.bcLeft.type = iconfig["boundary conditions"]["left"]["type"].as<string>();
  if (ideck.bcLeft.type=="moving wall" || ideck.bcLeft.type=="inlet") {
    ideck.bcLeft.velocity[0] = iconfig["boundary conditions"]["left"]["velocity"]["u"].as<double>();
    ideck.bcLeft.velocity[1] = iconfig["boundary conditions"]["left"]["velocity"]["v"].as<double>();
  }
  ideck.bcRight.type = iconfig["boundary conditions"]["right"]["type"].as<string>();
  if (ideck.bcRight.type=="moving wall" || ideck.bcRight.type=="inlet") {
    ideck.bcRight.velocity[0] = iconfig["boundary conditions"]["right"]["velocity"]["u"].as<double>();
    ideck.bcRight.velocity[1] = iconfig["boundary conditions"]["right"]["velocity"]["v"].as<double>();
  }
  ideck.bcBottom.type = iconfig["boundary conditions"]["bottom"]["type"].as<string>();
  if (ideck.bcBottom.type=="moving wall" || ideck.bcBottom.type=="inlet") {
    ideck.bcBottom.velocity[0] = iconfig["boundary conditions"]["bottom"]["velocity"]["u"].as<double>();
    ideck.bcBottom.velocity[1] = iconfig["boundary conditions"]["bottom"]["velocity"]["v"].as<double>();
  }
  ideck.bcTop.type = iconfig["boundary conditions"]["top"]["type"].as<string>();
  if (ideck.bcTop.type=="moving wall" || ideck.bcTop.type=="inlet") {
    ideck.bcTop.velocity[0] = iconfig["boundary conditions"]["top"]["velocity"]["u"].as<double>();
    ideck.bcTop.velocity[1] = iconfig["boundary conditions"]["top"]["velocity"]["v"].as<double>();
  }
  ideck.uinit = iconfig["initial conditions"]["u"].as<double>();
  ideck.vinit = iconfig["initial conditions"]["v"].as<double>();
  ideck.pinit = iconfig["initial conditions"]["p"].as<double>();
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
