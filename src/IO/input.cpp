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

ConfigData ConfigData::fromYAMLConfig(const YAML::Node& config) {
  ConfigData ideck;

  // domain settings
  ideck.lx = config["domain"]["lengths"]["x"].as<double>();
  ideck.ly = config["domain"]["lengths"]["y"].as<double>();
  ideck.nx = config["domain"]["dimensions"]["x"].as<int>();
  ideck.ny = config["domain"]["dimensions"]["y"].as<int>();
  // solver settings
  ideck.iter    = config["solver"]["iterations"].as<int>();
  ideck.fvflag  = config["output"]["flowviz"]["enabled"].as<bool>();
  ideck.fvfreq  = config["output"]["flowviz"]["frequency"].as<int>();
  ideck.resflag = config["output"]["residuals"]["enabled"].as<bool>();
  ideck.resfreq = config["output"]["residuals"]["frequency"].as<int>();
  ideck.resFile = config["output"]["residuals"]["file"].as<string>();
  ideck.pMethod = config["solver"]["pressure solver"]["method"].as<string>();
  if (ideck.pMethod == "SOR") {
    if (config["solver"]["pressure solver"]["omega"]) {
      ideck.sorOmega = config["solver"]["pressure solver"]["omega"].as<double>(); 
    } else {
      throw runtime_error("ERROR: INPUT DECK: SOR solver requires \"omega\".");
    }
  }
  // convergence criteria
  ideck.toler = config["convergence"]["residual"].as<double>();
  ideck.cfli  = config["dynamic CFL"]["cfli"].as<double>();
  ideck.cflf  = config["dynamic CFL"]["cflf"].as<double>();
  // IO parameters
  ideck.foutDir = config["output"]["flowviz"]["directory"].as<string>();
  ideck.ghost = config["output"]["flowviz"]["ghost cells"].as<bool>();
  // BC parameters
  ideck.bcLeft.type = config["boundary conditions"]["left"]["type"].as<string>();
  if (ideck.bcLeft.type=="moving wall" || ideck.bcLeft.type=="inlet") {
    ideck.bcLeft.velocity[0] = config["boundary conditions"]["left"]["velocity"]["u"].as<double>();
    ideck.bcLeft.velocity[1] = config["boundary conditions"]["left"]["velocity"]["v"].as<double>();
  }
  ideck.bcRight.type = config["boundary conditions"]["right"]["type"].as<string>();
  if (ideck.bcRight.type=="moving wall" || ideck.bcRight.type=="inlet") {
    ideck.bcRight.velocity[0] = config["boundary conditions"]["right"]["velocity"]["u"].as<double>();
    ideck.bcRight.velocity[1] = config["boundary conditions"]["right"]["velocity"]["v"].as<double>();
  }
  ideck.bcBottom.type = config["boundary conditions"]["bottom"]["type"].as<string>();
  if (ideck.bcBottom.type=="moving wall" || ideck.bcBottom.type=="inlet") {
    ideck.bcBottom.velocity[0] = config["boundary conditions"]["bottom"]["velocity"]["u"].as<double>();
    ideck.bcBottom.velocity[1] = config["boundary conditions"]["bottom"]["velocity"]["v"].as<double>();
  }
  ideck.bcTop.type = config["boundary conditions"]["top"]["type"].as<string>();
  if (ideck.bcTop.type=="moving wall" || ideck.bcTop.type=="inlet") {
    ideck.bcTop.velocity[0] = config["boundary conditions"]["top"]["velocity"]["u"].as<double>();
    ideck.bcTop.velocity[1] = config["boundary conditions"]["top"]["velocity"]["v"].as<double>();
  }
  ideck.uinit = config["initial conditions"]["u"].as<double>();
  ideck.vinit = config["initial conditions"]["v"].as<double>();
  return ideck;
}

ConfigData parseInputDeck(const string& inFile){

  // ... create struct for the yaml node
  try {
    YAML::Node config = YAML::LoadFile(inFile);
    return ConfigData::fromYAMLConfig(config);
  } catch (const std::exception& e) {
    std::cerr << "Error loading YAML input deck: " << e.what() << std::endl;
    throw;
  }

} // end function parseInputDeck

} // end namespace IO
