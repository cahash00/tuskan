/**
 * @file inputDeck.cpp
 * @brief input deck reader library. Going to use YAML-cpp to parse the input 
 *        file.
 * @author Caleb Hash
 * @date 2025-02-22
 */

#include <input.h>
#include <yaml-cpp/yaml.h>
#include <pprint.hpp>
#include <argparse/argparse.hpp>

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
  // BC parameters
  ideck.bcLeft = config["boundary conditions"]["left"].as<string>();
  ideck.bcRight = config["boundary conditions"]["right"].as<string>();
  ideck.bcBottom = config["boundary conditions"]["bottom"].as<string>();
  ideck.bcTop = config["boundary conditions"]["top"].as<string>();
  ideck.uinit = config["initial conditions"]["u"].as<double>();
  ideck.vinit = config["initial conditions"]["v"].as<double>();
  return ideck;
}

ConfigData parseInputDeck(const string& inFile){
  // power up the printer
  pprint::PrettyPrinter printer;

  // ... create struct for the yaml node
  try {
    YAML::Node config = YAML::LoadFile(inFile);
    return ConfigData::fromYAMLConfig(config);
  } catch (const std::exception& e) {
    std::cerr << "Error loading YAML input deck: " << e.what() << std::endl;
    throw;
  }

} // end function parseInputDeck

/**
 * Arugment parser
 */
void getUserInput(int argc, char* argv[], argparse::ArgumentParser& program) {
  // add arguments to the program
  program.add_argument("-i","--input")
         .required()
         .default_value("master.yaml")
         .help("Input deck for the calculation.");
  // try to parse the arguments - error out if it fails
  try {
    program.parse_args(argc,argv);
  } catch (const exception& err) {
    cerr << err.what() << endl;
    cerr << program << endl;
    throw runtime_error("ERROR: Could not process CLI inputs.");
  }
} // end function getUserInput

} // end namespace IO
