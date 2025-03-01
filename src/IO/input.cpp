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


using namespace std;

YAML::Node parse_user(const string& inFile){
  // power up the printer
  pprint::PrettyPrinter printer;

  // load up the master
  YAML::Node config = YAML::LoadFile(inFile);

  //check if it loaded properly
  if (!config) {
    throw runtime_error("Failed to load the YAML file: "+inFile);
  }
  return config;
}

/**
 * Arugment parser
 */
argparse::ArgumentParser getUserInput(int argc, char* argv) {
  argparse::ArgumentParser program("TUSKAN","0.0.0");
  program.add_argument("-i","--input")
    .required()
    .help("Input deck for the calculation.");
  try {
    program.parse_args(argc,argv);
  } catch (const exception& err) {
    cerr << err.what() << endl;
    cerr << program << endl;
    throw runtime_error("ERROR: Could not process CLI inputs.");
  }
  return program;
}
