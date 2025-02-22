/**
 * @file inputDeck.cpp
 * @brief input deck reader library. Going to use YAML-cpp to parse the input 
 *        file.
 * @author Caleb Hash
 * @date 2025-02-22
 */

#include <inputYAML.h>
#include "yaml-cpp/yaml.h"
#include <pprint.hpp>

YAML::Node parse_user(const std::string& inFile){
  // power up the printer
  pprint::PrettyPrinter printer;

  // load up the master
  YAML::Node config = YAML::LoadFile(inFile);

  //check if it loaded properly
  if (!config) {
    throw std::runtime_error("Failed to load the YAML file: "+inFile);
  }

  // print some stuff for testing
  printer.print(config["iterations"]);
  return config;
}
