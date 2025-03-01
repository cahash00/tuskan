#ifndef INPUT_H
#define INPUT_H

#include <yaml-cpp/yaml.h>
#include <argparse/argparse.hpp>

using namespace std;

YAML::Node parse_user(const string& inFile);

argparse::ArgumentParser getUserInput(int argc, char* argv)

#endif // INPUT_H
