#ifndef INPUT_H
#define INPUT_H

#include <yaml-cpp/yaml.h>
#include <argparse/argparse.hpp>

using namespace std;

extern YAML::Node config;

void parse_user(const string& inFile);

void getUserInput(int argc, char* argv[],
                  argparse::ArgumentParser& program);

#endif // INPUT_H
