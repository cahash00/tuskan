/**
 * @file args.cpp
 * @brief input deck reader library. Going to use YAML-cpp to parse the input 
 *        file.
 * @author Caleb Hash
 * @date April 06, 2025
 */

#include <args.h>
#include <argparse/argparse.hpp>
using namespace std;
namespace IO {

/**
 * Arugment parser
 */
void getUserInput(int argc, char* argv[], argparse::ArgumentParser& program) {
  // add arguments to the program
  program.add_argument("-i","--input")
         .required()
         .default_value("in2d.yaml")
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
