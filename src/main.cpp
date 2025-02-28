/**
 * @file main.cpp
 * @author Caleb Hash
 * @date 2025-02-20
 * @brief Two-Phase Flow Solver - HW4
 * @details first part of the final project - HW4
 */

#include <mesh.h>
#include <output.h>               // output library
#include <inputYAML.h>            // input deck reader
#include <matar.h>                // MATAR headers
#include <Kokkos_Core.hpp>        // Kokkos for initialization
#include <yaml-cpp/yaml.h>        // yaml-parser 
#include <argparse/argparse.hpp>  // argparse
#include <pprint.hpp>             // pretty printer 

using namespace std;
using namespace mtr;

int main(int argc, char* argv[]){
  // Fire up the printer
  pprint::PrettyPrinter printer;

  argparse::ArgumentParser program = getUserInput(argc,argv);

  /**
   *                           MAIN PROGRAM START
   */

  // assign the input file that was given
  auto inFile = program.get<string>("-i");

  // automatically initialize and finalize Kokkos with the main program
  Kokkos::ScopeGuard kokkos_guard(argc, argv); 
  YAML::Node config;
  
  // generate the domain
  double dx,dy;
  double lx = 2.0;
  double ly = 1.0;
  int nx = 10;
  int ny = 10;

  // call input file parser
  config = parse_user(inFile);
  // initialize the Kokkos matrices
  FMatrix<double> xc(nx,ny),yc(nx,ny),xn(nx+1,ny+1),yn(nx+1,ny+1);

  // Call the mesher
  mesher2D(lx,ly,nx,ny,xc,yc,xn,yn,dx,dy);

  // output grid file
  vtk_output_2D(nx,ny,xn,yn);

  return 0;
}
