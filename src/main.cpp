/**
 * @file main.cpp
 * @author Caleb Hash
 * @date 2025-02-20
 * @brief Two-Phase Flow Solver - HW4
 * @details first part of the final project - HW4
 */

#include <mesh.h>
#include <output.h>               // output library
#include <input.h>            // input deck reader
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
  
  // call input file parser
  YAML::Node config = parse_user(inFile);
  
  // generate the domain
  double dx,dy,dz;
  double lx = config["domain"]["lengths"]["x"].as<double>();
  double ly = config["domain"]["lengths"]["y"].as<double>();
  double lz = config["domain"]["lengths"]["z"].as<double>();
  int nx = config["domain"]["dimensions"]["x"].as<double>();
  int ny = config["domain"]["dimensions"]["y"].as<double>();
  int nz = config["domain"]["dimensions"]["z"].as<double>();

  // need to call getDomainIndices() here

  // get the TOTAL number of cells here so we can allocate properly

  // initialize the Kokkos matrices
  FMatrix<double> xc(nx,ny,nz),yc(nx,ny,nz),zc(nx,ny,nz),
                  xn(nx+1,ny+1,nz+1),yn(nx+1,ny+1,nz+1),zn(nx+1,ny+1,nz+1);
  
  // call mesh generator
  mesher3D(lx,ly,lz,nx,ny,nz,xc,yc,zc,xn,yn,zn,dx,dy,dz);

  // output grid file
  vtk_output_2D(nx,ny,xn,yn);

  return 0;
}
