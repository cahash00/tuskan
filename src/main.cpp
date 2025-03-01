/**
 * @file main.cpp
 * @author Caleb Hash
 * @date 2025-02-20
 * @brief Two-Phase Flow Solver - HW4
 * @details first part of the final project - HW4
 */

#include <mesh.h>
#include <output.h>               // output library
#include <input.h>                // input deck reader
#include <matar.h>                // MATAR headers
#include <Kokkos_Core.hpp>        // Kokkos for initialization
#include <yaml-cpp/yaml.h>        // yaml-parser 
#include <argparse/argparse.hpp>  // argparse
#include <pprint.hpp>             // pretty printer 
#include <params.h>

using namespace std;
using namespace mtr;

int main(int argc, char* argv[]){
  // Fire up the printer
  pprint::PrettyPrinter printer;

  argparse::ArgumentParser program("TUSKAN","0.0.0");
  getUserInput(argc,argv,program);

  /***************************************************************************
   *                           MAIN PROGRAM START                            *
   **************************************************************************/

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

  int istr,iend,jstr,jend,kstr,kend;
  getDomainIndices(nx,ny,nz,istr,iend,jstr,jend,kstr,kend);

  // get the TOTAL number of cells here so we can allocate properly
  vector<int> ndims(3,0);
  ndims[0] = nx+nghosts*2;
  ndims[1] = ny+nghosts*2;
  ndims[2] = nz+nghosts*2;
  printer.print("dimensions");
  printer.print(nx,ny,nz);

  // initialize the Kokkos matrices
  FMatrix<double> xc(ndims[0],ndims[1],ndims[2]),
                  yc(ndims[0],ndims[1],ndims[2]),
                  zc(ndims[0],ndims[1],ndims[2]),
                  xn(ndims[0]+1,ndims[1]+1,ndims[2]+1),
                  yn(ndims[0]+1,ndims[1]+1,ndims[2]+1),
                  zn(ndims[0]+1,ndims[1]+1,ndims[2]+1);
  
  // call mesh generator
  mesher3D(lx,ly,lz,nx,ny,nz,xc,yc,zc,xn,yn,zn,dx,dy,dz);

  // output grid file
  vtk_output_3D(nx,ny,nz,xn,yn,zn);

  return 0;
}
