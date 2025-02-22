/**
 * @file main.cpp
 * @author Caleb Hash
 * @date 2025-02-20
 * @brief Two-Phase Flow Solver - HW4
 * @details first part of the final project - HW4
 */

#include <output.h>
#include <inputYAML.h>
#include <matar.h>
#include <Kokkos_Core.hpp>
#include <yaml-cpp/yaml.h>
#include <argparse/argparse.hpp>
#include <pprint.hpp>

using namespace std;
using namespace mtr;

/**
 * @brief a 2D mesh builder
 * @param[in] lx length in the x-direction
 * @param[in] ly length in the y-direction
 * @param[in] nx number of cells in the x-direction
 * @param[in] ny number of cells in the y-direction
 * @param[inout] xc x coordinates for cell centers
 * @param[inout] yc y coordinates for cell centers
 * @param[inout] xn x coordinates for nodes
 * @param[inout] yn y coordinates for nodes
 * @param[inout] dx cell width in the x-direction
 * @param[inout] dy cell width in the y-direction
 */
void mesher2D(const double& lx, const double& ly, 
              const int& nx, const int& ny,
              FMatrixKokkos<double>& xc, FMatrixKokkos<double>& yc,
              FMatrixKokkos<double>& xn, FMatrixKokkos<double>& yn,
              double& dx, double& dy){
  // get dx and dy
  dx = lx/(double)nx;
  dy = ly/(double)ny;

  // calculate center values
  int i,j;
  DO2D(j,1,ny,i,1,nx,{
    xc(i,j) = (i-1)*dx + dx*0.5;
    yc(i,j) = (j-1)*dy + dy*0.5;
  });
  DO2D(j,1,ny+1,i,1,nx+1,{
      xn(i,j) = (i-1) * dx;
      yn(i,j) = (j-1) * dy;
  });
}

int main(int argc, char* argv[]){
  // Fire up the printer
  pprint::PrettyPrinter printer;

  /**
   * command line arguments
   * uses argparse to parse out command line inputs.
   */
  argparse::ArgumentParser program("TUSKAN","0.0.0");
  program.add_argument("-i","--input")
    .required()
    .help("Input deck for the calculation.");
  try {
    program.parse_args(argc,argv);
  } catch (const exception& err) {
    cerr << err.what() << endl;
    cerr << program << endl;
    return 1;
  }
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
  FMatrixKokkos<double> xc(nx,ny),yc(nx,ny),xn(nx+1,ny+1),yn(nx+1,ny+1);

  // Call the mesher
  mesher2D(lx,ly,nx,ny,xc,yc,xn,yn,dx,dy);

  // output grid file
  vtk_output_2D(nx,ny,xn,yn);

  return 0;
}
