/**
 * @file output.cpp
 * @brief output library for the solver. Writes in legacy vtk format for 
 *        structured grids
 * @todo add support for solution array eventually.
 * @author Caleb Hash
 * @date 2025-02-22
 */

#include <input.h>
#include <output.h>
#include <matar.h>
#include <pprint.hpp>
#include <string>
#include <format>
#include <iostream>
#include <params.h>
#include <input.h>

using namespace mtr;
using namespace std;

/**
 * @brief VTK output for a 2D mesh
 * @param[in] nx number of cells in y-direction
 * @param[in] ny number of cells in y-direction
 * @param[in] xn matrix of x nodal values
 * @param[in] yn matrix of y nodal values
 */
void vtk_output_2D(const int& ii,const FMatrix<double>& q,
                   const FMatrix<double>& xn,
                   const FMatrix<double>& yn) {
  // power on the printer
  pprint::PrettyPrinter printer;
  FILE *out;
  std::ostringstream foutss;
  foutss << config["output"]["directory"] << "/flow." << std::setw(4) 
    << std::setfill('0') << ii << ".vts";
  string fout = foutss.str();
  int nx = config["domain"]["dimensions"]["x"].as<int>();
  int ny = config["domain"]["dimensions"]["y"].as<int>();
  int nz = config["domain"]["dimensions"]["z"].as<int>();
  out = fopen(fout.c_str(),"w");
  fprintf(out,"# vtk DataFile Version 3.0\n");
  fprintf(out,"TUSKAN 2D Flow Solution File.\n");
  fprintf(out,"ASCII\n");
  fprintf(out,"DATASET STRUCTURED_GRID\n");
  fprintf(out,"DIMENSIONS %d %d %d\n",nx+1,ny+1,1);
  fprintf(out,"POINTS %d float\n",(nx+1)*(ny+1));
  DO_LOOP(j,jstr,jend,{
    DO_LOOP(i,istr,iend,{
      fprintf(out,"%f %f %f\n",xn(i,j,kstr),yn(i,j,kstr),0.0);
    });
  });
  fprintf(out,"POINT_DATA %d\n",(nx+1)*(ny+1));
  fprintf(out,"SCALARS u float 1\n");
  fprintf(out,"LOOKUP_TABLE default\n");
  DO_LOOP(j,jstr,jend,{
    DO_LOOP(i,istr,iend,{
      fprintf(out,"%f\n",q(1,i,j,kstr));
    });
  });
  fclose(out);
}

/**
 * @brief VTK output for a 3D mesh
 * @param[in] nx number of cells in x-direction
 * @param[in] ny number of cells in y-direction
 * @param[in] nz number of cells in z-direction
 * @param[in] xn matrix of x nodal values
 * @param[in] yn matrix of y nodal values
 * @param[in] zn matrix of z nodal values
 */
void vtk_output_3D(const int ii, const FMatrix<double>& q,
                   const FMatrix<double>& xn,
                   const FMatrix<double>& yn,
                   const FMatrix<double>& zn) {
  // power on the printer
  pprint::PrettyPrinter printer;
  FILE *out;
  std::ostringstream foutss;
  foutss << config["output"]["directory"] << "/flow." << std::setw(4) 
         << std::setfill('0') << ii << ".vts";
  string fout = foutss.str();

  // get sizing
  int nx = config["domain"]["dimensions"]["x"].as<int>();
  int ny = config["domain"]["dimensions"]["y"].as<int>();
  int nz = config["domain"]["dimensions"]["z"].as<int>();

  out = fopen(fout.c_str(),"w");
  fprintf(out,"# vtk DataFile Version 3.0\n");
  fprintf(out,"TUSKAN 3D Flow Solution File.\n");
  fprintf(out,"ASCII\n");
  fprintf(out,"DATASET STRUCTURED_GRID\n");
  fprintf(out,"DIMENSIONS %d %d %d\n",nx+1,ny+1,nz+1);
  fprintf(out,"POINTS %d float\n",(nx+1)*(ny+1)*(nz+1));
  DO_LOOP(k,kstr,kend,{
    DO_LOOP(j,jstr,jend,{
      DO_LOOP(i,istr,iend,{
        fprintf(out,"%f %f %f\n",xn(i,j,k),yn(i,j,k),zn(i,j,k));
      });
    });
  });
  fprintf(out,"POINT_DATA %d\n",(nx+1)*(ny+1)*(nz+1));
  fprintf(out,"SCALARS u float 1\n");
  fprintf(out,"LOOKUP_TABLE default\n");
  DO_LOOP(k,kstr,kend,{
    DO_LOOP(j,jstr,jend,{
      DO_LOOP(i,istr,iend,{
        fprintf(out,"%f\n",q(1,i,j,k));
      });
    });
  });
  fclose(out);
}
