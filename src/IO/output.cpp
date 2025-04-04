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

using namespace std;

namespace IO {

/**
 * @brief VTK output for a 2D mesh
 * @param[in] nx number of cells in y-direction
 * @param[in] ny number of cells in y-direction
 * @param[in] xn matrix of x nodal values
 * @param[in] yn matrix of y nodal values
 */
static void vtk_output_2D_node_IMPL(const string& caseName,
                               const string& foutDir,
                               const mtr::FMatrix<double>& u,
                               const mtr::FMatrix<double>& v,
                               const mtr::FMatrix<double>& p,
                               const mtr::FMatrix<double>& xn,
                               const mtr::FMatrix<double>& yn) {
  // power on the printer
  pprint::PrettyPrinter printer;
  FILE *out;
  std::ostringstream foutss;
  foutss << foutDir << "/flow." << caseName << ".vtk";
  string fout = foutss.str();
  out = fopen(fout.c_str(),"w");
  fprintf(out,"# vtk DataFile Version 3.0\n");
  fprintf(out,"TUSKAN 2D Flow Solution File.\n");
  fprintf(out,"ASCII\n");
  fprintf(out,"DATASET STRUCTURED_GRID\n");
  fprintf(out,"DIMENSIONS %d %d %d\n",nx+2,ny+2,1);
  fprintf(out,"POINTS %d float\n",(nx+2)*(ny+2));
  for (int j = jstr-1; j <= jend; j++) {
    for (int i = istr-1; i <= iend; i++) {
      fprintf(out,"%f %f %f\n",xn(i,j),yn(i,j),0.0);
    }
  }
  fprintf(out,"POINT_DATA %d\n",(nx+2)*(ny+2));
  fprintf(out,"SCALARS u float 1\n");
  fprintf(out,"LOOKUP_TABLE default\n");
  for (int j = jstr-1; j <= jend; j++) {
    for (int i = istr-1; i <= iend; i++) {
      fprintf(out,"%f\n",u(i,j));
    }
  }
  fprintf(out,"SCALARS v float 1\n");
  fprintf(out,"LOOKUP_TABLE default\n");
  for (int j = jstr-1; j <= jend; j++) {
    for (int i = istr-1; i <= iend; i++) {
      fprintf(out,"%f\n",v(i,j));
    }
  }
  fprintf(out,"SCALARS p float 1\n");
  fprintf(out,"LOOKUP_TABLE default\n");
  for (int j = jstr-1; j <= jend; j++) {
    for (int i = istr-1; i <= iend; i++) {
      fprintf(out,"%f\n",p(i,j));
    }
  }
  fclose(out);
}
/**
 * overloaded vtk node output functions
 */
void vtk_output_2D_node(const int& ii,
                   const string& foutDir,
                   const mtr::FMatrix<double>& u,
                   const mtr::FMatrix<double>& v,
                   const mtr::FMatrix<double>& p,
                   const mtr::FMatrix<double>& xn,
                   const mtr::FMatrix<double>& yn) {
  std::ostringstream foutss;
  foutss << setw(5) << std::setfill('0') << ii;
  string caseName = foutss.str();
  vtk_output_2D_node_IMPL(caseName,foutDir,u,v,p,xn,yn);
}
void vtk_output_2D_node(const string& caseName,
                   const string& foutDir,
                   const mtr::FMatrix<double>& u,
                   const mtr::FMatrix<double>& v,
                   const mtr::FMatrix<double>& p,
                   const mtr::FMatrix<double>& xn,
                   const mtr::FMatrix<double>& yn) {
  vtk_output_2D_node_IMPL(caseName,foutDir,u,v,p,xn,yn);
}

} // end namepsace IO
