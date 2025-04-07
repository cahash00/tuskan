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

void check_directories(const string& foutDir) {
  if (fs::exists(foutDir)) {
    // clear it
    for (const auto& entry : fs::directory_iterator(foutDir)) {
      fs::remove_all(entry.path());
    }
    IO::logger->info("  {} has been cleared.",foutDir);
  } else {
    fs::create_directories(foutDir);
    IO::logger->info("  {} has been created.",foutDir);
  }
} // end check_directories

/**
 * @brief VTK output for a 2D mesh
 * @param[in] nx number of cells in y-direction
 * @param[in] ny number of cells in y-direction
 * @param[in] xn matrix of x nodal values
 * @param[in] yn matrix of y nodal values
 */
static void vtk_output_2D_node_IMPL(const string& caseName,
                                    const string& foutDir,
                                    const bool ghost,
                                    const mtr::FMatrix<double>& xc,
                                    const mtr::FMatrix<double>& yc,
                                    const mtr::FMatrix<double>& u,
                                    const mtr::FMatrix<double>& v,
                                    const mtr::FMatrix<double>& p,
                                    const mtr::FMatrix<double>& phi) {
  /**
   * get cell-centered solution variables
   * pressure is already cell-centered
   */
  int xstr,xend,ystr,yend,dnx,dny;
  mtr::FMatrix<double> uc(p.dims(1),p.dims(2));
  mtr::FMatrix<double> vc(p.dims(1),p.dims(2));
  for (int j = jstr-1; j <= jend; j++) {
    for (int i = istr-1; i <= iend; i++) {
      uc(i,j) = 0.5*(u(i,j) + u(i+1,j));
    }
  }
  for (int j = jstr-1; j <= jend; j++) {
    for (int i = istr-1; i <= iend; i++) {
      vc(i,j) = 0.5*(v(i,j)+v(i,j+1));
    }
  }
  
  if (ghost) {
    xstr = istr;
    xend = iend-1;
    ystr = jstr;
    yend = jend-1;
    dnx = nx;
    dny = ny;
  } else {
    xstr = istr;
    xend = iend-1;
    ystr = jstr;
    yend = jend-1;
    dnx = nx;
    dny = ny;
  }
  

  /**
   * file output
   */
  FILE *out;
  std::ostringstream foutss;
  foutss << foutDir << "/flow." << caseName << ".vtk";
  string fout = foutss.str();
  out = fopen(fout.c_str(),"w");
  fprintf(out,"# vtk DataFile Version 3.0\n");
  fprintf(out,"TUSKAN 2D Flow Solution File.\n");
  fprintf(out,"ASCII\n");
  fprintf(out,"DATASET STRUCTURED_GRID\n");
  fprintf(out,"DIMENSIONS %d %d %d\n",dnx,dny,1);
  fprintf(out,"POINTS %d float\n",(dnx)*(dny));
  for (int j = ystr; j <= yend; j++) {
    for (int i = xstr; i <= xend; i++) {
      fprintf(out,"%f %f %f\n",xc(i,j),yc(i,j),0.0);
    }
  }
  fprintf(out,"POINT_DATA %d\n",(dnx)*(dny));
  fprintf(out,"VECTORS velocity float \n");
  for (int j = ystr; j <= yend; j++) {
    for (int i = xstr; i <= xend; i++) {
      fprintf(out,"%f %f %f\n",uc(i,j),vc(i,j),0.0);
    }
  }
  fprintf(out,"SCALARS p float 1\n");
  fprintf(out,"LOOKUP_TABLE default\n");
  for (int j = ystr; j <= yend; j++) {
    for (int i = xstr; i <= xend; i++) {
      fprintf(out,"%f\n",p(i,j));
    }
  }
  fprintf(out,"SCALARS phi float 1\n");
  fprintf(out,"LOOKUP_TABLE default\n");
  for (int j = ystr; j <= yend; j++) {
    for (int i = xstr; i <= xend; i++) {
      fprintf(out,"%f\n",phi(i,j));
    }
  }
  fclose(out);
}
/**
 * overloaded vtk node output functions
 */
void vtk_output_2D_node(const int& ii,
                   const string& foutDir,
                   const bool ghost,
                   const mtr::FMatrix<double>& xc,
                   const mtr::FMatrix<double>& yc,
                   const mtr::FMatrix<double>& u,
                   const mtr::FMatrix<double>& v,
                   const mtr::FMatrix<double>& p,
                   const mtr::FMatrix<double>& phi) {
  std::ostringstream foutss;
  foutss << setw(5) << std::setfill('0') << ii;
  string caseName = foutss.str();
  vtk_output_2D_node_IMPL(caseName,foutDir,ghost,xc,yc,u,v,p,phi);
}
void vtk_output_2D_node(const string& caseName,
                   const string& foutDir,
                   const bool ghost,
                   const mtr::FMatrix<double>& xc,
                   const mtr::FMatrix<double>& yc,
                   const mtr::FMatrix<double>& u,
                   const mtr::FMatrix<double>& v,
                   const mtr::FMatrix<double>& p,
                   const mtr::FMatrix<double>& phi) {
  vtk_output_2D_node_IMPL(caseName,foutDir,ghost,xc,yc,u,v,p,phi);
}

} // end namepsace IO
