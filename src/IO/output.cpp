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

using namespace mtr;
using namespace std;

/**
 * @brief VTK output for a 2D mesh
 * @param[in] nx number of cells in y-direction
 * @param[in] ny number of cells in y-direction
 * @param[in] xn matrix of x nodal values
 * @param[in] yn matrix of y nodal values
 */
void vtk_output_2D(const int& nx, const int& ny,
                   const FMatrix<double>& xn,
                   const FMatrix<double>& yn) {
  // power on the printer
  pprint::PrettyPrinter printer;
  FILE *out;
  out = fopen("test.vts","w");
  fprintf(out,"# vtk DataFile Version 3.0\n");
  fprintf(out,"Insert Any Information Here.\n");
  fprintf(out,"ASCII\n");
  fprintf(out,"DATASET STRUCTURED_GRID\n");
  fprintf(out,"DIMENSIONS %d %d %d\n",nx+1,ny+1,1);
  fprintf(out,"POINTS %d float\n",(nx+1)*(ny+1));
  DO_LOOP(j,1,ny+1,{
    DO_LOOP(i,1,nx+1,{
      fprintf(out,"%f %f %f\n",xn(i,j),yn(i,j),0.0);
    });
  });
  fprintf(out,"POINT_DATA %d\n",(nx+1)*(ny+1));
  fprintf(out,"SCALARS scalar_name float 1\n");
  fprintf(out,"LOOKUP_TABLE default\n");
  DO_LOOP(j,1,ny+1,{
    DO_LOOP(i,1,nx+1,{
      fprintf(out,"%f\n",0.0);
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
void vtk_output_3D(const int ii,const int& nx, const int& ny,const int& nz,
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

  out = fopen(fout.c_str(),"w");
  fprintf(out,"# vtk DataFile Version 3.0\n");
  fprintf(out,"Insert Any Information Here.\n");
  fprintf(out,"ASCII\n");
  fprintf(out,"DATASET STRUCTURED_GRID\n");
  fprintf(out,"DIMENSIONS %d %d %d\n",nx+1,ny+1,nz+1);
  fprintf(out,"POINTS %d float\n",(nx+1)*(ny+1)*(nz+1));
  DO_LOOP(k,1,nz+1,{
    DO_LOOP(j,1,ny+1,{
      DO_LOOP(i,1,nx+1,{
        fprintf(out,"%f %f %f\n",xn(i,j,k),yn(i,j,k),zn(i,j,k));
      });
    });
  });
  fprintf(out,"POINT_DATA %d\n",(nx+1)*(ny+1)*(nz+1));
  fprintf(out,"SCALARS scalar_name float 1\n");
  fprintf(out,"LOOKUP_TABLE default\n");
  DO_LOOP(k,1,nz+1,{
    DO_LOOP(j,1,ny+1,{
      DO_LOOP(i,1,nx+1,{
        fprintf(out,"%f\n",0.0);
      });
    });
  });
  fclose(out);
}

void splashScreen() {
  std::cout << "           __" << std::endl;
  std::cout << ".-.__      \\ .-.  ___  __" << std::endl;
  std::cout << "|_|  '--.-.-(   \\/\\;;\\_\\.-._______.-." << std::endl;
  std::cout << "(-)___     \\ \\ .-\\ \\;;\\(   \\       \\ \\" << std::endl;
  std::cout << " Y    '---._\\_((Q)) \\;;\\ .-\\     __(_)" << std::endl;
  std::cout << " I           __'-' / .--.((Q))---'    \\," << std::endl;
  std::cout << " I     ___.-:    \\|  |   \\'-'_          \\" << std::endl;
  std::cout << " A  .-'      \\ .-.\\   \\   \\ \\ '--.__     '\\" << std::endl;
  std::cout << " |  |____.----((Q))\\   \\__|--\\_      \\     '" << std::endl;
  std::cout << "    ( )        '-'  \\_  :  \\-' '--.___\\" << std::endl;
  std::cout << "     Y                \\  \\  \\       \\(_)" << std::endl;
  std::cout << "     I                 \\  \\  \\         \\," << std::endl;
  std::cout << "     I                  \\  \\  \\          \\" << std::endl;
  std::cout << "     A                   \\  \\  \\          '\\" << std::endl;
  std::cout << "     |                    \\  \\__|           '" << std::endl;
  std::cout << "                           \\_:.  \\" << std::endl;
  std::cout << "                             \\ \\  \\" << std::endl;
  std::cout << "                              \\ \\  \\" << std::endl;
  std::cout << "                               \\_\\_|" << std::endl;
}
