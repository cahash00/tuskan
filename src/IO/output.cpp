/**
 * @file output.cpp
 * @brief output library for the solver. Writes in legacy vtk format for 
 *        structured grids
 * @todo add support for solution array eventually.
 * @author Caleb Hash
 * @date 2025-02-22
 */

#include <output.h>
#include <matar.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkStructuredGrid.h>
#include <vtkXMLStructuredGridWriter.h>
#include <pprint.hpp>

using namespace mtr;

/**
 * @brief VTK output for a 2D mesh
 * @param[in] nx number of cells in y-direction
 * @param[in] ny number of cells in y-direction
 * @param[in] xn matrix of x nodal values
 * @param[in] yn matrix of y nodal values
 */
void vtk_output_2D(const int& nx, const int& ny,
                   const FMatrixKokkos<double>& xn,
                   const FMatrixKokkos<double>& yn) {
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
void vtk_output_3D(const int& nx, const int& ny,const int& nz,
                   const FMatrixKokkos<double>& xn,
                   const FMatrixKokkos<double>& yn,
                   const FMatrixKokkos<double>& zn) {
  // power on the printer
  pprint::PrettyPrinter printer;
  FILE *out;
  out = fopen("test.vts","w");
  fprintf(out,"# vtk DataFile Version 3.0\n");
  fprintf(out,"Insert Any Information Here.\n");
  fprintf(out,"ASCII\n");
  fprintf(out,"DATASET STRUCTURED_GRID\n");
  fprintf(out,"DIMENSIONS %d %d %d\n",nx+1,ny+1,1);
  fprintf(out,"POINTS %d float\n",(nx+1)*(ny+1)*(nz+1));
  DO_LOOP(k,1,nz+1,{
    DO_LOOP(j,1,ny+1,{
      DO_LOOP(i,1,nx+1,{
        fprintf(out,"%f %f %f\n",xn(i,j),yn(i,j),zn(i,j));
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
