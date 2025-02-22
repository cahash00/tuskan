#ifndef OUTPUT_H
#define OUTPUT_H

#include <matar.h>
using namespace mtr;

void vtk_output_2D(const int& nx, const int& ny,
                const FMatrixKokkos<double>& xn,
                const FMatrixKokkos<double>& yn);

void vtk_output_3D(const int& nx, const int& ny,const int& nz,
                   const FMatrixKokkos<double>& xn,
                   const FMatrixKokkos<double>& yn,
                   const FMatrixKokkos<double>& zn);

#endif // OUTPUT_H
