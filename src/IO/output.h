#ifndef OUTPUT_H
#define OUTPUT_H

#include <matar.h>

using namespace mtr;

void vtk_output_2D(const int& nx, const int& ny,
    const FMatrix<double>& xn,
    const FMatrix<double>& yn);

void vtk_output_3D(const int ii,const int& nx, const int& ny,const int& nz,
    const FMatrix<double>& xn,
    const FMatrix<double>& yn,
    const FMatrix<double>& zn);

void splashScreen();
#endif // OUTPUT_H
