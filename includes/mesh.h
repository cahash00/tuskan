#ifndef MESH_H
#define MESH_H

#include <matar.h>
using namespace std;
void mesher2D(const double& lx, const double& ly, 
              const int& nx, const int& ny,
              mtr::FMatrixKokkos<double>& xc, mtr::FMatrixKokkos<double>& yc,
              mtr::FMatrixKokkos<double>& xn, mtr::FMatrixKokkos<double>& yn,
              double& dx, double& dy);

#endif // MESH_H
