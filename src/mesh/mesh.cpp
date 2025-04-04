/**
 * @file mesh.cpp
 * @brief Implementation of the meshing functions
 * @author Caleb Hash
 * @date 2025-02-28
 */

#include <mesh.h>   // mesh library
#include <matar.h>  // MATAR headers
#include <params.h>


namespace mesh {
/**
 * Domain generation in 2D
 */
void mesher2D(const double& lx, const double& ly, 
              mtr::FMatrix<double>& xc, mtr::FMatrix<double>& yc,
              mtr::FMatrix<double>& xn, mtr::FMatrix<double>& yn,
              double& dx, double& dy){
  // get dx and dy
  dx = lx/static_cast<double>(nx);
  dy = ly/static_cast<double>(ny);
  
  // calculate cell center values
  int i,j;
  for (int j = jstr-nghosts; j <= jend+nghosts; j++) {
    for (int i = istr-nghosts; i <= iend+nghosts; i++) {
      xc(i,j) = (i-(1+nghosts))*dx + dx*0.5;
      yc(i,j) = (j-(1+nghosts))*dy + dy*0.5;
    }
  }

  // calculate the cell nodal points
  for (int j = jstr-nghosts; j <= jend+nghosts; j++) {
    for (int i = istr-nghosts; i <= iend+nghosts; i++) {
      xn(i,j) = (i-(1+nghosts)) * dx;
      yn(i,j) = (j-(1+nghosts)) * dy;
    }
  }
} // end mesher2D

} // end namespace mesh
