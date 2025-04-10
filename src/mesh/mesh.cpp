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
              mtr::FMatrix<double>& xu, mtr::FMatrix<double>& yu,
              mtr::FMatrix<double>& xv, mtr::FMatrix<double>& yv,
              double& dx, double& dy){
  // get dx and dy
  dx = lx/static_cast<double>(nx);
  dy = ly/static_cast<double>(ny);
  
  // calculate cell center values
  for (int j = jstr; j <= jend+1; j++) {
    for (int i = istr; i <= iend+1; i++) {
      xc(i,j) = (i-(istr))*dx + dx*0.5;
      yc(i,j) = (j-(jstr))*dy + dy*0.5;
    }
  }

  // calculate the cell nodal points for u and v
  for (int j = jstr-1; j <= jend+1; j++) {
    for (int i = istr-1; i <= iend+1; i++) {
      xu(i,j) = (i-(istr-1)) * dx;
      yu(i,j) = (j-(jstr)) * dy + 0.5*dy;
    }
  }
  for (int j = jstr-1; j <= jend+1; j++) {
    for (int i = istr-1; i <= iend+1; i++) {
      xv(i,j) = (i-(istr)) * dx + 0.5*dx;
      yv(i,j) = (j-(jstr-1)) * dy;
    }
  }
} // end mesher2D

} // end namespace mesh
