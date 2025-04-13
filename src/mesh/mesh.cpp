/**
 * @file mesh.cpp
 * @brief Implementation of the meshing functions
 * @author Caleb Hash
 * @date 2025-02-28
 */

#include <mesh.h>   // mesh library
#include <matar.h>  // MATAR headers
#include <params.h>
#include <types.h>
#include <communication.h>


namespace mesh {
/**
 * Domain generation in 2D
 */
void mesher2D(fmat<double>& xc, fmat<double>& yc,
              fmat<double>& xn, fmat<double>& yn,
              const double& dx, const double& dy){
  // calculate cell center values
  for (int j = jstr-nghosts; j <= jend; j++) {
    for (int i = istr-nghosts; i <= iend; i++) {
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
