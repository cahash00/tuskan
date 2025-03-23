/**
 * @file mesh.cpp
 * @brief Implementation of the meshing functions
 * @author Caleb Hash
 * @date 2025-02-28
 */

#include <mesh.h>   // mesh library
#include <matar.h>  // MATAR headers
#include <params.h>

using namespace std;
using namespace mtr;

/**********************
 * 2D Meshing Routine *
 **********************/
void mesher2D(const double& lx, const double& ly, 
              const int& nx,const int& ny,
              FMatrix<double>& xc, FMatrix<double>& yc,
              FMatrix<double>& xn, FMatrix<double>& yn,
              double& dx, double& dy){
  // get dx and dy
  dx = lx/static_cast<double>(nx);
  dy = ly/static_cast<double>(ny);
  
  // calculate cell center values
  int i,j;
  DO_ALL(j,jstr-nghosts,jend+nghosts,
         i,istr-nghosts,iend+nghosts,{
    xc(i,j) = (i-(1+nghosts))*dx + dx*0.5;
    yc(i,j) = (j-(1+nghosts))*dy + dy*0.5;
  });

  // calculate the cell nodal points
  DO_ALL(j,jstr-nghosts,jend+nghosts,
         i,istr-nghosts,iend+nghosts,{
    xn(i,j) = (i-(1+nghosts)) * dx;
    yn(i,j) = (j-(1+nghosts)) * dy;
  });
}
