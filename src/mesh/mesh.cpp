
#include <mesh.h>   // mesh library
#include <matar.h>  // MATAR headers

using namespace std;
using namespace mtr;
/**
 * @brief a 2D mesh builder
 * @param[in] lx length in the x-direction
 * @param[in] ly length in the y-direction
 * @param[in] nx number of cells in the x-direction
 * @param[in] ny number of cells in the y-direction
 * @param[inout] xc x coordinates for cell centers
 * @param[inout] yc y coordinates for cell centers
 * @param[inout] xn x coordinates for nodes
 * @param[inout] yn y coordinates for nodes
 * @param[inout] dx cell width in the x-direction
 * @param[inout] dy cell width in the y-direction
 */
void mesher2D(const double& lx, const double& ly, 
              const int& nx, const int& ny,
              FMatrixKokkos<double>& xc, FMatrixKokkos<double>& yc,
              FMatrixKokkos<double>& xn, FMatrixKokkos<double>& yn,
              double& dx, double& dy){
  // get dx and dy
  dx = lx/(double)nx;
  dy = ly/(double)ny;

  // calculate center values
  int i,j;
  DO2D(j,1,ny,i,1,nx,{
      xc(i,j) = (i-1)*dx + dx*0.5;
      yc(i,j) = (j-1)*dy + dy*0.5;
      });
  DO2D(j,1,ny+1,i,1,nx+1,{
      xn(i,j) = (i-1) * dx;
      yn(i,j) = (j-1) * dy;
      });
}
