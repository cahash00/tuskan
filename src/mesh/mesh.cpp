/**
 * @file mesh.cpp
 * @brief Implementation of the meshing functions
 * @author Caleb Hash
 * @date 2025-02-28
 */

#include <mesh.h>   // mesh library
#include <matar.h>  // MATAR headers

using namespace std;
using namespace mtr;

/**********************
 * 2D Meshing Routine *
 **********************/
void mesher2D(const double& lx, const double& ly, 
              const int& nx, const int& ny,
              FMatrix<double>& xc, FMatrix<double>& yc,
              FMatrix<double>& xn, FMatrix<double>& yn,
              double& dx, double& dy){
  // get dx and dy
  dx = lx/(double)nx;
  dy = ly/(double)ny;
  
  // add number of ghost cells to the domain
  int nghosts = 2; // 2 ghosts to support 2nd-order spatial discretization
  nx = nx+nghosts*2;
  ny = ny+nghosts*2;

  // calculate cell center values
  int i,j;
  DO2D(j,1,ny,i,1,nx,{
    xc(i,j) = (i-(1+nghosts))*dx + dx*0.5;
    yc(i,j) = (j-(1+nghosts))*dy + dy*0.5;
  });

  // calculate the cell nodal points
  DO2D(j,1,ny+1,i,1,nx+1,{
    xn(i,j) = (i-(1+nghosts)) * dx;
    yn(i,j) = (j-(1+nghosts)) * dy;
  });
}

/**********************
 * 3D Meshing Routine *
 **********************/
void mesher3D(const double& lx,const double& ly,const double& lz, 
              const int& nx,const int& ny,const int& nz,
              FMatrix<double>& xc,FMatrix<double>& yc,FMatrix<double>& zc,
              FMatrix<double>& xn,FMatrix<double>& yn,FMatrix<double>& zn,
              double& dx,double& dy,double& dz){
  // get dx and dy
  dx = lx/(double)nx;
  dy = ly/(double)ny;
  dz = lz/(double)nz;

  // add number of ghost cells to the domain
  int nghosts = 2; // 2 ghosts to support 2nd-order spatial discretization
  nx = nx+nghosts*2;
  ny = ny+nghosts*2;
  nz = nz+nghosts*2;

  // calculate cell center values
  int i,j,k;
  DO3D(k,1,nz,j,1,ny,i,1,nx,{
    xc(i,j,k) = (i-(1+nghosts)*dx) + dx*0.5;
    yc(i,j,k) = (j-(1+nghosts)*dy) + dy*0.5;
    zc(i,j,k) = (k-(1+nghosts)*dz) + dz*0.5;
  });

  // calculate the cell nodal points
  DO3D(k,1,nz+1j,1,ny+1,i,1,nx+1,{
    xn(i,j,k) = (i-(1+nghosts)) * dx;
    yn(i,j,k) = (j-(1+nghosts)) * dy;
    zn(i,j,k) = (k-(1+nghosts)) * dz;
  });
}
