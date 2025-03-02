/**
 * @file mesh.h
 * @brief Function declarations for all things mesh related. 
 * @author Caleb Hash
 * @date 2025-02-28
 */

#ifndef MESH_H
#define MESH_H

#include <matar.h>

using namespace mtr;
/**
 * @brief 2D cartesian mesh builder
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
              FMatrix<double>& xc, FMatrix<double>& yc,
              FMatrix<double>& xn, FMatrix<double>& yn,
              double& dx, double& dy);

/**
 * @brief 3D cartesian mesh builder
 * @param[in] lx length in the x-direction
 * @param[in] ly length in the y-direction
 * @param[in] lz length in the y-direction
 * @param[in] nx number of cells in the x-direction
 * @param[in] ny number of cells in the y-direction
 * @param[in] nz number of cells in the y-direction
 * @param[inout] xc x coordinates for cell centers
 * @param[inout] yc y coordinates for cell centers
 * @param[inout] zc y coordinates for cell centers
 * @param[inout] xn x coordinates for nodes
 * @param[inout] yn y coordinates for nodes
 * @param[inout] zn y coordinates for nodes
 * @param[inout] dx cell width in the x-direction
 * @param[inout] dy cell width in the y-direction
 * @param[inout] dz cell width in the y-direction
 */
void mesher3D(const double& lx,const double& ly,const double& lz, 
              int& nx,int& ny,int& nz,
              FMatrix<double>& xc,FMatrix<double>& yc,FMatrix<double>& zc,
              FMatrix<double>& xn,FMatrix<double>& yn,FMatrix<double>& zn,
              double& dx,double& dy,double& dz);

#endif // MESH_H
