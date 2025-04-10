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
namespace mesh {
void mesher2D(const double& lx, const double& ly, 
              FMatrix<double>& xc, FMatrix<double>& yc,
              FMatrix<double>& xu, FMatrix<double>& yu,
              FMatrix<double>& xv, FMatrix<double>& yv,
              double& dx, double& dy);
} // end namespace mesh
#endif // MESH_H
