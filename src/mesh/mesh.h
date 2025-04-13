/**
 * @file mesh.h
 * @brief Function declarations for all things mesh related. 
 * @author Caleb Hash
 * @date 2025-02-28
 */

#ifndef MESH_H
#define MESH_H

#include <types.h>
#include <communication.h>

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
void mesher2D(fmat<double>& xc, fmat<double>& yc,
              fmat<double>& xn, fmat<double>& yn,
              const double& dx, const double& dy);
} // end namespace mesh
#endif // MESH_H
