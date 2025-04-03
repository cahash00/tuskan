/**
 * @file BCs.h
 * @brief Boundary Condiitions file
 * @author Caleb Hash
 * @date 2025-02-25
 */
#ifndef BCS_H
#define BCS_H

#include <cmath>
#include <matar.h>
#include <params.h>

using namespace std;
using namespace mtr;

namespace BC {
/******************************************************************************/
void bc_noslip(FMatrix<double>& u);
void bc_periodic(FMatrix<double>& u);
/******************************************************************************/
} // end namespace BC
#endif // BCS_H
