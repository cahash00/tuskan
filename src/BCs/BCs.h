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

/******************************************************************************/
void bc_noslip(FMatrix<double>& u);
void bc_periodic(FMatrix<double>& u);
/******************************************************************************/
#endif // BCS_H
