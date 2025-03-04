/**
 * @file solver.cpp
 * @brief Solver code file
 * @author Caleb Hash
 * @date 2025-02-25
 */
#include <cmath>
#include <matar.h>
#include <params.h>

using namespace std;
using namespace mtr;

/******************************************************************************/
void bc_noslip(FMatrix<double>& q) {
  DO_LOOP(i,istr-nghosts,iend+nghosts,{
    q(1,i,jstr-1) = -q(1,i,jstr);
    q(1,i,jend+1) = -q(1,i,jend);
  });
}
/******************************************************************************/
void bc_periodic(FMatrix<double>& q) {
  DO_LOOP(j,jstr-nghosts,jend+nghosts,{
    q(1,istr-1,j) = q(1,iend,j);
    q(1,iend+1,j) = q(1,istr,j);
  });
}
/******************************************************************************/
