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
void bc_noslip(FMatrix<double>& u) {
  for (int j = jstr; j <= jend; j++) {
    for (int i = istr; i <= iend; i++) {
      u(i,jstr-1) = -u(i,jstr);
      u(i,jend+1) = -u(i,jend);
    }
  }
}
/******************************************************************************/
void bc_periodic(FMatrix<double>& u) {
  for (int j = jstr-nghosts; j <= jend+nghosts; j++) {
    u(istr-1,j) = u(iend,j);
    u(iend+1,j) = u(istr,j);
  }
}
/******************************************************************************/
