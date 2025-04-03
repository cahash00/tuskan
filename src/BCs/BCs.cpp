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

namespace BC {
/******************************************************************************/
void set_BCs() {

} // end set_BCs
/******************************************************************************/
void bc_noslip(FMatrix<double>& u) {
  // ... x-direction noslip
  // for (int i = istr; i <= iend; i++) {
  //   u(i,jstr-1) = -u(i,jstr);
  //   u(i,jend+1) = -u(i,jend);
  // }
  // ... y-direction noslip
  for (int j = jstr; j <= jend; j++) {
    u(istr-1,j) = -u(istr,j);
    u(iend+1,j) = -u(iend,j);
  }
} // end bc_noslip
/******************************************************************************/
void bc_periodic(FMatrix<double>& u) {
  // ... x-direction periodic
  // for (int j = jstr-nghosts; j <= jend+nghosts; j++) {
  //   u(istr-1,j) = u(iend,j);
  //   u(iend+1,j) = u(istr,j);
  // }
  // ... y-direction periodic
  for (int i = istr-nghosts; i <= iend+nghosts; i++) {
    u(i,jstr-1) = u(i,jend);
    u(i,jend+1) = u(i,jstr);
  }
} // end bc_periodic
} // end namespace BC
/******************************************************************************/
