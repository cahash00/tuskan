#ifndef PRESSURE_H
#define PRESSURE_H

#include <matar.h>
#include <BCs.h>

namespace psolve {

void SOR(const double& omega,
         mtr::FMatrix<double>& p,
         mtr::FMatrix<double>& ustar,
         mtr::FMatrix<double>& vstar,
         const double& dx,
         const double& dy,
         const double& dt,
         mtr::FMatrix<double>& rho,
         BC::bcTags bcTags);

}
#endif // end PRESSURE_H
