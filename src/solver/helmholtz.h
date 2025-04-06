#ifndef HELMHOLTZ_H
#define HELMHOLTZ_H

#include <matar.h>
#include <BCs.h>

namespace SOR {

void helmholtz(const double& omega,
               mtr::FMatrix<double>& ustar,
               mtr::FMatrix<double>& rhs,
               BC::bcTags bcTags,
               const double& dx,
               const double& dy,
               const double& dt,
               const double& nu);

void helmholtz_eigen(const double& omega,
    mtr::FMatrix<double>& ustar,
    mtr::FMatrix<double>& rhs,
    BC::bcTags bcTags,
    const double& dx,
    const double& dy,
    const double& dt,
    const double& nu);
} // end namespace SOR
#endif // end PRESSURE_H
