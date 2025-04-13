#ifndef PRESSURE_H
#define PRESSURE_H

#include <matar.h>
#include <BCs.h>
#include <types.h>

namespace psolve {

void SOR(const double& omega,
         fmat<double>& p,
         fmat<double>& ustar,
         fmat<double>& vstar,
         const double& dx,
         const double& dy,
         const double& dt,
         fmat<double>& rho,
         BC::bcTags bcTags);

}
#endif // end PRESSURE_H
