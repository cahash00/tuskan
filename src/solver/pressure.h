#ifndef PRESSURE_H
#define PRESSURE_H

#include <matar.h>

namespace SOR {

void pressure(const double& omega,
              mtr::FMatrix<double>& p,
              mtr::FMatrix<double>& ustar,
              mtr::FMatrix<double>& vstar,
              const double& dx,
              const double& dy,
              const double& dt,
              const double& rho);

}
#endif // end PRESSURE_H
