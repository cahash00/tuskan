#ifndef PRESSURE_H
#define PRESSURE_H

#include <matar.h>

namespace psolve {

void Jacobi(mtr::FMatrix<double>& p,
            mtr::FMatrix<double>& ustar,
            mtr::FMatrix<double>& vstar,
            const double& dx,
            const double& dy,
            const double& dt,
            const double& rho,
            int nx,
            int ny);

}
#endif // end PRESSURE_H
