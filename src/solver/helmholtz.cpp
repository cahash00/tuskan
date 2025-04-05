#include <cmath>
#include <matar.h>
#include <params.h>
#include <solver.h>
#include <BCs.h>
#include <pprint.hpp>

using namespace std;

namespace SOR {

void helmholtz(const double& omega,
               mtr::FMatrix<double>& ustar,
               mtr::FMatrix<double>& rhs,
               BC::bcTags bcTags,
               const double& dx,
               const double& dy,
               const double& dt,
               const double& nu) {
  
  // ... set values
  const double dx2 = dx*dx;
  const double dy2 = dy*dy;
  const double coeff = 2.0*nu/dt;
  const double denom = 2.0*(1.0/dx2+1.0/dy2) + coeff;

  // ... initialize pressure array
  mtr::FMatrix<double> u1(ustar.dims(1),ustar.dims(2));
  mtr::FMatrix<double> u2(ustar.dims(1),ustar.dims(2));

  // ... start the jacobi iterations
  // u1.set_values(0.001);
  // u2.set_values(0.001);

  // ... set initial guess to previous solution
  for (int j = jstr; j <= jend; j++) {
    for (int i = istr; i <= iend; i++) {
      u1(i,j) = ustar(i,j);
      u2(i,j) = ustar(i,j);
    }
  }
  

  // ... set pressure gradient boundary conditions
  int niter = 2000;
  for (int n = 0; n <= niter; n++) {
    BC::update_BCs(bcTags,u1);
    
    // ... loop over the domain
    for (int j = jstr; j <= jend; j++) {
      for (int i = istr; i <= iend; i++) {
        // ... variable dx and dy
        double term1 = (u1(i+1,j) + u2(i-1,j))/dx2 
                     + (u1(i,j+1) + u2(i,j-1))/dy2;

        u2(i,j) = (1.0-omega)*u1(i,j) + omega/denom*(term1 - rhs(i,j));
      } // end i-loop
    } // end j-loop

    // update the boundary conditions
    BC::update_BCs(bcTags,u2);

    double res = L2NORM(u1,u2);

    for (int j = jstr; j <= jend; j++) {
      for (int i = istr; i <= iend; i++) {
        u1(i,j) = u2(i,j);
      }
    }

    if (res < 2.0e-6 && n > 100) {
      break;
    } else if (n == niter) {
      cout << res << endl;
      // throw runtime_error("ERROR: SOR solver did not converge.");
    }
  } // end Jacobi loop
  
  // ... update the pressure with the SOR solution to pass back
  for (int j = jstr; j <= jend; j++) {
    for (int i = istr; i <= iend; i++) {
      ustar(i,j) = u2(i,j);
    }
  }

} // end helmholtz

} // end namespace SOR
