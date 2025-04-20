
#include <cmath>
#include <matar.h>
#include <params.h>
#include <solver.h>
#include <BCs.h>
#include <pprint.hpp>

using namespace std;

namespace psolve {

void SOR(const double& omega,
         mtr::FMatrix<double>& p,
         mtr::FMatrix<double>& ustar,
         mtr::FMatrix<double>& vstar,
         const double& dx,
         const double& dy,
         const double& dt,
         mtr::FMatrix<double>& rho,
         BC::bcTags bcTags) {
  
  // ... initialize pressure array
  const double dx2 = dx*dx;
  const double dy2 = dy*dy;
  // const double coeff = 2.0*(1.0/dx2+1.0/dy2);

  mtr::FMatrix<double> p1(p.dims(1),p.dims(2));
  mtr::FMatrix<double> p2(p.dims(1),p.dims(2));


  // ... set pressure gradient boundary conditions
  BC::update_BCs(bcTags,ustar,vstar,p);
  // ... start the jacobi iterations
  int jiter = 5000;
  for (int j = jstr-1; j <= jend+1; j++) {
    for (int i = istr-1; i <= iend+1; i++) {
      p1(i,j) = p(i,j);
      p2(i,j) = p(i,j);
    }
  }

  // ... set pressure gradient boundary conditions
  for (int n = 0; n <= jiter; n++) {
    BC::update_BCs(bcTags,ustar,vstar,p1);
    
    // ... loop over the domain
    for (int j = jstr; j <= jend; j++) {
      for (int i = istr; i <= iend; i++) {
        // Compute average densities at cell faces
        // - this is using harmonic averaging
        double rho_e = 2.0 * rho(i,j) * rho(i+1,j) / (rho(i,j) + rho(i+1,j));
        double rho_w = 2.0 * rho(i,j) * rho(i-1,j) / (rho(i,j) + rho(i-1,j));
        double rho_n = 2.0 * rho(i,j) * rho(i,j+1) / (rho(i,j) + rho(i,j+1));
        double rho_s = 2.0 * rho(i,j) * rho(i,j-1) / (rho(i,j) + rho(i,j-1));
        
        // ... variable dx and dy
        double term2 = 1.0 / dt * ((ustar(i,j)-ustar(i-1,j))/dx
                                 + (vstar(i,j)-vstar(i,j-1))/dy);

        double term1 = (1.0 / (rho_e * dx2)) * p1(i+1, j)
                     + (1.0 / (rho_w * dx2)) * p2(i-1, j)
                     + (1.0 / (rho_n * dy2)) * p1(i, j+1)
                     + (1.0 / (rho_s * dy2)) * p2(i, j-1);
        double coeff = (1.0 / (rho_e * dx2))
                     + (1.0 / (rho_w * dx2))
                     + (1.0 / (rho_n * dy2))
                     + (1.0 / (rho_s * dy2));

        p2(i,j) = (1.0-omega)*p1(i,j) 
                + omega/coeff*(term1 - term2);

      } // end i-loop
    } // end j-loop
    BC::update_BCs(bcTags,ustar,vstar,p2);

    for (int j = jstr-1; j <= jend+1; j++) {
      for (int i = istr-1; i <= iend+1; i++) {
        p1(i,j) = p2(i,j);
      }
    }


    double res = L2NORM(p1,p2);

    if (res < 2.0e-6 && n > 20) {
      break;
    } else if (n == jiter) {
      cout << res << endl;
      // throw runtime_error("ERROR: SOR solver did not converge.");
    }
  } // end Jacobi loop
  for (int j = jstr-1; j <= jend+1; j++) {
    for (int i = istr-1; i <= iend+1; i++) {
      p(i,j) = p2(i,j);
    }
  }

}

}
