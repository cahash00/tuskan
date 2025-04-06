
#include <cmath>
#include <matar.h>
#include <params.h>
#include <solver.h>
#include <BCs.h>
#include <pprint.hpp>

using namespace std;

namespace psolve {

void update_wall(mtr::FMatrix<double>& p) {
  for (int i = istr-1; i <= iend+1; i++) {
    // p(i,jend) = p(i,jend-1);
    // p(i,jstr-1) = p(i,jstr);
    p(i,jstr-1) = 1.0;
    p(i,jend) = 1.0 + -0.3*0.00066*(ny);
    p(i,jend+1) = 1.0 + -0.3*0.00066*(ny);
  }
  for (int j = jstr-1; j <= jend+1; j++) {
    p(istr-1,j) = p(istr,j);
    p(iend,j) = p(iend-1,j);
  }
}
void SOR(const double& omega,
         mtr::FMatrix<double>& p,
         mtr::FMatrix<double>& ustar,
         mtr::FMatrix<double>& vstar,
         const double& dx,
         const double& dy,
         const double& dt,
         mtr::FMatrix<double>& rho) {
  pprint::PrettyPrinter printer;
  // ... initialize pressure array
  const double dx2 = dx*dx;
  const double dy2 = dy*dy;
  const double coeff = 2.0*(1.0/dx2+1.0/dy2);

  mtr::FMatrix<double> p1(p.dims(1),p.dims(2));
  mtr::FMatrix<double> p2(p.dims(1),p.dims(2));

  // ... start the jacobi iterations
  int jiter = 2000;
  for (int j = jstr-1; j <= jend+1; j++) {
    for (int i = istr-1; i <= iend+1; i++) {
      p1(i,j) = p(i,j);
      p2(i,j) = p(i,j);
    }
  }

  // ... set pressure gradient boundary conditions
  update_wall(p1);
  update_wall(p2);

  // ... set pressure gradient boundary conditions
  for (int n = 0; n <= jiter; n++) {
    update_wall(p1);
    
    // ... loop over the domain
    for (int j = jstr; j <= jend-1; j++) {
      for (int i = istr; i <= iend-1; i++) {
        // ... variable dx and dy
        double term2 = rho(i,j) / dt * ((ustar(i+1,j)-ustar(i,j))/dx
                                      + (vstar(i,j+1)-vstar(i,j))/dy);

        double term1 = (p1(i+1,j) + p2(i-1,j))/dx2 
                     + (p1(i,j+1) + p2(i,j-1))/dy2;

        p2(i,j) = (1.0-omega)*p1(i,j) 
                + omega/coeff*(term1 - term2);

      } // end i-loop
    } // end j-loop
    update_wall(p2);

    for (int j = jstr-1; j <= jend; j++) {
      for (int i = istr-1; i <= iend; i++) {
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
  for (int j = jstr-1; j <= jend; j++) {
    for (int i = istr-1; i <= iend; i++) {
      p(i,j) = p2(i,j);
    }
  }

}

}
