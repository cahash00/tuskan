
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
         const double& rho) {
  pprint::PrettyPrinter printer;
  // ... initialize pressure array
  double dx2 = dx*dx;
  double dy2 = dy*dy;

  mtr::FMatrix<double> p1(p.dims(1),p.dims(2));
  mtr::FMatrix<double> p2(p.dims(1),p.dims(2));

  // ... start the jacobi iterations
  int jiter = 5000;
  p1.set_values(0.99);
  p2.set_values(0.99);

  // ... set pressure gradient boundary conditions
  for (int j = jstr-1; j <= jend+1; j++) {
    p1(istr-1,j) = 1.0;
    p1(iend,j) = 1.0 + -0.3*dx*(nx);
    p2(istr-1,j) = 1.0;
    p2(iend,j) = 1.0 + -0.3*dx*(nx);
  }
  for (int i = istr-1; i <= iend+1; i++) {
    p1(i,jend) = p1(i,jend-1);
    p1(i,jstr-1) = p1(i,jstr);
    p2(i,jend) = p2(i,jend-1);
    p2(i,jstr-1) = p2(i,jstr);
  }

  // ... set pressure gradient boundary conditions
  for (int n = 0; n <= jiter; n++) {
    
    // ... loop over the domain
    for (int j = jstr; j <= jend; j++) {
      for (int i = istr; i <= iend; i++) {
        // ... variable dx and dy
        double term2 = rho / dt * ((ustar(i+1,j)-ustar(i,j))/dx
                                 + (vstar(i,j+1)-vstar(i,j))/dy);

        double term1 = (p1(i+1,j) + p2(i-1,j))*dx2 
                     + (p1(i,j+1) + p2(i,j-1))*dy2;

        p2(i,j) = (1.0-omega)*p1(i,j) + omega/(2.0*(dy2+dx2))*(term1 - dx2*dy2*term2);
      } // end i-loop
    } // end j-loop
    // ... set pressure gradient boundary conditions
    for (int j = jstr; j <= jend; j++) {
      p2(istr,j) = 1.0;
      p2(iend,j) = 1.0 + -0.3*dx*(nx);
    }
    for (int i = istr; i <= iend; i++) {
      p2(i,jend) = p2(i,jend-1);
      p2(i,jstr-1) = p2(i,jstr);
    }

    for (int j = jstr; j <= jend; j++) {
      for (int i = istr; i <= iend; i++) {
        p1(i,j) = p2(i,j);
      }
    }

    // FIXME :: better to set p2 pressure BCs here? 

    double res = L2NORM(p1,p2,nx*ny);

    if (res < 2.0e-6 && n > 100) {
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
