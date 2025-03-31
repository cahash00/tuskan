
#include <cmath>
#include <matar.h>
#include <params.h>
#include <solver.h>
#include <BCs.h>
#include <pprint.hpp>

using namespace std;

namespace psolve {

void Jacobi(mtr::FMatrix<double>& p,
            mtr::FMatrix<double>& ustar,
            mtr::FMatrix<double>& vstar,
            const double& dx,
            const double& dy,
            const double& dt,
            const double& rho,
            int nx,
            int ny) {
  pprint::PrettyPrinter printer;
  // ... initialize pressure array
  double dx2 = dx*dx;
  double dy2 = dy*dy;
  
  mtr::FMatrix<double> p1(p.dims(1),p.dims(2));
  mtr::FMatrix<double> p2(p.dims(1),p.dims(2));
  
  // ... start the jacobi iterations
  int jiter = 5000;
  p1.set_values(0.99);
  
  for (int n = 0; n <= jiter; n++) {
    for (int i = istr; i <= iend; i++) {
      p1(i,jend) = p1(i,jend-1);
      p1(i,jstr-1) = p1(i,jstr);
    }
    // ... loop over the domain
    for (int j = jstr; j <= jend-1; j++) {
      p1(istr-1,j) = 1.0;
      p1(iend,j) = 1.0 + -0.3*dx*(nx);
      
      for (int i = istr; i <= iend-1; i++) {
        // ... variable dx and dy
        double term2 = rho / dt * ((ustar(i+1,j)-ustar(i,j))/dx
                                 + (vstar(i,j+1)-vstar(i,j))/dy);
        double term1 = (p1(i+1,j) + p1(i-1,j))*dy2 
                     + (p1(i,j+1) + p1(i,j-1))*dx2;
        p2(i,j) = 1.0/(2.0*(dy2+dx2))*(term1 - dx2*dy2*term2);
      } // end i-loopp
      p2(istr-1,j) = 1.0;
      p2(iend,j) = 1.0 + -0.3*dx*(nx);
    } // end j-loop
    for (int i = istr; i <= iend; i++) {
      p2(i,jend) = p2(i,jend-1);
      p2(i,jstr-1) = p2(i,jstr);
    }
    
    // double res = L2NORM(p1,p2,nx*ny);
    double res = 0.0;
    for (int j = jstr-1; j <= jend; j++) {
      for (int i = istr-1; i <= iend; i++) {
        double d = p2(i,j) - p1(i,j);
        res += d;
      }
    }
    
    if (res < 5.0e-4 && n > 100) {
      break;
    } else if (n == jiter) {
      cout << res << endl;
      throw runtime_error("ERROR: Jacobi solver did not converge.");
    }
    // update the solution

    for (int j = jstr-nghosts; j <= jend; j++) {
      for (int i = istr-nghosts; i <= iend; i++) {
        p1(i,j) = p2(i,j);
      }
    }

  } // end Jacobi loop
  for (int j = jstr-nghosts; j <= jend; j++) {
    for (int i = istr-nghosts; i <= iend; i++) {
      p(i,j) = p2(i,j);
    }
  }
}

void SOR(mtr::FMatrix<double>& p,
         mtr::FMatrix<double>& ustar,
         mtr::FMatrix<double>& vstar,
         const double& dx,
         const double& dy,
         const double& dt,
         const double& rho,
         int nx,
         int ny) {
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
  double omega = 1.7;

  for (int n = 0; n <= jiter; n++) {

    // ... set pressure gradient boundary conditions
    for (int j = jstr; j <= jend; j++) {
      for (int i = istr; i <= iend; i++) {
        // set pressure gradient for x-direction
        p2(istr-1,j) = 1.0;
        p2(iend+1,j) = 1.0 + -0.3*dx*(nx);
        p1(i,jend+1) = p1(i,jend);
        p1(i,jstr-1) = p1(i,jstr);
        // set pressure gradient for y-direction
        // p1(i,jstr-1) = 1.0;
        // p1(i,jend+1) = 1.0 + -0.3*dx*(nx);
        // p1(iend+1,j) = p1(iend,j);
        // p1(istr-1,j) = p1(istr,j);
        
      }
    }
    
    // ... loop over the domain
    for (int j = jstr; j <= jend; j++) {
      p1(istr-1,j) = 1.0;
      p1(iend+1,j) = 1.0 + -0.3*dx*(nx);

      for (int i = istr; i <= iend; i++) {
        // ... variable dx and dy
        double term2 = rho / dt * ((ustar(i+1,j)-ustar(i,j))/dx
                                 + (vstar(i,j+1)-vstar(i,j))/dy);

        double term1 = (p1(i+1,j) + p2(i-1,j))*dy2 
                     + (p1(i,j+1) + p2(i,j-1))*dx2;

        p2(i,j) = (1.0-omega)*p1(i,j) + omega/(2.0*(dy2+dx2))*(term1 - dx2*dy2*term2);
      } // end i-loop
    } // end j-loop
    for (int i = istr-1; i <= iend+1; i++) {
      p2(i,jend+1) = p2(i,jend);
      p2(i,jstr-1) = p2(i,jstr);
    }

    // double res = L2NORM(p1,p2,nx*ny);
    double res = 0.0;
    for (int j = jstr; j <= jend; j++) {
      for (int i = istr; i <= iend; i++) {
        double d = p2(i,j) - p1(i,j);
        res += d;
      }
    }

    if (res < 5.0e-6 && n > 100) {
      break;
    } else if (n == jiter) {
      cout << res << endl;
      throw runtime_error("ERROR: SOR solver did not converge.");
    }
    // update the solution

    for (int j = jstr-nghosts; j <= jend+nghosts; j++) {
      for (int i = istr-nghosts; i <= iend-nghosts; i++) {
        p1(i,j) = p2(i,j);
      }
    }

  } // end Jacobi loop
  for (int j = jstr-nghosts; j <= jend+nghosts; j++) {
    for (int i = istr-nghosts; i <= iend+nghosts; i++) {
      p(i,j) = p2(i,j);
    }
  }

}

}
