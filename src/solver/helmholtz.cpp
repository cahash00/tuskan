#include <cmath>
#include <matar.h>
#include <params.h>
#include <solver.h>
#include <BCs.h>
#include <pprint.hpp>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

using namespace std;

namespace SOR {

void helmholtz_eigen(const double& omega,
    mtr::FMatrix<double>& ustar,
    mtr::FMatrix<double>& rhs,
    BC::bcTags bcTags,
    const double& dx,
    const double& dy,
    const double& dt,
    const double& nu) {
    const double dx2 = dx*dx;
    const double dy2 = dy*dy;
    const double coeff = 2.0*nu/dt;
    const double denom = 2.0*(1.0/dx2+1.0/dy2) + coeff;
    int nx1 = nx+3;
    int ny1 = ny+3;
    Eigen::SparseMatrix<double> A(nx1 * ny1, nx1 * ny1);
    // Diagonal entries (the "1 + 2alpha" part)
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(nx * ny * 5); // Pre-allocate space

    for (int j = jstr-1; j <= jend+1; j++) {
      for (int i = istr-1; i <= iend+1; i++) {
        int idx = (j-1)*ny1+(i-1);

        // Main diagonal: (1 + 2alpha(dx^-2 + dy^-2))
        tripletList.push_back(Eigen::Triplet<double>(idx, idx, coeff));

        // Neighbors: (alpha * (dx^-2 or dy^-2))
        if (i < nx - 1) {  // right neighbor
          tripletList.push_back(Eigen::Triplet<double>(idx, idx+1, -coeff/dx2));
        }
        if (i > 0) {  // left neighbor
          tripletList.push_back(Eigen::Triplet<double>(idx, idx -1, -coeff/dx2));
        }
        if (j < ny - 1) {  // upper neighbor
          tripletList.push_back(Eigen::Triplet<double>(idx, idx+nx, -coeff/dy2));
        }
        if (j > 0) {  // lower neighbor
          tripletList.push_back(Eigen::Triplet<double>(idx, idx-nx, -coeff/dy2));
        }
      }
    }
    // Assemble the matrix
    A.setFromTriplets(tripletList.begin(), tripletList.end());

    // Create the right-hand side vector b (forcing term, source term)
    Eigen::VectorXd b(nx * ny);
    // Solve the system
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    auto test = solver.solve(b);
    for (int j = jstr-1; j <= jend+1; j++) {
      for (int i = istr-1; i <= iend+1; i++) {
        ustar(i,j) = test((j-1)*ny1+(i-1));
      }
    }
    
  
} // end helmholtz_eigen

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
    BC::update_BCs(bcTags,u2);
    
    // ... loop over the domain
    for (int j = jstr; j <= jend; j++) {
      for (int i = istr; i <= iend; i++) {
        // ... variable dx and dy
        double term1 = (u1(i+1,j) + u2(i-1,j))/dx2 
                     + (u1(i,j+1) + u2(i,j-1))/dy2;

        u2(i,j) = (1.0-omega)*u1(i,j) + omega/denom*(term1 - rhs(i,j));
      } // end i-loop
    } // end j-loop
    // cout << u2(20,10) << endl;

    // update the boundary conditions
    BC::update_BCs(bcTags,u2);

    double res = L2NORM(u1,u2);

    for (int j = jstr; j <= jend; j++) {
      for (int i = istr; i <= iend; i++) {
        u1(i,j) = u2(i,j);
      }
    }

    if (res < 2.0e-6 && n > 10) {
      cout << n << endl;
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
