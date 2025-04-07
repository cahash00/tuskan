

#include <matar.h>
#include <cmath>
#include <params.h>
#include <vector>
#include <BCs.h>

namespace levset {

void get_phi(mtr::FMatrix<double>& phi,
             const mtr::FMatrix<double>& xc,
             const mtr::FMatrix<double>& yc,
             const double& xd,
             const double& yd,
             const double& r_drop) {
  for (int j = jstr-1; j <= jend+1; j++) {
    for (int i = istr-1; i <= iend+1; i++) {
      phi(i,j) = sqrt(pow((xc(i,j)-xd),2) + pow((yc(i,j)-yd),2))-r_drop;
    }
  }
  
} // end get_phi

void heaviside(const double& M,
               const double& h,
               const mtr::FMatrix<double>& phi,
               mtr::FMatrix<double>& heavi) {
  const double Mh = M*h;
  for (int j = jstr-1; j <= jend+1; j++) {
    for (int i = istr-1; i <= iend+1; i++) {
      if (phi(i,j) < -Mh) {
        heavi(i,j) = 0.0;
      } else if (phi(i,j) > Mh) {
        heavi(i,j) = 1.0;
      } else {
        heavi(i,j) = 0.5*(1.0 + phi(i,j)/Mh + 1.0/M_PI*sin(M_PI*phi(i,j)/Mh));
      }
    }
  }
} // end marker

/******************************************************************************/
std::vector<double> get_phiHalf(const int i,
                                const int j,
                                const mtr::FMatrix<double>& u,
                                const mtr::FMatrix<double>& v,
                                const mtr::FMatrix<double>& phi) {
  std::vector<double> phiHalf(2,0.0);
  
  /**
     * i-direction: calculate D+ and D- (x-direction)
     */
    double dp1 = phi(i+1,j) - phi(i,j);   // D+ for phi(i,j)
    double dm1 = phi(i,j) - phi(i-1,j);   // D- for phi(i,j)

    // get interface values using simple 1st order CD
    if (u(i,j) > 0.0) {
        phiHalf[0] = phi(i,j) + 0.5*dm1;  // Backward difference
    }
    else if (u(i,j) < 0.0) {
        phiHalf[0] = phi(i+1,j) - 0.5*dp1;  // Forward difference
    }
    else {
        phiHalf[0] = phi(i,j);  // If velocity is zero, take current value
    }

    /**
     * j-direction: calculate D+ and D- (y-direction)
     */
    double dp1_y = phi(i,j+1) - phi(i,j);   // D+ for phi(i,j)
    double dm1_y = phi(i,j) - phi(i,j-1);   // D- for phi(i,j)

    // get interface values using simple 1st order CD
    if (v(i,j) > 0.0) {
        phiHalf[1] = phi(i,j) + 0.5*dm1_y;  // Backward difference
    }
    else if (v(i,j) < 0.0) {
        phiHalf[1] = phi(i,j+1) - 0.5*dp1_y;  // Forward difference
    }
    else {
        phiHalf[1] = phi(i,j);  // If velocity is zero, take current value
    }

    return phiHalf;

} // end get_phiHalf
/******************************************************************************/
double get_Lphi(const int i,
                const int j,
                const mtr::FMatrix<double>& u, 
                const mtr::FMatrix<double>& v, 
                const double dx,
                const double dy,
                mtr::FMatrix<double>& phi) {
  // ... get the cell-centered velocities
  const double ucell = 0.5*(u(i-1,j)+u(i,j));
  const double vcell = 0.5*(v(i,j-1)+v(i,j));

  // ... get the cell face phi values based on switch function
  std::vector<double> phiHalf = get_phiHalf(i,j,u,v,phi);
  
  // ... Calculate the spatial derivative using phiHalf values
  const double dphidx = (phiHalf[0] - phi(i, j)) / dx;
  const double dphidy = (phiHalf[1] - phi(i, j)) / dy;

  // ... Advection term calculation
  const double Lphi = - ucell * dphidx - vcell * dphidy;
  return Lphi;
}
/******************************************************************************/
void weno(const BC::bcTags bcTags,
          const double dx,
          const double dy,
          const double dt,
          const mtr::FMatrix<double>& u,
          const mtr::FMatrix<double>& v,
          mtr::FMatrix<double>& phi) {
  // ... initialize phi arrays
  mtr::FMatrix<double> phiStar(phi.dims(1),phi.dims(2));
  mtr::FMatrix<double> phi2(phi.dims(1),phi.dims(2));

  BC::update_BCs_phi(bcTags,phi);
  // ... calculate phi*
  for (int j = jstr; j <= jend-1; j++) {
    for (int i = istr; i <= iend-1; i++) {
      double Lphi = get_Lphi(i,j,u,v,dx,dy,phi);
      phiStar(i,j) = phi(i,j) + dt*Lphi;
    }
  }
  BC::update_BCs_phi(bcTags,phiStar);

  // ... correct phi*
  for (int j = jstr; j <= jend-1; j++) {
    for (int i = istr; i <= iend-1; i++) {
      double LphiN = get_Lphi(i,j,u,v,dx,dy,phi);
      double LphiStar = get_Lphi(i,j,u,v,dx,dy,phiStar);
      phi2(i,j) = phi(i,j) + dt*0.5*(LphiN+LphiStar);
    }
  }
  BC::update_BCs_phi(bcTags,phi2);

  // ... update the phi solution
  for (int j = jstr-1; j <= jend; j++) {
    for (int i = istr-1; i <= iend; i++) {
      phi(i,j) = phi2(i,j);
    }
  }
} // end weno
double surfaceTension(const int i, const int j,
                      const double dx, const double dy,
                      const mtr::FMatrix<double>& phi) {
    // Compute gradients of phi using central differences
    double dphidx = (phi(i+1, j) - phi(i-1, j)) / (2.0 * dx);
    double dphidy = (phi(i, j+1) - phi(i, j-1)) / (2.0 * dy);

    // Compute magnitude of the gradient
    double grad_mag = std::sqrt(dphidx * dphidx + dphidy * dphidy + 1e-12);

    // Normalized components of the gradient
    double nx = dphidx / grad_mag;
    double ny = dphidy / grad_mag;

    // Compute second derivatives using central differences
    double d2phidx2 = (phi(i+1, j) - 2.0 * phi(i, j) + phi(i-1, j)) / (dx * dx);
    double d2phidy2 = (phi(i, j+1) - 2.0 * phi(i, j) + phi(i, j-1)) / (dy * dy);
    double d2phidxdy = (phi(i+1, j+1) - phi(i+1, j-1) 
                      - phi(i-1, j+1) + phi(i-1, j-1)) / (4.0 * dx * dy);

    // Divergence of normalized gradient
    double curvature = (d2phidx2 * ny * ny - 2.0 * d2phidxdy * nx * ny 
                      + d2phidy2 * nx * nx) / grad_mag;

    return curvature;
}

void reinitialize(BC::bcTags bcTags,
                  const double& dx,
                  const double& dy,
                  const double& dtau,
                  const int isteps,
                  mtr::FMatrix<double>& phi) {
  const double eps=1.0e-8;
  mtr::FMatrix<double> sign0(phi.dims(1),phi.dims(2));
  mtr::FMatrix<double> phi0(phi.dims(1),phi.dims(2));
  mtr::FMatrix<double> phi2(phi.dims(1),phi.dims(2));
  for (int j = jstr-1; j <= jend+1; j++) {
    for (int i = istr-1; i <= iend+1; i++) {
      phi2(i,j) = phi(i,j);
      phi0(i,j) = phi(i,j);
      double gphi0 = sqrt(pow((phi(i+1,j) - phi(i-1,j)) / (2.0 * dx),2)
                        + pow((phi(i,j+1) - phi(i,j-1)) / (2.0 * dy),2));
      sign0(i,j) = phi(i,j) / sqrt(phi(i,j)*phi(i,j)+gphi0*gphi0*dx*dx+eps*eps);
    }
  }
  

  // ... update the phi matrix
  for (int n = 0; n < isteps; n++) {
    for (int j = jstr; j <= jend-1; j++) {
      for (int i = istr; i <= iend-1; i++) {
        double a = (phi(i,j)-phi(i-1,j))/dx;
        double b = (phi(i+1,j)-phi(i,j))/dx;
        double c = (phi(i,j)-phi(i,j-1))/dy;
        double d = (phi(i,j+1)-phi(i,j))/dy;
        double ap = max(a,0.0);
        double am = min(a,0.0);
        double bp = max(b,0.0);
        double bm = min(b,0.0);
        double cp = max(c,0.0);
        double cm = min(c,0.0);
        double dp = max(d,0.0);
        double dm = min(d,0.0);
        double sp = max(sign0(i,j),0.0);
        double sm = min(sign0(i,j),0.0);
        double H = sp*(sqrt(max(ap*ap,bm*bm) + max(cp*cp,dm*dm))-1.0)
                 + sm*(sqrt(max(am*am,bp*bp) + max(cm*cm,dp*dp))-1.0);
        //> subcell fix
        //  Russo and Smereka (2000)
        if (phi0(i+1,j)*phi0(i,j) < 0.0 ||
            phi0(i-1,j)*phi0(i,j) < 0.0 ||
            phi0(i,j+1)*phi0(i,j) < 0.0 || 
            phi0(i,j-1)*phi0(i,j) < 0.0) {
          double D = 2.0*dx*phi0(i,j) / sqrt( pow(phi0(i+1,j)-phi0(i-1,j),2) 
                                            + pow(phi0(i,j+1)-phi0(i,j-1),2) );
          phi2(i,j) = phi(i,j) - dtau*(sign0(i,j)*abs(phi(i,j))-D);
        } else {
          phi2(i,j) = phi(i,j) - dtau*H;
        }
      }
    }
    // BC::update_BCs_phi(bcTags,phi2);
    for (int j = jstr; j <= jend-1; j++) {
      for (int i = istr; i <= iend-1; i++) {
        phi(i,j) = phi2(i,j);
      }
    }
  }
} // end reinitialize

} // end namespace levset
