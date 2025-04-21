

#include <matar.h>
#include <cmath>
#include <params.h>
#include <vector>
#include <BCs.h>

namespace levset {

/******************************************************************************/
void get_phi(mtr::FMatrix<double>& phi,
             const mtr::FMatrix<double>& xc,
             const mtr::FMatrix<double>& yc,
             const double& xd,
             const double& yd,
             const double& r_drop) {
  for (int j = jstr-1; j <= jend; j++) {
    for (int i = istr-1; i <= iend; i++) {
      phi(i,j) = sqrt(pow((xc(i,j)-xd),2) + pow((yc(i,j)-yd),2))-r_drop;
      double phi_water = xc(i,j)-3.5;
      double phi_air = -yc(i,j)+1.0+0.1*0.5*cos(2*M_PI*xc(i,j)/0.5);
      phi(i,j) = max(phi_water,phi_air);
      phi(i,j) = phi_air;
    }
  }
  
} // end get_phi

/******************************************************************************/
void heaviside(const double& M,
               const double& h,
               const mtr::FMatrix<double>& phi,
               mtr::FMatrix<double>& heavi) {
  const double Mh = M*h;
  for (int j = jstr-1; j <= jend; j++) {
    for (int i = istr-1; i <= iend; i++) {
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
void advecPhi(const BC::bcTags bcTags,
          const double dx,
          const double dy,
          const double dt,
          const mtr::FMatrix<double>& u,
          const mtr::FMatrix<double>& v,
          mtr::FMatrix<double>& phi) {
  // ... initialize phi arrays
  mtr::FMatrix<double> phiStar(phi.dims(1),phi.dims(2));
  mtr::FMatrix<double> phi2(phi.dims(1),phi.dims(2));

  BC::update_BCs_phi(bcTags,dx,dy,phi);
  // ... simple first-order Euler time discretization
  for (int j = jstr; j <= jend-1; j++) {
    for (int i = istr; i <= iend-1; i++) {
      double ucell=0.0,vcell=0.0;
      ucell = (u(i+1,j)+u(i,j))*0.5;
      vcell = (v(i,j+1)+v(i,j))*0.5;

      double phix =0.0,phiy=0.0;
      if (ucell > 0.0) {
        phix = (phi(i,j)-phi(i-1,j))/dx;
      } else if (ucell < 0.0) {
        phix = (phi(i+1,j)-phi(i,j))/dx;
      } else {
        phix = 0.0;
      }
      if (vcell > 0.0) {
        phiy = (phi(i,j)-phi(i,j-1))/dy;
      } else if (vcell < 0.0) {
        phiy = (phi(i,j+1)-phi(i,j))/dy;
      } else {
        phiy = 0.0;
      }
      
      phi2(i,j) = phi(i,j) - dt*(ucell*phix + vcell*phiy);
    }
  }
  
  BC::update_BCs_phi(bcTags,dx,dy,phi2);

  // ... update the phi solution
  for (int j = jstr; j <= jend; j++) {
    for (int i = istr; i <= iend; i++) {
      phi(i,j) = phi2(i,j);
    }
  }
} // end advecPhi

/******************************************************************************/
void surfaceTension(mtr::FMatrix<double>& Fx,
                    mtr::FMatrix<double>& Fy,
                    mtr::FMatrix<double>& phi,
                    mtr::FMatrix<double>& kappa,
                    const double Mh,
                    const double sigma,
                    const double dx,
                    const double dy) {
  std::vector<double> fst(2,0.0);
  double eps=1e-12;

  for (int j = jstr; j <= jend-1; j++) {
    for (int i = istr; i <= iend-1; i++) {
      double delta = 0.0;
      double phixx,phiyy,phixy,ny1,nx1;
      if (abs(phi(i,j)) < Mh) {
        // Compute gradients of phi using central differences
        double phix = (phi(i+1,j) - phi(i-1,j)) / (2.0*dx);
        double phiy = (phi(i,j+1) - phi(i,j-1)) / (2.0*dy);

        // Compute magnitude of the gradient
        double grad_mag = sqrt(phix*phix + phiy*phiy + eps);

        // Normalized components of the gradient
        nx1 = phix / grad_mag;
        ny1 = phiy / grad_mag;

        // Compute second derivatives using central differences
        phixx = (phi(i+1,j) - 2.0*phi(i,j) + phi(i-1,j))/(dx*dx);
        phiyy = (phi(i,j+1) - 2.0*phi(i,j) + phi(i,j-1))/(dy*dy);
        phixy = (phi(i+1,j+1) - phi(i+1,j-1)
                          - phi(i-1,j+1) + phi(i-1,j-1)) / (4.0*dx*dy);

        kappa(i,j) = (phixx*ny1*ny1 - 2.0*phixy*nx1*ny1
            + phiyy*nx1*nx1) / grad_mag;
        delta = 1.0/(2.0*Mh)*(1.0+cos(M_PI*phi(i,j)/Mh));
      } else {
        kappa(i,j) = 0.0;
        delta = 0.0;
      }
      Fx(i,j) = sigma*kappa(i,j)*delta*nx1;
      Fy(i,j) = sigma*kappa(i,j)*delta*ny1;
    }
  }
}

/******************************************************************************/
void reinitialize(BC::bcTags bcTags,
                  const double& dx,
                  const double& dy,
                  const double& dtau,
                  const int isteps,
                  mtr::FMatrix<double>& phi) {
  const double eps = 1e-12;
  const double Mh = 3.0 * dx;

  mtr::FMatrix<double> phi0(phi.dims(1), phi.dims(2));
  mtr::FMatrix<double> phi2(phi.dims(1), phi.dims(2));
  mtr::FMatrix<double> sign0(phi.dims(1),phi.dims(2));

  for (int j = jstr; j <= jend-1; ++j) {
    for (int i = istr; i <= iend-1; ++i) {
      phi0(i,j) = phi(i,j);
    }
  }

  BC::update_BCs_phi(bcTags, dx, dy, phi0);

  for (int j = jstr; j <= jend-1; ++j) {
    for (int i = istr; i <= iend-1; ++i) {
      double dphix = (phi0(i+1,j) - phi0(i-1,j)) / (2.0 * dx);
      double dphiy = (phi0(i,j+1) - phi0(i,j-1)) / (2.0 * dy);
      double denom = sqrt(phi0(i,j)*phi0(i,j) + (dphix*dphix + dphiy*dphiy)*dx*dx + eps);
      sign0(i,j) = phi0(i,j) / denom;
    }
  }

  for (int n = 0; n < isteps; ++n) {
    BC::update_BCs_phi(bcTags, dx, dy, phi);

    for (int j = jstr; j <= jend-1; ++j) {
      for (int i = istr; i <= iend-1; ++i) {
        // ... Compute Godunov upwind gradients
        double a = (phi(i,j) - phi(i-1,j)) / dx;
        double b = (phi(i+1,j) - phi(i,j)) / dx;
        double c = (phi(i,j) - phi(i,j-1)) / dy;
        double d_ = (phi(i,j+1) - phi(i,j)) / dy;

        double ap = std::max(a, 0.0), am = std::min(a, 0.0);
        double bp = std::max(b, 0.0), bm = std::min(b, 0.0);
        double cp = std::max(c, 0.0), cm = std::min(c, 0.0);
        double dp = std::max(d_, 0.0), dm = std::min(d_, 0.0);

        double sp = std::max(sign0(i,j), 0.0);
        double sm = std::min(sign0(i,j), 0.0);

        double grad_plus = sqrt(std::max(ap*ap, bm*bm) + std::max(cp*cp, dm*dm));
        double grad_minus = sqrt(std::max(am*am, bp*bp) + std::max(cm*cm, dp*dp));
        double H = sp * (grad_plus - 1.0) + sm * (grad_minus - 1.0);

        // ... Subcell fix: Russo & Smereka (2000)
        bool crosses = (phi0(i+1,j)*phi0(i,j) < 0.0) ||
                       (phi0(i-1,j)*phi0(i,j) < 0.0) ||
                       (phi0(i,j+1)*phi0(i,j) < 0.0) ||
                       (phi0(i,j-1)*phi0(i,j) < 0.0);

        if (crosses) {
          double dphix0 = (phi0(i+1,j) - phi0(i-1,j)) / (2.0 * dx);
          double dphiy0 = (phi0(i,j+1) - phi0(i,j-1)) / (2.0 * dy);
          double D = 2.0 * dx * phi0(i,j) /
                     std::sqrt(dphix0*dphix0 + dphiy0*dphiy0 + eps);
          phi2(i,j) = phi(i,j) - dtau * (sign0(i,j) * std::abs(phi(i,j)) - D);
        } else {
          phi2(i,j) = phi(i,j) - dtau * H;
        }
      }
    }

    // Copy phi2 into phi in-place
    for (int j = jstr; j <= jend; ++j)
      for (int i = istr; i <= iend; ++i)
        phi(i,j) = phi2(i,j);
  }

  BC::update_BCs_phi(bcTags, dx, dy, phi);
}

/******************************************************************************/
double getVolume(mtr::FMatrix<double>& heavi,
                 const double dx,
                 const double dy) {
  double vol = 0.0;
  for (int j = jstr; j <= jend-1; j++) {
    for (int i = istr; i <= iend-1; i++) {
      vol = vol + (1.0-heavi(i,j))*dx*dy;
    }
  }
  return vol;
}

/******************************************************************************/
double getLength(mtr::FMatrix<double>& phi,
                 const double dx,
                 const double dy,
                 const double Mh) {
  double len=0.0,delta=0.0;
  for (int j = jstr; j <= jend-1; j++) {
    for (int i = istr; i <= iend-1; i++) {
      if (abs(phi(i,j)) < Mh) {
        delta = 1.0/(2.0*Mh)*(1+cos(M_PI*phi(i,j)/Mh));
      } else {
        delta = 0.0;
      }
      len = len + delta*dx*dy;
    }
  }
  return len;
}
  
/******************************************************************************/
void volumeCorrection(mtr::FMatrix<double>& phi,
                      const double Mh,
                      const double V0,
                      const double Vn,
                      const double Ln) {
  double dV = (V0-Vn)/Ln;
  double delta = 0.0;
  for (int j = jstr; j <= jend-1; j++) {
    for (int i = istr; i <= iend-1; i++) {
      if (abs(phi(i,j)) < Mh) {
        phi(i,j) = phi(i,j) - dV;
      }
    }
  }
}
  

} // end namespace levset
