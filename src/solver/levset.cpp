

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
  for (int j = jstr-1; j <= jend; j++) {
    for (int i = istr-1; i <= iend; i++) {
      phi(i,j) = sqrt(pow((xc(i,j)-xd),2) + pow((yc(i,j)-yd),2))-r_drop;
      double phi_water = xc(i,j)-1.0;
      double phi_air = yc(i,j)-1.0;
      phi(i,j) = max(phi_water,phi_air);
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
  for (int j = jstr; j <= jend; j++) {
    for (int i = istr; i <= iend; i++) {
      double ucell=0.0,vcell=0.0;
      ucell = (u(i,j)+u(i-1,j))*0.5;
      vcell = (v(i,j)+v(i,j-1))*0.5;

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

void reinitialize(BC::bcTags bcTags,
                  const double& dx,
                  const double& dy,
                  const double& dtau,
                  const int isteps,
                  mtr::FMatrix<double>& phi) {
  const double eps=1.0e-12;
  mtr::FMatrix<double> sign0(phi.dims(1),phi.dims(2));
  mtr::FMatrix<double> phi0(phi.dims(1),phi.dims(2));
  mtr::FMatrix<double> phi2(phi.dims(1),phi.dims(2));
  BC::update_BCs_phi(bcTags,dx,dy,phi0);
  for (int j = jstr; j <= jend-1; j++) {
    for (int i = istr; i <= iend-1; i++) {
      phi2(i,j) = phi(i,j);
      phi0(i,j) = phi(i,j);
      double gphi0 = sqrt(pow((phi(i+1,j) - phi(i-1,j)) / (2.0 * dx),2)
                        + pow((phi(i,j+1) - phi(i,j-1)) / (2.0 * dy),2));
      sign0(i,j) = phi(i,j) / sqrt(phi(i,j)*phi(i,j)+gphi0*gphi0*dx*dx+eps*eps);
    }
  }
  BC::update_BCs_phi(bcTags,dx,dy,phi0);
  BC::update_BCs_phi(bcTags,dx,dy,phi2);
  

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
        //  Vastly improves the ability to reduce the amount of volume lost
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
    
    BC::update_BCs_phi(bcTags,dx,dy,phi2);
    for (int j = jstr; j <= jend-1; j++) {
      for (int i = istr; i <= iend-1; i++) {
        phi(i,j) = phi2(i,j);
      }
    }
  }
} // end reinitialize

double getVolume(mtr::FMatrix<double>& heavi,
                 const double dx,
                 const double dy) {
  double vol = 0.0;
  for (int j = jstr; j <= jend; j++) {
    for (int i = istr; i <= iend; i++) {
      vol = vol + (1.0-heavi(i,j))*dx*dy;
    }
  }
  return vol;
}

double getLength(mtr::FMatrix<double>& phi,
                 const double dx,
                 const double dy,
                 const double Mh) {
  double len=0.0,delta=0.0;
  for (int j = jstr; j <= jend; j++) {
    for (int i = istr; i <= iend; i++) {
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
  
  

void volumeCorrection(mtr::FMatrix<double>& phi,
                      const double Mh,
                      const double V0,
                      const double Vn,
                      const double Ln) {
  double dV = (V0-Vn)/Ln;
  double delta = 0.0;
  for (int j = jstr; j <= jend; j++) {
    for (int i = istr; i <= iend; i++) {
      if (abs(phi(i,j)) < Mh) {
        phi(i,j) = phi(i,j) - dV;
      }
    }
  }
}
  

} // end namespace levset
