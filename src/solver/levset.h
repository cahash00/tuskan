#ifndef LEVSET_H
#define LEVSET_H

#include <matar.h>
#include <BCs.h>

namespace levset {


void get_phi(mtr::FMatrix<double>& phi,
             const mtr::FMatrix<double>& xc,
             const mtr::FMatrix<double>& yc,
             const double& xd,
             const double& yd,
             const double& r_drop);

/******************************************************************************/
void heaviside(const double& M,
               const double& h,
               const mtr::FMatrix<double>& phi,
               mtr::FMatrix<double>& heavi);

/******************************************************************************/
void weno(BC::bcTags bcTags,
          const double dx,
          const double dy,
          const double dt,
          const mtr::FMatrix<double>& u,
          const mtr::FMatrix<double>& v,
          mtr::FMatrix<double>& phi);

/******************************************************************************/
void reinitialize(BC::bcTags bcTags,
                  const double& dx,
                  const double& dy,
                  const double& dtau,
                  const int isteps,
                  mtr::FMatrix<double>& phi);

  /******************************************************************************/
double curvature(const int i, const int j,
                 const double dx, const double dy,
                 const mtr::FMatrix<double>& phi);


} // end namespace levset

#endif // end LEVSET_H
