#ifndef LEVSET_H
#define LEVSET_H

#include <matar.h>
#include <BCs.h>
#include <types.h>

namespace levset {


void get_phi(fmat<double>& phi,
             const fmat<double>& xc,
             const fmat<double>& yc,
             const double& xd,
             const double& yd,
             const double& r_drop);

/******************************************************************************/
void heaviside(const double& M,
               const double& h,
               const fmat<double>& phi,
               fmat<double>& heavi);

/******************************************************************************/
void advecPhi(BC::bcTags bcTags,
              const double dx,
              const double dy,
              const double dt,
              const fmat<double>& u,
              const fmat<double>& v,
              fmat<double>& phi);

/******************************************************************************/
void reinitialize(BC::bcTags bcTags,
                  const double& dx,
                  const double& dy,
                  const double& dtau,
                  const int isteps,
                  fmat<double>& phi);

/******************************************************************************/
void surfaceTension(fmat<double>& Fx,
                    fmat<double>& Fy,
                    fmat<double>& phi,
                    fmat<double>& kappa,
                    const double Mh,
                    const double sigma,
                    const double dx,
                    const double dy);

/******************************************************************************/
double getLength(fmat<double>& phi,
                 const double dx,
                 const double dy,
                 const double Mh);

/******************************************************************************/
void volumeCorrection(fmat<double>& phi,
                      const double Mh,
                      const double V0,
                      const double Vn,
                      const double Ln);

/******************************************************************************/
double getVolume(fmat<double>& heavi,
                 const double dx,
                 const double dy);

} // end namespace levset

#endif // end LEVSET_H
