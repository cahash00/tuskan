#ifndef SOLVER_H
#define SOLVER_H

#include <types.h>
#include <input.h>
#include <matar.h>

using namespace std;

/*******************************************************************************
 * @brief get the advection contribution
 ******************************************************************************/
double getAdvecU(const int& i,
                 const int& j,
                 const double rdx,
                 const double rdy,
                 const fmat<double>& u,
                 const fmat<double>& v);
double getAdvecV(const int& i,
                 const int& j,
                 const double rdx,
                 const double rdy,
                 const fmat<double>& u,
                 const fmat<double>& v);

/*******************************************************************************
 * @brief get the diffusion contribution
 ******************************************************************************/
double getDiffU(const int& i,
                const int& j,
                const double rdx,
                const double rdy,
                const fmat<double>& nu,
                const fmat<double>& u,
                const fmat<double>& v);
double getDiffV(const int& i,
                const int& j,
                const double rdx,
                const double rdy,
                const fmat<double>& nu,
                const fmat<double>& u,
                const fmat<double>& v);

/*******************************************************************************
 * @brief get the L2NORM between two arrays
 ******************************************************************************/
double L2NORM(const fmat<double>& m1, 
              const fmat<double>& m2);

/*******************************************************************************
 * @brief get the global minimum timestep
 ******************************************************************************/
double get_min_dt(const double& cfl, 
                  const double& dx,
                  const double& dy,
                  const fmat<double>& u,
                  const fmat<double>& v,
                  const double& nu);

/*******************************************************************************
 * @brief initialize the domain and arrays
 ******************************************************************************/
void initialize_solution(IO::ConfigData& config,
                         fmat<double>& u, 
                         fmat<double>& v,
                         fmat<double>& u2,
                         fmat<double>& v2,
                         fmat<double>& u_old,
                         fmat<double>& v_old,
                         fmat<double>& ustar,
                         fmat<double>& vstar,
                         fmat<double>& p);
double dynamic_cfl(const int ii,
                   double cfl,
                   double ires,
                   double res0,
                   double res1,
                   double resmax,
                   double cfli,
                   double cflf,
                   const int nx,
                   const int ny);
#endif
