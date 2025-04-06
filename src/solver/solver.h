#ifndef SOLVER_H
#define SOLVER_H

#include <matar.h>

using namespace std;

/*******************************************************************************
 * @brief get the advection contribution
 ******************************************************************************/
double getAdvecU(const int& i,
                 const int& j,
                 const double rdx,
                 const double rdy,
                 mtr::FMatrix<double>& u,
                 mtr::FMatrix<double>& v);
double getAdvecV(const int& i,
                 const int& j,
                 const double rdx,
                 const double rdy,
                 mtr::FMatrix<double>& u,
                 mtr::FMatrix<double>& v);

/*******************************************************************************
 * @brief get the diffusion contribution
 ******************************************************************************/
double getDiffU(const int& i,
                const int& j,
                const double rdx,
                const double rdy,
                mtr::FMatrix<double>& u,
                mtr::FMatrix<double>& v);
double getDiffV(const int& i,
                const int& j,
                const double rdx,
                const double rdy,
                mtr::FMatrix<double>& u,
                mtr::FMatrix<double>& v);

/*******************************************************************************
 * @brief get the L2NORM between two arrays
 ******************************************************************************/
double L2NORM(mtr::FMatrix<double>& m1, 
              mtr::FMatrix<double>& m2);

/*******************************************************************************
 * @brief get the global minimum timestep
 ******************************************************************************/
double get_min_dt(const double& cfl, 
                  const double& dx,
                  const double& dy,
                  mtr::FMatrix<double>& u,
                  mtr::FMatrix<double>& v,
                  const double& nu);

/*******************************************************************************
 * @brief initialize the domain and arrays
 ******************************************************************************/
void initialize_solution(const double& uinit,
                         const double& vinit,
                         mtr::FMatrix<double>& u, 
                         mtr::FMatrix<double>& v,
                         mtr::FMatrix<double>& u2,
                         mtr::FMatrix<double>& v2,
                         mtr::FMatrix<double>& ustar,
                         mtr::FMatrix<double>& vstar,
                         mtr::FMatrix<double>& p);
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
