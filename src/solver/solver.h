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
/******************************************************************************/
double getAdvecV(const int& i,
                 const int& j,
                 const double rdx,
                 const double rdy,
                 mtr::FMatrix<double>& u,
                 mtr::FMatrix<double>& v);

/*******************************************************************************
 * @brief get the diffusion contribution
 ******************************************************************************/
double getDiffu(const int& i,
                const int& j,
                const double rdx,
                const double rdy,
                mtr::FMatrix<double>& u,
                mtr::FMatrix<double>& v);
/******************************************************************************/
double getDiffv(const int& i,
                const int& j,
                const double rdx,
                const double rdy,
                mtr::FMatrix<double>& u,
                mtr::FMatrix<double>& v);

/*******************************************************************************
 * @brief get the L2NORM between two arrays
 ******************************************************************************/
double L2NORM(mtr::FMatrix<double>& uexact, 
              mtr::FMatrix<double>& u, 
              const int& N);

/*******************************************************************************
 * @brief get the global minimum timestep
 ******************************************************************************/
double get_min_dt(const double& cfl, 
                  const double& dx,
                  mtr::FMatrix<double>& u);

/*******************************************************************************
 * @brief initialize the domain and arrays
 ******************************************************************************/
void initialize_solution(mtr::FMatrix<double>& u, 
                         mtr::FMatrix<double>& v,
                         mtr::FMatrix<double>& u2,
                         mtr::FMatrix<double>& v2,
                         mtr::FMatrix<double>& p,
                         mtr::FMatrix<double>& ustar,
                         mtr::FMatrix<double>& vstar);

#endif
