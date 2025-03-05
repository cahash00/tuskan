/**
 * @file solver.cpp
 * @brief Solver code file
 * @author Caleb Hash
 * @date 2025-02-25
 */
#include <cmath>
#include <matar.h>
#include <params.h>

using namespace std;

/******************************************************************************/
double getAdvec(const int& i,const int& j,
                const double rdx,const double rdy,
                mtr::FMatrix<double>& q) {
  double advec;
  advec = rdx*pow((q(1,i+1,j)+q(1,i,j))*0.5,2)
        - rdx*pow((q(1,i,j)+q(1,i-1,j))*0.5,2);
        // + rdy*(q(1,i,j)+q(1,i,j+1))*0.5 * (q(2,i,j)+q(2,i+1,j))*0.5 
        // - rdy*(q(1,i,j)+q(1,i,j-1))*0.5 * (q(2,i+1,j-1)+q(2,i,j-1))*0.5;
  return advec;
}

/******************************************************************************/
double getDiffu(const int& i,const int& j,
                const double rdx,const double rdy,
                mtr::FMatrix<double>& q) {
  double diffu;
  diffu = ( q(1,i+1,j) - 2.0*q(1,i,j) + q(1,i-1,j))*rdx*rdx + 
          ( q(1,i,j+1) - 2.0*q(1,i,j) + q(1,i,j-1))*rdy*rdx;
  return diffu;
}

/******************************************************************************/
double L2NORM(mtr::FMatrix<double>& m1, mtr::FMatrix<double>& m2, 
              const int& N) {
  double l2norm = {0.0};
  DO2D(j,jstr,jend,i,istr,iend,{
    l2norm = l2norm + pow(m1(1,i,j) - m2(1,i,j),2);
  });
  l2norm = sqrt(l2norm / static_cast<double>(N));
  return l2norm;
}
/******************************************************************************/
double get_min_dt(const double& cfl, const double& dx,mtr::FMatrix<double>& q) {
  double minval = 1e5;
  double dt = {0.0};
  DO_LOOP(j,jstr-nghosts,jend+nghosts,{
    DO_LOOP(i,istr-nghosts,iend+nghosts,{
      dt = min(minval,cfl*dx/abs(q(1,i,j)));
    });
  });
  return dt;
}
/******************************************************************************/
void update_solution(mtr::FMatrix<double>& q, mtr::FMatrix<double> q2) {
  DO_LOOP(j,jstr,jend,{
    DO_LOOP(i,istr,iend,{
      q(1,i,j) = q2(1,i,j);
    });
  });
}
/******************************************************************************/
void initialize_solution(mtr::FMatrix<double>& q, mtr::FMatrix<double> q2) {
  DO_LOOP(j,jstr-nghosts,jend+nghosts,{
    DO_LOOP(i,istr-nghosts,iend+nghosts,{
      q(1,i,j) = 0.007; // average u
      q(2,i,j) = 0.0; // zero y-velocity
      q2(1,i,j) = 0.0;
      q2(2,i,j) = 0.0;
    });
  });
}
