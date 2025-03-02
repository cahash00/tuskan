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
  advec = rdx * ( pow((q(1,i+1,j) + q(1,i,j))*0.5,2) - 
                  pow((q(1,i,j)   + q(1,i-1,j))*0.5,2) ) 
        + rdy * ((q(1,i,j) + q(1,i,j+1)) * (q(2,i,j)    +q(2,i+1,j))*0.25 -
                 (q(1,i,j) + q(1,i,j-1)) * (q(2,i+1,j-1)+q(2,i,j-1))*0.25);
  return advec;
}

/******************************************************************************/
double getDiffu(const int& i,const int& j,
                const double rdx,const double rdy,
                mtr::FMatrix<double>& q) {
  double diffu;
  diffu = ( q(1,i+1,j) + 
            q(1,i-1,j) + 
            q(1,i,j+1) + 
            q(1,i,j-1) - 
            4.0*q(1,i,j) ) * rdx*rdx;
  return diffu;
}

/******************************************************************************/
double L2NORM(mtr::FMatrix<double>& qexact, mtr::FMatrix<double>& q, 
              const int& N) {
  double l2norm = {0.0};
  DO2D(j,jstr,jend,i,istr,iend,{
    l2norm = l2norm + pow(qexact(1,i,j) - q(1,i,j),2);
  });
  l2norm = sqrt(l2norm / N);
  return l2norm;
}
