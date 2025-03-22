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
                mtr::FMatrix<double>& u,
                mtr::FMatrix<double>& v) {
  double advec;
  advec = rdx*pow((u(i+1,j)+u(i,j))*0.5,2)
        - rdx*pow((u(i,j)+u(i-1,j))*0.5,2)
        + rdy*(u(i,j)+u(i,j+1))*0.5 * (v(i,j)+v(i+1,j))*0.5 
        - rdy*(u(i,j)+u(i,j-1))*0.5 * (v(i+1,j-1)+v(i,j-1))*0.5;
  return advec;
}

/******************************************************************************/
double getDiffu(const int& i,
                const int& j,
                const double rdx,
                const double rdy,
                mtr::FMatrix<double>& u,
                mtr::FMatrix<double>& v) {
  double diffu;

  // Second derivatives (viscous terms)
  double d2udx2 = (u(i+1,j) - 2.0*u(i,j) + u(i-1,j)) * rdx * rdx;
  double d2udy2 = (u(i,j+1) - 2.0*u(i,j) + u(i,j-1)) * rdy * rdy;

  // Cross-diffusion term
  double d2uvdxdy = (v(i+1,j+1) - v(i-1,j+1) 
                   - v(i+1,j-1) + v(i-1,j-1)) 
                  * 0.25 * rdx * rdy;

  diffu = d2udx2 + d2udy2 + d2uvdxdy;

  return diffu;
}
/******************************************************************************/
double getDiffv(const int& i, 
                const int& j,
                const double rdx, 
                const double rdy,
                mtr::FMatrix<double>& u,
                mtr::FMatrix<double>& v) {
  double diffv;

  // Second derivatives (viscous terms)
  double d2vdx2 = (v(i+1,j) - 2.0*v(i,j) + v(i-1,j)) * rdx * rdx;
  double d2vdy2 = (v(i,j+1) - 2.0*v(i,j) + v(i,j-1)) * rdy * rdy;

  // Cross-diffusion term
  double d2uvdxdy = (u(i+1,j+1) - u(i-1,j+1) 
                   - u(i+1,j-1) + u(i-1,j-1)) 
                  * 0.25 * rdx * rdy;

  diffv = d2vdx2 + d2vdy2 + d2uvdxdy;

  return diffv;
}

/******************************************************************************/
double L2NORM(mtr::FMatrix<double>& m1, 
              mtr::FMatrix<double>& m2, 
              const int& N) {
  // assert that the matrices must be equal in rank
  assert(m1.order() == m2.order());
  
  double l2norm = 0.0;
  double lsum = 0.0;
  double total = 0.0;
  for (int j = jstr; j <= jend; j++) {
    for (int i = istr; i <= iend; i++) {
      double d = m2(i,j) - m1(i,j);
      l2norm += d*d;
      total += 1.0;
    }
  }
  l2norm = sqrt(l2norm / total);
  return l2norm;
}
/******************************************************************************/
double get_min_dt(const double& cfl, 
                  const double& dx,
                  mtr::FMatrix<double>& u) {
  double minval = 1e5;
  double dt = 1e5;
  DO_ALL(j,jstr-nghosts,jend+nghosts,
         i,istr-nghosts,iend+nghosts,{
    dt = min(dt,cfl*dx/abs(u(i,j)));
  });
  return dt;
}
/******************************************************************************/
void initialize_solution(mtr::FMatrix<double>& u,
                         mtr::FMatrix<double>& v,
                         mtr::FMatrix<double>& u2,
                         mtr::FMatrix<double>& p,
                         mtr::FMatrix<double>& ustar,
                         mtr::FMatrix<double>& vstar) {
  u.set_values(0.007);
  v.set_values(0.0);
  v.set_values(0.0);
  u2.set_values(0.0);
  p.set_values(0.0);
  ustar.set_values(0.0);
  vstar.set_values(0.0);
}
