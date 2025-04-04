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
double getAdvecU(const int& i,const int& j,
                const double rdx,const double rdy,
                mtr::FMatrix<double>& u,
                mtr::FMatrix<double>& v) {
  double advec;
  advec = rdx*pow((u(i+1,j)+u(i,j))*0.5,2)
        - rdx*pow((u(i,j)+u(i-1,j))*0.5,2)
        + rdy*(u(i,j)+u(i,j+1))*0.5 * (v(i-1,j+1)+v(i,j+1))*0.5 
        - rdy*(u(i,j)+u(i,j-1))*0.5 * (v(i,j)+v(i-1,j))*0.5;
  return advec;
}
/******************************************************************************/
double getAdvecV(const int& i,const int& j,
                 const double rdx,const double rdy,
                 mtr::FMatrix<double>& u,
                 mtr::FMatrix<double>& v) {
  double advec;
  advec = rdy*pow((v(i,j)+v(i,j+1))*0.5,2)
        - rdy*pow((v(i,j)+v(i,j-1))*0.5,2)
        + rdx*(u(i+1,j)+u(i+1,j-1))*0.5 * (v(i,j)+v(i+1,j))*0.5
        - rdx*(u(i,j)+u(i,j-1))*0.5 * (v(i,j)+v(i-1,j))*0.5;
  return advec;
}
/******************************************************************************/
double getDiffU(const int& i,
                const int& j,
                const double rdx,
                const double rdy,
                mtr::FMatrix<double>& u,
                mtr::FMatrix<double>& v) {
  double diffu;
  diffu = ( u(i+1,j) - 2.0*u(i,j) + u(i-1,j))*rdx*rdx + 
          ( u(i,j+1) - 2.0*u(i,j) + u(i,j-1))*rdy*rdy;
  return diffu;
}
/******************************************************************************/
double getDiffV(const int& i,
                const int& j,
                const double rdx,
                const double rdy,
                mtr::FMatrix<double>& u,
                mtr::FMatrix<double>& v) {
  double diffv;
  diffv = ( v(i+1,j) - 2.0*v(i,j) + v(i-1,j))*rdx*rdx + 
          ( v(i,j+1) - 2.0*v(i,j) + v(i,j-1))*rdy*rdy;
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
                  mtr::FMatrix<double>& u,
                  mtr::FMatrix<double>& v) {
  double minval = 1e5;
  double dt = 1e5;
  for (int j = jstr-nghosts; j <= jend+nghosts; j++) {
    for (int i = istr-nghosts; i <= iend+nghosts; i++) {
      double vmag = sqrt(u(i,j)*u(i,j) + v(i,j)*v(i,j));
      dt = min(dt,cfl*dx/vmag);
    }
  }
  return dt;
}
/******************************************************************************/
void initialize_solution(mtr::FMatrix<double>& u,
                         mtr::FMatrix<double>& v,
                         mtr::FMatrix<double>& u2,
                         mtr::FMatrix<double>& v2,
                         mtr::FMatrix<double>& ustar,
                         mtr::FMatrix<double>& vstar,
                         mtr::FMatrix<double>& p) {
  u.set_values(0.000); // average u
  v.set_values(0.007);
  u2.set_values(0.0);
  v2.set_values(0.0);
  ustar.set_values(0.0);
  vstar.set_values(0.0);
  p.set_values(0.0);
}
/******************************************************************************/
double dynamic_cfl(const int ii,
                   mtr::FMatrix<double>& u,
                   mtr::FMatrix<double>& u2,
                   mtr::FMatrix<double>& v,
                   mtr::FMatrix<double>& v2,
                   double cfl,
                   double ires,
                   double resmax,
                   double cfli,
                   double cflf,
                   const int nx,
                   const int ny) {
  double res0,res1,cfl0; 
  if (ii > 0) res1 = ires;
  double cflb = cfl; // store current cfl
  ires = L2NORM(u,u2,nx*ny);
  resmax = max(resmax,ires);
  if (ii==0) {
    res0 = ires;
    res1 = ires;
  }
  if (ires == resmax) cfl0 = cfl; // if res is higher, keep
  if (ires < res1 && ires < res0) {
    cfl = cfl0*resmax/ires; // if res is lower, increase CFL
  }
  cfl = max(cfl,cflb);
  cfl = min(cflf,max(cfl,cfli));
  return cfl;
} // end dynamic_cfl
