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
  double u_n,u_s,u_e,u_w;
  double v_n,v_s,v_e,v_w;
  // u components
  u_n = 0.5*(u(i,j+1) + u(i,j));
  u_s = 0.5*(u(i,j-1) + u(i,j));
  u_e = 0.5*(u(i+1,j) + u(i,j));
  u_w = 0.5*(u(i-1,j) + u(i,j));
  // v components
  v_n = 0.5*(v(i,j)   + v(i+1,j));
  v_s = 0.5*(v(i,j-1) + v(i+1,j-1));
  v_e = 0.5*(v(i,j)   + v(i,j-1));
  v_w = 0.5*(v(i+1,j) + v(i+1,j-1));
  
  advec = rdx*(u_e*u_e-u_w*u_w) + rdy*(u_n*v_n-u_s*v_s);
  return advec;
}
// double getAdvecU(const int& i,const int& j,
//     const double rdx,const double rdy,
//     mtr::FMatrix<double>& u,
//     mtr::FMatrix<double>& v) {
//   double advec;
//   double ue2, uw2, unvn, usvs;

//   // ∂(u²)/∂x using 1st-order upwind
//   double ue = 0.5*(u(i+1,j) + u(i,j));
//   double uw = 0.5*(u(i-1,j) + u(i,j));
//   if (ue >= 0)
//     ue2 = u(i,j) * ue;
//   else
//     ue2 = u(i+1,j) * ue;
//   if (uw >= 0)
//     uw2 = u(i-1,j) * uw;
//   else
//     uw2 = u(i,j) * uw;

//   // ∂(uv)/∂y using 1st-order upwind
//   double vn = 0.5*(v(i,j)   + v(i+1,j));     // v at (i+1/2,j+1/2)
//   double vs = 0.5*(v(i,j-1) + v(i+1,j-1));   // v at (i+1/2,j-1/2)
//   double un = 0.5*(u(i,j+1) + u(i,j));
//   double us = 0.5*(u(i,j-1) + u(i,j));
//   if (vn >= 0)
//     unvn = vn * u(i,j);
//   else
//     unvn = vn * u(i,j+1);
//   if (vs >= 0)
//     usvs = vs * u(i,j-1);
//   else
//     usvs = vs * u(i,j);

//   advec = rdx*(ue2 - uw2) + rdy*(unvn - usvs);
//   return advec;
// }
/******************************************************************************/
double getAdvecV(const int& i,const int& j,
                 const double rdx,const double rdy,
                 mtr::FMatrix<double>& u,
                 mtr::FMatrix<double>& v) {
  double advec;
  double u_n,u_s,u_e,u_w;
  double v_n,v_s,v_e,v_w;
  // v components
  v_n = 0.5*(v(i,j) + v(i,j+1));
  v_s = 0.5*(v(i,j) + v(i,j-1));
  v_e = 0.5*(v(i,j) + v(i+1,j));
  v_w = 0.5*(v(i,j) + v(i-1,j));
  // u components
  u_n = 0.5*(u(i-1,j+1) + u(i,j+1));
  u_s = 0.5*(u(i-1,j) + u(i,j));
  u_e = 0.5*(u(i,j+1) + u(i,j));
  u_w = 0.5*(u(i-1,j+1) + u(i-1,j));

  advec = rdy*( v_n*v_n-v_s*v_s ) + rdx*( v_e*u_e-v_w*u_w );
  return advec;
}
// double getAdvecV(const int& i,const int& j,
//                         const double rdx,const double rdy,
//                         mtr::FMatrix<double>& u,
//                         mtr::FMatrix<double>& v) {
//   double advec;
//   double vn2, vs2, veue, vwuw;

//   // ∂(v²)/∂y with upwind
//   double vn = 0.5*(v(i,j+1) + v(i,j));
//   double vs = 0.5*(v(i,j-1) + v(i,j));
//   if (vn >= 0)
//     vn2 = v(i,j) * vn;
//   else
//     vn2 = v(i,j+1) * vn;
//   if (vs >= 0)
//     vs2 = v(i,j-1) * vs;
//   else
//     vs2 = v(i,j) * vs;

//   // ∂(uv)/∂x with upwind
//   double ue = 0.5*(u(i,j+1) + u(i,j));
//   double uw = 0.5*(u(i-1,j+1) + u(i-1,j));
//   double ve = 0.5*(v(i,j) + v(i+1,j));
//   double vw = 0.5*(v(i,j) + v(i-1,j));
//   if (ve >= 0)
//     veue = ve * u(i-1,j+1);
//   else
//     veue = ve * u(i,j+1);
//   if (vw >= 0)
//     vwuw = vw * u(i-2,j+1);  // watch indexing — may need bounds check!
//   else
//     vwuw = vw * u(i-1,j+1);

//   advec = rdy*(vn2 - vs2) + rdx*(veue - vwuw);
//   return advec;
// }
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
              mtr::FMatrix<double>& m2) {
  // assert that the matrices must be equal in rank
  assert(m1.order() == m2.order());
  
  double l2norm = 0.0;
  double lsum = 0.0;
  double total = 0.0;
  for (int j = jstr; j <= jend; j++) {
    for (int i = istr; i <= iend; i++) {
      double d = m2(i,j) - m1(i,j) + 1.0e-15;
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
                  const double& dy,
                  mtr::FMatrix<double>& u,
                  mtr::FMatrix<double>& v) {
  double dt=1.0e5;
  for (int j = jstr; j <= jend; j++) {
    for (int i = istr; i <= iend; i++) {
      dt = min(abs(cfl*dy/(v(i,j)+1.0e-15)),abs(cfl*dx/(u(i,j)+1.0e-15)));
    }
  }
  return dt;
}
/******************************************************************************/
void initialize_solution(const double& uinit,
                         const double& vinit,
                         mtr::FMatrix<double>& u,
                         mtr::FMatrix<double>& v,
                         mtr::FMatrix<double>& u_old,
                         mtr::FMatrix<double>& v_old,
                         mtr::FMatrix<double>& u2,
                         mtr::FMatrix<double>& v2,
                         mtr::FMatrix<double>& ustar,
                         mtr::FMatrix<double>& vstar,
                         mtr::FMatrix<double>& p) {
  u.set_values(uinit);
  v.set_values(vinit);
  u_old.set_values(uinit);
  v_old.set_values(vinit);
  u2.set_values(uinit);
  v2.set_values(vinit);
  ustar.set_values(uinit);
  vstar.set_values(vinit);
  p.set_values(1.0);
}
/******************************************************************************/
double dynamic_cfl(const int ii,
                   double cfl,
                   double ires,
                   double res0,
                   double res1,
                   double resmax,
                   double cfli,
                   double cflf,
                   const int nx,
                   const int ny) {
  double cfl0,cflnew; 
  if (ii > 0) res1 = ires;
  resmax = max(resmax,ires);
  if (ii==0) {
    res0 = ires;
    res1 = ires;
  }
  if (ires == resmax) cfl0 = cfl; // if res is higher, keep
  if (ires < res1 && ires < res0) {
    cflnew = cfl0*resmax/ires; // if res is lower, increase CFL
  }
  cflnew = max(cflnew,cfl);
  cflnew = min(cflf,max(cflnew,cfli));
  return cflnew;
} // end dynamic_cfl
