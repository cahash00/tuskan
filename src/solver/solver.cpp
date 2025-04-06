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
  double advec = 0.0;
  double ue2   = 0.0;
  double uw2   = 0.0;
  double unvn  = 0.0;
  double usvs  = 0.0;

  // du2/dx using 1st-order upwind
  double ue = 0.5*(u(i+1,j) + u(i,j));
  double uw = 0.5*(u(i-1,j) + u(i,j));
  if (ue >= 0) {
    ue2 = u(i,j) * ue;
  } else {
    ue2 = u(i+1,j) * ue;
  }
  if (uw >= 0) {
    uw2 = u(i-1,j) * uw;
  } else {
    uw2 = u(i,j) * uw;
  }

  // duv/dy using 1st-order upwind
  double vn = 0.5*(v(i,j)   + v(i+1,j));
  double vs = 0.5*(v(i,j-1) + v(i+1,j-1));
  if (vn >= 0) {
    unvn = vn * u(i,j);
  } else {
    unvn = vn * u(i,j+1);
  }
  if (vs >= 0) {
    usvs = vs * u(i,j-1);
  } else {
    usvs = vs * u(i,j);
  }

  advec = rdx*(ue2 - uw2) + rdy*(unvn - usvs);
  return advec;
}
/******************************************************************************/
double getAdvecV(const int& i,const int& j,
                 const double rdx,const double rdy,
                 mtr::FMatrix<double>& u,
                 mtr::FMatrix<double>& v) {
  double advec = 0.0;
  double vn2   = 0.0;
  double vs2   = 0.0;
  double veue  = 0.0;
  double vwuw  = 0.0;

  // dv2/dy using 1st-order upwind
  double vn = 0.5*(v(i,j+1) + v(i,j));
  double vs = 0.5*(v(i,j-1) + v(i,j));
  if (vn >= 0) {
    vn2 = v(i,j) * vn;
  } else {
    vn2 = v(i,j+1) * vn;
  }
  if (vs >= 0) {
    vs2 = v(i,j-1) * vs;
  } else {
    vs2 = v(i,j) * vs;
  }

  // duv/dx using 1st-order upwind
  double ve = 0.5*(v(i,j) + v(i+1,j));
  double vw = 0.5*(v(i,j) + v(i-1,j));
  if (ve >= 0) {
    veue = ve * u(i-1,j+1);
  } else {
    veue = ve * u(i,j+1);
  }
  if (vw >= 0) {
    vwuw = vw * u(i-2,j+1);
  } else {
    vwuw = vw * u(i-1,j+1);
  }

  advec = rdy*(vn2 - vs2) + rdx*(veue - vwuw);
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
              mtr::FMatrix<double>& m2) {
  // assert that the matrices must be equal in rank
  assert(m1.order() == m2.order());
  
  double l2norm = 0.0;
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
                  mtr::FMatrix<double>& v,
                  const double& nu) {
  // double umax,vmax = 0.0;
  double umax = 0.0;
  double vmax = 0.0;
  for (int j = jstr-1; j <= jend+1; j++) {
    for (int i = istr-1; i <= iend+1; i++) {
      umax = max(umax,abs(u(i,j)));
      vmax = max(vmax,abs(v(i,j)));
    }
  }
  // find timestep limit due to advection
  double dt1 = 1e5;
  if (umax > 0.0 && vmax > 0.0) {
    for (int j = jstr-1; j <= jend+1; j++) {
      for (int i = istr-1; i <= iend+1; i++) {
        dt1 = cfl*min(dx/umax,dy/vmax);
      }
    }
  }
  double dt2 = 0.5*dx*dx*dy*dy/(nu*(dx*dx+dy*dy));
  double dt = min(dt1,dt2);
  
  return dt;
}
/******************************************************************************/
void initialize_solution(const double& uinit,
                         const double& vinit,
                         mtr::FMatrix<double>& u,
                         mtr::FMatrix<double>& v,
                         mtr::FMatrix<double>& u2,
                         mtr::FMatrix<double>& v2,
                         mtr::FMatrix<double>& ustar,
                         mtr::FMatrix<double>& vstar,
                         mtr::FMatrix<double>& p) {
  u.set_values(uinit); // average u
  v.set_values(vinit);
  u2.set_values(uinit);
  v2.set_values(vinit);
  ustar.set_values(uinit);
  vstar.set_values(vinit);
  p.set_values(0.0);
}
/******************************************************************************/
