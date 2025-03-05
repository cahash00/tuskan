#ifndef SOLVER_H
#define SOLVER_H

#include <matar.h>

using namespace std;

double getAdvec(const int& i,const int& j,
                const double rdx,const double rdy,
                mtr::FMatrix<double>& q);

double getDiffu(const int& i,const int& j,
                const double rdx,const double rdy,
                mtr::FMatrix<double>& q);

double L2NORM(mtr::FMatrix<double>& qexact, mtr::FMatrix<double>& q, 
              const int& N);

double get_min_dt(const double& cfl, const double& dx,FMatrix<double>& q);
void update_solution(mtr::FMatrix<double>& q, mtr::FMatrix<double> q2);
void initialize_solution(mtr::FMatrix<double>& q, mtr::FMatrix<double> q2);
#endif
