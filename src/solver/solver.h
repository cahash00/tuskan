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

#endif
