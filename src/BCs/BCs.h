/**
 * @file BCs.h
 * @brief Boundary Condiitions file
 * @author Caleb Hash
 * @date 2025-02-25
 */
#ifndef BCS_H
#define BCS_H

#include <cmath>
#include <matar.h>
#include <params.h>
#include <vector>
#include <input.h>

namespace BC {


struct boundarySpecs {
  mtr::FMatrix<int> bvals;
  std::vector<double> vel;
  double pressure;
  double rho;
  double nu;
  // constructor
  boundarySpecs(int dim);
};
struct bcTags {
  boundarySpecs Left;
  boundarySpecs Right;
  boundarySpecs Bottom;
  boundarySpecs Top;
  bcTags(int nx, int ny);
};

bcTags tag_BCs(IO::ConfigData config,
               mtr::FMatrix<double>& xn,
               const int nx,
               const int ny);

void update_BCs(bcTags tags, 
                mtr::FMatrix<double>& u,
                mtr::FMatrix<double>& v,
                mtr::FMatrix<double>& p);
void update_BCs_phi(bcTags tags, 
                    const double dx,
                    const double dy,
                    mtr::FMatrix<double>& phi);
void update_BCs_rho(bcTags tags, 
    mtr::FMatrix<double>& phi);
void update_BCs_nu(bcTags tags, 
    mtr::FMatrix<double>& phi);


} // end namespace BC
#endif // BCS_H
