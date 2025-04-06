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

namespace BC {


struct bcTags {
  mtr::FMatrix<int> Left;
  mtr::FMatrix<int> Right;
  mtr::FMatrix<int> Bottom;
  mtr::FMatrix<int> Top;
  bcTags(int nx, int ny);
};

bcTags tag_BCs(const std::vector<std::string>& bc_tags,
               const int nx,
               const int ny);

void update_BCs(bcTags tags, 
    mtr::FMatrix<double>& u,
    mtr::FMatrix<double>& v,
    mtr::FMatrix<double>& p);


} // end namespace BC
#endif // BCS_H
