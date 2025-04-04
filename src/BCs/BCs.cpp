/**
 * @file solver.cpp
 * @brief Solver code file
 * @author Caleb Hash
 * @date 2025-02-25
 */
#include <cmath>
#include <matar.h>
#include <params.h>
#include <vector>
#include <map>

namespace BC {


struct bcTags {
  mtr::FMatrix<int> Left;
  mtr::FMatrix<int> Right;
  mtr::FMatrix<int> Bottom;
  mtr::FMatrix<int> Top;
  bcTags(int nx, int ny) :
    Left(ny),Right(ny),Bottom(nx),Top(nx) {}
};

/**
 * @brief tag_BCs
 * @detail goes over the domain and tags the edges for boundary conditions
 */
bcTags tag_BCs(const std::vector<std::string>& bc_tags,
               const int nx,
               const int ny) {
  std::map<std::string,int> bcTypes = {
    {"noslip",0},
    {"periodic",7},
    {"inlet",1},
    {"outlet",2}
  };
  bcTags tags(nx,ny);
  tags.Left.set_values(bcTypes[bc_tags[0]]);
  tags.Right.set_values(bcTypes[bc_tags[1]]);
  tags.Bottom.set_values(bcTypes[bc_tags[2]]);
  tags.Top.set_values(bcTypes[bc_tags[3]]);

  // check to make sure periodics align
  if (bc_tags[0]=="periodic") {
    assert(bc_tags[0]==bc_tags[1]);
  } else if (bc_tags[2]=="periodic") {
    assert(bc_tags[2]==bc_tags[3]);
  }
  return tags;
}

void update_BCs(bcTags tags, mtr::FMatrix<double>& q) {
  // update the left boundary
  for (int j = jstr-nghosts; j <= jend+nghosts; j++) {
    if (tags.Left(j)==0) {
      q(istr-1,j) = -q(istr,j);
    } else if (tags.Left(j)==7) {
      q(istr-1,j) = q(iend,j);
    }
    if (tags.Right(j)==0) {
      q(iend+1,j) = -q(iend,j);
    } else if (tags.Right(j)==7) {
      q(iend+1,j) = q(istr,j);
    }
  }

  for (int i = istr-nghosts; i <= iend+nghosts; i++) {
    if (tags.Bottom(i)==0) {
      q(i,jstr-1) = -q(i,jstr);
    } else if (tags.Bottom(i)==7) {
      q(i,jstr-1) = q(i,jend);
    }
    if (tags.Top(i)==0) {
      q(i,jend+1) = -q(i,jend);
    } else if (tags.Top(i)==7) {
      q(i,jend+1) = q(i,jstr);
    }
  }
  
} // end update_BCs


} // end namespace BC
