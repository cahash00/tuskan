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

void update_BCs(bcTags tags, 
                mtr::FMatrix<double>& u,
                mtr::FMatrix<double>& v,
                mtr::FMatrix<double>& p) {
  double ub=0.0;
  // update the left boundary
  for (int j = jstr-nghosts; j <= jend+nghosts; j++) {
    if (tags.Left(j)==0) {
      // noslip 
      u(istr-1,j) = ub-u(istr,j);
      v(istr-1,j) = ub-v(istr,j);
      p(istr-1,j) = p(istr,j);
    } else if (tags.Left(j)==7) {
      // periodic
      u(istr-1,j) = u(iend,j);
      v(istr-1,j) = v(iend,j);
    }
    if (tags.Right(j)==0) {
      // noslip 
      u(iend+1,j) = ub-u(iend,j);
      v(iend+1,j) = ub-v(iend,j);
      p(iend,j) = p(iend-1,j);
    } else if (tags.Right(j)==7) {
      // periodic
      u(iend+1,j) = u(istr,j);
      v(iend+1,j) = v(istr,j);
    }
  }

  for (int i = istr-nghosts; i <= iend+nghosts; i++) {
    if (tags.Bottom(i)==0) {
      // noslip 
      u(i,jstr-1) = 0.0-u(i,jstr);
      v(i,jstr-1) = 0.0-v(i,jstr);
      p(i,jstr-1) = p(i,jstr);
    } else if (tags.Bottom(i)==7) {
      // periodic
      u(i,jstr-1) = u(i,jend);
      v(i,jstr-1) = v(i,jend);
    }
    if (tags.Top(i)==0) {
      // noslip 
      double uwall=0.007;
      u(i,jend+1) = 2.0*uwall-u(i,jend);
      v(i,jend+1) = 0.0-v(i,jend);
      p(i,jend) = p(i,jend-1);
    } else if (tags.Top(i)==7) {
      // periodic
      u(i,jend+1) = u(i,jstr);
      v(i,jend+1) = v(i,jstr);
    }
  }
  
} // end update_BCs


} // end namespace BC
