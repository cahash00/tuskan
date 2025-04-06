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
#include <input.h>

namespace BC {

struct boundarySpecs {
  mtr::FMatrix<int> bvals;
  std::vector<double> vel = std::vector<double>(2,0.0);
  double pressure = 0.0;
  // constructor
  explicit boundarySpecs(int dim) : bvals(dim) {}
};

struct bcTags {
  boundarySpecs Left;
  boundarySpecs Right;
  boundarySpecs Bottom;
  boundarySpecs Top;
  bcTags(int nx, int ny) :
    Left(ny),Right(ny),Bottom(nx),Top(nx) {}
};

/**
 * @brief tag_BCs
 * @detail goes over the domain and tags the edges for boundary conditions
 */
bcTags tag_BCs(IO::ConfigData config,
               const int nx,
               const int ny) {
  std::map<std::string,int> bcTypes = {
    {"noslip wall",0},
    {"moving wall",1},
    {"inlet",2},
    {"outlet",3},
    {"periodic",7}
  };
  bcTags tags(nx,ny);
  tags.Left.bvals.set_values(bcTypes[config.bcLeft.type]);
  tags.Right.bvals.set_values(bcTypes[config.bcRight.type]);
  tags.Bottom.bvals.set_values(bcTypes[config.bcBottom.type]);
  tags.Top.bvals.set_values(bcTypes[config.bcTop.type]);
  for (int i = 0; i < 2; i++) {
    tags.Left.vel[i] = config.bcLeft.velocity[i];
    tags.Right.vel[i] = config.bcRight.velocity[i];
    tags.Bottom.vel[i] = config.bcBottom.velocity[i];
    tags.Top.vel[i] = config.bcTop.velocity[i];
  }

  // check to make sure periodics align
  if (config.bcLeft.type=="periodic") {
    assert(config.bcLeft.type==config.bcRight.type);
  } else if (config.bcTop.type=="periodic") {
    assert(config.bcBottom.type==config.bcTop.type);
  }
  return tags;
}

void update_BCs(bcTags tags, 
                mtr::FMatrix<double>& u,
                mtr::FMatrix<double>& v,
                mtr::FMatrix<double>& p) {
  // update the left + right boundaries
  for (int j = jstr-nghosts; j <= jend+nghosts; j++) {
    if (tags.Left.bvals(j)==0) {
      // noslip 
      u(istr-1,j) = -u(istr,j);
      v(istr-1,j) = -v(istr,j);
      p(istr-1,j) = p(istr,j);
    } else if (tags.Left.bvals(j)==1) {
      // moving wall
      u(istr-1,j) = 2.0*tags.Left.vel[0]-u(istr,j);
      v(istr-1,j) = 2.0*tags.Left.vel[1]-v(istr,j);
      p(istr-1,j) = p(istr,j);
    } else if (tags.Left.bvals(j)==7) {
      // periodic
      u(istr-1,j) = u(iend,j);
      v(istr-1,j) = v(iend,j);
    }
    if (tags.Right.bvals(j)==0) {
      // noslip 
      u(iend+1,j) = -u(iend,j);
      v(iend+1,j) = -v(iend,j);
      p(iend,j) = p(iend-1,j);
    } else if (tags.Right.bvals(j)==1) {
      // moving wall
      u(iend+1,j) = 2.0*tags.Right.vel[0]-u(iend,j);
      v(iend+1,j) = 2.0*tags.Right.vel[1]-v(iend,j);
      p(iend,j) = p(iend-1,j);
    } else if (tags.Right.bvals(j)==7) {
      // periodic
      u(iend+1,j) = u(istr,j);
      v(iend+1,j) = v(istr,j);
    }
  }

  // update top and bottom boundaries
  for (int i = istr-nghosts; i <= iend+nghosts; i++) {
    if (tags.Bottom.bvals(i)==0) {
      // noslip wall 
      u(i,jstr-1) = -u(i,jstr);
      v(i,jstr-1) = -v(i,jstr);
      p(i,jstr-1) = p(i,jstr);
    } else if (tags.Bottom.bvals(i)==1) {
      // moving wall
      u(i,jstr-1) = -u(i,jstr);
      v(i,jstr-1) = -v(i,jstr);
      p(i,jstr-1) = p(i,jstr);
    } else if (tags.Bottom.bvals(i)==7) {
      // periodic
      u(i,jstr-1) = u(i,jend);
      v(i,jstr-1) = v(i,jend);
    }
    if (tags.Top.bvals(i)==0) {
      // noslip 
      u(i,jend+1) = -u(i,jend);
      v(i,jend+1) = -v(i,jend);
      p(i,jend) = p(i,jend-1);
    } else if (tags.Top.bvals(i)==1) {
      // moving wall
      u(i,jend+1) = 2.0*tags.Top.vel[0]-u(i,jend);
      v(i,jend+1) = 2.0*tags.Top.vel[1]-v(i,jend);
      p(i,jend) = p(i,jend-1);
    } else if (tags.Top.bvals(i)==7) {
      // periodic
      u(i,jend+1) = u(i,jstr);
      v(i,jend+1) = v(i,jstr);
    }
  }
  
} // end update_BCs


} // end namespace BC
