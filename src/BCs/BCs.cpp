/**
 * @file solver.cpp
 * @brief Solver code file
 * @author Caleb Hash
 */
#include <cmath>
#include <matar.h>
#include <params.h>
#include <vector>
#include <map>
#include <input.h>
#include <types.h>

namespace BC {

struct boundarySpecs {
  fmat<int> bvals;
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
  tags.Top.pressure = config.bcTop.pressure;
  tags.Bottom.pressure = config.bcBottom.pressure;
  tags.Right.pressure = config.bcRight.pressure;
  tags.Left.pressure = config.bcLeft.pressure;

  // check to make sure periodics align
  if (config.bcLeft.type=="periodic") {
    assert(config.bcLeft.type==config.bcRight.type);
  } else if (config.bcTop.type=="periodic") {
    assert(config.bcBottom.type==config.bcTop.type);
  }
  return tags;
}

void update_BCs(bcTags tags, 
                fmat<double>& u,
                fmat<double>& v,
                fmat<double>& p) {
  // update the left + right boundaries
  for (int j = jstr-nghosts; j <= jend+nghosts; j++) {
    /**
     * LEFT BOUNDARY
     */
    if (tags.Left.bvals(j)==0) {
      // noslip 
      u(istr,j) = 0.0;
      v(istr-1,j) = -v(istr,j);
      p(istr-1,j) = p(istr,j);
    } else if (tags.Left.bvals(j)==1) {
      // moving wall
      u(istr-1,j) = tags.Left.vel[0];
      v(istr-1,j) = 2.0*tags.Left.vel[1]-v(istr,j);
      p(istr-1,j) = p(istr,j);
    } else if (tags.Left.bvals(j)==7) {
      // periodic
      u(istr-1,j) = u(iend-1,j);
      v(istr-1,j) = v(iend,j);
      p(istr-1,j) = tags.Left.pressure;
    }
    /**
     * RIGHT BOUNDARY
     */
    if (tags.Right.bvals(j)==0) {
      // noslip 
      u(iend,j) = 0.0;
      v(iend+1,j) = -v(iend,j);
      p(iend,j) = p(iend-1,j);
    } else if (tags.Right.bvals(j)==1) {
      // moving wall
      u(iend,j) = tags.Right.vel[0];
      v(iend+1,j) = 2.0*tags.Right.vel[1]-v(iend,j);
      p(iend,j) = p(iend-1,j);
    } else 
    if (tags.Right.bvals(j)==7) {
      // periodic
      u(iend,j) = u(istr,j);
      v(iend+1,j) = v(istr,j);
      p(iend,j) = tags.Right.pressure;
    }
  }

  // update top and bottom boundaries
  for (int i = istr-nghosts; i <= iend+nghosts; i++) {
    /**
     * BOTTOM BOUNDARY
     */
    if (tags.Bottom.bvals(i)==0) {
      // noslip wall 
      u(i,jstr-1) = -u(i,jstr);
      v(i,jstr) = 0.0;
      p(i,jstr-1) = p(i,jstr);
    } else if (tags.Bottom.bvals(i)==1) {
      // moving wall
      u(i,jstr-1) = 2.0*tags.Bottom.vel[0]-u(i,jstr);
      v(i,jstr-1) = tags.Bottom.vel[1];
      p(i,jstr-1) = p(i,jstr);
    } else if (tags.Bottom.bvals(i)==7) {
      // periodic
      u(i,jstr-1) = u(i,jend);
      v(i,jstr-1) = v(i,jend-1);
      p(i,jstr-1) = tags.Bottom.pressure;
    }
    /**
     * TOP BOUNDARY
     */
    if (tags.Top.bvals(i)==0) {
      // noslip 
      u(i,jend+1) = -u(i,jend);
      v(i,jend) = 0.0;
      p(i,jend) = p(i,jend-1);
    } else if (tags.Top.bvals(i)==1) {
      // moving wall
      u(i,jend) = 2.0*tags.Top.vel[0]-u(i,jend);
      v(i,jend-1) = tags.Top.vel[1];
      p(i,jend) = p(i,jend-1);
    } else if (tags.Top.bvals(i)==7) {
      // periodic
      u(i,jend+1) = u(i,jstr);
      v(i,jend) = v(i,jstr);
      p(i,jend) = tags.Top.pressure;
    }
  }

  
} // end update_BCs

void update_BCs_phi(bcTags tags, 
                    fmat<double>& phi) {
  // update the left + right boundaries
  for (int j = jstr-nghosts; j <= jend+nghosts; j++) {
    /**
     * LEFT BOUNDARY
     */
    if (tags.Left.bvals(j)==0) {
      // noslip 
      phi(istr-1,j) = phi(istr,j);
    } else if (tags.Left.bvals(j)==1) {
      // moving wall
      phi(istr-1,j) = phi(istr,j);
    } else if (tags.Left.bvals(j)==7) {
      // periodic
      phi(istr-1,j) = phi(iend-1,j);
    }
    /**
     * RIGHT BOUNDARY
     */
    if (tags.Right.bvals(j)==0) {
      // noslip 
      phi(iend,j) = phi(iend-1,j);
    } else if (tags.Right.bvals(j)==1) {
      // moving wall
      phi(iend,j) = phi(iend-1,j);
    } else 
      if (tags.Right.bvals(j)==7) {
        // periodic
        phi(iend,j) = phi(istr,j);
      }
  }

  // update top and bottom boundaries
  for (int i = istr-nghosts; i <= iend+nghosts; i++) {
    /**
     * BOTTOM BOUNDARY
     */
    if (tags.Bottom.bvals(i)==0) {
      // noslip wall 
      phi(i,jstr-1) = phi(i,jstr);
    } else if (tags.Bottom.bvals(i)==1) {
      // moving wall
      phi(i,jstr-1) = phi(i,jstr);
    } else if (tags.Bottom.bvals(i)==7) {
      // periodic
      phi(i,jstr-1) = phi(i,jend-1);
    }
    /**
     * TOP BOUNDARY
     */
    if (tags.Top.bvals(i)==0) {
      // noslip 
      phi(i,jend) = phi(i,jend-1);
    } else if (tags.Top.bvals(i)==1) {
      // moving wall
      phi(i,jend) = phi(i,jend-1);
    } else if (tags.Top.bvals(i)==7) {
      // periodic
      phi(i,jend) = phi(i,jstr);
    }
  }
} // end update_BCs


} // end namespace BC
