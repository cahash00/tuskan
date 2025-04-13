/**
 * @file BCs.h
 * @brief Boundary Condiitions file
 * @author Caleb Hash
 */
#ifndef BCS_H
#define BCS_H

#include <cmath>
#include <matar.h>
#include <params.h>
#include <vector>
#include <input.h>
#include <types.h>

namespace BC {


struct boundarySpecs {
  fmat<int> bvals;
  std::vector<double> vel;
  double pressure;
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

bcTags tag_BCs(IO::ConfigData config);

void update_BCs(bcTags tags, 
                fmat<double>& u,
                fmat<double>& v,
                fmat<double>& p);
void update_BCs_phi(bcTags tags, 
                    fmat<double>& phi);


} // end namespace BC
#endif // BCS_H
