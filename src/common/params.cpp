
#include <params.h>
#include <inputYAML.h>
#include <yaml-cpp/yaml.cpp>
#include <vector>
#include <params.h>

const int nghosts = 2;

void geDomainIndices(int nx,int ny,int nz,
                     int istr,int iend,
                     int jstr,int jend,
                     int kstr,int kend) {
  istr = nghosts;
  iend = nx;
  jstr = nghosts;
  jend = ny;
  kstr = nghosts;
  kend = nz;
}
