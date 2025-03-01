
#include <params.h>

using namespace std;

extern const int nghosts = 2;

void getDomainIndices(int nx,int ny,int nz,
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
