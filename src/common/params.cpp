
#include <params.h>

using namespace std;

const int nghosts = 2;
int istr,iend,jstr,jend,kstr,kend;

void getDomainIndices(int nx,int ny,int nz) {
  istr = nghosts+1;
  iend = nx+istr;
  jstr = nghosts+1;
  jend = ny+jstr;
  kstr = nghosts+1;
  kend = nz+kstr;
}
