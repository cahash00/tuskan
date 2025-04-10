
#include <params.h>

using namespace std;

const int nghosts = 1;
int istr,iend,jstr,jend,nx,ny;

void getDomainIndices(int inx, int iny) {
  istr = 2;
  iend = inx+1;
  jstr = 2;
  jend = iny+1;
  nx = inx;
  ny = iny;
}
