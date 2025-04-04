
#include <params.h>

using namespace std;

const int nghosts = 1;
int istr,iend,jstr,jend,nx,ny;

void getDomainIndices(int inx, int iny) {
  istr = nghosts+1;
  iend = inx+istr;
  jstr = nghosts+1;
  jend = iny+jstr;
  nx = inx;
  ny = iny;
}
