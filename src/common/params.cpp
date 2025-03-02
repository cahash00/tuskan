
#include <params.h>

using namespace std;

const int nghosts = 1;
int istr,iend,jstr,jend;

void getDomainIndices(int nx,int ny) {
  istr = nghosts+1;
  iend = nx+istr;
  jstr = nghosts+1;
  jend = ny+jstr;
}
