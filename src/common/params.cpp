
#include <params.h>
#include <mpi.h>
#include <communication.h>

using namespace std;

const int nghosts = 1;
int istr,iend,jstr,jend,nx,ny;

void getDomainIndices(comm::Decomp& d) {
  istr = nghosts+d.i0;
  iend = d.nxl+istr;
  jstr = nghosts+d.j0;
  jend = d.nyl+jstr;
}
