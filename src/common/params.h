#ifndef PARAMS_H
#define PARAMS_H

using namespace std;
extern const int nghosts;
extern int istr,iend,jstr,jend,kstr,kend;

void getDomainIndices(int nx,int ny,int nz);
  
#endif // PARAMS_H
