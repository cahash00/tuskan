#ifndef PARAMS_H
#define PARAMS_H

using namespace std;
extern const int nghosts;
extern int istr,iend,jstr,jend;
extern int nx,ny;

void getDomainIndices(int inx,int iny);
  
#endif // PARAMS_H
