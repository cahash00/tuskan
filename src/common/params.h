#ifndef PARAMS_H
#define PARAMS_H

using namespace std;
extern const int nghosts;

void getDomainIndices(int nx,int ny,int nz,
                     int istr,int iend,
                     int jstr,int jend,
                     int kstr,int kend);
  
#endif // PARAMS_H
