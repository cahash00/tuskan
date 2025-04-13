#ifndef PARAMS_H
#define PARAMS_H

#include <mpi.h>
#include <communication.h>

using namespace std;
extern const int nghosts;
extern int istr,iend,jstr,jend;

void getDomainIndices();
  
#endif // PARAMS_H
