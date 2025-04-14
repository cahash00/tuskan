#ifndef COMMUNICATION_H
#define COMMUNICATION_H

#include <mpi.h>
#include <matar.h>
#include <types.h>

namespace comm {

struct Decomp {
  int rank;
  int size;
  int px,py;
  int iproc,jproc;
  int nxl,nyl;
  int i0,j0;
  int nr,nl,nb,nt; // neighbor procs
  MPI_Comm cart_comm;
};

/**
 * Domain decomposition and partitioning initializer
 */
Decomp decomposer(const int nxglob,
                  const int nyglob);
void print_decomposer_info(const Decomp &d);
void exchange_ghost_cells(const std::string face,
                          fmat<double>& var);

} // end namespace comm
extern comm::Decomp d;
#endif // end COMMUNICATION_H
