#include <mpi.h>
#include <iostream>
#include <communication.h>
#include <types.h>

comm::Decomp d;


namespace comm {

/**
 * domain decomposition and partitioning initializer
 */
Decomp decomposer(const int nxglob,
                  const int nyglob) {

  MPI_Comm_rank(MPI_COMM_WORLD,&d.rank);
  MPI_Comm_size(MPI_COMM_WORLD,&d.size);

  int dims[2] = {0,0};
  int periods[2] = {0,0};
  MPI_Dims_create(d.size,2,dims);

  d.px = dims[0];
  d.py = dims[1];

  MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,1,&d.cart_comm);

  int coords[2];
  MPI_Cart_coords(d.cart_comm,d.rank,2,coords);

  // assign which proc this
  d.iproc = coords[0];
  d.jproc = coords[1]; 

  // get number of cells that will live on this proc
  d.nxl = nxglob / d.px;
  d.nyl = nyglob / d.py;

  // get the global index of the starting cell on this proc
  d.i0 = d.iproc*d.nxl+1;
  d.j0 = d.jproc*d.nyl+1;
  
  // identify neighbor ranks for ghost exchanges
  MPI_Cart_shift(d.cart_comm, 0, 1, &d.nl, &d.nr);  // X-direction
  MPI_Cart_shift(d.cart_comm, 1, 1, &d.nb, &d.nt);  // Y-direction

  return d;
} // end decomposer

/******************************************************************************/
void print_decomposer_info(const Decomp &d) {
    std::cout << "[Rank " << d.rank << "] px = " << d.px << ", py = " << d.py << std::endl;
    std::cout << "[Rank " << d.rank << "] iproc = " << d.iproc << ", jproc = " << d.jproc << std::endl;
    std::cout << "[Rank " << d.rank << "] nx_local = " << d.nxl << ", ny_local = " << d.nyl << std::endl;
    std::cout << "[Rank " << d.rank << "] i0 = " << d.i0 << ", j0 = " << d.j0 << std::endl;
    std::cout << "[Rank " << d.rank << "] neighbors (L,R,B,T): "
              << d.nl << ", " << d.nr << ", "
              << d.nb << ", " << d.nt << std::endl;
} // end print_decomposer_info

/******************************************************************************/
void exchange_ghost_cells(const std::string face,
                          fmat<double>& var) {
  MPI_Request reqs[8];
  MPI_Status stats[8];

  // ... LEFT BOUNDARY
  if (d.nl != MPI_PROC_NULL) {
    for (int j = 1; j <= d.nyl+1; j++) {
      MPI_Sendrecv(&var(d.i0  ,j),1,MPI_DOUBLE,d.nl,0,
                   &var(d.i0-1,j),1,MPI_DOUBLE,d.nl,0,
                   d.cart_comm,MPI_STATUS_IGNORE);
    }
  } 
  
  if (d.nr != MPI_PROC_NULL) {
    for (int j = 1; j <= d.nyl+1; j++) {
      MPI_Sendrecv(&var(d.nxl+1,j),1,MPI_DOUBLE,d.nr,0,
                   &var(d.nxl  ,j),1,MPI_DOUBLE,d.nr,0,
                   d.cart_comm,MPI_STATUS_IGNORE);
    }
  } 
  if (d.nb != MPI_PROC_NULL) {
    for (int i = 1; i <= d.nxl+1; i++) {
      MPI_Sendrecv(&var(i,d.j0-1),1,MPI_DOUBLE,d.nb,1,
                   &var(i,d.j0  ),1,MPI_DOUBLE,d.nb,1,
                   d.cart_comm,MPI_STATUS_IGNORE);
    }
  }
  if (d.nt != MPI_PROC_NULL) {
    for (int i = 1; i <= d.nxl+1; i++) {
      MPI_Sendrecv(&var(i,d.nyl  ),1,MPI_DOUBLE,d.nt,1,
                   &var(i,d.nyl+1),1,MPI_DOUBLE,d.nt,1,
                   d.cart_comm,MPI_STATUS_IGNORE);
    }
  } 
  MPI_Barrier(MPI_COMM_WORLD);  // ensure all ranks have finished setup/output
}
    



} // end namespace comm
