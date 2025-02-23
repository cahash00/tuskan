# tuskan
Two-phase Unsteady Solver using Kokkos for Applied Navier-stokes equations 

## Current / Planned Features
C++ project using Kokkos + MATAR
- [ ] MPI support
- [x] YAML input file
- [ ] Structured mesh support
  - [x] 2D rectilinear mesh generation given a $\Delta x$ and $\Delta y$
  - [x] 3D rectilinear mesh generation, given $\Delta x$, $\Delta y$, and $\Delta z$
  - [ ] Multiblock Mesh support - GridPro
- [ ] Pressure solver
