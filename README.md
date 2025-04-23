# tuskan
Two-phase Unsteady Solver using Kokkos for Applied Navier-stokes equations 

C++ project using Kokkos + MATAR. I would love to do this in Fortran, but to
get good practice with Kokkos + MATAR (which I use for my research) I will do
this in C++ using Kokkos and MATAR. Fortran is better though. Using Kokkos +
MATAR should give me the opportunity to use GPU at some point. 

## Current / Planned Features
- [ ] MPI support
- [x] YAML input file
- [ ] Structured mesh support
  - [x] 2D rectilinear mesh generation given a $\Delta x$ and $\Delta y$
  - [x] 3D rectilinear mesh generation, given $\Delta x$, $\Delta y$, and $\Delta z$
  - [ ] Multiblock Mesh support - GridPro
- [x] Pressure solver
- [x] Levelset implementation
- [x] Phi redistancing
- [x] Volume Correction
- [x] Surface tension
