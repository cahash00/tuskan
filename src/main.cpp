/**
 * @file main.cpp
 * @author Caleb Hash
 * @date April 13, 2025
 * @brief Two-Phase Flow Solver
 * @details First part of the final project
 */

#include <mesh.h>
#include <output.h>
#include <args.h>
#include <input.h>
#include <matar.h>
#include <Kokkos_Core.hpp>
#include <argparse/argparse.hpp>
#include <pprint.hpp>
#include <params.h>
#include <generalUtils.h>
#include <solver.h>
#include <levset.h>
#include <pressure.h>
#include <logger.h>
#include <restart.h>
#include <cmath>
#include <BCs.h>
#include <fstream>
#include <string>
#include <types.h>
#include <mpi.h>
#include <communication.h>

using namespace std;

/**
 * Main program
 */
int main(int argc, char* argv[]){
  int rank=0,size=0;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  // ... create the logger
  IO::init_logger();
  /***************************************************************************
   *                           MAIN PROGRAM START                            *
   **************************************************************************/
  // ... Parse out the user input form command line
  argparse::ArgumentParser program("TUSKAN","0.0.0");
  IO::getUserInput(argc,argv,program);
  // assign the input file that was given
  auto inFile = program.get<string>("-i");

  // ... parse the YAML input deck
  IO::ConfigData config = IO::parseInputDeck(inFile);
  // check to see if the output directory exists
  IO::check_directories(config.fv.dir);
  
  // ... initialize and finalize Kokkos with the main scope
  if (rank==0) if(rank==0) IO::logger->info("Initializing Kokkos");
  Kokkos::ScopeGuard kokkos_guard(argc, argv); 
  // Kokkos::initialize(argc,argv);
  if (rank==0) if(rank==0) IO::logger->info("  done");

  // ... call decomposer
  comm::Decomp d = comm::decomposer(config.mesh.nx,config.mesh.ny);
  comm::print_decomposer_info(d);
  // ... get domain stats
  getDomainIndices();

  // ... initialize the MATAR matrices for the domain
  // use local number of cells rather than the global ncells
  const int isize = d.nxl+2;
  const int jsize = d.nyl+2;
  fmat<double> xc(isize,jsize),yc(isize,jsize),
               xn(isize+1,jsize+1),yn(isize+1,jsize+1);
  fmat<double> uc(isize,jsize),vc(isize,jsize);
  fmat<double> p(isize+1,jsize+1);
  fmat<double> phi(isize+1,jsize+1);
  fmat<double> heavi(isize+1,jsize+1);
  fmat<double> rho(isize+1,jsize+1);
  fmat<double> kappa(isize+1,jsize+1);
  fmat<double> nu(isize+1,jsize+1);
  fmat<double> ustar(isize+1,jsize+1);
  fmat<double> vstar(isize+1,jsize+1);
  fmat<double> u(isize+1,jsize+1);
  fmat<double> v(isize+1,jsize+1);
  fmat<double> u_old(isize+1,jsize+1);
  fmat<double> v_old(isize+1,jsize+1);
  fmat<double> u2(isize+1,jsize+1);
  fmat<double> v2(isize+1,jsize+1);
  fmat<double> Fx(isize+1,jsize+1);
  fmat<double> Fy(isize+1,jsize+1);
  heavi.set_values(0.0);
  phi.set_values(0.0);
  
  /***********************
   * call mesh generator *
   ***********************/
  Timer timer;
  timer.start();
  if(rank==0) IO::logger->info("Generating 2D mesh");
  double dx = config.mesh.lx/static_cast<double>(config.mesh.nx);
  double dy = config.mesh.ly/static_cast<double>(config.mesh.ny);
  mesh::mesher2D(xc,yc,xn,yn,dx,dy);
  timer.stop();
  if(rank==0) IO::logger->info("  done ({} seconds)",timer.time());
  if(rank==0) IO::logger->info("Tagging boundaries");
  BC::bcTags bcTags = BC::tag_BCs(config);
  if(rank==0) IO::logger->info("  done");

  /******************
   * initialization *
   ******************/
  ofstream logFile(config.res.file, ios::out);
  if (!logFile) {
      std::cerr << "Error opening log file!\n";
      return -1;
  }
  timer.start();
  if(rank==0) IO::logger->info("Initializing the domain");
  double rdx  = 1.0/dx;  // reciprocal of dx
  double rdy  = 1.0/dy;  // reciprocal of dx
  if(rank==0) IO::logger->info("  dx: {},dy: {}",dx,dy);

  // ... initialize domain
  initialize_solution(config,u,v,u2,v2,u_old,v_old,ustar,vstar,p);
  if (config.restart.load) {
    restart::load("u.rank"+std::to_string(rank)+".converged",u);
    restart::load("v.rank"+std::to_string(rank)+".converged",v);
    restart::load("p.rank"+std::to_string(rank)+".converged",p);
  }
  BC::update_BCs(bcTags,u,v,p);

  // ... initialize phi
  levset::get_phi(phi,xc,yc,config.drop.x,config.drop.y,config.drop.r);

  // ... use phi to initialize rho and nu
  const double rhol = config.iliq.rho;
  const double rhog = config.igas.rho;
  const double nul = config.iliq.mu / rhol;
  const double nug = config.igas.mu / rhog;
  const double sigma = config.drop.sigma;
  const double Mh = config.drop.M*min(dx,dy);
  levset::heaviside(config.drop.M,min(dx,dy),phi,heavi);
  
  // get initial volume to use for volume correction
  double V0 = levset::getVolume(heavi,dx,dy);

  // initialize density and nu fields
  for (int j = jstr-1; j <= jend; j++) {
    for (int i = istr-1; i <= iend; i++) {
      if (config.drop.enabled==false) {
        heavi(i,j) = 1.0;
      }
      rho(i,j) = rhog*heavi(i,j) + rhol*(1.0-heavi(i,j));
      nu(i,j)  = nug*heavi(i,j) + nul*(1.0-heavi(i,j));
      u(i,j) = 0.0;
      p(i,j) = 1.0;
    }
  }
  levset::surfaceTension(Fx,Fy,phi,kappa,Mh,sigma,dx,dy);
  IO::vtk_output_2D("00000",config.fv.dir,config.fv.ghost,
                    xn,yn,u,v,"p",p,"rho",rho,"nu",nu,"phi",phi,
                    "kappa",kappa,"Fx",Fx,"Fy",Fy,"heavi",heavi);
  
  timer.stop();
  if(rank==0) IO::logger->info("  done ({} seconds)",timer.time());
  
  /********************
   * main solver loop *
   *******************/
  // ... initialize doubles
  double ires,res0,res1,cfl0,resmax = 0.0;
  double cfl = config.solver.cfli;
  int finalIter = 0;
  const double dtau = 0.5*min(dx,dy);

  // ... start solver & timer
  timer.start();
  if(rank==0) IO::logger->info("Starting Main Solver");
  cout << rank << " " << istr << endl;
  cout << rank << " " << iend << endl;
  cout << rank << " " << jstr << endl;
  cout << rank << " " << jend << endl;

  MPI_Barrier(MPI_COMM_WORLD);  // ensure all ranks have finished setup/output
  for (int ii = 1; ii <= config.solver.iter; ii++) {
    // ... update boundary conditions
    BC::update_BCs(bcTags,u,v,p);
    levset::surfaceTension(Fx,Fy,phi,kappa,Mh,sigma,dx,dy);
    MPI_Barrier(MPI_COMM_WORLD);  // ensure all ranks have finished setup/output

    // ... get the minimum dt in the domain for current iteration
    double dt = get_min_dt(cfl,dx,dy,u,v,max(nug,nul));


    // ... loop over domain for predictor step
    for (int j = jstr; j <= jend; j++) {
      for (int i = istr; i <= iend; i++) {
        std::vector<double> advec(2,0.0);
        std::vector<double> advec_old(2,0.0);
        std::vector<double> diffu(2,0.0);
        std::vector<double> ab2(2,0.0);

        // get the advection term
        advec[0] = getAdvecU(i,j,rdx,rdy,u,v);
        advec[1] = getAdvecV(i,j,rdx,rdy,u,v);
        
        // get the diffusion term
        diffu[0] = getDiffU(i,j,rdx,rdy,nu,u,v);
        diffu[1] = getDiffV(i,j,rdx,rdy,nu,u,v);

        // AB2 method for convection
        advec_old[0] = getAdvecU(i,j,rdx,rdy,u_old,v_old);
        advec_old[1] = getAdvecV(i,j,rdx,rdy,u_old,v_old);
        ab2[0] = 1.5*advec[0]-0.5*advec_old[0];
        ab2[1] = 1.5*advec[1]-0.5*advec_old[1];

        // predictor step - explicit
        double rhoi = (rho(i,j)+rho(i-1,j))*0.5;
        double rhoj = (rho(i,j)+rho(i,j-1))*0.5;
        double fx = (Fx(i,j)+Fx(i-1,j))*0.5;
        double fy = (Fy(i,j)+Fy(i,j-1))*0.5;
        ustar(i,j) = u(i,j) + dt*(-ab2[0] + diffu[0]-fx/rhoi);
        vstar(i,j) = v(i,j) + dt*(-ab2[1] + diffu[1]-fy/rhoj);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);  // ensure all ranks have finished setup/output
    cout << "after advection" << endl;

    // ... solve for the pressure
    psolve::SOR(config.solver.omega,p,ustar,vstar,dx,dy,dt,rho,bcTags);
    BC::update_BCs(bcTags,ustar,vstar,p);
    MPI_Barrier(MPI_COMM_WORLD);  // ensure all ranks have finished setup/output

    
    // ... apply the pressure correctior
    for (int j = jstr; j <= jend; j++) {
      for (int i = istr; i <= iend; i++) {
        double rhoi = (rho(i,j)+rho(i-1,j))*0.5;
        double rhoj = (rho(i,j)+rho(i,j-1))*0.5;
        double dpdx = (p(i,j) - p(i-1,j)) / (dx);
        double dpdy = (p(i,j) - p(i,j-1)) / (dy);
        u2(i,j) = ustar(i,j) - 1.0/rhoi*dt*dpdx;
        v2(i,j) = vstar(i,j) - 1.0/rhoj*dt*dpdy;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);  // ensure all ranks have finished setup/output
    cout << "after pressure" << endl;
    BC::update_BCs_phi(bcTags,rho);
    cout << "after BC" << endl;
    BC::update_BCs_phi(bcTags,nu);
    cout << "after BC" << endl;
    BC::update_BCs(bcTags,u2,v2,p);
    cout << "after BC" << endl;
    MPI_Barrier(MPI_COMM_WORLD);  // ensure all ranks have finished setup/output
    cout << "after BC" << endl;

    // ... solve advection eq for phi
    double Vn = 0.0;
    if (config.drop.enabled==true) {
      levset::advecPhi(bcTags,dx,dy,dt,u2,v2,phi);
      if (config.levset.reinit) {
        levset::reinitialize(bcTags,dx,dy,dtau,config.levset.ireinit,phi);
      }
      levset::heaviside(config.drop.M,min(dx,dy),phi,heavi);
      Vn = levset::getVolume(heavi,dx,dy);
      double Ln = levset::getLength(phi,dx,dy,Mh);
      levset::volumeCorrection(phi,Mh,V0,Vn,Ln);
    }
    MPI_Barrier(MPI_COMM_WORLD);  // ensure all ranks have finished setup/output
    cout << "after phi advection" << endl;
    
    // ... update density and viscosity values
    for (int j = jstr-1; j <= jend; j++) {
      for (int i = istr-1; i <= iend; i++) {
        rho(i,j) = rhog*heavi(i,j) + rhol*(1.0-heavi(i,j));
        nu(i,j)  = nug*heavi(i,j) + nul*(1.0-heavi(i,j));
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);  // ensure all ranks have finished setup/output

    // ... output intermediate flowviz
    if (config.fv.enabled) {
      if (ii % config.fv.freq == 0) {
        std::ostringstream foutss;
        foutss << setw(5) << std::setfill('0') << ii << "." << rank;
        string caseName = foutss.str();
        if (config.fv.mode=="center") {
          IO::getCellCenter(u,v,uc,vc);
          IO::vtk_output_2D(caseName,config.fv.dir,false,xc,yc,uc,vc,
                            "p",p,"rho",rho,"nu",nu,"phi",phi,"kappa",kappa,
                            "Fx",Fx,"Fy",Fy,"heavi",heavi);
        } else {
          IO::vtk_output_2D(caseName,config.fv.dir,config.fv.ghost,xn,yn,u,v,
                            "p",p,"rho",rho,"nu",nu,"phi",phi,"kappa",kappa,
                            "Fx",Fx,"Fy",Fy,"heavi",heavi);
        }
      }
    } 
    MPI_Barrier(MPI_COMM_WORLD);  // ensure all ranks have finished setup/output
    cout << "after flowviz" << endl;

    // ... Dynamic CFL
    if (ii > 1) 
      res1 = ires;
    double cflb = cfl;
    ires = max(L2NORM(u,u2),L2NORM(v,v2));
    resmax = max(resmax,ires);
    if (ii==1) {
      res0 = ires;
      res1 = ires;
    }
    if (ires == resmax) 
      cfl0 = cfl;
    if (ires < res1 && ires < res0) 
      cfl = cfl0*resmax/ires;
    cfl = max(cfl,cflb);
    cfl = min(config.solver.cflf,max(cfl,config.solver.cfli)); 

    MPI_Barrier(MPI_COMM_WORLD);  // ensure all ranks have finished setup/output
    // ... store the previous timestep and update
    for (int j = jstr-1; j <= jend; j++) {
      for (int i = istr-1; i <= iend; i++) {
        u_old(i,j) = u(i,j);
        v_old(i,j) = v(i,j);
        u(i,j) = u2(i,j);
        v(i,j) = v2(i,j);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);  // ensure all ranks have finished setup/output
    
    // ... calculate residuals
    logFile << ii << " " << ires << "\n";
    if (config.res.enabled) {
      if (ii % config.res.freq == 0) {
        if(rank==0) IO::logger->info("  iter {:04}, cfl: {:4e},dt: {:4e}, res: {:4e}, dV: {:4e}",ii,cfl,dt,ires/res0,V0-Vn);
        logFile.flush();
      }
    }

    // ... check steady state convergence
    if (ires/res0 < config.solver.toler && ii > 1000) {
      finalIter=ii;
      break;
    } else if (ires/res0 > 1.0e20) {
      break;
    }
  } // end of ii-loop

  timer.stop();
  if(rank==0) IO::logger->info("  done ({} seconds)",timer.time());
  if(rank==0) IO::logger->info("Average time / iteration: {} seconds",
                   timer.time()/static_cast<double>(finalIter));

  // ... output section
  if (config.fv.enabled) {
    if(rank==0) IO::logger->info("Outputting final flow solution");
    IO::vtk_output_2D("final",config.fv.dir,false,xc,yc,u,v,
                      "p",p,"rho",rho,"nu",nu,"phi",phi,"heavi",heavi);
  } else {
    IO::logger->warn("Output was disabled.");
  }

  // ... output restart file
  if (config.restart.save) {
    restart::save("u.rank"+std::to_string(rank)+".restart",u);
    restart::save("v.rank"+std::to_string(rank)+".restart",v);
    restart::save("p.rank"+std::to_string(rank)+".restart",p);
  }


  // ... Final run summary output here
  if(rank==0) IO::logger->info("Done!");
  MPI_Finalize();

  return 0;
}
