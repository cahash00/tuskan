/**
 * @file main.cpp
 * @author Caleb Hash
 * @date March 02, 2025
 * @brief Two-Phase Flow Solver - HW4
 * @details first part of the final project - HW4
 */

#include <mesh.h>
#include <output.h>               // output library
#include <args.h>
#include <input.h>                // input deck reader
#include <matar.h>                // MATAR headers
#include <Kokkos_Core.hpp>        // Kokkos for initialization
#include <argparse/argparse.hpp>  // argparse
#include <pprint.hpp>             // pretty printer 
#include <params.h>
#include <generalUtils.h>
#include <solver.h>
#include <pressure.h>
#include <logger.h>
#include <cmath>
#include <BCs.h>
#include <fstream>
#include <string>

using namespace std;
typedef mtr::FMatrix<double> fmatD;
typedef mtr::FMatrix<int> fmatI;

/**
 * Main program
 */
int main(int argc, char* argv[]){
  // ... Fire up the printer
  pprint::PrettyPrinter printer;

  // ... create the logger
  IO::init_logger();

  // ... Parse out the user input form command line
  argparse::ArgumentParser program("TUSKAN","0.0.0");
  IO::getUserInput(argc,argv,program);
  // assign the input file that was given
  auto inFile = program.get<string>("-i");


  /***************************************************************************
   *                           MAIN PROGRAM START                            *
   **************************************************************************/

  // ... initialize and finalize Kokkos with the main scope
  IO::logger->info("Initializing Kokkos");
  Kokkos::ScopeGuard kokkos_guard(argc, argv); 
  IO::logger->info("  done");
  
  // ... parse the YAML input deck
  IO::ConfigData config = IO::parseInputDeck(inFile);
  // check to see if the output directory exists
  IO::check_directories(config.foutDir);

  
  // ... get domain stats
  getDomainIndices(config.nx,config.ny);
  vector<int> ndims(2,0);
  ndims[0] = nx+nghosts*2;
  ndims[1] = ny+nghosts*2;

  // ... initialize the MATAR matrices for the domain
  fmatD xc(ndims[0],ndims[1]),yc(ndims[0],ndims[1]),
        xn(ndims[0]+1,ndims[1]+1),yn(ndims[0]+1,ndims[1]+1);
  fmatD p(ndims[0]+1,ndims[1]+1);
  fmatD phi(ndims[0],ndims[1]);
  fmatD mark(ndims[0],ndims[1]);
  fmatD rho(ndims[0]+1,ndims[1]+1);
  fmatD ustar(ndims[0]+1,ndims[1]+1);
  fmatD vstar(ndims[0]+1,ndims[1]+1);
  fmatD u(ndims[0]+1,ndims[1]+1);
  fmatD v(ndims[0]+1,ndims[1]+1);
  fmatD u_old(ndims[0]+1,ndims[1]+1);
  fmatD v_old(ndims[0]+1,ndims[1]+1);
  fmatD u2(ndims[0]+1,ndims[1]+1);
  fmatD v2(ndims[0]+1,ndims[1]+1);
  fmatD uexact(ndims[0]+1,ndims[1]+1);
  
  // ... call mesh generator
  Timer timer;
  timer.start();
  IO::logger->info("Generating 2D mesh");
  double dx,dy = 0.0;
  mesh::mesher2D(config.lx,config.ly,xc,yc,xn,yn,dx,dy);
  timer.stop();
  IO::logger->info("  done ({} seconds)",timer.time());
  IO::logger->info("Tagging boundaries");
  BC::bcTags bcTags = BC::tag_BCs(config,u.dims(1),u.dims(2));
  IO::logger->info("  done");

  // ... initialization
  ofstream logFile(config.resFile, ios::out);
  if (!logFile) {
      std::cerr << "Error opening log file!\n";
      return -1;
  }
  timer.start();
  IO::logger->info("Initializing the domain");
  rho.set_values(1.0e3); // density
  double nu   = 1.0e-6;  // kinematic viscosity m^2/s
  double rdx  = 1.0/dx;  // reciprocal of dx
  double rdy  = 1.0/dy;  // reciprocal of dx
  IO::logger->info("  dx: {},dy: {}",dx,dy);

  // initialize domain
  initialize_solution(config,
                      u,v,
                      u2,v2,
                      u_old,v_old,
                      ustar,vstar,
                      p);
  BC::update_BCs(bcTags,u,v,p);
  
  timer.stop();
  IO::logger->info("  done ({} seconds)",timer.time());
  
  /********************
   * main solver loop *
   *******************/
  // ... initialize doubles
  double ires,res0,res1,cfl0,resmax = 0.0;
  double cfl = config.cfli;
  int finalIter = 0;

  // ... start solver & timer
  timer.start();
  IO::logger->info("Starting Main Solver");
  for (int ii = 0; ii < config.iter; ii++) {
    // ... update boundary conditions
    BC::update_BCs(bcTags,u,v,p);

    // ... get the minimum dt in the domain for current iteration
    double dt = get_min_dt(cfl,dx,dy,u,v,nu);

    // ... loop over domain for predictor step
    for (int j = jstr; j <= jend; j++) {
      for (int i = istr; i <= iend; i++) {
        std::vector<double> advec(2,0.0);
        std::vector<double> advec_old(2,0.0);
        std::vector<double> diffu(2,0.0);
        std::vector<double> ab2(2,0.0);
        // get the advection term
        advec[0] = getAdvecU(i,j,rdx,rdy,u,v);
        // advec[1] = getAdvecV(i,j,rdx,rdy,u,v);
        // get the diffusion term
        diffu[0] = getDiffU(i,j,rdx,rdy,u,v);
        // diffu[1] = getDiffV(i,j,rdx,rdy,u,v);
        // AB2 method for convection
        advec_old[0] = getAdvecU(i,j,rdx,rdy,u_old,v_old);
        // advec_old[1] = getAdvecV(i,j,rdx,rdy,u_old,v_old);
        ab2[0] = 1.5*advec[0]-0.5*advec_old[0];
        // ab2[1] = 1.5*advec[1]-0.5*advec_old[1];
        // predictor step - explicit
        ustar(i,j) = u(i,j) + dt*(-ab2[0] + nu*diffu[0]);
        // vstar(i,j) = v(i,j) + dt*(-ab2[1] + nu*diffu[1]);
      }
    }
    BC::update_BCs(bcTags,ustar,vstar,p);

    // ... solve for the pressure
    psolve::SOR(config.sorOmega,p,ustar,vstar,dx,dy,dt,rho,bcTags);
    
    // ... apply the pressure correctior
    for (int j = jstr; j <= jend; j++) {
      for (int i = istr; i <= iend; i++) {
        double dpdx = (p(i,j) - p(i-1,j)) / (dx);
        double dpdy = (p(i,j) - p(i,j-1)) / (dy);
        u2(i,j) = ustar(i,j) - 1.0/rho(i,j)*dt*dpdx;
        // v2(i,j) = vstar(i,j) - 1.0/rho(i,j)*dt*dpdy;
      }
    }
    BC::update_BCs(bcTags,u2,v2,p);

    // ... output intermediate flowviz
    if (config.fvflag) {
      if (ii % config.fvfreq == 0) {
        IO::vtk_output_2D_node(ii,config.foutDir,config.ghost,xn,yn,u,v,p,phi);
      }
    } 

    // ... Dyanmic CFL
    if (ii > 0) 
      res1 = ires;
    double cflb = cfl;
    ires = max(L2NORM(u,u2),L2NORM(v,v2));
    resmax = max(resmax,ires);
    if (ii==0) {
      res0 = ires;
      res1 = ires;
    }
    if (ires == resmax) 
      cfl0 = cfl; // if res is higher, keep
    if (ires < res1 && ires < res0) 
      cfl = cfl0*resmax/ires;
    cfl = max(cfl,cflb);
    cfl = min(config.cflf,max(cfl,config.cfli)); 

    // ... store the previous timestep
    for (int j = jstr-1; j <= jend+1; j++) {
      for (int i = istr-1; i <= iend+1; i++) {
        u_old(i,j) = u(i,j);
        // v_old(i,j) = v(i,j);
      }
    }
    
    // ... update the u array with the updated solution array
    for (int j = jstr; j <= jend; j++) {
      for (int i = istr; i <= iend; i++) {
        u(i,j) = u2(i,j);
        // v(i,j) = v2(i,j);
      }
    }
    
    // ... calculate residuals
    logFile << ii << " " << ires << "\n";
    if (config.resflag) {
      if (ii % config.resfreq == 0) {
        IO::logger->info("  iter {:04}, cfl: {:5e},dt: {:5e}, res: {:5e}",ii,cfl,dt,ires/res0);
        logFile.flush();
      }
    }

    // ... check steady state convergence
    if (ires/res0 < config.toler && ii > 1000) {
      finalIter=ii;
      break;
    } else if (ires/res0 > 1.0e20) {
      break;
    }
  } // end of ii-loop
  timer.stop();
  IO::logger->info("  done ({} seconds)",timer.time());
  IO::logger->info("Average time / iteration: {} seconds",
                   timer.time()/static_cast<double>(finalIter));

  /**
   * output section
   */
  if (config.fvflag) {
    IO::logger->info("Outputting final flow solution");
    IO::vtk_output_2D_node(string("final"),config.foutDir,config.ghost,xn,yn,u,v,p,phi);
  } else {
    IO::logger->warn("Output was disabled.");
  }

  /**
   * Final run summary output here
   */
  IO::logger->info("Done!");

  return 0;
}
