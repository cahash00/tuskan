/**
 * @file main.cpp
 * @author Caleb Hash
 * @date March 02, 2025
 * @brief Two-Phase Flow Solver - HW4
 * @details first part of the final project - HW4
 */

#include <mesh.h>
#include <output.h>               // output library
#include <input.h>                // input deck reader
#include <matar.h>                // MATAR headers
#include <Kokkos_Core.hpp>        // Kokkos for initialization
#include <yaml-cpp/yaml.h>        // yaml-parser 
#include <argparse/argparse.hpp>  // argparse
#include <pprint.hpp>             // pretty printer 
#include <params.h>
#include <generalUtils.h>
#include <solver.h>
#include <pressure.h>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/pattern_formatter.h>
#include <chrono>
#include <memory>
#include <cmath>
#include <BCs.h>
#include <fstream>
#include <string>

using namespace std;
using namespace mtr;

/**
 * Main program
 */
int main(int argc, char* argv[]){
  // Fire up the printer
  pprint::PrettyPrinter printer;

  /**
   * create the logger
   */
  auto logger = spdlog::stdout_color_mt("console");
  // Apply a custom formatter with relative timestamps
  auto formatter = std::make_unique<spdlog::pattern_formatter>();
  formatter->add_flag<customSPDLOG>('R'); // Use %R for relative time
  formatter->set_pattern("%R [%^%l%$] %v"); // Example: [12.345] [INFO] Message
  logger->set_formatter(std::move(formatter));
  // Set the logger as default
  spdlog::set_default_logger(logger);

  /**
   * Parse out the user input form command line
   */
  argparse::ArgumentParser program("TUSKAN","0.0.0");
  getUserInput(argc,argv,program);
  // assign the input file that was given
  auto inFile = program.get<string>("-i");


  /***************************************************************************
   *                           MAIN PROGRAM START                            *
   **************************************************************************/

  // ... initialize and finalize Kokkos with the main scope
  spdlog::info("Initializing Kokkos");
  Kokkos::ScopeGuard kokkos_guard(argc, argv); 
  spdlog::info(" done");
  
  // ... call input file parser
  parse_user(inFile);
  
  /**
   * parse the yaml file
   */
  double dx,dy = {0.0};
  // domain settings
  double lx = config["domain"]["lengths"]["x"].as<double>();
  double ly = config["domain"]["lengths"]["y"].as<double>();
  int nx    = config["domain"]["dimensions"]["x"].as<int>();
  int ny    = config["domain"]["dimensions"]["y"].as<int>();
  // solver settings
  int iter     = config["solver"]["iterations"].as<int>();
  double cfl   = config["solver"]["CFL"].as<double>();
  bool fvflag  = config["output"]["flowviz"]["enabled"].as<bool>();
  int fvfreq   = config["output"]["flowviz"]["frequency"].as<int>();
  bool resflag = config["output"]["residuals"]["enabled"].as<bool>();
  int resfreq  = config["output"]["residuals"]["frequency"].as<int>();
  string resFile  = config["output"]["residuals"]["file"].as<string>();
  string pMethod = config["solver"]["pressure solver"]["method"].as<string>();
  int pIter = config["solver"]["pressure solver"]["iterations"].as<int>();
  if (pMethod == "SOR") {
    double sorWeight = config["solver"]["pressure solver"]["SOR weight"].as<double>(); 
  }
  // convergence criteria
  double toler   = config["convergence"]["residual"].as<double>();
  double cfli    = config["dynamic CFL"]["cfli"].as<double>();
  double cflf    = config["dynamic CFL"]["cflf"].as<double>();
  bool dcfl      = config["dynamic CFL"]["enabled"].as<bool>();

  printer.print(config["case name"]);

  // ... get domain stats
  getDomainIndices(nx,ny);
  vector<int> ndims(2,0);
  ndims[0] = nx+nghosts*2;
  ndims[1] = ny+nghosts*2;

  // ... initialize the MATAR matrices for the domain
  FMatrix<double> xc(ndims[0],ndims[1]),
                  yc(ndims[0],ndims[1]),
                  xn(ndims[0]+1,ndims[1]+1),
                  yn(ndims[0]+1,ndims[1]+1);
  FMatrix<double> p(ndims[0]+1,ndims[1]+1);
  FMatrix<double> ustar(ndims[0]+1,ndims[1]+1);
  FMatrix<double> vstar(ndims[0]+1,ndims[1]+1);
  FMatrix<double> u(ndims[0]+1,ndims[1]+1);
  FMatrix<double> v(ndims[0]+1,ndims[1]+1);
  FMatrix<double> u2(ndims[0]+1,ndims[1]+1);
  FMatrix<double> v2(ndims[0]+1,ndims[1]+1);
  FMatrix<double> uexact(ndims[0]+1,ndims[1]+1);
  
  // ... call mesh generator
  Timer timer;
  timer.start();
  spdlog::info("Generating 2D mesh");
  mesher2D(lx,ly,nx,ny,xc,yc,xn,yn,dx,dy);
  timer.stop();
  spdlog::info("  done ({} seconds)",timer.time());

  // ... initialization
  ofstream logFile(resFile, ios::out);
  if (!logFile) {
      std::cerr << "Error opening log file!\n";
      return -1;
  }
  timer.start();
  spdlog::info("Initializing the domain");
  double rho  = 1.0e3;   // density
  double rrho = 1.0e-3;  // reciprocal of density
  double nu   = 1.0e-6;    // kinematic viscosity m^2/s
  double mu   = nu*rho;  // dynamic viscosity
  double rdx  = 1.0/dx;  // reciprocal of dx
  double rdy  = 1.0/dy;  // reciprocal of dx
  double dt   = {0};     // initialize dt
  double dpdx = -0.3;    // analytical solution for dpdx
  spdlog::info("  dx: {},dy: {}",dx,dy);

  // initialize domain and calculate exact solution
  initialize_solution(u,v,u2,ustar,vstar,p);
  timer.stop();
  spdlog::info("  done ({} seconds)",timer.time());
  
  /**********
   * main solver loop
   **********/
  double ires,res0,res1,cfl0,resmax = {0.0};
  timer.start();
  spdlog::info("Starting Main Solver");
  for (int ii = 0; ii < iter; ii++) {

    // ... update boundary conditions
    bc_noslip(u);
    bc_periodic(u);

    // ... get the minimum dt in the domain for current iteration
    dt = get_min_dt(cfl,dx,u);

    // ... loop over domain for predictor step
    for (int j = jstr; j <= jend; j++) {
      for (int i = istr; i <= iend; i++) {
        // get the advection term
        double advec = getAdvec(i,j,rdx,rdy,u,v);
        // get the diffusion term
        double diffu = getDiffu(i,j,rdx,rdy,u);
        // predictor step
        ustar(i,j) = u(i,j) + dt*(-advec + nu*diffu);
      }
    }

    // ... solve for the pressure correction term
    // set the ghost cells for the ustar
    bc_noslip(ustar);
    bc_periodic(ustar);
    // p.set_values(0.0);
    if (pMethod == "Jacobi") {
      // psolve::Jacobi(p,ustar,vstar,dx,dy,dt,rho,nx,ny);
      psolve::SOR(p,ustar,vstar,dx,dy,dt,rho,nx,ny);
    }
    
    // ... apply the pressure correctior
    for (int j = jstr; j <= jend; j++) {
      for (int i = istr; i <= iend; i++) {
        dpdx = (p(i,j) - p(i-1,j)) / (dx);
        u2(i,j) = ustar(i,j) - rrho*dt*dpdx;
      }
    }

    // ... output intermediate flowviz
    if (fvflag) {
      if (ii % fvfreq == 0) {
        vtk_output_2D(ii,u,xn,yn);
      }
    } 
    
    // ... Dyanmic CFL
    if (dcfl) {
      if (ii > 0) res1 = ires;
      double cflb = cfl; // store current cfl
      ires = L2NORM(u,u2,nx*ny);
      resmax = max(resmax,ires);
      if (ii==0) {
        res0 = ires;
        res1 = ires;
      }
      if (ires == resmax) cfl0 = cfl; // if res is higher, keep
      if (ires < res1 && ires < res0) {
        cfl = cfl0*resmax/ires; // if res is lower, increase CFL
      }
      cfl = max(cfl,cflb);
      cfl = min(cflf,max(cfl,cfli));
    } else {
      ires = L2NORM(u,u2,nx*ny);
      if (ii==0) res0 = ires;
    }
    // ... update the u array with the updated solution array
    for (int j = jstr; j <= jend; j++) {
      for (int i = istr; i <= iend; i++) {
        u(i,j) = u2(i,j);
      }
    }
    
    // ... calculate residuals
    logFile << ii << " " << ires << "\n";
    if (resflag) {
      if (ii % resfreq == 0) {
        spdlog::info("  iter {:04}, cfl: {:5e}, res: {:5e}",ii,cfl,ires/res0);
        logFile.flush();
      }
    }

    // exit if converged
    if (ires/res0 < toler && ii > 1000) {
      iter=ii;
      break;
    }
  } // end of ii-loop
  timer.stop();
  spdlog::info("  done ({} seconds)",timer.time());
  spdlog::info("Average time / iteration: {} seconds",
               timer.time()/static_cast<double>(iter));

  /**
   * output section
   */
  if (fvflag) {
    spdlog::info("Outputting final flow solution");
    vtk_output_2D(string("final"),u,xn,yn);
    // output the values along the channel
    ofstream fout("compare.dat",ios::out);
    DO_LOOP(j,jstr-nghosts,jend+nghosts,{
      DO_LOOP(i,istr-nghosts,iend+nghosts,{
        uexact(i,j) = 1.0/(2.0*mu)*-0.3*(yc(i,j)*yc(i,j)-ly*yc(i,j));
      });
    });
    DO_LOOP(j,jstr-1,jend+1,{
      fout << yc(nx/2,j) << " " << u(nx/2,j) 
           << " " << uexact(nx/2,j) << endl;
    });
    spdlog::info("Outputting exact flow solution");
    vtk_output_2D(string("exact"),uexact,xn,yn);
  } else {
    spdlog::warn("Output was disabled.");
  }

  /**
   * Final run summary output here
   */
  spdlog::info("Done!");

  return 0;
}
