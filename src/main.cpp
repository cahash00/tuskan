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
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/pattern_formatter.h>
#include <chrono>
#include <memory>

using namespace std;
using namespace mtr;

// Custom formatter for relative timestamps
class RelativeTimeFormatter : public spdlog::custom_flag_formatter {
private:
    std::chrono::steady_clock::time_point start_time;

public:
    RelativeTimeFormatter() : start_time(std::chrono::steady_clock::now()) {}

    void format(const spdlog::details::log_msg&,
                const std::tm&,
                spdlog::memory_buf_t& dest) override {
        auto now = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - start_time);
        
        // Correct usage: append formatted string to dest
        fmt::format_to(std::back_inserter(dest), "[{:03}.{:03}]", duration.count() / 1000, duration.count() % 1000);
    }

    std::unique_ptr<custom_flag_formatter> clone() const override {
        return std::make_unique<RelativeTimeFormatter>();
    }
};

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
  formatter->add_flag<RelativeTimeFormatter>('R'); // Use %R for relative time
  formatter->set_pattern("%R [%^%l%$] %v"); // Example: [12.345] [INFO] Message
  logger->set_formatter(std::move(formatter));
  // Set the logger as default
  spdlog::set_default_logger(logger);

  /**
   * Parse out the user input form command line
   */
  argparse::ArgumentParser program("TUSKAN","0.0.0");
  getUserInput(argc,argv,program);

  /***************************************************************************
   *                           MAIN PROGRAM START                            *
   **************************************************************************/
  // splashScreen();

  // assign the input file that was given
  auto inFile = program.get<string>("-i");

  // automatically initialize and finalize Kokkos with the main program
  Kokkos::ScopeGuard kokkos_guard(argc, argv); 
  spdlog::info("Kokkos has been initialized");
  
  // call input file parser
  parse_user(inFile);
  
  // generate the domain
  double dx,dy,dz;
  double lx = config["domain"]["lengths"]["x"].as<double>();
  double ly = config["domain"]["lengths"]["y"].as<double>();
  double lz = config["domain"]["lengths"]["z"].as<double>();
  int nx = config["domain"]["dimensions"]["x"].as<int>();
  int ny = config["domain"]["dimensions"]["y"].as<int>();
  int nz = config["domain"]["dimensions"]["z"].as<int>();
  int iter = config["solver"]["iterations"].as<int>();
  double cfl = config["solver"]["CFL"].as<double>();
  double toler = config["solver"]["tolerance"].as<double>();

  getDomainIndices(nx,ny,nz);

  // get the TOTAL number of cells here so we can allocate properly
  vector<int> ndims(3,0);
  ndims[0] = nx+nghosts*2;
  ndims[1] = ny+nghosts*2;
  ndims[2] = nz+nghosts*2;

  // initialize the Kokkos matrices for the domain
  FMatrix<double> xc(ndims[0],ndims[1],ndims[2]),
                  yc(ndims[0],ndims[1],ndims[2]),
                  zc(ndims[0],ndims[1],ndims[2]),
                  xn(ndims[0]+1,ndims[1]+1,ndims[2]+1),
                  yn(ndims[0]+1,ndims[1]+1,ndims[2]+1),
                  zn(ndims[0]+1,ndims[1]+1,ndims[2]+1);

  // initialize flow solution variables - staggered mesh method
  // q(1)  rho
  // q(2)  u
  // q(3)  v
  // q(4)  w
  // q(5)  temperature
  // q(6)  pressure
  FMatrix<double> q(6,ndims[0]+1,ndims[1]+1,ndims[2]+1);
  FMatrix<double> q2(6,ndims[0]+1,ndims[1]+1,ndims[2]+1);
  FMatrix<double> uexact(1,ndims[0]+1,ndims[1]+1,ndims[2]+1);
  int ii,jj,kk;

  
  // call mesh generator
  Timer timeMe;
  timeMe.start();
  if (nz==1) {
    spdlog::info("Generating 2D mesh");
    mesher3D(lx,ly,lz,nx,ny,nz,xc,yc,zc,xn,yn,zn,dx,dy,dz);
  } else {
    spdlog::info("Generating 3D mesh");
    mesher3D(lx,ly,lz,nx,ny,nz,xc,yc,zc,xn,yn,zn,dx,dy,dz);
  }
  printer.print(xn(5,5,3));
  timeMe.stop();
  spdlog::info("  Done ({} seconds)",timeMe.time());

  // ... initialization
  timeMe.start();
  spdlog::info("Initializing the domain");
  double rho = 1.0e3;  // density
  double rrho = 1.0e3; // reciprocal of density
  double nu = 1e-6;    // kinematic viscosity m^2/s
  double mu = nu*rho;  // dynamic viscosity
  double rdx = 1.0/dx; // reciprocal of dx
  double rdy = 1.0/dy; // reciprocal of dx
  double dt = 0.01;

  DO_ALL(k,kstr-nghosts,kend+nghosts,
         j,jstr-nghosts,jend+nghosts,
         i,istr-nghosts,iend+nghosts,{
    q(1,i,j,k) = 10.0*ly; // average u
    q(2,i,j,k) = 0.0; // zero y-velocity
    q(3,i,j,k) = 0.0; // zero z-velocity
    uexact(1,i,j,k) = -60.0/ly/ly*(yn(i,j,k)*yn(i,j,k)-ly*yn(i,j,k)); 
  });
  vtk_output_3D(777,uexact,xn,yn,zn);
  q2 = 0.0;
  q2 = 0.0;
  q2 = 0.0;
  timeMe.stop();
  spdlog::info("  Done ({} seconds)",timeMe.time());
  
  /**
   * main solver loop
   * Will loop over and solve until L2Norm is met
   */
  timeMe.start();
  spdlog::info("Starting Main Solver");
  for (int ii = 0; ii < iter; ii++) {
    // enforce boundary conditions
    DO_ALL(k,kstr-nghosts,kend+nghosts,
           j,jstr-nghosts,jend+nghosts,{
      q(1,istr-1,j,k) = q(1,iend,j,k);
      q(1,iend+1,j,k) = q(1,istr,j,k);
    });
    DO_ALL(i,istr-nghosts,iend+nghosts,
           k,kstr-nghosts,kend+nghosts,{
      q(1,i,jstr-1,k) = -q(1,i,jstr,k);
      q(1,i,jend+1,k) = -q(1,i,jend,k);

    });

    // enforce boundary conditions
    DO_LOOP(k,kstr,kend,{
      DO_LOOP(j,jstr,jend,{
        DO_LOOP(i,istr,iend,{
          // calculate the timestep based on CFL
          dt = 1.0e-6;
          // get the pressure gradient term
          double dpdx = -120.0*nu*rho/(ly*ly);

          // get the advection term
          double advec = getAdvec(i,j,k,rdx,rdy,q);
          // get the diffusion term
          double diffu = getDiffu(i,j,k,rdx,rdy,q);
          // get u*
          double ustar = q(1,i,j,k) + dt*(-advec + nu*diffu);
          
          // get the u^n+1 values
          q2(1,i,j,k) = ustar - rrho*dt * dpdx;
          // printer.print(dpdx,advec,diffu,ustar,q2(1,i,j,k));
        }); // end of i-loop
      }); // end of j-loop
    }); // end of k-loop

    // output intermediate flowviz
    if (config["output"]["enabled"].as<bool>()) {
      if ((ii % config["output"]["frequency"].as<int>()) == 0) {
        if (nz==1) {
          spdlog::info("Outputting 2D flow solution");
          vtk_output_2D(ii,q,xn,yn);
        } else {
          spdlog::info("Outputting 3D flow solution");
          vtk_output_3D(ii,q,xn,yn,zn);
        }
      }
    } 
    // update the q array with the updated solution array
    q = q2;
    double l2norm = L2NORM(uexact,q,nx*ny*nz);
    spdlog::info("  iter {:04}, res: {}",ii,l2norm);
    if (l2norm<=toler) {
      iter==ii;
      break;
    }
  } // end of ii-loop
  timeMe.stop();
  spdlog::info("  Done ({} seconds)",timeMe.time());
  spdlog::info("Average time / iteration: {} seconds",
               timeMe.time()/static_cast<double>(iter));

  /**
   * output section
   */
  if (config["output"]["enabled"].as<bool>()) {
    if (nz==1) {
      spdlog::info("Outputting 2D flow solution");
      vtk_output_2D(7777,uexact,xn,yn);
    } else {
      spdlog::info("Outputting 3D flow solution");
      vtk_output_3D(7778,q,xn,yn,zn);
    }
    spdlog::info("Output has finished.");
  } else {
    spdlog::warn("Output was disabled.");
  }

  /**
   * Final run summary output here
   */
  spdlog::info("Done!");

  return 0;
}
