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
#include <cmath>
#include <BCs.h>

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
  double dx,dy;
  double lx = config["domain"]["lengths"]["x"].as<double>();
  double ly = config["domain"]["lengths"]["y"].as<double>();
  int nx = config["domain"]["dimensions"]["x"].as<int>();
  int ny = config["domain"]["dimensions"]["y"].as<int>();
  int iter = config["solver"]["iterations"].as<int>();
  double cfl = config["solver"]["CFL"].as<double>();
  double toler = config["solver"]["tolerance"].as<double>();
  bool fvflag = config["output"]["flowviz"]["enabled"].as<bool>();
  int fvfreq = config["output"]["flowviz"]["frequency"].as<int>();
  bool resflag = config["output"]["residuals"]["enabled"].as<bool>();
  int resfreq = config["output"]["residuals"]["frequency"].as<int>();

  getDomainIndices(nx,ny);

  // get the TOTAL number of cells here so we can allocate properly
  vector<int> ndims(2,0);
  ndims[0] = nx+nghosts*2;
  ndims[1] = ny+nghosts*2;

  // initialize the Kokkos matrices for the domain
  FMatrix<double> xc(ndims[0],ndims[1]),
                  yc(ndims[0],ndims[1]),
                  xn(ndims[0]+1,ndims[1]+1),
                  yn(ndims[0]+1,ndims[1]+1);

  // initialize flow solution variables - staggered mesh method
  // q(1)  rho
  // q(2)  u
  // q(3)  v
  // q(4)  w
  // q(5)  temperature
  // q(6)  pressure
  FMatrix<double> q(2,ndims[0]+1,ndims[1]+1);
  FMatrix<double> q2(2,ndims[0]+1,ndims[1]+1);
  FMatrix<double> uexact(1,ndims[0]+1,ndims[1]+1);
  
  // call mesh generator
  Timer timeMe;
  timeMe.start();
  spdlog::info("Generating 2D mesh");
  mesher2D(lx,ly,nx,ny,xc,yc,xn,yn,dx,dy);
  timeMe.stop();
  spdlog::info("  Done ({} seconds)",timeMe.time());

  // ... initialization
  timeMe.start();
  spdlog::info("Initializing the domain");
  double rho = 1.0e3;  // density
  double rrho = 1.0e-3; // reciprocal of density
  double nu = 1e-6;    // kinematic viscosity m^2/s
  double mu = nu*rho;  // dynamic viscosity
  double rdx = 1.0/dx; // reciprocal of dx
  double rdy = 1.0/dy; // reciprocal of dx
  spdlog::info("  dx: {},dy: {}",dx,dy);
  double dt = {0};
  double dpdx = -0.3;

  // initialize domain and calculate exact solution
  DO_LOOP(j,jstr-nghosts,jend+nghosts,{
    DO_LOOP(i,istr-nghosts,iend+nghosts,{
      printer.print(i,j);
      q(1,i,j) = 0.0002; // average u
      q(2,i,j) = 0.0; // zero y-velocity
      q2(1,i,j) = 0.0;
      q2(2,i,j) = 0.0;
      uexact(1,i,j) = 1.0/(2.0*mu)*dpdx*(yn(i,j)*yn(i,j)-ly*yn(i,j)); 
    });
  });
  vtk_output_2D(777,uexact,xn,yn);
  printer.print("Here after vtk out");

  timeMe.stop();
  spdlog::info("  Done ({} seconds)",timeMe.time());
  
  /**
   * main solver loop
   * Will loop over and solve until L2Norm is met
   */
  timeMe.start();
  spdlog::info("Starting Main Solver");
  for (int ii = 0; ii < iter; ii++) {

    // ... update boundary conditions
    bc_noslip(q);
    bc_periodic(q);

    // ... get the minimum dt in the domain for current iteration
    double minval = 1e5;
    DO_LOOP(i,istr,iend,{
      DO_LOOP(j,jstr,jend,{
        dt = min(minval,cfl*dx/abs(q(1,i,j)));
      });
    });

    // ... loop over domain
    DO_LOOP(j,jstr,jend,{
      DO_LOOP(i,istr,iend,{
        // get the advection term
        double advec = getAdvec(i,j,rdx,rdy,q);
        // get the diffusion term
        double diffu = getDiffu(i,j,rdx,rdy,q);
        // Predictor step
        double ustar = q(1,i,j) + dt*(-advec + nu*diffu);
        // get the u^n+1 values
        q2(1,i,j) = ustar - rrho*dt*dpdx;
      }); // end of i-loop
    }); // end of j-loop

    // ... output intermediate flowviz
    if (fvflag) {
      if (ii % fvfreq == 0) {
        vtk_output_2D(ii,q,xn,yn);
      }
    } 
    // ... update the q array with the updated solution array
    double resid = L2NORM(q,q2,nx*ny);
    DO_LOOP(j,jstr-nghosts,jend+nghosts,{
      DO_LOOP(i,istr-nghosts,iend+nghosts,{
        q(1,i,j) = q2(1,i,j);
      }); // end of i-loop
    }); // end of j-loop

    // ... calculate residuals
    double l2norm = L2NORM(uexact,q,nx*ny);
    if (resflag) {
      printer.print("here after resflag");
      if (ii % resfreq == 0) {
        spdlog::info("  iter {:04}, comp: {:5e}, res: {:5e}",ii,l2norm,resid);
      }
    }

    // exit if converged
    if (resid < toler || l2norm < toler) {
      iter=ii;
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
  if (fvflag) {
    spdlog::info("Outputting final flow solution");
    string fout = "final";
    vtk_output_2D(fout,q,xn,yn);
  } else {
    spdlog::warn("Output was disabled.");
  }

  /**
   * Final run summary output here
   */
  spdlog::info("Done!");

  return 0;
}
