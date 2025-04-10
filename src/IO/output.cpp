/**
 * @file output.cpp
 * @brief output library for the solver. Writes in legacy vtk format for 
 *        structured grids
 * @todo add support for solution array eventually.
 * @author Caleb Hash
 * @date 2025-02-22
 */

#include <input.h>
#include <output.h>
#include <matar.h>
#include <pprint.hpp>
#include <string>
#include <format>
#include <iostream>
#include <params.h>
#include <input.h>
#include <fstream>

using namespace std;

namespace IO {

void check_directories(const string& foutDir) {
  if (fs::exists(foutDir)) {
    // clear it
    for (const auto& entry : fs::directory_iterator(foutDir)) {
      fs::remove_all(entry.path());
    }
    IO::logger->info("  {} has been cleared.",foutDir);
  } else {
    fs::create_directories(foutDir);
    IO::logger->info("  {} has been created.",foutDir);
  }
} // end check_directories

void getCellCenter(mtr::FMatrix<double>& u,
    mtr::FMatrix<double>& v,
    mtr::FMatrix<double>& uc,
    mtr::FMatrix<double>& vc) {
  for (int j = jstr; j <= jend; j++) {
    for (int i = istr; i <= iend; i++) {
      uc(i,j) = 0.5*(u(i,j) + u(i-1,j));
    }
  }
  for (int j = jstr; j <= jend; j++) {
    for (int i = istr; i <= iend; i++) {
      vc(i,j) = 0.5*(v(i,j)+v(i,j-1));
    }
  }
}

} // end namepsace IO
