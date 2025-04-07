/**
 * @file restart.cpp
 * @brief restart library for the solver. Writes in legacy vtk format for 
 *        structured grids
 * @author Caleb Hash
 * @date April 07, 2025
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

namespace restart {
/******************************************************************************/
// Function to save the matrix to a binary file
void save(const std::string& filename, 
          const mtr::FMatrix<double>& u) {
  std::ofstream file(filename, std::ios::binary);
  if (!file.is_open()) {
    std::cerr << "Error: Could not open file " << filename << " for writing.\n";
    return;
  }

  for (int j = jstr; j <= jend; j++) {
    for (int i = istr; i <= iend; i++) {
      file << u(i,j) << " ";
    }
    file << "\n";
  }
  file.close();
  std::cout << "Matrix saved to " << filename << std::endl;
}

// Function to load the matrix from a binary file
void load(const std::string& filename,
          const mtr::FMatrix<double>& u) {
  std::ifstream file(filename, std::ios::binary);
  if (!file.is_open()) {
    std::cerr << "Error: Could not open file " << filename << " for reading.\n";
  }

  for (int j = jstr; j <= jend; j++) {
    for (int i = istr; i <= iend; i++) {
      file >> u(i,j);
    }
  }
  
  file.close();
  std::cout << "Matrix loaded from " << filename << std::endl;
}

} // end namepsace IO
