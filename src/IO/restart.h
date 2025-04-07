#ifndef RESTART_H
#define RESTART_H
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
          const mtr::FMatrix<double>& u);

// Function to load the matrix from a binary file
void load(const std::string& filename,
          const mtr::FMatrix<double>& u);

} // end namepsace restart

#endif // end RESTART_H
