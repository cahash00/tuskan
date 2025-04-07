#ifndef OUTPUT_H
#define OUTPUT_H

#include <matar.h>
#include <iostream>
#include <filesystem>
#include <string>
#include <logger.h>

using namespace std;
namespace fs = std::filesystem;

namespace IO {
/******************************************************************************/
void check_directories(const string& foutDir);

void vtk_output_2D_node(const int& ii,
                        const string& foutDir,
                        const bool ghost,
                        const mtr::FMatrix<double>& xc,
                        const mtr::FMatrix<double>& yc,
                        const mtr::FMatrix<double>& u,
                        const mtr::FMatrix<double>& v,
                        const mtr::FMatrix<double>& p,
                        const mtr::FMatrix<double>& phi,
                        const mtr::FMatrix<double>& mark);
void vtk_output_2D_node(const string& caseName, 
                        const string& foutDir,
                        const bool ghost,
                        const mtr::FMatrix<double>& xc,
                        const mtr::FMatrix<double>& yc,
                        const mtr::FMatrix<double>& u,
                        const mtr::FMatrix<double>& v,
                        const mtr::FMatrix<double>& p,
                        const mtr::FMatrix<double>& phi,
                        const mtr::FMatrix<double>& mark);
/******************************************************************************/
} // end namespace IO
#endif // OUTPUT_H
