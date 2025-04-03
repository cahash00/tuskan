#ifndef OUTPUT_H
#define OUTPUT_H

#include <matar.h>

using namespace std;

namespace IO {
/******************************************************************************/
void vtk_output_2D(const int& ii,
                   const string& foutDir,
                   const mtr::FMatrix<double>& u,
                   const mtr::FMatrix<double>& v,
                   const mtr::FMatrix<double>& p,
                   const int& nx,
                   const int& ny,
                   const mtr::FMatrix<double>& xn,
                   const mtr::FMatrix<double>& yn);
void vtk_output_2D(const string& caseName, 
                   const string& foutDir,
                   const mtr::FMatrix<double>& u,
                   const mtr::FMatrix<double>& v,
                   const mtr::FMatrix<double>& p,
                   const int& nx,
                   const int& ny,
                   const mtr::FMatrix<double>& xn,
                   const mtr::FMatrix<double>& yn);
/******************************************************************************/
} // end namespace IO
#endif // OUTPUT_H
