#ifndef OUTPUT_H
#define OUTPUT_H

#include <matar.h>

using namespace std;

/******************************************************************************/
void vtk_output_2D(const int& ii,
                   const string& foutDir,
                   const mtr::FMatrix<double>& q,
                   const int& nx,
                   const int& ny,
                   const mtr::FMatrix<double>& xn,
                   const mtr::FMatrix<double>& yn);
void vtk_output_2D(const string& caseName, 
                   const string& foutDir,
                   const mtr::FMatrix<double>& q,
                   const int& nx,
                   const int& ny,
                   const mtr::FMatrix<double>& xn,
                   const mtr::FMatrix<double>& yn);
/******************************************************************************/
#endif // OUTPUT_H
