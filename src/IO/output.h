#ifndef OUTPUT_H
#define OUTPUT_H

#include <matar.h>

using namespace mtr;
using namespace std;

/******************************************************************************/
void vtk_output_2D(const int& ii, const FMatrix<double>& q,
    const FMatrix<double>& xn,
    const FMatrix<double>& yn);
void vtk_output_2D(const string& caseName, const FMatrix<double>& q,
    const FMatrix<double>& xn,
    const FMatrix<double>& yn);
/******************************************************************************/
#endif // OUTPUT_H
