#ifndef GRIDPRO_H
#define GRIDPRO_H

#include <iostream>
#include <matar.h>

namespace gridpro {

int get_nblk(const std::string gridfile);

mtr::FMatrix<int> get_dims(const std::string gridfile);

void reader(const std::string filename,
            const int& nblk,
            mtr::FMatrix<double>& xc,
            mtr::FMatrix<double>& yc);

} // end namespace gridpro

#endif // end GRIDPRO_H
