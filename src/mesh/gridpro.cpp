
#include <gridpro.h>
#include <matar.h>
#include <iostream>
#include <fstream>
#include <pprint.hpp>
#include <params.h>

namespace gridpro {

int get_nblk(const std::string gridfile) {
  // ... open the connectivity file
  std::ifstream connfile(gridfile+".conn");
  std::string content;
  int nblk = 0;
  std::getline(connfile,content);
  nblk = std::stoi(content.substr(0,2));
  connfile.close();
  return nblk;
} // end function get_dims

mtr::FMatrix<int> get_dims(const std::string gridfile) {

}

void reader(const std::string gridfile,
            const int& nblk,
            mtr::FMatrix<double>& xc,
            mtr::FMatrix<double>& yc) {
  // fire up the printer
  pprint::PrettyPrinter printer;

  // ... open the grid file
  std::ifstream ifn(gridfile);
  int i = 0;
  int iblk,jblk;
  for (int n = 0; n < nblk; n++) {
    std::string content;
    std::getline(ifn,content);
    std::istringstream line(content);

    // get the bumber of cells for the current block
    if (content.length() < 10) {
      line >> iblk >> jblk;
      for (int i = 0; i < iblk; i++) {
        for (int j = 0; j < jblk; j++) {
          std::getline(ifn,content);
          std::istringstream line(content);
          line >> xc(i,j) >> yc(i,j);
          printer.print(xc(i,j),yc(i,j));
        }
      }
    }
  }
  ifn.close();
}

} // end namespace gridpro
