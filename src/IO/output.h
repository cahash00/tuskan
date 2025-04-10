#ifndef OUTPUT_H
#define OUTPUT_H

#include <matar.h>
#include <iostream>
#include <filesystem>
#include <string>
#include <logger.h>
#include <params.h>

using namespace std;
namespace fs = std::filesystem;

namespace IO {
/******************************************************************************/
void check_directories(const string& foutDir);

/******************************************************************************/
// save an array to a restart file
void save_restart(const std::string& filename, 
    const mtr::FMatrix<double>& u);

/******************************************************************************/
// load an array from a restart file
void load_restart(const std::string& filename,
    const mtr::FMatrix<double>& u);

/******************************************************************************/
inline void add_array(FILE* out,const bool ghost) {};

template <typename T, typename... Args>
void add_array(FILE* out,const bool ghost,const string& label,T&& array, Args&&... args) {
  double xstr,ystr,xend,yend;
  if (ghost) {
    xstr = istr-1;
    xend = iend+1;
    ystr = jstr-1;
    yend = jend+1;
  } else {
    xstr = istr;
    xend = iend;
    ystr = jstr;
    yend = jend;
  }

  fprintf(out,"SCALARS %s float 1\n",label.c_str());
  fprintf(out,"LOOKUP_TABLE default\n");
  for (int j = ystr; j <= yend; j++) {
    for (int i = xstr; i <= xend; i++) {
      fprintf(out,"%f\n",array(i,j));
    }
  }

  add_array(out,ghost,forward<Args>(args)...);
};

template <typename... Args>
  void vtk_output_2D(const string& caseName,
      const string& foutDir,
      const bool ghost,
      const mtr::FMatrix<double>& xc,
      const mtr::FMatrix<double>& yc,
      const mtr::FMatrix<double>& u,
      const mtr::FMatrix<double>& v,
      Args&&... args) {
    /**
     * get cell-centered solution variables
     * pressure is already cell-centered
     */
    int xstr,xend,ystr,yend,dnx,dny;

    if (ghost) {
      xstr = istr-1;
      xend = iend+1;
      ystr = jstr-1;
      yend = jend+1;
      dnx = nx+2;
      dny = ny+2;
    } else {
      xstr = istr;
      xend = iend;
      ystr = jstr;
      yend = jend;
      dnx = nx;
      dny = ny;
    }


    // ... file output
    FILE *out;
    std::ostringstream foutss;
    foutss << foutDir << "/flow." << caseName << ".vtk";
    string fout = foutss.str();
    out = fopen(fout.c_str(),"w");
    fprintf(out,"# vtk DataFile Version 3.0\n");
    fprintf(out,"TUSKAN 2D Flow Solution File.\n");
    fprintf(out,"ASCII\n");
    fprintf(out,"DATASET STRUCTURED_GRID\n");
    fprintf(out,"DIMENSIONS %d %d %d\n",dnx,dny,1);
    fprintf(out,"POINTS %d float\n",(dnx)*(dny));
    for (int j = ystr; j <= yend; j++) {
      for (int i = xstr; i <= xend; i++) {
        fprintf(out,"%f %f %f\n",xc(i,j),yc(i,j),0.0);
      }
    }
    fprintf(out,"POINT_DATA %d\n",(dnx)*(dny));
    fprintf(out,"VECTORS velocity float \n");
    for (int j = ystr; j <= yend; j++) {
      for (int i = xstr; i <= xend; i++) {
        fprintf(out,"%f %f %f\n",u(i,j),v(i,j),0.0);
      }
    }
    add_array(out,ghost,forward<Args>(args)...);
    fclose(out);
};

void getCellCenter(mtr::FMatrix<double>& u,
                   mtr::FMatrix<double>& v,
                   mtr::FMatrix<double>& uc,
                   mtr::FMatrix<double>& vc);

/******************************************************************************/
} // end namespace IO
#endif // OUTPUT_H
