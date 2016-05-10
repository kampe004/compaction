#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstring>

#include "logging.h"
#include "dynamicmodel.h"
#include "modelstate.h"
#include "meteo.h"
#include "constants.h"
#include "surfacedensity.h"

namespace DSM{ 

ModelState::ModelState(Meteo& meteo, SurfaceDensity& surf, DynamicModel& dm) : _meteo(meteo), _dm(dm), _surf(surf) {
   // initialize grid with a very tiny top layer
   Layer layer;
   layer.T      = _meteo.surfaceTemperature();
   layer.dens   = 200.;
   layer.dz     = 0.00001;
   layer.d      = fs_dend;
   layer.s      = fs_sphere;
   layer.gs     = fs_gs;
   _grid.push_back(layer);
}

double ModelState::getDepthOfDensity(const double dens) {
   double diffi, diffip1;
   double dzi, dzip1;
   int tidx = gridsize()-1; // index of top layer
   if (_grid[tidx].dens >= dens) return _grid[tidx].dz/2; // node depth of top layer

   double tot_depth = _grid[tidx].dz; // interface depth
   for (int i = tidx-1; i >= 0; i--) {
      if (_grid[i].dens >= dens ){
         diffi = std::abs(_grid[i].dens - dens); // delta 1
         diffip1 = std::abs(_grid[i+1].dens - dens); // delta 2
         
         dzi = tot_depth + _grid[i].dz / 2; // node depth
         dzip1 = tot_depth - _grid[i+1].dz / 2; // node depth

         return diffi/(diffi+diffip1) * dzip1 + diffip1/(diffi+diffip1) * dzi; // linear interpolation
      }
      tot_depth = tot_depth + _grid[i].dz;
   }
   return -1.0; // dens not found
}

bool ModelState::hasReachedDensity(const double dens){
   return this->getDepthOfDensity(dens) > 0.; 
}

double ModelState::getZ550() {
   return getDepthOfDensity(550.);
}

double ModelState::getZ830() {
   return getDepthOfDensity(830.); }

double ModelState::totalDepth() {
   double tot_depth = 0;
   for (int i = gridsize()-1; i >= 0; i--) {
      tot_depth = tot_depth + _grid[i].dz;
   };
   return tot_depth;
}

void ModelState::printSummary() {
   double tot_mass = 0;
   double tot_depth = 0;
   for (int i = gridsize()-1; i >= 0; i--) {
      tot_mass += _grid[i].dz * _grid[i].dens;
      tot_depth = tot_depth + _grid[i].dz;
   }
   logger << "Total mass = " << tot_mass << std::endl;
   logger << "Total depth = " << tot_depth << std::endl;

   double z550 = getZ550();
   if (z550 > 0.){
      logger << "z550 depth = " << z550 << std::endl;
   }
   double z830 = getZ830();
   if (z830 > 0.){
      logger << "z830 depth = " << z830 << std::endl;
   }
}

void ModelState::writeModelState() {
   const int nFields = 8;
   static const int PREC = 4; // precision
   const char * names[nFields];
   int k = 0;
   names[k++] = "node depth[m] ";
   names[k++] = "dz[m] ";
   names[k++] = "mass[kg] ";
   names[k++] = "rho[kg/m3] ";
   names[k++] = "T[K] ";
   names[k++] = "grain_gs[m] ";
   names[k++] = "grain_d ";
   names[k++] = "grain_s ";
   int fldlen[nFields];
   for (int i = 0; i<nFields; i++) {
      fldlen[i] = std::max((int)strlen(names[i]), PREC+3); 
   }

   std::ofstream fout;
   fout.open("layers.txt");

   /* write header */
   fout << gridsize() << std::endl;
   for (int i = 0; i<nFields; i++) {
      fout << std::setw(fldlen[i]) << std::left << names[i];
   }
   fout << std::endl;

   /* write values */
   double node_depth = 0;
   fout << std::setprecision(PREC);
   for (int i = gridsize()-1; i >= 0; i--) {
      node_depth += _grid[i].dz * 0.5;
      k = 0;
      fout << std::setw(fldlen[k++]) << std::left << node_depth;
      fout << std::setw(fldlen[k++]) << std::left << _grid[i].dz;
      fout << std::setw(fldlen[k++]) << std::left << _grid[i].dz * _grid[i].dens;
      fout << std::setw(fldlen[k++]) << std::left << _grid[i].dens;
      fout << std::setw(fldlen[k++]) << std::left << _grid[i].T;
      fout << std::setw(fldlen[k++]) << std::left << _grid[i].gs;
      fout << std::setw(fldlen[k++]) << std::left << _grid[i].d;
      fout << std::setw(fldlen[k++]) << std::left << _grid[i].s;
      fout << std::endl;
      node_depth += _grid[i].dz * 0.5;
   }
   fout.close();

   std::ofstream f_z550;
   f_z550.open("z550.txt");
   f_z550 << getZ550() << std::endl ;
   f_z550.close();

   std::ofstream f_z830;
   f_z830.open("z830.txt");
   f_z830 << getZ830() << std::endl ;
   f_z830.close();
}

double ModelState::maxTemp() {
   double retval = _grid[0].T;
   for (int i = 1; i < gridsize(); i++){
     retval = std::max(retval,_grid[i].T);
   }
   return retval;
}

double ModelState::minTemp() {
   double retval = _grid[0].T;
   for (int i = 1; i < gridsize(); i++){
     retval = std::min(retval,_grid[i].T);
   }
   return retval;
}

double ModelState::maxDens() {
   double retval = _grid[0].dens;
   for (int i = 1; i < gridsize(); i++){
     retval = std::max(retval,_grid[i].dens);
   }
   return retval;
}

void ModelState::combineLayers(Layer& lay1, Layer& lay2) {
   /* the resulting layer is stored in 'a' */
   const double m1 = lay1.dz * lay1.dens;
   const double m2 = lay2.dz * lay2.dens;
   lay1.d   = (lay1.d * m1 + lay2.d * m2)/(m1+m2);
   lay1.s   = (lay1.s * m1 + lay2.s * m2)/(m1+m2);
   lay1.gs  = (lay1.gs* m1 + lay2.gs* m2)/(m1+m2);
   lay1.dz  = lay1.dz + lay2.dz;
   lay1.dens = (m1+m2)/lay1.dz;
}

} // namespace
