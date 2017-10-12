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
   layer.dz     = 0.1;
   layer.dnd    = fs_dnd;
   layer.sph    = fs_sph;
   layer.gs     = fs_gs;
   _grid.push_back(layer);

   // ice grid, 100 x 10 cm layers = 10 m
   layer.T      = 250.;
   layer.dens   = rho_i;
   layer.dz     = 0.1;
   double tot_dz = 0.0;
   for (int i = 0; i< 12; i++){
      _icegrid.push_back(Layer(layer)); // use copy constructor to create unique layers
      tot_dz += layer.dz;
      layer.dz *= 1.5;
   }
   std::cout << "created " << _icegrid.size() << " ice layers, combined depth = " << tot_dz << " m" << std::endl;
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
      logger << "final z550 depth = " << z550 << std::endl;
   }
   double z830 = getZ830();
   if (z830 > 0.){
      logger << "final z830 depth = " << z830 << std::endl;
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
      fldlen[i] = std::max((int)strlen(names[i])+1, PREC+5); 
   }

   std::ofstream fout;
   fout.open("layers.txt");

   /* write header */
   fout << gridsize() << std::endl;
   for (int i = 0; i<nFields; i++) {
      fout << std::setw(fldlen[i]) << std::left << names[i] << " ";
   }
   fout << std::endl;

   /* write values */
   double node_depth = 0;
   fout << std::setprecision(PREC);
   for (int i = gridsize()-1; i >= 0; i--) {
      node_depth += _grid[i].dz * 0.5;
      k = 0;
      fout << std::setw(fldlen[k++]) << std::left << node_depth << " ";
      fout << std::setw(fldlen[k++]) << std::left << _grid[i].dz << " ";
      fout << std::setw(fldlen[k++]) << std::left << _grid[i].dz * _grid[i].dens << " ";
      fout << std::setw(fldlen[k++]) << std::left << _grid[i].dens << " ";
      fout << std::setw(fldlen[k++]) << std::left << _grid[i].T << " ";
      fout << std::setw(fldlen[k++]) << std::left << _grid[i].gs << " ";
      fout << std::setw(fldlen[k++]) << std::left << _grid[i].dnd << " ";
      fout << std::setw(fldlen[k++]) << std::left << _grid[i].sph << " ";
      fout << std::endl;
      node_depth += _grid[i].dz * 0.5;
   }
   fout.close();

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
   /* the merged layer is stored in 'lay1' 
      Caller needs to delete lay2 afterwards!! */
   const double m1 = lay1.dz * lay1.dens;
   const double m2 = lay2.dz * lay2.dens;
   const double heat1 = m1 * cpice * lay1.T;
   const double heat2 = m2 * cpice * lay2.T;
   lay1.dnd = (lay1.dnd * m1 + lay2.dnd * m2)/(m1+m2);
   lay1.sph = (lay1.sph * m1 + lay2.sph * m2)/(m1+m2);
   lay1.gs  = (lay1.gs* m1 + lay2.gs* m2)/(m1+m2);
   lay1.dz  = lay1.dz + lay2.dz;
   lay1.dens = (m1+m2)/lay1.dz;
   lay1.T   = (heat1+heat2)/(cpice * lay1.dz * lay1.dens); // Heat is conserved
}

} // namespace
