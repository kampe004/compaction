#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

#include "logging.h"
#include "dynamicmodel.h"
#include "modelstate.h"
#include "meteo.h"
#include "surfacedensity.h"

namespace DSM{ 

ModelState::ModelState(Meteo& meteo, SurfaceDensity& surf, DynamicModel& dm) : _meteo(meteo), _dm(dm), _surf(surf) {
   // initialize grid with a very tiny top layer
   Layer layer;
   layer.T      = _meteo.surfaceTemperature();
   //layer.dens   = _dm.meteo->surfaceDensity();
   layer.dz     = 0.00001;
   _grid.push_back(layer);
}

double ModelState::getDepthOfDensity(double dens) {
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

bool ModelState::hasReachedDensity(double dens){
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
   std::ofstream f_layers;
   f_layers.open("layers.csv");
   f_layers << gridsize() << std::endl;
   for (int i = gridsize()-1; i >= 0; i--) {
      f_layers << _grid[i].dz << "," << _grid[i].dz * _grid[i].dens<< "," << _grid[i].dens << "," << _grid[i].T << std::endl;
   }
   f_layers.close();

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

} // namespace
