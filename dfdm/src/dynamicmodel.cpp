
#include <iostream>
#include <memory>
#include <ctime>

#include "dynamicmodel.h"
#include "config.h"
#include "constants.h"
#include "logging.h"
#include "meteo.h"

namespace DSM{ 

DynamicModel::DynamicModel(){
   logger << "DynamicModel()" << std::endl;
   _nt = 0;
   const char * option_name = "general:dt";
   _dt = (long) config.getInt(option_name, true, 1, (int)1e9, 0);
}

void DynamicModel::run(){
   //Snowpack sp; 
   std::unique_ptr<Meteo> meteo = instantiate_meteo(*this);
   //std::unique_ptr<FsdParam> fsd = instantiate_fsd(_dict)
   //std::unique_ptr<CompactionParam> compact;
   //std::unique_ptr<MetamorphismParam> metamorph;
   
   meteo->surfaceTemperature();
   long dt_per_year = sec_in_year / _dt;
   logger << "INFO: dt_per_year = " << dt_per_year << std::endl;

//    std::clock_t start;
//    int kYear = 0;
//    while(! core->hasReachedDensity(850.) && kYear < settings.max_year && core->totalDepth() < settings.max_depth) {
// //   for (int year = 0; year < nyears; year++){
//       start = clock();
//       for(int tstep = 0; tstep < dt_per_year; tstep++)
//          core->runTimeStep(dt);
//       double elapsed = ((double)(clock() - start)) / CLOCKS_PER_SEC;
//       logger  << "year=" << kYear
//             << ", Tmin=" << core->minTemp() 
//             << ", Tmax=" << core->maxTemp() 
//             << ", rho_max=" << core->maxDens()
//             << ", sec/year=" << elapsed 
//             << ", year/hour=" << 3600./elapsed << std::endl;
//       kYear++; 
//    }
//    //core->printIceCoreSummary();
//    //core->writeFiles();
}


// void ModelState::runTimeStep(long dt) {
//    if (!is_initialized){
//       logger << "ERROR: core is not initialized, run init() first!" << std::endl;
//       std::abort();
//    }
//    accumulate(dt);
//    compact(dt);
// 
//    // check minimum thickness, start with second layer (index N-2)
//    int i = grid.size()-2;
//    double combined_mass;
//    while (i >= 0 ) { // Use while instead of for-loop as size may change during loop
//       if (grid[i].dz < dzmin && gridsize() > 1){
//          //logger << "DEBUG: merging" << std::endl;
//          // Layer is too thin, combine with neighbor
//          if (i == 0  || grid[i-1].dz > grid[i+1].dz  ) {
//             /* CASE 1: Bottom layer can only be combined with the one above */
//             /* CASE 2: Determine smallest neighbour */
//             combined_mass = grid[i+1].dens * grid[i+1].dz + grid[i].dens * grid[i].dz;
//             grid[i+1].dz   = grid[i+1].dz + grid[i].dz;
//             grid[i+1].dens  = combined_mass / grid[i+1].dz;
//          } else {
//             combined_mass = grid[i-1].dens * grid[i-1].dz + grid[i].dens * grid[i].dz;
//             grid[i-1].dz   = grid[i-1].dz + grid[i].dz;
//             grid[i-1].dens  = combined_mass / grid[i-1].dz;
//          }
//          grid.erase(grid.begin() + i); // NOTE: this is very inefficient
//       }
//       i--;
//    }
// 
//    // remove bottom layers that have attained ice density
//    //while(grid.front().dens > 900.) 
//    //   grid.erase(grid.begin());
//    if (grid.front().dens > 900.) {
//    for (int i = 0; i < gridsize(); i++) {
//       logger << "i = " << i << " dens = " << grid[i].dens << "\n";
//       } 
//    }
// 
//    if (settings.have_diffusion && gridsize() > 1)
//       heatDiffusion(dt);
// 
//    // update internal time variable
//    current_time = current_time + dt;
// }
// 
// void ModelState::accumulate(long dt) {
//    double dens = this->surfaceDensity(current_time);
//    double acc = this->accumulationRate(current_time) * 1e-3 * dt; // convert from [mm/s] to [m]
//    double dz = (rho_w/dens) * acc;
// 
//    Layer& top = grid.back();  // top layer
//    top.dens = top.dz/(top.dz+dz)*top.dens + dz/(top.dz+dz)*dens;
//    top.dz = top.dz + dz;
// 
//    if (top.dz > dzmax ) {
//       // Split in unequal parts 4:1
//       Layer layer;
//       layer.dens  = top.dens;
//       layer.dz   = top.dz/5;
//       layer.T    = top.T;
// 
//       top.dz = top.dz * 4/5; 
// 
//       grid.push_back(layer);
//    }
// }
// 
// void ModelState::compact(long dt) {
//    // enum class DensificationMethod {Ligtenberg2011, Anderson1976, Barnola1991, Spencer2001, BarnolaSpencer}
//    switch (settings.dm) {
//       case DensificationMethod::Ligtenberg2011 : compactLigtenberg2011(dt); break;
//       case DensificationMethod::Anderson1976   : compactAnderson1976(dt); break;
//       case DensificationMethod::Barnola1991   : compactBarnola1991(dt); break;
//       case DensificationMethod::Spencer2001   : compactSpencer2001(dt); break;
// //      case DensificationMethod::BarnolaSpencer : compactBarnolaSpencer(dt); break;
//       case DensificationMethod::HerronLangway  : compactHerronLangway(dt); break;
//       default:
//          logger << "ERROR: programmer error: unknown densification method" << std::endl;
//          std::abort();
//    }
// }
// 

} // namespace
