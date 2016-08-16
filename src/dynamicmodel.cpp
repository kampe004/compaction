
#include <iostream>
#include <memory>
#include <ctime>
#include <cmath>

#include "dynamicmodel.h"
#include "history.h"
#include "modelstate.h"
#include "config.h"
#include "constants.h"
#include "logging.h"
#include "meteo.h"
#include "surfacedensity.h"
#include "metamorphism.h"
#include "compaction.h"
#include "heatsolver.h"
#include "utils.h"

namespace DSM{ 

DynamicModel::DynamicModel() : _nt(0) { }

void DynamicModel::run(){
   /* MAIN DRIVER ROUTINE
      initializes submodels, runs timeloop and writes output */
   const char * option_name;

   /* read time step */
   option_name = "general:dt";
   _dt = (long) config.getInt(option_name, true, 1, (int)1e9, -1);
   const long dt_per_year = sec_in_year / _dt;
   logger << "INFO: dt_per_year = " << dt_per_year << std::endl;

   /* history stuff */
   option_name = "general:history_offset";
   const int hist_offset = (long) config.getInt(option_name, false, 0, (int)1e9, 1);

   /* read stopping criteria */
   option_name = "stopping_criteria:which_stop";
   const int which_stop = config.getInt(option_name, false, 0, 5, 0);
   bool stop_year = false;
   bool stop_dens = false;
   bool stop_depth = false;
   switch (which_stop) {
      case 0   :  stop_year = true; break;
      case 1   :  stop_dens = true; break;
      case 2   :  stop_depth = true; break;
      case 3   :  stop_year = true; stop_dens = true; break;
      case 4   :  stop_year = true; stop_depth = true; break;
      case 5   :  stop_year = true; stop_dens = true; stop_depth = true; break;
      default  : 
         logger << "ERROR: unknown value: " << which_stop << " for config option " << option_name << std::endl;
         std::abort();
   }
   int max_years;
   double max_dens, max_depth;
   if (stop_year) {
      option_name = "stopping_criteria:years";
      max_years = config.getInt(option_name, true, 1, (int)1e9, -1);
      logger << "INFO: max_years = " << max_years << std::endl;
   }
   if (stop_dens) {
      option_name = "stopping_criteria:dens";
      max_dens = config.getDouble(option_name, true, 1., 1000., -1.);
      logger << "INFO: max_dens = " << max_dens << std::endl;
   }
   if (stop_depth) {
      option_name = "stopping_criteria:depth";
      max_depth = config.getDouble(option_name, true, 1., 1000., -1.);
      logger << "INFO: max_depth = " << max_depth << std::endl;
   }

   /* setup submodels */
   std::unique_ptr<Meteo> meteo           = instantiate_meteo(*this);
   std::unique_ptr<SurfaceDensity> surf   = instantiate_surfacedensity(*meteo);
   ModelState mstate(*meteo, *surf, *this);
   std::unique_ptr<Metamorphism> mm       = instantiate_metamorphism(mstate, *this);
   std::unique_ptr<Compaction> comp       = instantiate_compaction(mstate, *this);
   std::unique_ptr<HeatSolver> heat       = instantiate_heatsolver(mstate, *this);
   History history(mstate);
   
   /* print some meteo statistics */
   logger << "INFO: initial surface temperature is: " << meteo->surfaceTemperature() << std::endl;
   logger << "INFO: initial surface wind is: " << meteo->surfaceWind() << std::endl;
   logger << "INFO: initial accumulation rate is: " << meteo->accumulationRate() << std::endl;

   /* start of the main time loop */
   std::clock_t start;
   int kYear = 0;
   while( 
         (stop_year ? kYear < max_years : true)
      && (stop_dens ? !mstate.hasReachedDensity(max_dens) : true)
      && (stop_depth ? mstate.totalDepth() < max_depth : true) ){

      start = clock();
      for(int tstep = 0; tstep < dt_per_year; tstep++) {
         runTimeStep(mstate, *mm, *comp, *heat);
         if (kYear >= hist_offset)
            history.update();
      }
      double elapsed = ((double)(clock() - start)) / CLOCKS_PER_SEC;

      logger  << "year=" << kYear
            << ", Tmin=" << mstate.minTemp() 
            << ", Tmax=" << mstate.maxTemp() 
            << ", rho_max=" << mstate.maxDens()
            << ", sec/year=" << elapsed 
            << ", year/hour=" << 3600./elapsed << std::endl;
      kYear++; 
   }
   mstate.printSummary();
   mstate.writeModelState();
   history.writeHistory();
}

void DynamicModel::runTimeStep(ModelState& mstate, Metamorphism& mm, Compaction& comp, HeatSolver& heat) {
   accumulate(mstate);
   mm.metamorphism();
   comp.compaction();
   doGridChecks(mstate);
   heat.heatdiffusion();
   _nt++;
}

void DynamicModel::doGridChecks(ModelState& mstate){
   /* PERFORMS A NUMBER OF CHECKS ON THE CONSISTENCY AND EFFICIENCY OF THE GRID */

   Grid& grid = mstate.getGrid();
   int i;

   /* check minimum thickness (top layer, index N-1) */
   i = grid.size()-1;
   if (grid[i].dz < dzmin_top && mstate.gridsize() > 1){
      // Layer is too thin, combine with neighbor
      mstate.combineLayers(grid[i-1], grid[i]);
      grid.pop_back();
   }

   /* check minimum thickness, start with second layer (index N-2) */
   i = grid.size()-2;
   while (i >= 0 ) { // Use while instead of for-loop as size may change during loop
      if (grid[i].dz < dzmin && mstate.gridsize() > 1){
         // Layer is too thin, combine with neighbor
         if (i == 0  || grid[i-1].dz > grid[i+1].dz  ) {
            /* CASE 1: Bottom layer can only be combined with the one above */
            /* CASE 2: Determine smallest neighbour */
            mstate.combineLayers(grid[i+1], grid[i]);
         } else {
            mstate.combineLayers(grid[i-1], grid[i]);
         }
         grid.erase(grid.begin() + i); // NOTE: this is very inefficient
      }
      i--;
   }


   /* check maximum layer thickness */
   Layer& top = grid.back();
   while (grid.back().dz > dzmax ) {
      /* Split in unequal parts 4:1 
         all vars are identical except thickness */
      top = grid.back();
      Layer layer = top;
      layer.dz    = top.dz/5.;
      top.dz      = top.dz * 4./5; 
      grid.push_back(layer);
   }

   /* remove bottom layers that have attained ice density */
   while(!grid.empty() && grid.front().dens > 900.) 
      grid.erase(grid.begin());

   if (grid.empty()) {
      logger << "ERROR: all layers turned to ice and were removed! This probably means that something is wrong with the compaction routine or the accumulation flux" << std::endl;
      std::abort();
   }
}

void DynamicModel::accumulate(ModelState& mstate) {
   const double acc = mstate.getMeteo().accumulationRate();
   if ( doubles_equal(acc, 0.0)) {return;}
   if ( acc < 0.) {
      accumulateNegative(mstate);
   } else {
      accumulatePositive(mstate);
   }
}

void DynamicModel::accumulateNegative(ModelState& mstate) {
   Grid& grid = mstate.getGrid();
   const double evap = std::abs(mstate.getMeteo().accumulationRate() * 1e-3 * _dt);
   while (grid.back().dz * grid.back().dens / rho_w < evap) {
      if (grid.size() < 2) {
         logger << "WARNING: not enough mass present for sublimation / evaporation" << std::endl;
         return;
      }
      logger << "grid size = " << grid.size() << std::endl;
      logger << "DEBUG combining [ " << (grid.end()[-2]).dz * (grid.end()[-2]).dens / rho_w << ", "
             << grid.back().dz * grid.back().dens / rho_w << " ] mm W.E. to meet " << evap << " mm W.E. of subl." << std::endl;
      mstate.combineLayers(grid.end()[-2], grid.back());
      grid.pop_back();
   }
   grid.back().dz -= evap * (rho_w/grid.back().dens);
//   logger << "DEBUG removing " << evap * (rho_w/grid.back().dens) << " m from top layer" << std::endl;
   if (grid.back().dz <= 0.) {
      logger << "ERROR: negative layer thickness after evaporation (programmer error)" << std::endl;
      std::abort();
   }
   return;
}

void DynamicModel::accumulatePositive(ModelState& mstate) {
   Grid& grid = mstate.getGrid();
   const double dens = mstate.getSurf().density();
   const double acc = mstate.getMeteo().accumulationRate() * 1e-3 * _dt; // convert from [mm/s] to [m]
   const double dz = (rho_w/dens) * acc;
   Layer& top = grid.back();

   // grain parameters are mass weighted
   const double U = mstate.getMeteo().surfaceWind();
   const double dfall = std::min(std::max(1.29-0.17*U, 0.20), 1.);
   const double sfall = std::min(std::max(0.08*U + 0.38, 0.5), 0.9);
   const double new_snow = dz*dens;
   const double old_snow = top.dz*top.dz;

#ifdef DEBUG
   const double dold = top.dnd;
   const double sold = top.sph;
#endif
   
   top.dnd    = (top.dnd*old_snow + dfall*new_snow)/(new_snow+old_snow);
   top.sph    = (top.sph*old_snow + sfall*new_snow)/(new_snow+old_snow);
   top.gs   = (top.gs*old_snow + fs_gs*new_snow)/(new_snow+old_snow);

#ifdef DEBUG
   if (top.dnd > 1. || top.dnd < 0) {
      logger << "DEBUG top.dnd = " << top.dnd << ", U = " << U << std::endl;
      logger << "DEBUG dfall = " << dfall << ", sfall = " << sfall << std::endl;
      logger << "DEBUG dold = " << dold << ", sold = " << sold << std::endl;
      logger << "DEBUG new_snow = " << new_snow << ", old_snow = " << old_snow << std::endl;
      std::abort();
   }
#endif

   top.dens = top.dz/(top.dz+dz)*top.dens + dz/(top.dz+dz)*dens;
   top.dz = top.dz + dz;

}

} // namespace
