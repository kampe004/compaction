
#include <iostream>
#include <memory>
#include <ctime>

#include "dynamicmodel.h"
#include "modelstate.h"
#include "config.h"
#include "constants.h"
#include "logging.h"
#include "meteo.h"
#include "surfacedensity.h"
#include "metamorphism.h"
#include "compaction.h"

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

   /* read heat diffusion option */
   option_name = "physics:heat_diffusion";
   const int heat_diffusion = config.getInt(option_name, false, 0, 1, 1);
   switch (heat_diffusion) {
      case 0   :  _has_heat = false; break;
      case 1   :  _has_heat = true; break;
   }

   /* setup submodels */
   std::unique_ptr<Meteo> meteo = instantiate_meteo(*this);
   std::unique_ptr<SurfaceDensity> surf = instantiate_surfacedensity(*meteo);
   ModelState mstate(*meteo, *surf, *this);
   std::unique_ptr<Metamorphism> mm = instantiate_metamorphism(mstate, *this);
   std::unique_ptr<Compaction> comp = instantiate_compaction(mstate, *this);
   
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
         runTimeStep(mstate, *mm, *comp);
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
}

void DynamicModel::runTimeStep(ModelState& mstate, Metamorphism& mm, Compaction& comp) {
   accumulate(mstate);
   mm.metamorphism();
   comp.compaction();
   doGridChecks(mstate);

   if (_has_heat && mstate.gridsize() > 1)
      heatDiffusionShallow(mstate);

   _nt++;
}

void DynamicModel::doGridChecks(ModelState& mstate){
   /* PERFORMS A NUMBER OF CHECKS ON THE CONSISTENCY AND EFFICIENCY OF THE GRID */

   /* check minimum thickness, start with second layer (index N-2) */
   Grid& grid = mstate.getGrid();
   int i = grid.size()-2;
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

   /* remove bottom layers that have attained ice density */
   while(!grid.empty() && grid.front().dens > 900.) 
      grid.erase(grid.begin());

   if (grid.empty()) {
      logger << "ERROR: all layers turned to ice and were removed! This probably means that something is wrong with the compaction routine or the accumulation flux" << std::endl;
      std::abort();
   }
   //if (grid.front().dens > 900.) {
   //   logger << "WARNING: ice density reached! " << std::endl;
   //   for (int i = 0; i < mstate.gridsize(); i++) 
   //      logger << "i = " << i << " dens = " << grid[i].dens << std::endl;
   //   }
   }
}

void DynamicModel::accumulate(ModelState& mstate) {
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
   top.d    = (top.d*old_snow + dfall*new_snow)/(new_snow+old_snow);
   top.s    = (top.s*old_snow + sfall*new_snow)/(new_snow+old_snow);
   top.gs   = (top.gs*old_snow + fs_gs*new_snow)/(new_snow+old_snow);

   top.dens = top.dz/(top.dz+dz)*top.dens + dz/(top.dz+dz)*dens;
   top.dz = top.dz + dz;

   if (top.dz > dzmax ) {
      // Split in unequal parts 4:1
      Layer layer = top;
      layer.dz    = top.dz/5.;
      top.dz      = top.dz * 4./5; 
      grid.push_back(layer);
   }
}

void DynamicModel::heatDiffusion(ModelState& mstate) {
   /* Uses forward Euler method for diffusion equation*/ 
   Grid& grid = mstate.getGrid();
   int Np = mstate.gridsize();
   double tc;
   double td[Np];
   double T_new[Np];

   for (int i = 0; i < mstate.gridsize(); i++) {
      tc = tcair + (7.75e-5*grid[i].dens + 1.105e-6*pow(grid[i].dens,2))*(tcice - tcair); // thermal conductivity [W/m/K]
      td[i] = tc/(grid[i].dens * cp); // thermal diffusivity [m2/s]

      double stab_crit = td[i] * _dt / pow(grid[i].dz,2);
      if (stab_crit > 0.5) {
          logger << "WARNING: stability criterium r < 0.5 violated, r = " << stab_crit << std::endl;
      }
   }

   /* TOP: Dirichlet b.c. */
   T_new[Np-1] = mstate.getMeteo().surfaceTemperature(); // top boundary equal to surface temperature

   /* BOTTOM: Neumann b.c. */
   T_new[0] = td[0]*_dt/grid[0].dz * (grid[1].T-grid[0].T) + grid[0].T;

   for (int i = 1; i < Np-1; i++) {
      T_new[i] = td[i]*_dt/pow(grid[i].dz,2)*(grid[i-1].T - 2 * grid[i].T + grid[i+1].T) + grid[i].T; // # Forward in Time Centered in Space
   }

   /* Update temperature */
   for (int i = 0; i < Np; i++) {
      grid[i].T = T_new[i];
   }
}

void DynamicModel::heatDiffusionShallow(ModelState& mstate) {
   /* Uses forward Euler method for diffusion equation
      only solves for the shallow part (top 10 metres) 
      assume constant temperature below */ 
   static const double solve_lim = 10.0;

   Grid& grid = mstate.getGrid();
   int Np = mstate.gridsize();
   double tc;
   double td;
   double T_new[Np];

   /* TOP: Dirichlet b.c. */
   T_new[Np-1] = mstate.getMeteo().surfaceTemperature(); // top boundary equal to surface temperature
   double total_depth = grid[Np-1].dz;

   int i;
   for (i = Np-2; i > 0; i--) {
      tc = tcair + (7.75e-5*grid[i].dens + 1.105e-6*pow(grid[i].dens,2))*(tcice - tcair); // thermal conductivity [W/m/K]
      td = tc/(grid[i].dens * cp); // thermal diffusivity [m2/s]

      double stab_crit = td * _dt / pow(grid[i].dz,2);
      if (stab_crit > 0.5) {
          logger << "WARNING: stability criterium r < 0.5 violated, r = " << stab_crit << std::endl;
      }
      T_new[i] = td*_dt/pow(grid[i].dz,2)*(grid[i-1].T - 2 * grid[i].T + grid[i+1].T) + grid[i].T; // # Forward in Time Centered in Space

      total_depth += grid[i].dz;
      if (total_depth > solve_lim) {
         break;
      }
   }

   if (i == 0) {
      /* BOTTOM: Neumann b.c. */
      T_new[0] = td*_dt/grid[0].dz * (grid[1].T-grid[0].T) + grid[0].T;
   }

   /* Update temperature */
   for (int k = i; k < Np; k++) {
      grid[k].T = T_new[k];
   }
}

} // namespace
