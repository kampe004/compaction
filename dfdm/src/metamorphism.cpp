
#include "constants.h"
#include "config.h"
#include "utils.h"
#include "logging.h"
#include "metamorphism.h"
#include "modelstate.h"
#include "dynamicmodel.h"

namespace DSM{ 

std::unique_ptr<Metamorphism> instantiate_metamorphism(ModelState& mstate, DynamicModel& dm){
   const char * option_name = "physics:which_metamorphism";
   int which_metamorphism = config.getInt(option_name, false, 0, 2, 1);

   switch (which_metamorphism) {
      case 0   : return { std::make_unique<MetamorphismNone>(mstate, dm) };
      case 1   : return { std::make_unique<MetamorphismAnderson1976>(mstate, dm) };
      case 2   : return { std::make_unique<MetamorphismCROCUS>(mstate, dm) };
      default:
         logger << "ERROR: unknown value: " << which_metamorphism << " for config option " << option_name << std::endl;
         std::abort();
   }
}

Metamorphism::Metamorphism(ModelState& mstate, DynamicModel& dm) : _mstate(mstate), _dm(dm) { } 

MetamorphismNone::MetamorphismNone(ModelState& mstate, DynamicModel& dm) : Metamorphism(mstate, dm) { 
   logger << "MetamorphismNone()" << std::endl; 
} 

MetamorphismAnderson1976::MetamorphismAnderson1976(ModelState& mstate, DynamicModel& dm) : Metamorphism(mstate, dm) { 
   logger << "MetamorphismAnderson1976()" << std::endl; 
} 

MetamorphismCROCUS::MetamorphismCROCUS(ModelState& mstate, DynamicModel& dm) : Metamorphism(mstate, dm) { 
   logger << "MetamorphismNoneCROCUS" << std::endl; 
} 

void MetamorphismNone::metamorphism() { 
   return; 
}

void MetamorphismAnderson1976::metamorphism() {
   /* Compaction due to destructive metamorphism (Anderson 1976) 
      implementation based on CLM 4.5 Tech Note */
   Grid& grid = _mstate.getGrid();
   static const double Tf = T0;
   static const double c3 = 2.777e-6;
   static const double c4 = 0.04;

   double layer_mass;
   
   for (int i = grid.size()-1; i >= 0; i--) {
      layer_mass = grid[i].dz * grid[i].dens;
      double c1 = grid[i].dens < 100. ? 1.0 : exp(-0.046*(grid[i].dens-100.));
      double c2 = 1 ; // assuming no liquid! 
      double cr1 = -c3*c2*c1*exp(-c4*(Tf-grid[i].T)); // compaction due to destructive metamorphism
      //double cr = cr1 + cr2;
      grid[i].dz = grid[i].dz * (1+(cr1*_dm.getDt()));
      grid[i].dens = (layer_mass/grid[i].dz); // # mass conservation
   }
}

void MetamorphismCROCUS::metamorphism() {
   /* Empirical laws for dry snow metamorphism. From: Vionnet et al, 2012: table 1.
      TODO: Tgrad > 15 is ommitted here (no depth hoar formation)
      */
   const double deltat = (double)_dm.getDt() / (double) sec_in_day;

   Grid& grid = _mstate.getGrid();
   int Np = grid.size();
   if (Np > 1) return; // can't compute gradients on single snow layer

   double Tgrad;
   for (int i = 0; i < Np; i++) {
      if (i == 0) {
         // SNOW BOTTOM, one-sided up
         Tgrad = std::abs((grid[1].T-grid[0].T)/(0.5*grid[0].dz+0.5*grid[1].dz));
      } else if (i == Np-1) {
         // SNOW TOP, one-sided down
         Tgrad = std::abs((grid[Np-1].T-grid[Np-2].T)/(0.5*grid[Np-1].dz+0.5*grid[Np-2].dz));
      } else {
         // central difference
         Tgrad = std::abs(grid[i+1].T-grid[i-1].T)/(0.5*grid[i+1].dz+grid[i].dz+grid[i-1].dz);
      }

      if (doubles_equal(grid[i].d, 0.0)) {
         // Non-dendritic snow
         if (Tgrad <= 5.) {
            grid[i].s = grid[i].s + deltat * 1e9 * exp(-6000./grid[i].T);
         } else {
            grid[i].s = grid[i].s + deltat * -2.e8 * exp(-6000./grid[i].T) * pow(Tgrad,0.4);
         }
      } else {
         // Dendritic snow
         if (Tgrad <= 5.) {
            grid[i].d = grid[i].d + deltat * -2.e8 * exp(-6000./grid[i].T);
            grid[i].s = grid[i].s + deltat * 1e9 * exp(-6000./grid[i].T);
         } else {
            grid[i].d = grid[i].d + deltat * -2.e8 * exp(-6000./grid[i].T) * pow(Tgrad,0.4);
            grid[i].s = grid[i].s + deltat * -2.e8 * exp(-6000./grid[i].T) * pow(Tgrad,0.4);
         }
      }
   }
}

} // namespace
