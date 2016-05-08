
#include "constants.h"
#include "config.h"
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
   
   for (int i = _mstate.gridsize()-1; i >= 0; i--) {
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
   return;
}

} // namespace
