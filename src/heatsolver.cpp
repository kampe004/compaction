
#include <math.h>
#include <vector>

#include "constants.h"
#include "config.h"
#include "utils.h"
#include "logging.h"
#include "heatsolver.h"
#include "modelstate.h"
#include "dynamicmodel.h"
#include "meteo.h"

namespace DSM{ 

std::unique_ptr<HeatSolver> instantiate_heatsolver(ModelState& mstate, DynamicModel& dm){
   const char * option_name = "physics:which_diffusion";
   const int which_diffusion = config.getInt(option_name, false, 0, 2, 2);

   switch (which_diffusion) {
      case 0   : return { std::make_unique<HeatSolverNone>(mstate, dm) };
      case 1   : return { std::make_unique<HeatSolverExplicit>(mstate, dm) };
      case 2   : return { std::make_unique<HeatSolverImplicit>(mstate, dm) };
      default:
         logger << "ERROR: unknown value: " << which_diffusion << " for config option " << option_name << std::endl;
         std::abort();
   }
}

HeatSolver::HeatSolver(ModelState& mstate, DynamicModel& dm) : _mstate(mstate), _dm(dm) { } 

HeatSolverNone::HeatSolverNone(ModelState& mstate, DynamicModel& dm) : HeatSolver(mstate, dm) { 
   logger << "HeatSolverNone()" << std::endl; 
} 

HeatSolverExplicit::HeatSolverExplicit(ModelState& mstate, DynamicModel& dm) : HeatSolver(mstate, dm) { 
   logger << "HeatSolverExplicit()" << std::endl; 
} 

HeatSolverImplicit::HeatSolverImplicit(ModelState& mstate, DynamicModel& dm) : HeatSolver(mstate, dm) { 
   logger << "HeatSolverImplicit()" << std::endl; 
} 

void HeatSolverNone::heatdiffusion() { 
   return; 
}

void HeatSolverExplicit::heatdiffusion() {
   /* Forward Euler timestepping 
      In our model, layers are counted bottom-up, which means that grid index 0 is the 
      bottom layer and layer Np-1 is the top layer. 

      Here we follow this approach for consistency, counting layer i from 0 to Np-1.
      Therefore, for a given layer i, the layer i+1 sits above it and i-1 below it.
      We want the difference zi(i+1)-zi(i) to be positive, so zi is counted as positive from the 
      bottom interface (ice or soil), which means it is actually height not depth.

      Interfaces are regarded as located *below* the node with the same index, i.e.
      interface 0 is the bottom interface and interface Np is the surface.
   */
   Grid& grid = _mstate.getGrid();
   const int Np = grid.size();
   const double dt = (double)_dm.getDt(); // time in seconds
   const double Tskin = (double)_mstate.getMeteo().surfaceTemperature();

   double zi[Np];       // node location [m]
   double zii[Np+1];    // interface location (interface sits below node of the same index)
   double tc[Np];       // thermal conductivity at nodes [W/m/K]
   double tci[Np+1];    // thermal conductivity at interfaces [W/m/K]
   double Hfi[Np+1];    // heat flux at interface [W/m2]
   double Tgrad;        // temperature gradient
   double cpv;          // volumetric heat capacity [J/m3/K]

   double tot_depth = 0.0;
   for (int i = 0; i<Np; i++) {
      zii[i]      = tot_depth; 
      tot_depth   = tot_depth + grid[i].dz;
      zi[i]       = tot_depth - 0.5*grid[i].dz;
      /* Thermal conductivity as in Jordan (1991) */
      tc[i] = tcair + (7.75e-5*grid[i].dens + 1.105e-6*pow(grid[i].dens,2))*(tcice - tcair); 
   }
   zii[Np] = tot_depth;

   Hfi[0] = 0.0; // bottom boundary condition
   tci[0] = -1.0;
   for (int i = 1; i<Np; i++) {
      /* interpolation of thermal conductivity to the interface, 
         as CLM45 Tech note eq 6.11 
         Note: interfaces are located *below* the node with the same index */
      tci[i]   = tc[i]*tc[i-1]*(zi[i]-zi[i-1]) / (tc[i]*(zii[i]-zi[i-1]) + tc[i-1]*(zi[i]-zii[i]));
      /* heat flux (positively downwards) */
      Tgrad    = (grid[i].T - grid[i-1].T) / (zi[i]-zi[i-1]);
      Hfi[i]    = tci[i] * Tgrad; 
   }
   Hfi[Np] = tc[Np-1] * (Tskin-grid[Np-1].T) / (zii[Np] - zi[Np-1]);

   for (int i = 0; i<Np; i++) {
      cpv         = grid[i].dens * cpice;
      grid[i].T   = grid[i].T + dt / (cpv*grid[i].dz) * (Hfi[i+1] - Hfi[i]);
   }
}

void HeatSolverImplicit::heatdiffusion() {
   /* Crank-Nicholson implicit timestepping 
      See notes for the explicit solver above. 
      Here a system of equations is built that is solved with LAPACK. 
      The system is tri-diagonal. 
   */
   Grid& grid = _mstate.getGrid();
   const int Np = grid.size();
   if (Np < 2) return;
   const double dt = (double)_dm.getDt(); // time in seconds
   const double Tskin = (double)_mstate.getMeteo().surfaceTemperature();
   static const double alpha = 0.5; // Crank-Nicholson factor

   double zi[Np];       // node location [m]
   double zii[Np+1];    // interface location (interface sits below node of the same index)
   double tc[Np];       // thermal conductivity at nodes [W/m/K]
   double tci[Np+1];    // thermal conductivity at interfaces [W/m/K]
   double Hfi[Np+1];    // heat flux at interface [W/m2]
   double Tgrad;        // temperature gradient
   double cpv[Np];      // volumetric heat capacity [J/m3/K]

   double tot_depth = 0.0;
   for (int i = 0; i<Np; i++) {
      zii[i]      = tot_depth; 
      tot_depth   = tot_depth + grid[i].dz;
      zi[i]       = tot_depth - 0.5*grid[i].dz;
      /* Thermal conductivity as in Jordan (1991) */
      tc[i] = tcair + (7.75e-5*grid[i].dens + 1.105e-6*pow(grid[i].dens,2))*(tcice - tcair); 
   }
   zii[Np] = tot_depth;

   for (int i = 1; i<Np; i++) {
      /* interpolation of thermal conductivity to the interface, 
         as CLM45 Tech note eq 6.11 
         Note: interfaces are located *below* the node with the same index */
      tci[i]   = tc[i]*tc[i-1]*(zi[i]-zi[i-1]) / (tc[i]*(zii[i]-zi[i-1]) + tc[i-1]*(zi[i]-zii[i]));
      /* heat flux (positively downwards) */
      Tgrad    = (grid[i].T - grid[i-1].T) / (zi[i]-zi[i-1]);
      Hfi[i]    = tci[i] * Tgrad; 
   }
   Hfi[0]   = 0.0; // bottom boundary condition
   Hfi[Np]  = tc[Np-1] * (Tskin-grid[Np-1].T) / (zii[Np] - zi[Np-1]); // top boundary condition
   tci[0]   = -999.0;
   tci[Np]  = -999.0;

   /* helper variables for LAPACK call */
   int dim = Np;
   int nrhs = 1;
   int ldb = dim;
   int info;

   std::vector<double> diag, rhs, dl, du;
   diag.resize(Np);
   rhs.resize(Np);
   dl.resize(Np-1);
   du.resize(Np-1);

   /* RIGHT HAND SIDE */
   for (int i = 0; i<Np; i++) {
      cpv[i]         = grid[i].dens * cpice; 
      rhs[i] = grid[i].T * (cpv[i] * grid[i].dz)/dt + 
         alpha * (Hfi[i+1]-Hfi[i]);
   }
   /* DIAGONAL */
   for (int i = 1; i<Np-1; i++) {
      diag[i] = (cpv[i]*grid[i].dz)/dt + 
         (1-alpha) * (tci[i+1]/(zi[i+1]-zi[i]) + tci[i]/(zi[i]-zi[i-1]) );
   }
   /* BOUNDARY CONDITIONS */
   /* at the bottom, the flux vanishes */
   diag[0]     = (cpv[0]*grid[0].dz)/dt + 
         (1-alpha) * ( tci[1]/(zi[1]-zi[0]) ); 
   /* at the top, Dirichlet boundary condition adds constant element to RHS */
   diag[Np-1] =  (cpv[Np-1]*grid[Np-1].dz)/dt + 
         (1-alpha) * ( tc[Np-1]/(zii[Np]-zi[Np-1]) + tci[Np-1]/(zi[Np-1]-zi[Np-2]) ); 
   rhs[Np-1] += (1-alpha)*tc[Np-1]*Tskin/(zii[Np]-zi[Np-1]);
   /* LOWER DIAGONAL */
   for (int i = 1; i<Np; i++) {
      dl[i-1] = - (1-alpha) * tci[i]/(zi[i]-zi[i-1]);
   }
   /* UPPER DIAGONAL */
   for (int i = 0; i<Np-1; i++) {
      du[i] = - (1-alpha) * tci[i+1]/(zi[i+1]-zi[i]);
   }
   /* SOLVE */
   dgtsv_(&dim, &nrhs, & *dl.begin(), & *diag.begin(), & *du.begin(), & *rhs.begin(), &ldb, &info);
   if (info != 0) {
      logger << "WARNING: Lapack routine DGTSV returned with info = " << info << std::endl;
   }
   for (int i = 0; i<Np; i++) {
      grid[i].T = rhs[i];
      if (grid[i].T > 280.) {
         logger << "ERROR: unrealistic temperature i = " << i << ", T[i] = " << grid[i].T << std::endl;
         logger << "ERROR: Probably this is caused by a high temperature gradient through a thin layer (too small dx /  too big dt)" << std::endl;
         std::abort();
      }
   }
}

} // namespace
