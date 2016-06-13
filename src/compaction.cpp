
#include <math.h>

#include "constants.h"
#include "config.h"
#include "logging.h"
#include "compaction.h"
#include "modelstate.h"
#include "dynamicmodel.h"
#include "meteo.h"
#include "utils.h"

namespace DSM{ 

std::unique_ptr<Compaction> instantiate_compaction(ModelState& mstate, DynamicModel& dm){
   const char * option_name = "physics:which_compaction";
   int which_compaction = config.getInt(option_name, false, 0, 6, 1);

   switch (which_compaction) {
      case 0   : return { std::make_unique<CompactionNone>(mstate, dm) };
      case 1   : return { std::make_unique<CompactionHerronLangway>(mstate, dm) };
      case 2   : return { std::make_unique<CompactionAnderson>(mstate, dm) };
      case 3   : return { std::make_unique<CompactionBarnolaPimienta>(mstate, dm) };
      case 4   : return { std::make_unique<CompactionLigtenberg>(mstate, dm) };
      case 5   : return { std::make_unique<CompactionCROCUS>(mstate, dm) };
      case 6   : return { std::make_unique<CompactionCROCUSDriftOnly>(mstate, dm) };
      default:
         logger << "ERROR: unknown value: " << which_compaction << " for config option " << option_name << std::endl;
         std::abort();
   }
}

Compaction::Compaction(ModelState& mstate, DynamicModel& dm) : _mstate(mstate), _dm(dm) { }

CompactionNone::CompactionNone(ModelState& mstate, DynamicModel& dm) : Compaction(mstate, dm) {
   logger << "CompactionNone()" << std::endl; 
}

CompactionHerronLangway::CompactionHerronLangway(ModelState& mstate, DynamicModel& dm) : Compaction(mstate, dm) {
   logger << "CompactionHerronLangway()" << std::endl; 
}

CompactionAnderson::CompactionAnderson(ModelState& mstate, DynamicModel& dm) : Compaction(mstate, dm) { 
   logger << "CompactionAnderson()" << std::endl; 
}

CompactionBarnolaPimienta::CompactionBarnolaPimienta(ModelState& mstate, DynamicModel& dm) : Compaction(mstate, dm) { 
   logger << "CompactionBarnolaPimienta()" << std::endl; 
}

CompactionLigtenberg::CompactionLigtenberg(ModelState& mstate, DynamicModel& dm) : Compaction(mstate, dm) { 
   logger << "CompactionLigtenberg()" << std::endl; 
}

CompactionCROCUS::CompactionCROCUS(ModelState& mstate, DynamicModel& dm) : Compaction(mstate, dm) { 
   logger << "CompactionCROCUS()" << std::endl; 
}

CompactionCROCUSDriftOnly::CompactionCROCUSDriftOnly(ModelState& mstate, DynamicModel& dm) : CompactionCROCUS(mstate, dm) { 
   logger << "CompactionCROCUSDriftOnly()" << std::endl; 
}

void CompactionLigtenberg::compaction() {
   /*  Densification as equation [Ar10T] from Ligtenberg2011
      with additional scaling parameters MO for Antarctica */
   Grid& grid = _mstate.getGrid();
   const double dt = (double)_dm.getDt();
   static const double E_c=60.e3; //#  (kJ/mol)
   static const double E_g=42.4e3; //# (kJ/mol)  
   const double acc_year = _mstate.getMeteo().annualAcc(); 
   const double Tsmean = _mstate.getMeteo().annualTskin();
   double MO, C;
   double layer_mass;
   for (int i = 0; i < grid.size(); i++) {
      layer_mass = grid[i].dz * grid[i].dens;
      if (grid[i].dens <= 550.){
         MO = 1.435 - 0.151 * log(acc_year);
         C=  0.07;
      } else {
         MO =  2.366 - 0.293 * log(acc_year);
         C =  0.03;
      }
      grid[i].dens = grid[i].dens + (dt/sec_in_year)*MO*C*acc_year*g*(rho_i-grid[i].dens)*exp(-E_c/(R*grid[i].T)+E_g/(R*Tsmean));
      grid[i].dz = layer_mass/grid[i].dens; // mass conservation
   }
}

void CompactionAnderson::compaction() {
   /* Compaction due to destructive metamorphism (Anderson 1976) and overburden (Anderson 1976)
      Formulas taken from CLM 4.5 Tech Note */
   Grid& grid = _mstate.getGrid();
   const double dt = (double)_dm.getDt();

   static const double eta0 = 9e5;
   static const double c5 = 0.08;
   static const double c6 = 0.023;

   double overburden = 0;
   double layer_mass;
   
   for (int i = grid.size()-1; i >= 0; i--) {
      layer_mass = grid[i].dz * grid[i].dens;
      const double P = overburden + 0.5*layer_mass;
      const double eta = eta0 * exp(c5*(T0-grid[i].T) + c6*grid[i].dens); //# viscocity, assuming no liquid
      const double cr2 = -P/eta;

      //double cr = cr1 + cr2;
      grid[i].dz = grid[i].dz * (1+(cr2*dt));
      grid[i].dens = (layer_mass/grid[i].dz); // # mass conservation

      overburden = overburden + layer_mass;
   }
}

void CompactionHerronLangway::compaction() {
   /* Densification as Herron & Langway 1980 
      empirical model from the analysis of several Antarctic 
      and Greenland ice core density profiles  */
   Grid& grid = _mstate.getGrid();
   const double acc_year = _mstate.getMeteo().annualAcc(); 
   const double dt = (double)_dm.getDt();

   double layer_mass;
   for (int i = grid.size()-1; i >= 0; i--) {
      layer_mass = grid[i].dz * grid[i].dens;
      if (grid[i].dens <= 550.){
         const double k0 = 11. * exp(-10160. / (R*grid[i].T));
         grid[i].dens = grid[i].dens + (dt/(double)sec_in_year) * k0 * acc_year * 1e-3 * (rho_i - grid[i].dens);
      } else { 
         const double k1 = 575. * exp(-21400. / (R*grid[i].T));
         grid[i].dens = grid[i].dens + (dt/(double)sec_in_year) * k1 * sqrt(acc_year*1e-3) * (rho_i - grid[i].dens);
      }
      grid[i].dz = layer_mass/grid[i].dens; // # mass conservation
   }
}

void CompactionBarnolaPimienta::compaction() {
   /* Densification as Barnola & Pimienta 1991 
      From the surface to rho = 550 kg/m3 the Herron & Langway expression
      is used, below 550 the Pimienta expression.  */
   Grid& grid = _mstate.getGrid();
   const double acc_year = _mstate.getMeteo().annualAcc(); 
   const double dt = (double)_dm.getDt();

   double overburden = 0;
   double layer_mass;
   for (int i = grid.size()-1; i >= 0; i--) {
      layer_mass = grid[i].dz * grid[i].dens;
      if (grid[i].dens <= 550.){
         // Herron & Langway
         const double k0 = 11. * exp(-10160./(R*grid[i].T));
         grid[i].dens = grid[i].dens + (dt/sec_in_year) * k0 * acc_year * 1e-3 * (rho_i - grid[i].dens);
      } else { 
         // Barnola & Pimienta 1991
         static const double A0 = 2.54e4; // in the original paper: [2.54e4 MPa^-3 s-1] 
         static const double Q = 60.e3; // activation energy 60e3 [J/mol] = 60 [kJ/mol]
         double ff;   
         if (grid[i].dens <= 800.){
            // f_e
            static const double alpha = -37.455;
            static const double beta = 99.743;
            static const double delta = -95.027;
            static const double gamma = 30.673;
            const double dens_g_cm3 = grid[i].dens * 1e-3; // convert from kg/m3 to g/cm3
            const double exponent = alpha*pow(dens_g_cm3,3) + beta*pow(dens_g_cm3,2) + delta*dens_g_cm3 + gamma;
            //logger << "DEBUG exponent = " << exponent << "\n";
            ff = pow(10, exponent);
         } else {
            // f_s
            ff = (3./16.) * (1-grid[i].dens/rho_i) / pow(1 - pow(1-grid[i].dens/rho_i, 1./3.),3);
         }
         const double P = (overburden + 0.5*layer_mass)*g*1e-6/1.0;  // 1 Pa = 1 kg / (m * s2). We have unit area, convert to MPa
         //P = P*rho_i/grid[i].dens; // convert to grain-load stress LvK
         grid[i].dens = grid[i].dens + dt * A0*exp(-Q/(R*grid[i].T)) * ff * pow(P,3.0) * grid[i].dens; 
      }
      grid[i].dz = layer_mass/grid[i].dens; // # mass conservation
      overburden = overburden + layer_mass;
      //logger << "DEBUG dens =  = " << grid[i].dens << "\n";
   }
}

void CompactionCROCUS::compaction() {
   /* Compaction due to overburden (Vionnet et al, 2012) 
      TODO: wet snow not implemented (f1 factor)
      TODO: depth hoar formation not implemented (f2 factor)
      */
   Grid& grid = _mstate.getGrid();
   const double dt = (double)_dm.getDt(); // time in seconds

   static const double eta0 = 7.62237e6; 
   static const double an = 0.1; 
   static const double bn = 0.023; 
   static const double cn = 1./250.;
   static const double f1f2 = 4.0;
   
   double overburden = 0;
   double layer_mass, P, eta, cr;
   
   for (int i = grid.size()-1; i >= 0; i--) {
      layer_mass = grid[i].dz * grid[i].dens;
      P = (overburden + 0.5*layer_mass)*g; // sigma_i
      eta = f1f2 * eta0 * grid[i].dens * cn * exp(an*(T0-grid[i].T) + bn*grid[i].dens); // viscocity, assuming no liquid
      cr = -P / eta;
      grid[i].dz = grid[i].dz * (1+(cr*dt));
      grid[i].dens = (layer_mass/grid[i].dz); // mass conservation
      overburden = overburden + layer_mass;
   }
   compactionWindDrift();
}

void CompactionCROCUS::compactionWindDrift(){
   /* Compaction due to wind drift */
   Grid& grid = _mstate.getGrid();
   const int Np = grid.size();
   if (Np < 2) return; // no metamorphism
   static const double rho_min = 50.; // minimum density [kg/m3]
   static const double rho_max = 350.; // maximum density [kg/m3]
   static const double tau_ref = (double)48 * 3600; // reference time [s]
   const double U = _mstate.getMeteo().surfaceWind();
   const double dt = (double)_dm.getDt(); // time step [s]
   double Frho;
   double MO; // mobility index
   double SI; // driftability index
   double tau; // time characteristic
   double zi = 0.; // pseudo-depth in m
   double gamma_drift;
   double u0, afac; // TUNING
   for (int i = grid.size()-1; i >= 0; i--) {
      Frho = 1.25 - 0.0042*(std::max(rho_min, grid[i].dens)-rho_min);
      if (doubles_equal(grid[i].dnd, 0.0)) {
         // Non-dendritic snow
         MO = 0.34 * (-0.583*grid[i].gs - 0.833*grid[i].sph + 0.833) + 0.66*Frho;
      } else { 
         // dendritic snow
         MO = 0.34 * (0.75*grid[i].dnd - 0.5*grid[i].sph + 0.5) + 0.66*Frho;
      }

      /* LvK: tuning
      static const double wind_fac = -0.085;
      u0 = log((-1.-MO)/-2.868) / -0.085;
      u0 = u0 - 2.;
      afac = (-1. - MO) / exp(wind_fac * u0 ); // based on exponent factor, determine pre-factor by setting SI = 0 to the same U-value
      SI = afac * exp(wind_fac*U) + 1 + MO; 
      */
      SI = -2.868 * exp(-0.085*U)+1.+MO;

      if (SI <= 0.0) break; // first non-mobile layer has been reached
      SI = std::min(SI, 3.25);
      gamma_drift = std::max(0., SI*exp(-zi/0.1));
      tau = tau_ref / gamma_drift;
      if (doubles_equal(grid[i].dnd, 0.0)) {
         // Non-dendritic snow
         grid[i].sph = grid[i].sph + dt * (1-grid[i].sph)/tau;
         grid[i].gs = grid[i].gs + dt * 5e-4/tau;
      } else { 
         // dendritic snow
         grid[i].dnd = grid[i].dnd + dt * grid[i].dnd/(2*tau);
         grid[i].sph = grid[i].sph + dt * (1-grid[i].sph)/tau;
      }
      // consistency checks
      grid[i].dnd = std::max(std::min(1.0, grid[i].dnd), 0.0);
      grid[i].sph = std::max(std::min(1.0, grid[i].sph), 0.0);
      grid[i].gs = std::max(std::min(fs_gs_max, grid[i].gs), 0.0);

      // density evolution
      //logger << "DEBUG dens = " << grid[i].dens << ", U = " << U << ", zi = " << zi << ", drho / 24h = " << 24*3600 * std::max(0.0, rho_max - grid[i].dens)/tau << std::endl;
      grid[i].dens = grid[i].dens + dt * std::max(0.0, rho_max - grid[i].dens)/tau;

      if (gamma_drift > 100.0) {
         logger << "warning: large gamma\n";
         logger << "i / Np = " << i << " / " << grid.size() << "\n";
         logger << "U = " << U << "\n";
         logger << "SI = " << SI << "\n";
         logger << "MO = " << MO << "\n";
         logger << "zi = " << zi << "\n";
         //logger << "dens_old = " << dens_old << "\n";
         //logger << "dens_new = " << grid[i].dens << "\n";
         //logger << "d= " << s_old << "\n";
         //logger << "s= " << d_old << "\n";
         //logger << "gs= " << gs_old << "\n";
      }
      logger.flush();  
      zi += grid[i].dz * (3.25 - SI);
 
   }
}

void CompactionCROCUSDriftOnly::compaction() {
   compactionWindDrift();
}


// void CompactionSpencer::compaction() {
//    /*  Densification as Spencer2001
//       Three stages are defined:
//          1) rho < 550.
//          2) rho > 550. but below close-off density (COD)
//          3) rho > COD
//       In the paper, COD is made dependent on temperature at the close-off density Tc (usually
//       quite close to T10m) using the following expression:
//          COD = [944.6 - 6.15e-2 * Tc - 1.52e-5 * Tc**2] / [0.959 + 6.59e-4 * Tc - 3.62e-8 * Tc**2] 
//      */
//    double overburden = 0;
//    double Peff;
// 
//    double Tc = grid[0].T; // Take temperature at bottom as Tc to avoid recalculation of Tc at every layer TODO: use layer temperature
//    double rho_cod = (944.6 - 6.15e-2 * Tc - 1.52e-5 * pow(Tc,2)) / (0.959 + 6.59e-4 * Tc - 3.62e-8 * pow(Tc,2));
// 
//    double C1,C2,C3,C4,C5;
//    double layer_mass;
// 
//    for (int i = grid.size()-1; i >= 0; i--) {
//       layer_mass = grid[i].dz * grid[i].dens;
//       Peff = (overburden + 0.5*layer_mass)*g;  // 1 Pa = 1 kg / (m * s2). We have unit area by which we divide (not shown)
//  //      std::cout << "DEBUG Peff = " << Peff << std::endl;
//       Peff = Peff*rho_i/grid[i].dens; // correction for relative density 
// //      std::cout << "DEBUG Peff = " << Peff << std::endl;
//       // TODO: correction for relative grain-contact area
//       
//       if (grid[i].dens < 550.) {
//          C1 = 3.38e9; // per annum
//          C2 = 46.8e3;  // converted from kJ to J
//          C3 = 0.000121;
//          C4 = -0.689;
//          C5 = -0.149;
//       } else if (grid[i].dens < rho_cod) {
//          C1 = 9.06e8;
//          C2 = 41.0e3;
//          C3 = 0.0856;
//          C4 = -1.05;
//          C5 = -0.0202;
//       } else {
//          C1 = 1.38e7;
//          C2 = 30.1e3;
//          C3 = 0.284;
//          C4 = -0.0734;
//          C5 = 0.00322;
//       }
//       
//       C1 = C1 /(double) sec_in_year; // per annum to per second
// //      std::cout << "DEBUG C1 = " << C1 << std::endl;
// //      std::cout << "DEBUG dens = " << grid[i].dens << std::endl;
//       grid[i].dens = grid[i].dens + (double)dt * grid[i].dens * 
//          C1 * exp(-C2/(R*grid[i].T)) *
//          (1.-grid[i].dens/rho_i) *
//          pow(1.-pow(1.-grid[i].dens/rho_i,C3),C4) * 
//          pow(Peff,C5);
// 
// //      std::cout << "DEBUG dens = " << grid[i].dens << std::endl;
//       grid[i].dz = layer_mass/grid[i].dens; // # mass conservation
//       overburden = overburden + layer_mass;
//    }
// }
//  

} // namespace
