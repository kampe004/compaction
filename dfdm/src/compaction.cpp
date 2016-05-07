#include <cmath>

#include "constants.h"
#include "config.h"
#include "logging.h"
#include "compaction.h"
#include "modelstate.h"
#include "dynamicmodel.h"
#include "meteo.h"

namespace DSM{ 

std::unique_ptr<Compaction> instantiate_compaction(ModelState& mstate, DynamicModel& dm){
   const char * option_name = "physics:which_compaction";
   int which_compaction = config.getInt(option_name, false, 0, 4, 0);

   switch (which_compaction) {
      case 0   : return { std::make_unique<CompactionHerronLangway>(mstate, dm) };
      case 1   : return { std::make_unique<CompactionAnderson>(mstate, dm) };
      case 2   : return { std::make_unique<CompactionBarnolaPimienta>(mstate, dm) };
      case 3   : return { std::make_unique<CompactionLigtenberg>(mstate, dm) };
      case 4   : return { std::make_unique<CompactionCROCUS>(mstate, dm) };
      default:
         logger << "ERROR: unknown value: " << which_compaction << " for config option " << option_name << std::endl;
         std::abort();
   }
}

Compaction::Compaction(ModelState& mstate, DynamicModel& dm) : _mstate(mstate), _dm(dm) { } 
CompactionHerronLangway::CompactionHerronLangway(ModelState& mstate, DynamicModel& dm) : Compaction(mstate, dm) { } 
CompactionAnderson::CompactionAnderson(ModelState& mstate, DynamicModel& dm) : Compaction(mstate, dm) { } 
CompactionBarnolaPimienta::CompactionBarnolaPimienta(ModelState& mstate, DynamicModel& dm) : Compaction(mstate, dm) { } 
CompactionLigtenberg::CompactionLigtenberg(ModelState& mstate, DynamicModel& dm) : Compaction(mstate, dm) { } 
CompactionCROCUS::CompactionCROCUS(ModelState& mstate, DynamicModel& dm) : Compaction(mstate, dm) { } 

void CompactionLigtenberg::compaction() {
   /*  Densification as equation [Ar10T] from Ligtenberg2011
      with additional scaling parameters MO for Antarctica */
   Grid& grid = _mstate.getGrid();
   const double dt = _dm.getDt();
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
      grid[i].dens = grid[i].dens + ((double)dt/sec_in_year)*MO*C*acc_year*g*(rho_i-grid[i].dens)*exp(-E_c/(R*grid[i].T)+E_g/(R*Tsmean));
      grid[i].dz = layer_mass/grid[i].dens; // mass conservation
   }
}

void CompactionAnderson::compaction() {
   /* Compaction due to destructive metamorphism (Anderson 1976) and overburden (Anderson 1976)
      Formulas taken from CLM 4.5 Tech Note */
   Grid& grid = _mstate.getGrid();
   const double dt = _dm.getDt();

   static const double Tf = T0;
   static const double eta0 = 9e5;
   static const double c5 = 0.08;
   static const double c6 = 0.023;

   double overburden = 0;
   double layer_mass;
   
   for (int i = grid.size()-1; i >= 0; i--) {
      layer_mass = grid[i].dz * grid[i].dens;
      const double P = overburden + 0.5*layer_mass;
      const double eta = eta0 * exp(c5*(Tf-grid[i].T) + c6*grid[i].dens); //# viscocity, assuming no liquid
      const double cr2 = -P/eta;

      //double cr = cr1 + cr2;
      grid[i].dz = grid[i].dz * (1+(cr2*dt));
      grid[i].dens = (layer_mass/grid[i].dz); // # mass conservation

      overburden = overburden + layer_mass;
   }
}

void CompactionHerronLangway::compaction() {
   /*  Densification as Herron & Langway 1980 
      emperical model from the analysis of several Antarctic 
      and Greenland ice core density profiles  */
   Grid& grid = _mstate.getGrid();
   const double acc_year = _mstate.getMeteo().annualAcc(); 
   const double dt = _dm.getDt();

   double layer_mass;
   for (int i = grid.size()-1; i >= 0; i--) {
      layer_mass = grid[i].dz * grid[i].dens;
      if (grid[i].dens <= 550.){
         const double k0 = 11. * exp(-10160. / (R*grid[i].T));
         grid[i].dens = grid[i].dens + ((double)dt/sec_in_year) * k0 * acc_year * 1e-3 * (rho_i - grid[i].dens);
      } else { 
         const double k1 = 575. * exp(-21400. / (R*grid[i].T));
         grid[i].dens = grid[i].dens + ((double)dt/sec_in_year) * k1 * sqrt(acc_year*1e-3) * (rho_i - grid[i].dens);
      }
      grid[i].dz = layer_mass/grid[i].dens; // # mass conservation
   }
}

void CompactionBarnolaPimienta::compaction() {
   /*  Densification as Barnola & Pimienta 1991 
      From the surface to rho = 550 kg/m3 the Herron & Langway expression
      is used, below 550 the Pimienta expression.  */
   Grid& grid = _mstate.getGrid();
   const double acc_year = _mstate.getMeteo().annualAcc(); 
   const double dt = _dm.getDt();

   double overburden = 0;
   double layer_mass;
   for (int i = grid.size()-1; i >= 0; i--) {
      layer_mass = grid[i].dz * grid[i].dens;
      if (grid[i].dens <= 550.){
         // Herron & Langway
         const double k0 = 11. * exp(-10160./(R*grid[i].T));
         grid[i].dens = grid[i].dens + ((double)dt/sec_in_year) * k0 * acc_year * 1e-3 * (rho_i - grid[i].dens);
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
         grid[i].dens = grid[i].dens + (double)dt * A0*exp(-Q/(R*grid[i].T)) * ff * pow(P,3.0) * grid[i].dens; 
      }
      grid[i].dz = layer_mass/grid[i].dens; // # mass conservation
      overburden = overburden + layer_mass;
      //logger << "DEBUG dens =  = " << grid[i].dens << "\n";
   }
}

void CompactionCROCUS::compaction() {
   return;
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
