/* --------------------------------------------------------------------------------

   DESCRIPTION
     Abstract base class for ice core classes
   
   DATE
     15-JAN-2016
   
   AUTHOR
     L.vankampenhout@uu.nl
   -------------------------------------------------------------------------------- */

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

#include "icecoresite.h"
#include "logging.h"

namespace DSM{ 

ModelState::ModelState(Settings& s) {
   settings = s;
}

void ModelState::init(){
   // initialize grid with a very tiny top layer
   Layer layer;
   layer.T      = this->surfaceTemperature(current_time);
   layer.dens   = this->surfaceDensity(current_time);
   layer.dz     = 0.00001;
   grid.push_back(layer);
   is_initialized = true;
}

long ModelState::getSecondsSinceStart(){
   return current_time;
}


void ModelState::compactLigtenberg2011(long dt) {
   /*  Densification as equation [Ar10T] from Ligtenberg2011
      with additional scaling parameters MO for Antarctica
       */
   static const double E_c=60.e3; //#  (kJ/mol)
   static const double E_g=42.4e3; //# (kJ/mol)  
   double acc_year = annualIntegratedAccumulation(); 
   double Tsmean = annualMeanSurfaceTemperature();
   double MO, C;
   double layer_mass;
   for (int i = 0; i < gridsize(); i++) {
      layer_mass = grid[i].dz * grid[i].dens;
      if (grid[i].dens <= 550.){
         MO = 1.435 - 0.151 * log(acc_year);
         C=  0.07;
      } else {
         MO =  2.366 - 0.293 * log(acc_year);
         C =  0.03;
      }
      grid[i].dens = grid[i].dens + ((double)dt/sec_in_year)*MO*C*acc_year*g*(rho_i-grid[i].dens)*exp(-E_c/(R*grid[i].T)+E_g/(R*Tsmean));
      grid[i].dz = layer_mass/grid[i].dens; // # mass conservation
   }
}

void ModelState::compactAnderson1976(long dt) {
   /* Compaction due to destructive metamorphism (Anderson 1976) and overburden (Anderson 1976) */
   /* Formulas taken from CLM 4.5 Tech Note */
   static const double Tf = T0;
   static const double c3 = 2.777e-6;
   static const double c4 = 0.04;

   double overburden = 0;
   double layer_mass;
   
   for (int i = gridsize()-1; i >= 0; i--) {
      layer_mass = grid[i].dz * grid[i].dens;
      double P = overburden + 0.5*layer_mass;

      double c1 = exp(-0.046*(grid[i].dens-100.)); // # assuming density above 100 kg/m3
      double c2 = 1 ; //# assuming no liquid
      double cr1 = -c3*c2*c1*exp(-c4*(Tf-grid[i].T)); // # compaction due to destructive metamorphism
      double eta = settings.eta0 * exp(settings.c5*(Tf-grid[i].T) + settings.c6*grid[i].dens); //# viscocity, assuming no liquid
      double cr2 = -P/eta;

      double cr = cr1 + cr2;
      grid[i].dz = grid[i].dz * (1+(cr*dt));
      grid[i].dens = (layer_mass/grid[i].dz); // # mass conservation

      overburden = overburden + layer_mass;
   }
}

void ModelState::compactHerronLangway(long dt) {
   /*  Densification as Herron & Langway 1980 
      emperical model from the analysis of several Antarctic 
      and Greenland ice core density profiles  
   */
   double layer_mass;
   for (int i = gridsize()-1; i >= 0; i--) {
      layer_mass = grid[i].dz * grid[i].dens;
      if (grid[i].dens <= 550.){
         double k0 = 11. * exp(-10160. / (R*grid[i].T));
         grid[i].dens = grid[i].dens + ((double)dt/sec_in_year) * k0 * annualIntegratedAccumulation() * 1e-3 * (rho_i - grid[i].dens);
      } else { 
         double k1 = 575. * exp(-21400. / (R*grid[i].T));
         grid[i].dens = grid[i].dens + ((double)dt/sec_in_year) * k1 * sqrt(annualIntegratedAccumulation()*1e-3) * (rho_i - grid[i].dens);
      }
      grid[i].dz = layer_mass/grid[i].dens; // # mass conservation
   }
}

void ModelState::compactBarnola1991(long dt) {
   /*  Densification as Barnola & Pimienta 1991 
      From the surface to rho = 550 kg/m3 the Herron & Langway expression
      is used, below 550 the Pimienta expression.
     */
   double overburden = 0;
   double layer_mass;
   for (int i = gridsize()-1; i >= 0; i--) {
      layer_mass = grid[i].dz * grid[i].dens;
      if (grid[i].dens <= 550.){
         // Herron & Langway
         double k0 = 11. * exp(-10160./(R*grid[i].T));
         grid[i].dens = grid[i].dens + ((double)dt/sec_in_year) * k0 * annualIntegratedAccumulation() * 1e-3 * (rho_i - grid[i].dens);
      } else { 
         // Barnola & Pimienta 1991
         double A0 = 2.54e4; // in the original paper: [2.54e4 MPa^-3 s-1] 
         double Q = 60.e3; // activation energy 60e3 [J/mol] = 60 [kJ/mol]
         double ff;   
         if (grid[i].dens <= 800.){
            // f_e
            double alpha = -37.455;
            double beta = 99.743;
            double delta = -95.027;
            double gamma = 30.673;
            double dens_g_cm3 = grid[i].dens * 1e-3; // convert from kg/m3 to g/cm3
            double exponent = alpha*pow(dens_g_cm3,3) + beta*pow(dens_g_cm3,2) + delta*dens_g_cm3 + gamma;
            //logger << "DEBUG exponent = " << exponent << "\n";
            ff = pow(10, exponent);
         } else {
            // f_s
            ff = (3./16.) * (1-grid[i].dens/rho_i) / pow(1 - pow(1-grid[i].dens/rho_i, 1./3.),3);
         }
         double P = (overburden + 0.5*layer_mass)*g*1e-6/1.0;  // 1 Pa = 1 kg / (m * s2). We have unit area, convert to MPa
         //P = P*rho_i/grid[i].dens; // convert to grain-load stress LvK
         grid[i].dens = grid[i].dens + (double)dt * A0*exp(-Q/(R*grid[i].T)) * ff * pow(P,3.0) * grid[i].dens; 
      }
      grid[i].dz = layer_mass/grid[i].dens; // # mass conservation
      overburden = overburden + layer_mass;
      //logger << "DEBUG dens =  = " << grid[i].dens << "\n";
   }
}

void ModelState::compactSpencer2001(long dt) {
   /*  Densification as Spencer2001
      Three stages are defined:
         1) rho < 550.
         2) rho > 550. but below close-off density (COD)
         3) rho > COD
      In the paper, COD is made dependent on temperature at the close-off density Tc (usually
      quite close to T10m) using the following expression:
         COD = [944.6 - 6.15e-2 * Tc - 1.52e-5 * Tc**2] / [0.959 + 6.59e-4 * Tc - 3.62e-8 * Tc**2] 
     */
   double overburden = 0;
   double Peff;

   double Tc = grid[0].T; // Take temperature at bottom as Tc to avoid recalculation of Tc at every layer TODO: use layer temperature
   double rho_cod = (944.6 - 6.15e-2 * Tc - 1.52e-5 * pow(Tc,2)) / (0.959 + 6.59e-4 * Tc - 3.62e-8 * pow(Tc,2));

   double C1,C2,C3,C4,C5;
   double layer_mass;

   for (int i = gridsize()-1; i >= 0; i--) {
      layer_mass = grid[i].dz * grid[i].dens;
      Peff = (overburden + 0.5*layer_mass)*g;  // 1 Pa = 1 kg / (m * s2). We have unit area by which we divide (not shown)
 //      std::cout << "DEBUG Peff = " << Peff << std::endl;
      Peff = Peff*rho_i/grid[i].dens; // correction for relative density 
//      std::cout << "DEBUG Peff = " << Peff << std::endl;
      // TODO: correction for relative grain-contact area
      
      if (grid[i].dens < 550.) {
         C1 = 3.38e9; // per annum
         C2 = 46.8e3;  // converted from kJ to J
         C3 = 0.000121;
         C4 = -0.689;
         C5 = -0.149;
      } else if (grid[i].dens < rho_cod) {
         C1 = 9.06e8;
         C2 = 41.0e3;
         C3 = 0.0856;
         C4 = -1.05;
         C5 = -0.0202;
      } else {
         C1 = 1.38e7;
         C2 = 30.1e3;
         C3 = 0.284;
         C4 = -0.0734;
         C5 = 0.00322;
      }
      
      C1 = C1 /(double) sec_in_year; // per annum to per second
//      std::cout << "DEBUG C1 = " << C1 << std::endl;
//      std::cout << "DEBUG dens = " << grid[i].dens << std::endl;
      grid[i].dens = grid[i].dens + (double)dt * grid[i].dens * 
         C1 * exp(-C2/(R*grid[i].T)) *
         (1.-grid[i].dens/rho_i) *
         pow(1.-pow(1.-grid[i].dens/rho_i,C3),C4) * 
         pow(Peff,C5);

//      std::cout << "DEBUG dens = " << grid[i].dens << std::endl;
      grid[i].dz = layer_mass/grid[i].dens; // # mass conservation
      overburden = overburden + layer_mass;
   }
}
 
void ModelState::heatDiffusion(long dt) {
   /* Uses forward Euler method for diffusion equation*/ 
   double T_new[gridsize()];
   double td[gridsize()];
   double tc;

   for (int i = 0; i < gridsize(); i++) {
      // # thermal conductivity 
      tc = tcair + (7.75e-5*grid[i].dens + 1.105e-6*pow(grid[i].dens,2))*(tcice - tcair);
      // # thermal diffusivity
      td[i] = tc/(grid[i].dens * cp);
   }

   // Dirichlet boundary condition at the top
   int Np = gridsize();
   T_new[Np-1] = this->surfaceTemperature(current_time); // # top boundary equal to surface temperature

   // Neumann boundary condition at the bottom
   T_new[0] = td[0]*dt/grid[0].dz * (grid[1].T-grid[0].T) + grid[0].T;

   for (int i = 1; i < gridsize()-1; i++) {
      //double stab_crit = td[i] * dt / pow(grid[i].dz,2);
      //if (stab_crit > 0.5) {
      //   logger << "ERROR: stability criterium violated, r = " << stab_crit << std::endl;
      //   std::abort();
      //} else {
         T_new[i] = td[i]*dt/pow(grid[i].dz,2)*(grid[i-1].T-2*grid[i].T+grid[i+1].T) + grid[i].T; // # Forward in Time Centered in Space
      //}
   }

   // Update temperatures
   for (int i = 0; i < gridsize(); i++) {
      grid[i].T = T_new[i];
   }
}

double ModelState::getDepthOfDensity(double dens) {
   double diffi, diffip1;
   double dzi, dzip1;
   int tidx = gridsize()-1; // index of top layer
   if (grid[tidx].dens >= dens) return grid[tidx].dz/2; // node depth of top layer

   double tot_depth = grid[tidx].dz; // interface depth
   for (int i = tidx-1; i >= 0; i--) {
      if (grid[i].dens >= dens ){
         diffi = std::abs(grid[i].dens - dens); // delta 1
         diffip1 = std::abs(grid[i+1].dens - dens); // delta 2
         
         dzi = tot_depth + grid[i].dz / 2; // node depth
         dzip1 = tot_depth - grid[i+1].dz / 2; // node depth

         return diffi/(diffi+diffip1) * dzip1 + diffip1/(diffi+diffip1) * dzi; // linear interpolation
      }
      tot_depth = tot_depth + grid[i].dz;
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
   return getDepthOfDensity(830.);
}

double ModelState::totalDepth() {
   double tot_depth = 0;
   for (int i = gridsize()-1; i >= 0; i--) {
      tot_depth = tot_depth + grid[i].dz;
   };
   return tot_depth;
}

void ModelState::printIceCoreSummary() {
   logger << "Core = "+this->toString() << std::endl;
   double tot_mass = 0;
   double tot_depth = 0;
   for (int i = gridsize()-1; i >= 0; i--) {
      tot_mass += grid[i].dz * grid[i].dens;
      tot_depth = tot_depth + grid[i].dz;
      //logger << "Layer " << gridsize()-i-1 << " : h=" << grid[i].dz << ", rho=" << grid[i].dens << ", T=" << grid[i].T << std::endl;
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

void ModelState::writeFiles() {
   std::ofstream f_layers;
   f_layers.open("layers.csv");
   f_layers << gridsize() << std::endl;
   for (int i = gridsize()-1; i >= 0; i--) {
      f_layers << grid[i].dz << "," << grid[i].dz * grid[i].dens<< "," << grid[i].dens << "," << grid[i].T << std::endl;
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

std::string ModelState::toString()
{ 
   std::ostringstream s; 
   s << settings.dm;
   if (settings.have_diffusion) {
      s << "_heatTrue";
   } else {
      s << "_heatFalse";
   }
   return s.str();
}

double ModelState::surfaceDensity(long time) {
   /* after Helsen (2008) */
   return -154.91+1.4266*(73.6+1.06*Ts_ann_mean+0.0669*acc_ann_mean+4.77*w10m_ann_mean);
   //return settings.rho_s;
}

}
