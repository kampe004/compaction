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

namespace Densification{ 

IceCoreSite::IceCoreSite(Settings& s) {
    settings = s;
}

void IceCoreSite::init(){
    // initialize grid with a very tiny top layer
    Layer layer;
    layer.T       = this->surfaceTemperature(current_time);
    layer.dens    = this->surfaceDensity(current_time);
    layer.dz      = 0.00001;
    layer.mass    = layer.dz * layer.dens; 
    grid.push_back(layer);
    is_initialized = true;
}

void IceCoreSite::runTimeStep(long dt) {
    if (!is_initialized){
        logger << "ERROR: core is not initialized, run init() first!" << std::endl;
        std::abort();
    }
    accumulate(dt);
    compact(dt);

    // check minimum thickness, start with second layer (index N-2)
    int i = grid.size()-2;
    while (i >= 0 ) { // Use while instead of for-loop as size may change during loop
        if (grid[i].dz < dzmin && gridsize() > 1){
            //logger << "DEBUG: merging" << std::endl;
            // Layer is too thin, combine with neighbor
            if (i == 0  || grid[i-1].dz > grid[i+1].dz  ) {
                /* CASE 1: Bottom layer can only be combined with the one above */
                /* CASE 2: Determine smallest neighbour */
                grid[i+1].mass  = grid[i+1].mass + grid[i].mass;
                grid[i+1].dz    = grid[i+1].dz + grid[i].dz;
                grid[i+1].dens  = grid[i+1].mass / grid[i+1].dz;
            } else {
                grid[i-1].mass  = grid[i-1].mass + grid[i].mass;
                grid[i-1].dz    = grid[i-1].dz + grid[i].dz;
                grid[i-1].dens  = grid[i-1].mass / grid[i-1].dz;
            }
            grid.erase(grid.begin() + i); // NOTE: this is very inefficient
        }
        i--;
    }

    // remove bottom layers that have attained ice density
    //while(grid.front().dens > 900.) 
    //    grid.erase(grid.begin());
    if (grid.front().dens > 900.) {
    for (int i = 0; i < gridsize(); i++) {
        logger << "i = " << i << " dens = " << grid[i].dens << "\n";
       } 
    }

    if (settings.have_diffusion && gridsize() > 1)
        heatDiffusion(dt);

    // update internal time variable
    current_time = current_time + dt;
}

long IceCoreSite::getSecondsSinceStart(){
    return current_time;
}

void IceCoreSite::accumulate(long dt) {
    double dens = this->surfaceDensity(current_time);
    double acc = this->accumulationRate(current_time) * 1e-3 * dt; // convert from [mm/s] to [m]
    double dz = (rho_w/dens) * acc;

    Layer& top = grid.back();  // top layer
    top.dens = top.dz/(top.dz+dz)*top.dens + dz/(top.dz+dz)*dens;
    top.dz = top.dz + dz;
    top.mass = top.dz * top.dens;

    if (top.dz > dzmax ) {
        // Split in unequal parts 4:1
        Layer layer;
        layer.dens  = top.dens;
        layer.mass  = top.mass/5;
        layer.dz    = top.dz/5;
        layer.T     = top.T;

        top.dz = top.dz * 4/5; 
        top.mass = top.mass * 4/5;

        grid.push_back(layer);
    }
}

void IceCoreSite::compact(long dt) {
    // enum class DensificationMethod {Ligtenberg2011, Anderson1976, Barnola1991, Spencer2001, BarnolaSpencer}
    switch (settings.dm) {
        case DensificationMethod::Ligtenberg2011 : compactLigtenberg2011(dt); break;
        case DensificationMethod::Anderson1976   : compactAnderson1976(dt); break;
        case DensificationMethod::Barnola1991    : compactBarnola1991(dt); break;
        case DensificationMethod::Spencer2001    : //compactSpencer2001(dt); break;
        case DensificationMethod::BarnolaSpencer : //compactBarnolaSpencer(dt); break;
        case DensificationMethod::HerronLangway  : compactHerronLangway(dt); break;
        default:
            logger << "ERROR: programmer error: unknown densification method" << std::endl;
            std::abort();
    }
}

void IceCoreSite::compactLigtenberg2011(long dt) {
    /*  Densification as equation [Ar10T] from Ligtenberg2011
        with additional scaling parameters MO for Antarctica
         */
    static const double E_c=60.e3; //#  (kJ/mol)
    static const double E_g=42.4e3; //# (kJ/mol)  
    double acc_year = annualIntegratedAccumulation(); 
    double Tsmean = annualMeanSurfaceTemperature();
    double MO, C;
    for (int i = 0; i < gridsize(); i++) {
        if (grid[i].dens <= 550.){
            MO = 1.435 - 0.151 * log(acc_year);
            C=  0.07;
        } else {
            MO =  2.366 - 0.293 * log(acc_year);
            C =  0.03;
        }
        grid[i].dens = grid[i].dens + ((double)dt/sec_in_year)*MO*C*acc_year*g*(rho_i-grid[i].dens)*exp(-E_c/(R*grid[i].T)+E_g/(R*Tsmean));
        grid[i].dz = grid[i].mass/grid[i].dens; // # mass conservation
    }
}

void IceCoreSite::compactAnderson1976(long dt) {
    /* Compaction due to destructive metamorphism (Anderson 1976) and overburden (Anderson 1976) */
    /* Formulas taken from CLM 4.5 Tech Note */
    static const double Tf = T0;
    static const double c3 = 2.777e-6;
    static const double c4 = 0.04;

    double overburden = 0;
    //for (int i = 0; i < gridsize(); i++) {
    for (int i = gridsize()-1; i >= 0; i--) {
        double P = overburden + 0.5*grid[i].mass;

        double c1 = exp(-0.046*(grid[i].dens-100.)); // # assuming density above 100 kg/m3
        double c2 = 1 ; //# assuming no liquid
        double cr1 = -c3*c2*c1*exp(-c4*(Tf-grid[i].T)); // # compaction due to destructive metamorphism
        double eta = settings.eta0 * exp(settings.c5*(Tf-grid[i].T) + settings.c6*grid[i].dens); //# viscocity, assuming no liquid
        double cr2 = -P/eta;

        double cr = cr1 + cr2;
        grid[i].dz = grid[i].dz * (1+(cr*dt));
        grid[i].dens = (grid[i].mass/grid[i].dz); // # mass conservation

        overburden = overburden + grid[i].mass;
    }
}

void IceCoreSite::compactHerronLangway(long dt) {
    /*  Densification as Herron & Langway 1980 
        emperical model from the analysis of several Antarctic 
        and Greenland ice core density profiles  
    */
    double overburden = 0;
    for (int i = gridsize()-1; i >= 0; i--) {
        if (grid[i].dens <= 550.){
            double k0 = 11. * exp(-10160. / (R*grid[i].T));
            grid[i].dens = grid[i].dens + ((double)dt/sec_in_year) * k0 * annualIntegratedAccumulation() * 1e-3 * (rho_i - grid[i].dens);
        } else { 
            double k1 = 575. * exp(-21400. / (R*grid[i].T));
            grid[i].dens = grid[i].dens + ((double)dt/sec_in_year) * k1 * sqrt(annualIntegratedAccumulation()*1e-3) * (rho_i - grid[i].dens);
        }
        grid[i].dz = grid[i].mass/grid[i].dens; // # mass conservation
    }
}

void IceCoreSite::compactBarnola1991(long dt) {
    /*  Densification as Barnola & Pimienta 1991 
        From the surface to rho = 550 kg/m3 the Herron & Langway expression
        is used, below 550 the Pimienta expression.
      */
    double overburden = 0;
    for (int i = gridsize()-1; i >= 0; i--) {
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
            double P = (overburden + 0.5*grid[i].mass)*g*1e-6/1.0;  // 1 Pa = 1 kg / (m * s2). We have unit area, convert to MPa
            grid[i].dens = grid[i].dens + (double)dt * A0*exp(-Q/(R*grid[i].T)) * ff * pow(P,3.0) * grid[i].dens; 
        }
        grid[i].dz = grid[i].mass/grid[i].dens; // # mass conservation
        overburden = overburden + grid[i].mass;
        //logger << "DEBUG dens =  = " << grid[i].dens << "\n";
    }
}

void IceCoreSite::heatDiffusion(long dt) {
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
        //    logger << "ERROR: stability criterium violated, r = " << stab_crit << std::endl;
        //    std::abort();
        //} else {
            T_new[i] = td[i]*dt/pow(grid[i].dz,2)*(grid[i-1].T-2*grid[i].T+grid[i+1].T) + grid[i].T; // # Forward in Time Centered in Space
        //}
    }

    // Update temperatures
    for (int i = 0; i < gridsize(); i++) {
        grid[i].T = T_new[i];
    }
}

double IceCoreSite::getDepthOfDensity(double dens) {
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

bool IceCoreSite::hasReachedDensity(double dens){
    return this->getDepthOfDensity(dens) > 0.; 
}

double IceCoreSite::getZ550() {
    return getDepthOfDensity(550.);
}

double IceCoreSite::getZ830() {
    return getDepthOfDensity(830.);
}

double IceCoreSite::totalDepth() {
    double tot_depth = 0;
    for (int i = gridsize()-1; i >= 0; i--) {
        tot_depth = tot_depth + grid[i].dz;
    };
    return tot_depth;
}

void IceCoreSite::printIceCoreSummary() {
    logger << "Core = "+this->toString() << std::endl;
    double tot_mass = 0;
    double tot_depth = 0;
    for (int i = gridsize()-1; i >= 0; i--) {
        tot_mass = tot_mass + grid[i].mass;
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

void IceCoreSite::writeFiles() {
    std::ofstream f_layers;
    f_layers.open("layers.csv");
    f_layers << gridsize() << std::endl;
    for (int i = gridsize()-1; i >= 0; i--) {
        f_layers << grid[i].dz << "," << grid[i].mass << "," << grid[i].dens << "," << grid[i].T << std::endl;
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

std::string IceCoreSite::toString()
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

}
