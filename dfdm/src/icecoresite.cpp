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

#include "icecoresite.h"
#include "logging.h"

namespace Densification{ 

IceCoreSite::IceCoreSite(Settings& settings) {
    densification_method = settings.dm;
    have_diffusion = settings.heat;
    eta0 = settings.eta0;
    c5 = settings.c5;
    c6 = settings.c6;
}

void IceCoreSite::init(){
    // initialize grid with a very tiny top layer
    Layer layer;
    layer.T       = this->surfaceTemperature();
    layer.dens    = this->surfaceDensity();
    layer.dz      = 0.001;
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

    // check minimum thickness
    int i = grid.size()-2;
    while (i >= 0 ) { // Use while instead of for-loop as size may change during loop
        if (grid[i].dz < dzmin && gridsize() > 1){
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
    while(grid.front().dens > 900.) 
        grid.erase(grid.begin());

    if (have_diffusion && gridsize() > 1)
        heatDiffusion(dt);

    // update internal time variable
    current_time = current_time + dt;
}

void IceCoreSite::accumulate(long dt) {
    double dens = this->surfaceDensity();
    double acc = this->accumulationRate() * 1e-3 * dt;
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
        grid.push_back(layer);

        top.dz = top.dz * 4/5; 
        top.mass = top.mass * 4/5;
    }
}

void IceCoreSite::compact(long dt) {
    if (densification_method == DensificationMethod::Ar10T) {
        compactAr10T(dt);
    } else if (densification_method == DensificationMethod::Overburden) {
        compactOverburden(dt);
    } else {
        logger << "ERROR: unknown densification method" << std::endl;
        std::abort();
    }
}

void IceCoreSite::compactAr10T(long dt) {
    /*  Densification as equation [Ar10T] from Ligtenberg2011 */
    static const double E_c=60.e3; //#  (kJ/mol)
    static const double E_g=42.4e3; //# (kJ/mol)  
    double acc_year = annualAccumulation(); 
    double Tsmean = annualSurfaceTemperature();
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

void IceCoreSite::compactOverburden(long dt) {
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
        double eta = eta0 * exp(c5*(Tf-grid[i].T) + c6*grid[i].dens); //# viscocity, assuming no liquid
        double cr2 = -P/eta;

        double cr = cr1 + cr2;
        grid[i].dz = grid[i].dz * (1+(cr*dt));
        grid[i].dens = (grid[i].mass/grid[i].dz); // # mass conservation

        overburden = overburden + grid[i].mass;
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
    T_new[Np-1] = this->surfaceTemperature(); // # top boundary equal to surface temperature

    // Neumann boundary condition at the bottom
    T_new[0] = td[0]*dt/grid[0].dz * (grid[1].T-grid[0].T) + grid[0].T;

    for (int i = 1; i < gridsize()-1; i++) {
        double stab_crit = td[i] * dt / pow(grid[i].dz,2);
        if (stab_crit > 0.5) {
            logger << "ERROR: stability criterium violated, r = " << stab_crit << std::endl;
            std::abort();
        } else {
            T_new[i] = td[i]*dt/pow(grid[i].dz,2)*(grid[i-1].T-2*grid[i].T+grid[i+1].T) + grid[i].T; // # Forward in Time Centered in Space
        }
    }

    // Update temperatures
    for (int i = 0; i < gridsize(); i++) {
        grid[i].T = T_new[i];
    }
}

double IceCoreSite::getDepthOfDensity(double dens) {
    double tot_depth = 0.0;
    for (int i = gridsize()-1; i >= 0; i--) {
        tot_depth = tot_depth + grid[i].dz;
        if (grid[i].dens >= dens ){
            // TODO: do some kind of interpolation here
            return tot_depth;
        }
    }
    return -1.0;
}

double IceCoreSite::getZ550() {
    return getDepthOfDensity(550.);
}

double IceCoreSite::getZ830() {
    return getDepthOfDensity(830.);
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
    if (densification_method == DensificationMethod::Ar10T) {
        s << "Ar10T"; 
    } else if (densification_method == DensificationMethod::Overburden) {
        s << "Overburden" << "_" << eta0 << "_" << c5 << "_" << c6;
    } else {
        s << "UNKNOWN"; 
    }
    if (have_diffusion) {
        s << "_heatTrue";
    } else {
        s << "_heatFalse";
    }
    return s.str();
}

}
