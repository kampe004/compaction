/* --------------------------------------------------------------------------------
   
   DESCRIPTION
      simple firn model with compaction and heat diffusion
   
   DATE
      15-JAN-2016
   
   AUTHOR
      L.vankampenhout@uu.nl
   -------------------------------------------------------------------------------- */

#include <iostream>
#include <ctime>

#include "idealizedcoresite.h"

using namespace Densification; 

int main(){
    /*  Idealized */
    IdealizedCoreSite core;

    /* eta = Overburden.eta0 * math.exp(Overburden.c5*(Overburden.Tf-T) + Overburden.c6*dens) # viscocity, assuming no liquid */
    core.c5 = 0.16;

    core.init(DensificationMethod::Ar10T, true);

    /* Time loop */
    long sec_in_3h = 3*3600;
    long dt = sec_in_3h;  // timestep of 3 hours in seconds
    long dt_per_year = ((long)sec_in_year)/dt;
    long nyears = 50;  // simulation duration in years
    long Nt = (nyears*365*24*60*60)/dt; // number of points on time axis

    printf("DEBUG Nt=%ld \n",Nt);
    printf("DEBUG dt_per_year=%ld \n",dt_per_year);

    std::cout << "INFO: starting simulation with nyears = " << nyears << std::endl;
    std::clock_t start;
    for (int year = 0; year < nyears; year++){
        std::cout << "INFO year = " << year << std::endl;
        start = clock();
        for(int tstep = 0; tstep < dt_per_year; tstep++){
            core.runTimeStep(dt);
        }

//         """ print some yearly diagnostics """ 
        std::cout << "core = " << core.toString() << ", Tmin = " << core.minTemp() << ", Tmax = " << core.maxTemp() 
            << ", rho_max = " << core.maxDens() << std::endl;

        double elapsed = ((double)(clock() - start)) / CLOCKS_PER_SEC;
        std::cout << "INFO: sec/year = " << elapsed << ", year/hour = " << 3600./elapsed << std::endl;
    }

    core.printIceCoreSummary();
    core.writeFiles();
}

