/* --------------------------------------------------------------------------------
   
   DESCRIPTION
      simple firn model with compaction and heat diffusion
   
   DATE
      15-JAN-2016
   
   AUTHOR
      L.vankampenhout@uu.nl
   -------------------------------------------------------------------------------- */

#include <cstdlib>
#include <iostream>
#include <ctime>

#include "iniparser.h"
#include "idealizedcoresite.h"

using namespace Densification; 

/* define data container that holds all user settings */
typedef struct settings settings; 
struct settings {
    /* overburden parameters */
    double eta0;
    double c5;
    double c6;
};

void readSettingsFromIniFile(char ininame[], settings s){
    dictionary* d = iniparser_load(ininame);
    s.eta0 = iniparser_getdouble(d, "overburden:eta0", -1.0);
    if (s.eta0 < 0.0) {
        std::cout << "ERROR: overburden:eta0 not defined in INI file or is negative" << std::endl;
        std::abort();
    } else {
        std::cout << "INFO: using overburden:eta0 = " << s.eta0 << std::endl;
    }
    s.c5 = iniparser_getdouble(d, "overburden:c5", -1.0);
    if (s.c5 < 0.0) {
        std::cout << "ERROR: overburden:c5 not defined in INI file or is negative" << std::endl;
        std::abort();
    } else {
        std::cout << "INFO: using overburden:c5 = " << s.c5 << std::endl;
    }
    s.c6 = iniparser_getdouble(d, "overburden:c6", -1.0);
    if (s.c6 < 0.0) {
        std::cout << "ERROR: overburden:c6 not defined in INI file or is negative" << std::endl;
        std::abort();
    } else {
        std::cout << "INFO: using overburden:c6 = " << s.c6 << std::endl;
    }
    iniparser_freedict(d);
}

int main(){
    /*  parse config file */
    char ininame[] = "settings.ini";
    settings my_settings;
    readSettingsFromIniFile(ininame, my_settings);

    IdealizedCoreSite core;

    /* configure run settings */
    core.c5 = my_settings.c5;
    core.eta0 = my_settings.eta0;
    core.c6 = my_settings.c6;

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

