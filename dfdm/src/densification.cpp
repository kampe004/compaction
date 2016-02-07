/* --------------------------------------------------------------------------------
   
   DESCRIPTION
      simple firn model with overburden compaction and heat diffusion
   
   DATE
      15-JAN-2016
   
   AUTHOR
      L.vankampenhout@uu.nl
   -------------------------------------------------------------------------------- */

#include <cstdlib>
#include <iostream>
#include <ctime>

#include "logging.h"
#include "settings.h"
#include "iniparser.h"
#include "idealizedcoresite.h"
#include "netcdfcoresite.h"

namespace Densification{
std::ofstream logger;
}

using namespace Densification; 

void readSettingsFromIniFile(char ininame[], Settings& s){
    dictionary* d = iniparser_load(ininame);

    s.heat = iniparser_getboolean(d, "general:heat", -1);
    const char * tmp = iniparser_getstring(d, "general:compaction", "NOT_SPECIFIED");
    if (std::string(tmp) == "overburden"){
        s.dm = DensificationMethod::Overburden;
    } else {
        s.dm = DensificationMethod::Ar10T;
    }
    s.eta0 = iniparser_getdouble(d, "overburden:eta0", -1.0);
    s.c5 = iniparser_getdouble(d, "overburden:c5", -1.0);
    s.c6 = iniparser_getdouble(d, "overburden:c6", -1.0);

    s.forcing_dt = iniparser_getint(d, "forcing:dt", -1);
    tmp = iniparser_getstring(d, "forcing:f_acc", "NOT_SPECIFIED");
    s.f_acc = std::string(tmp);
    tmp = iniparser_getstring(d, "forcing:f_wind10m", "NOT_SPECIFIED");
    s.f_wind10m = std::string(tmp);
    tmp = iniparser_getstring(d, "forcing:f_tskin", "NOT_SPECIFIED");
    s.f_tskin = std::string(tmp);

    logger << "INFO: general:compaction = " << tmp << std::endl;
    logger << "INFO: general:heat = " << s.heat << std::endl;
    logger << "INFO: overburden:eta0 = " << s.eta0 << std::endl;
    logger << "INFO: overburden:c5 = " << s.c5 << std::endl;
    logger << "INFO: overburden:c6 = " << s.c6 << std::endl;
    logger << "INFO: forcing:dt = " << s.forcing_dt << std::endl;
    logger << "INFO: forcing:f_acc = " << s.f_acc << std::endl;
    logger << "INFO: forcing:f_w10m = " << s.f_wind10m << std::endl;
    logger << "INFO: forcing:f_tskin = " << s.f_tskin << std::endl;


    iniparser_freedict(d);
}

int main(){
    logger.open("densification.log");
    char ininame[] = "settings.ini";
    Settings settings;
    readSettingsFromIniFile(ininame, settings);

    //IdealizedCoreSite core(settings);
    NetcdfCoreSite core(settings);
    core.init();

    /* Time loop */
    long sec_in_3h = 3*3600;
    long dt = sec_in_3h;  // timestep of 3 hours in seconds
    long dt_per_year = sec_in_year/dt;
    long nyears = 3;  // simulation duration in years

    logger << "INFO: dt=" << dt << std::endl;
    logger << "INFO: nyears=" << nyears << std::endl;
    logger << "INFO: dt_per_year=" << dt_per_year << std::endl;
    logger << "INFO: starting simulation of " << core.toString() << std::endl;

    std::clock_t start;
    for (int year = 0; year < nyears; year++){
        start = clock();
        for(int tstep = 0; tstep < dt_per_year; tstep++)
            core.runTimeStep(dt);
        double elapsed = ((double)(clock() - start)) / CLOCKS_PER_SEC;
        logger  << "year=" << year
                << ", Tmin=" << core.minTemp() 
                << ", Tmax=" << core.maxTemp() 
                << ", rho_max=" << core.maxDens()
                << ", sec/year=" << elapsed 
                << ", year/hour=" << 3600./elapsed << std::endl;
    }
    core.printIceCoreSummary();
    core.writeFiles();
    logger.close();
}

