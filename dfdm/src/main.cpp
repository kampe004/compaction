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
// need unique instance of logger in namespace
std::ofstream logger; 

// need unique implementation of DensificationMethod << operator in namespace
std::ostream& operator<<(std::ostream& os, DensificationMethod dm) {
   switch (dm) {
      case DensificationMethod::Ligtenberg2011   : os << "Ligtenberg"; break;
      case DensificationMethod::Anderson1976     : os << "Anderson"; break;
      case DensificationMethod::Barnola1991      : os << "Barnola"; break;
      case DensificationMethod::Spencer2001      : os << "Spencer"; break;
      case DensificationMethod::BarnolaSpencer   : os << "BarnolaSpencer"; break;
      case DensificationMethod::HerronLangway    : os << "HerronLangway"; break;
      default                            : os.setstate(std::ios_base::failbit);
   }
   return os;
}

std::ostream& operator<<(std::ostream& os, ForcingMethod dm) {
   switch (dm) {
      case ForcingMethod::Idealized       : os << "Ideal"; break;
      case ForcingMethod::Netcdf         : os << "NetCDF"; break;
      default                        : os.setstate(std::ios_base::failbit);
   }
   return os;
}
}

using namespace Densification; 

void readSettingsFromIniFile(char ininame[], Settings& s){
   dictionary* d = iniparser_load(ininame);

   /****** GENERAL ******/
   const char * tmp = iniparser_getstring(d, "general:compaction", "NOT_SPECIFIED");
   if (std::string(tmp) == "Ligtenberg2011"){ // enum class DensificationMethod {Ligtenberg2011, Anderson1976, Barnola1991, Spencer2001, BarnolaSpencer}
      s.dm = DensificationMethod::Ligtenberg2011;
   } else if (std::string(tmp) == "Anderson1976"){
      s.dm = DensificationMethod::Anderson1976;
   } else if (std::string(tmp) == "Barnola1991"){
      s.dm = DensificationMethod::Barnola1991;
   } else if (std::string(tmp) == "Spencer2001"){
      s.dm = DensificationMethod::Spencer2001;
   } else if (std::string(tmp) == "BarnolaSpencer"){
      s.dm = DensificationMethod::BarnolaSpencer;
   } else if (std::string(tmp) == "HerronLangway"){
      s.dm = DensificationMethod::HerronLangway;
   } else {
      logger << "ERROR: unknown general:compaction' value: '" << std::string(tmp) << "'" << std::endl;
      std::abort();
   }

   tmp = iniparser_getstring(d, "general:forcing", "NOT_SPECIFIED");
   if (std::string(tmp) == "netcdf"){ // enum class ForcingMethod {Idealized, Netcdf};
      s.fm = ForcingMethod::Netcdf;
   } else if (std::string(tmp) == "idealized"){
      s.fm = ForcingMethod::Idealized;
   } else {
      logger << "ERROR: unknown 'general:forcing' value: '" << std::string(tmp) << "'" << std::endl;
      logger << "Valid values are 'netcdf' and 'idealized'" << std::endl;
      std::abort();
   }

   s.have_diffusion = iniparser_getboolean(d, "general:heat", -1);
   s.max_depth = iniparser_getdouble(d, "general:max_depth", 200.0);
   s.max_year = iniparser_getint(d, "general:max_year", 4000);

   logger << "INFO: general:compaction = " << tmp << std::endl;
   logger << "INFO: general:heat = " << s.have_diffusion << std::endl;
   logger << "INFO: general:max_depth = " << s.max_depth << std::endl;
   logger << "INFO: general:max_year = " << s.max_year << std::endl;

   /****** DENSITY ******/
   s.rho_s = iniparser_getdouble(d, "density:rho_s", 100.0);
   logger << "INFO: density:rho_s = " << s.rho_s << std::endl;
   
   /****** OVERBURDEN ******/
   s.eta0 = iniparser_getdouble(d, "overburden:eta0", -1.0);
   s.c5 = iniparser_getdouble(d, "overburden:c5", -1.0);
   s.c6 = iniparser_getdouble(d, "overburden:c6", -1.0);
   logger << "INFO: overburden:eta0 = " << s.eta0 << std::endl;
   logger << "INFO: overburden:c5 = " << s.c5 << std::endl;
   logger << "INFO: overburden:c6 = " << s.c6 << std::endl;

   if (s.fm == ForcingMethod::Netcdf) {
      /****** NETCDF FORCING ******/
      s.forcing_dt = iniparser_getint(d, "forcing:dt", -1);
      tmp = iniparser_getstring(d, "forcing:f_acc", "NOT_SPECIFIED");
      s.f_acc = std::string(tmp);
      tmp = iniparser_getstring(d, "forcing:f_wind10m", "NOT_SPECIFIED");
      s.f_wind10m = std::string(tmp);
      tmp = iniparser_getstring(d, "forcing:f_tskin", "NOT_SPECIFIED");
      s.f_tskin = std::string(tmp);
      logger << "INFO: using Netcdf forcing" << std::endl;
      logger << "INFO: forcing:dt = " << s.forcing_dt << std::endl;
      logger << "INFO: forcing:f_acc = " << s.f_acc << std::endl;
      logger << "INFO: forcing:f_w10m = " << s.f_wind10m << std::endl;
      logger << "INFO: forcing:f_tskin = " << s.f_tskin << std::endl;
   } else if (s.fm == ForcingMethod::Idealized) {
      /****** IDEALIZED FORCING ******/
      s.T_mean = iniparser_getdouble(d, "forcing:Tmean", -1.0);
      s.Tamp = iniparser_getdouble(d, "forcing:Tamp", 10.);
      s.acc_mean = iniparser_getdouble(d, "forcing:acc", -1.0);
      s.v10m_mean = iniparser_getdouble(d, "forcing:v10m", -1.0);

      logger << "INFO: using Netcdf forcing" << std::endl;
      logger << "INFO: forcing:Tmean = " << s.T_mean << std::endl;
      logger << "INFO: forcing:Tamp = " << s.Tamp << std::endl;
      logger << "INFO: forcing:acc = " << s.acc_mean << std::endl;
      logger << "INFO: forcing:v10m = " << s.v10m_mean << std::endl;
   }

   iniparser_freedict(d);
}

int main(){
   logger.open("densification.log");
   logger.precision(16);
   char ininame[] = "settings.ini";
   Settings settings;
   readSettingsFromIniFile(ininame, settings);

   IceCoreSite *core;

   if (settings.fm == ForcingMethod::Netcdf) {
      core = new NetcdfCoreSite(settings);
   } else if (settings.fm == ForcingMethod::Idealized) {
      core = new IdealizedCoreSite(settings);
   } else {
      logger << "ERROR: value not implemented, ForcingMethod = " << settings.fm << std::endl;
      std::abort();
   }
   core->init();

   /* Time loop */
   long sec_in_3h = 3*3600;
   long dt = sec_in_3h;  // timestep of 3 hours in seconds
   long dt_per_year = sec_in_year/dt;

   logger << "INFO: dt=" << dt << std::endl;
//   logger << "INFO: nyears=" << nyears << std::endl;
   logger << "INFO: dt_per_year=" << dt_per_year << std::endl;
   logger << "INFO: starting simulation of " << core->toString() << std::endl;

   std::clock_t start;
   int kYear = 0;
   while(! core->hasReachedDensity(850.) && kYear < settings.max_year && core->totalDepth() < settings.max_depth) {
//   for (int year = 0; year < nyears; year++){
      start = clock();
      for(int tstep = 0; tstep < dt_per_year; tstep++)
         core->runTimeStep(dt);
      double elapsed = ((double)(clock() - start)) / CLOCKS_PER_SEC;
      logger  << "year=" << kYear
            << ", Tmin=" << core->minTemp() 
            << ", Tmax=" << core->maxTemp() 
            << ", rho_max=" << core->maxDens()
            << ", sec/year=" << elapsed 
            << ", year/hour=" << 3600./elapsed << std::endl;
      kYear++; 
   }
   core->printIceCoreSummary();
   core->writeFiles();
   logger.close();
}

