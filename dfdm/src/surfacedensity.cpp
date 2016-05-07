
#include "surfacedensity.h"
#include "config.h"
#include "meteo.h"
#include "logging.h"

namespace DSM{ 

std::unique_ptr<SurfaceDensity> instantiate_surfacedensity(Meteo& meteo){
   const char * option_name = "fresh_snow_density:which_fsd";
   int which_fsd = config.getInt(option_name, false, 0, 1, 1);

   switch (which_fsd) {
      case 0   : return { std::make_unique<SurfaceDensityConstant>(meteo) };
      case 1   : return { std::make_unique<SurfaceDensityHelsen2008>(meteo) };
      case 2   : return { std::make_unique<SurfaceDensityLenaerts2012>(meteo) };
      case 3   : return { std::make_unique<SurfaceDensityCROCUS>(meteo) };
      default:
         logger << "ERROR: unknown value: " << which_fsd << " for config option " << option_name << std::endl;
         std::abort();
   }
}

SurfaceDensity::SurfaceDensity(Meteo& meteo) : _meteo(meteo) { } 

SurfaceDensityConstant::SurfaceDensityConstant(Meteo& meteo) : SurfaceDensity(meteo) { } 

SurfaceDensityHelsen2008::SurfaceDensityHelsen2008(Meteo& meteo) : SurfaceDensity(meteo) { }

SurfaceDensityLenaerts2012::SurfaceDensityLenaerts2012(Meteo& meteo) : SurfaceDensity(meteo) { }

SurfaceDensityCROCUS::SurfaceDensityCROCUS(Meteo& meteo) : SurfaceDensity(meteo) { }

double SurfaceDensityConstant::density() {
   return _val;
}

double SurfaceDensityHelsen2008::density() {
   /* after Helsen (2008) */
   return -154.91+1.4266*(73.6+1.06*_meteo.annualTskin()+0.0669*_meteo.annualAcc()+4.77*_meteo.annualW10m());
}

double SurfaceDensityLenaerts2012::density() {
   /* Lenaerts et al. 2012 
      formula 11
      The multiple linear regression that relates fresh snow density to mean surface temperature (Tsfc,Acc) and 10 m wind speed (U10m,Acc) during accumulation
      */
   double a = 97.5;
   double b = 0.77; 
   double c = 4.49;
   return a + b*_meteo.surfaceTemperature() + c*_meteo.surfaceWind();
}

double SurfaceDensityCROCUS::density() {
   return 1.0;
}

} // namespace
