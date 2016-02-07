/* --------------------------------------------------------------------------------
   
   DESCRIPTION
      Idealized forcing
   
   DATE
      15-JAN-2016
   
   AUTHOR
      L.vankampenhout@uu.nl
   -------------------------------------------------------------------------------- */

#include <iostream>
#include <sstream>
#include <cmath>

#include "idealizedcoresite.h"

namespace Densification{ 

double IdealizedCoreSite::surfaceDensity(long time) {
    /* 
    Surface density after Helsen (2008) 
    
    time in seconds since start of the simulations
    */
    return -151.94+1.4266*(73.6+1.06*Tsmean+0.0669*acc+4.77*v_10m);
}

double IdealizedCoreSite::surfaceTemperature(long time) {
    /* """ Surface temperature (varies over time, cosine function) """ */
    double amp     = 10.0; // amplitude
    double period  = sec_in_year; // period of a year in sec
    double pi = 3.14;
    double T       = amp*cos(2*pi*(double)time/period) + Tsmean;
    return T;
}

double IdealizedCoreSite::accumulationRate(long time) {
    return acc/sec_in_year;
}

double IdealizedCoreSite::annualIntegratedAccumulation() {
    return acc;
}

double IdealizedCoreSite::annualMeanSurfaceTemperature() {
    return Tsmean;
}

std::string IdealizedCoreSite::toString() {
    std::ostringstream s; 
    s << "Ideal_" << super::toString();
    return s.str();
}

} // namespace

