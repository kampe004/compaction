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

double IdealizedCoreSite::surfaceDensity() {
    /* after Helsen (2008) */
    double rho = -151.94+1.4266*(73.6+1.06*Tsmean+0.0669*acc+4.77*v_10m);
    rho = 300.;
    return rho;
}

double IdealizedCoreSite::surfaceTemperature() {
    /* """ Surface temperature (varies over time, cosine function) """ */
    double amp     = 10.0; // amplitude
    double period  = sec_in_year; // period of a year in sec
    double pi = 3.14;
    double T       = amp*cos(2*pi*(double)current_time/period) + Tsmean;
    return T;
}

double IdealizedCoreSite::accumulationRate() {
    return acc/sec_in_year;
}

double IdealizedCoreSite::annualAccumulation() {
    return acc;
}

double IdealizedCoreSite::annualSurfaceTemperature() {
    return Tsmean;
}

// double IdealizedCoreSite::annualSurfaceDensity();
// {
//     return 
// }

std::string IdealizedCoreSite::toString() {
    std::ostringstream s; 
    s << "Ideal_" << super::toString();
    return s.str();
}

} // namespace

