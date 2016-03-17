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
#include "logging.h"

namespace Densification{ 

IdealizedCoreSite::IdealizedCoreSite(Settings& settings) : IceCoreSite(settings) {
    /* we take all parameters from settings */
    logger << "INFO: mean annual T = " << settings.T_mean << std::endl;
    logger << "INFO: mean annual Tamp = " << settings.Tamp << std::endl;
    logger << "INFO: mean annual acc = " << settings.acc_mean << std::endl;
    logger << "INFO: mean annual v10m = " << settings.v10m_mean<< std::endl;

    if (settings.T_mean < 0.) {
        logger << "ERROR: invalid value for T_mean = " << settings.T_mean << std::endl;
        std::abort();
    }
    if (settings.Tamp < 0.) {
        logger << "ERROR: invalid value for Tamp = " << settings.Tamp << std::endl;
        std::abort();
    }
    if (settings.acc_mean < 0.) {
        logger << "ERROR: invalid value for acc_mean = " << settings.acc_mean << std::endl;
        std::abort();
    }
    if (settings.v10m_mean < 0.) {
        logger << "ERROR: invalid value for v10m_mean = " << settings.v10m_mean << std::endl;
        std::abort();
    }
}

double IdealizedCoreSite::surfaceDensity(long time) {
    /* 
    Surface density after Helsen (2008) 
    time in seconds since start of the simulations
    */
    return -154.91+1.4266*(73.6+1.06*settings.T_mean+0.0669*settings.acc_mean+4.77*settings.v10m_mean);
    //return settings.rho_s;
}

double IdealizedCoreSite::surfaceTemperature(long time) {
    /* """ Surface temperature (varies over time, cosine function) """ */
    double period  = sec_in_year; // period of a year in sec
    double T       = settings.Tamp*cos(2*M_PI*(double)time/period) + settings.T_mean;
    return T;
}

double IdealizedCoreSite::accumulationRate(long time) {
    return settings.acc_mean/sec_in_year;
}

double IdealizedCoreSite::annualIntegratedAccumulation() {
    return settings.acc_mean;
}

double IdealizedCoreSite::annualMeanSurfaceTemperature() {
    return settings.T_mean;
}

std::string IdealizedCoreSite::toString() {
    std::ostringstream s; 
    s << "Ideal_" << super::toString();
    return s.str();
}

} // namespace

