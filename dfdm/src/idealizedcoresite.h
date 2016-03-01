#ifndef IDEALIZEDCORESITE_H
#define IDEALIZEDCORESITE_H

#include "icecoresite.h"

namespace Densification { 

class IdealizedCoreSite: public IceCoreSite {
    /// core with idealized meteorological forcing

public:
    IdealizedCoreSite(Settings& settings) : IceCoreSite(settings) {};
    std::string toString();

protected:
    double surfaceDensity(long time);
    double surfaceTemperature(long time);
    double accumulationRate(long time);
    double annualIntegratedAccumulation();
    double annualMeanSurfaceTemperature();
//    double annualSurfaceDensity();

private:
    typedef IceCoreSite super;

    const double acc=34.; // #Average Anual accumulation (mm/yr)
    const double v_10m=9.; //   ;#10m windspeed (m/s)
    const double Tsmean = 273.15 - 53; //  # annual mean surface temperature (K)
};

}

#endif
