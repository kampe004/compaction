#ifndef IDEALIZEDCORESITE_H
#define IDEALIZEDCORESITE_H

#include "icecoresite.h"

namespace Densification { 

class IdealizedCoreSite: public IceCoreSite {
    /// core with idealized meteorological forcing
private:
    typedef IceCoreSite super;

    const double acc=300.; // #Average Anual accumulation (mm/yr)
    const double v_10m=10.; //   ;#10m windspeed (m/s)
    const double Tsmean = 250.; //  # annual mean surface temperature (K)

protected:
    double surfaceDensity(long time);
    double surfaceTemperature(long time);
    double accumulationRate(long time);
    double annualIntegratedAccumulation();
    double annualMeanSurfaceTemperature();
//    double annualSurfaceDensity();

public:
    IdealizedCoreSite() {};
    IdealizedCoreSite(Settings& settings) : IceCoreSite(settings) {};
    std::string toString();
};

}

#endif
