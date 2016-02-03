#ifndef IDEALIZEDCORESITE_H
#define IDEALIZEDCORESITE_H

#include "icecoresite.h"

namespace Densification { 

class IdealizedCoreSite: public IceCoreSite {
    /// core with idealized forcing at the surface
private:
    typedef IceCoreSite super;

    const double acc=300.; // #Average Anual accumulation (mm/yr)
    const double v_10m=10.; //   ;#10m windspeed (m/s)
    const double Tsmean = 250.; //  # annual mean surface temperature (K)

public:
    IdealizedCoreSite() {};
    double surfaceDensity();
    double surfaceTemperature();
    double accumulationRate();
    double annualAccumulation();
    double annualSurfaceTemperature();
//    double annualSurfaceDensity();
    std::string toString();
};

}

#endif
