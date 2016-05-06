#ifndef IDEALIZEDCORESITE_H
#define IDEALIZEDCORESITE_H

#include "icecoresite.h"

namespace DSM{ 

class IdealizedCoreSite: public IceCoreSite {
   /// core with idealized meteorological forcing

public:
   IdealizedCoreSite(Settings& settings);
   std::string toString();

protected:
   double surfaceTemperature(long time);
   double accumulationRate(long time);
   double annualIntegratedAccumulation();
   double annualMeanSurfaceTemperature();
//   double annualSurfaceDensity();

private:
   typedef IceCoreSite super;
};

}

#endif
