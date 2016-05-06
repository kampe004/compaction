#ifndef NETCDFCORESITE_H
#define NETCDFCORESITE_H

#include <vector>
#include "icecoresite.h"

namespace DSM{ 

class NetcdfCoreSite: public IceCoreSite {
   /* meteorological forcing from NetCDF */

public:
   NetcdfCoreSite(Settings& settings);
   std::string toString();

protected:
   double surfaceTemperature(long time);
   double accumulationRate(long time);

private:
   typedef IceCoreSite super;
   long forcing_dt; // sampling interval in seconds in forcing file

   // raw data from NetCDF files
   std::vector<double> acc_all;
   std::vector<double> tskin_all;
   std::vector<double> w10m_all;

   void readForcing(Settings& settings);
};
}

#endif
