#ifndef NETCDFCORESITE_H
#define NETCDFCORESITE_H

#include <vector>
#include "icecoresite.h"

namespace Densification { 

class NetcdfCoreSite: public IceCoreSite {
    /// core with meteorological forcing from NetCDF files
private:
    typedef IceCoreSite super;

    long forcing_dt; // sampling interval in seconds in forcing file

    // mean annual statistics for use in Helsen2008
    double acc_ann_mean;
    double w10m_ann_mean;
    double Ts_ann_mean;

    // raw data from NetCDF files
    std::vector<double> acc_all;
    std::vector<double> tskin_all;
    std::vector<double> w10m_all;

    void readForcing(Settings& settings);

protected:
    double surfaceDensity(long time);
    double surfaceTemperature(long time);
    double accumulationRate(long time);
    double annualIntegratedAccumulation();
    double annualMeanSurfaceTemperature();
//    double annualSurfaceDensity();

public:
    NetcdfCoreSite();
    NetcdfCoreSite(Settings& settings);
    std::string toString();
};

}

#endif
