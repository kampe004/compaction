/* --------------------------------------------------------------------------------
   
   DESCRIPTION
     Ice core which reads meteorological forcing from NetCDF files

     Requires NetCDF-4 C++ library to be installed: 
      https://www.unidata.ucar.edu/software/netcdf/docs/cxx4/index.html
     Macports package 'netcdf-cxx4'
   
   DATE
     5-FEB-2016
   
   AUTHOR
     L.vankampenhout@uu.nl
   -------------------------------------------------------------------------------- */

#include <iostream>
#include <sstream>
#include <cmath>
#include <netcdf>
#include <iomanip>
#include <numeric>

#include "settings.h"
#include "logging.h"
#include "netcdfcoresite.h"

namespace Densification{ 

NetcdfCoreSite::NetcdfCoreSite(Settings& settings) : IceCoreSite(settings) {
   readForcing(settings);
}

void NetcdfCoreSite::readForcing(Settings& settings){
   /* find forcing files and read them into memory entirely */

   if (settings.forcing_dt < 0) {
      logger << "ERROR: forcing timestep dt not set!" << std::endl;
      std::abort();
   } else {
      forcing_dt = settings.forcing_dt;
   }

   logger << "INFO: opening " << settings.f_acc << std::endl;
   logger << "INFO: opening " << settings.f_wind10m << std::endl;
   logger << "INFO: opening " << settings.f_tskin << std::endl;

   netCDF::NcFile nc_file_acc(settings.f_acc, netCDF::NcFile::read);
   netCDF::NcFile nc_file_w10m(settings.f_wind10m, netCDF::NcFile::read);
   netCDF::NcFile nc_file_tskin(settings.f_tskin, netCDF::NcFile::read);

   netCDF::NcDim time  = nc_file_acc.getDim("time");
   netCDF::NcDim time2 = nc_file_w10m.getDim("time");
   netCDF::NcDim time3 = nc_file_tskin.getDim("time");

   int nt = time.getSize();
   logger << "INFO: number of timesteps in forcing file " << nt << std::endl;
   if (nt != time2.getSize() || nt != time3.getSize()) {
      logger << "ERROR: forcing files have different number of timesteps!\n";
      std::abort();
   }

   acc_all.resize(nt); 
   tskin_all.resize(nt); 
   w10m_all.resize(nt); 

   // TODO: dynamical variable names from INI file
   netCDF::NcVar var1 = nc_file_acc.getVar("acc");
   netCDF::NcVar var2 = nc_file_w10m.getVar("w10m");
   netCDF::NcVar var3 = nc_file_tskin.getVar("tskin");
   
   var1.getVar(&acc_all.front());
   var2.getVar(&w10m_all.front());
   var3.getVar(&tskin_all.front());

   logger << "INFO: succesfully read " << acc_all.size() << " timesteps from NetCDF files" << std::endl;

   long total_sec = (long)settings.forcing_dt * (long)(nt);
   double nyears = (double)total_sec / (double)sec_in_year;

   logger << "INFO: forcing file has " << nyears << " year(s) of data" << std::endl; 

   double acc_total = std::accumulate(acc_all.begin(), acc_all.end(), 0.0);
   double w10m_total = std::accumulate(w10m_all.begin(), w10m_all.end(), 0.0);
   double tskin_total = std::accumulate(tskin_all.begin(), tskin_all.end(), 0.0);

   acc_ann_mean = acc_total / nyears;
   w10m_ann_mean = w10m_total / nt; 
   Ts_ann_mean = tskin_total / nt; 

   logger << "INFO: mean annual acc = " << acc_ann_mean << std::endl;
   logger << "INFO: mean annual w10m = " << w10m_ann_mean << std::endl;
   logger << "INFO: mean annual Ts = " << Ts_ann_mean << std::endl;

   // TODO: do interpolation to finer timegrid here, in advance? 
}

double NetcdfCoreSite::surfaceDensity(long time) {
   /* after Helsen (2008) */
   return -154.91+1.4266*(73.6+1.06*Ts_ann_mean+0.0669*acc_ann_mean+4.77*w10m_ann_mean);
   //return settings.rho_s;
}

double NetcdfCoreSite::surfaceTemperature(long time) {
   int idx = time / forcing_dt; // index in forcing array
   idx = idx % tskin_all.size();
   return tskin_all[idx];
}

double NetcdfCoreSite::accumulationRate(long time) {
   int idx = time / forcing_dt; // index in forcing array
   idx = idx % acc_all.size();
   return acc_all[idx] / forcing_dt;
}

double NetcdfCoreSite::annualIntegratedAccumulation() {
   return acc_ann_mean;
}

double NetcdfCoreSite::annualMeanSurfaceTemperature() {
   return Ts_ann_mean;
}

std::string NetcdfCoreSite::toString() {
   std::ostringstream s; 
   s << "Netcdf_" << super::toString();
   return s.str();
}

} // namespace

