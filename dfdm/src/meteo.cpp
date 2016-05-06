#include <iostream>
#include <sstream>
#include <cmath>
#include <netcdf>
#include <iomanip>
#include <numeric>
#include <string>

#include "config.h"
#include "logging.h"
#include "constants.h"
#include "meteo.h"
#include "dynamicmodel.h"

namespace DSM{ 

std::unique_ptr<Meteo> instantiate_meteo(DynamicModel& dm){
   const char * option_name = "forcing:which_forcing";
   int which_forcing = config.getInt(option_name, true, 0, 1, 0);

   switch (which_forcing) {
      case 0   : return { std::make_unique<MeteoIdealized>(dm) };
      case 1   : return { std::make_unique<MeteoNetcdf>(dm) };
      default:
         logger << "ERROR: unknown value: " << which_forcing << " for config option " << option_name << std::endl;
         std::abort();
   }
}


/* class Meteo*/
Meteo::Meteo(DynamicModel& dm) : _dm(dm) {
   logger << "Meteo()" << std::endl;
}

/* class MeteoIdealized*/
MeteoIdealized::MeteoIdealized(DynamicModel& dm) : Meteo(dm){
   logger << "MeteoIdealized()" << std::endl;
   logger << "INFO: using idealized forcing" << std::endl;

   const char * option_name;
   option_name = "forcing:ideal_T_mean";
   _ideal_T_mean = config.getDouble(option_name, true, 0, 1000, 0);
   option_name = "forcing:ideal_T_amp";
   _ideal_T_amp = config.getDouble(option_name, true, 0, 1000, 0);
   option_name = "forcing:ideal_acc";
   _ideal_acc = config.getDouble(option_name, true, 0, 1000, 0);
   option_name = "forcing:ideal_w10m";
   _ideal_w10m = config.getDouble(option_name, true, 0, 1000, 0);
}

/* class MeteoNetcdf */
MeteoNetcdf::MeteoNetcdf(DynamicModel& dm) : Meteo(dm){
   logger << "MeteoNetcdf()" << std::endl;
   logger << "INFO: using Netcdf forcing" << std::endl;

   const char * option_name = "forcing:netcdf_dt";
   _dt_forcing = config.getInt(option_name, true, 0, (int)1e9, 0);
   if (_dm.getDt() != _dt_forcing) {
      logger << "ERROR: " << option_name << " [ " << _dt_forcing << " ] must equal model dt [ " << _dm.getDt() << " ]" << std::endl;
      std::abort();
   }
   logger << "INFO: " << option_name << " = " << _dt_forcing << std::endl;
   readForcing();
}


void MeteoNetcdf::readForcing(){
   /* Open NetCDF forcing files and read them into memory entirely 
      currently it is assumed that the files describe a single point, 
      thus they are dependent on time only. */

   std::string f_acc = config.getFilename("forcing:netcdf_acc", true, true, "");
   std::string f_w10m = config.getFilename("forcing:netcdf_w10m", true, true, "");
   std::string f_tskin = config.getFilename("forcing:netcdf_tskin", true, true, "");

   logger << "INFO: opening " << f_acc << std::endl;
   logger << "INFO: opening " << f_w10m << std::endl;
   logger << "INFO: opening " << f_tskin << std::endl;

   netCDF::NcFile nc_file_acc(f_acc, netCDF::NcFile::read);
   netCDF::NcFile nc_file_w10m(f_w10m, netCDF::NcFile::read);
   netCDF::NcFile nc_file_tskin(f_tskin, netCDF::NcFile::read);

   netCDF::NcDim time1 = nc_file_acc.getDim("time");
   netCDF::NcDim time2 = nc_file_w10m.getDim("time");
   netCDF::NcDim time3 = nc_file_tskin.getDim("time");

   unsigned nt = time1.getSize();
   logger << "INFO: number of timesteps in forcing file " << nt << std::endl;
   if (nt != time2.getSize() || nt != time3.getSize()) {
      logger << "ERROR: forcing files have different number of timesteps between them" << std::endl;
      logger << f_acc << " has nt = " << time1.getSize() << std::endl;
      logger << f_w10m << " has nt = " << time2.getSize() << std::endl;
      logger << f_tskin << " has nt = " << time3.getSize() << std::endl;
      std::abort();
   }

   _acc_all.resize(nt); 
   _tskin_all.resize(nt); 
   _w10m_all.resize(nt); 

   netCDF::NcVar var1 = nc_file_acc.getVar("acc");
   netCDF::NcVar var2 = nc_file_w10m.getVar("w10m");
   netCDF::NcVar var3 = nc_file_tskin.getVar("tskin");
   
   var1.getVar(&_acc_all.front());
   var2.getVar(&_w10m_all.front());
   var3.getVar(&_tskin_all.front());
   logger << "INFO: read " << _acc_all.size() << " timesteps from NetCDF files" << std::endl;

   long total_sec = (long)_dt_forcing * (long)(nt);
   double nyears = (double)total_sec / (double)sec_in_year;
   logger << "INFO: forcing file has " << nyears << " year(s) of data" << std::endl; 

   double acc_total = std::accumulate(_acc_all.begin(), _acc_all.end(), 0.0);
   double w10m_total = std::accumulate(_w10m_all.begin(), _w10m_all.end(), 0.0);
   double tskin_total = std::accumulate(_tskin_all.begin(), _tskin_all.end(), 0.0);

   _acc_ann_mean = acc_total / nyears;
   _w10m_ann_mean = w10m_total / nt; 
   _Ts_ann_mean = tskin_total / nt; 
   logger << "INFO: mean annual acc = "   << _acc_ann_mean << std::endl;
   logger << "INFO: mean annual w10m = "  << _w10m_ann_mean << std::endl;
   logger << "INFO: mean annual Ts = "    << _Ts_ann_mean << std::endl;
}


double MeteoNetcdf::surfaceTemperature() {
   /* return surface temperature for current time step
      currently assumes that dt = dt_forcing (this checked in the constructor) 
      */
   long nt = _dm.getCurrentTimeStep();
   int idx = (int) nt; // index in forcing array
   idx = idx % _tskin_all.size();
   return _tskin_all[idx];
}

double MeteoNetcdf::accumulationRate() {
   /* return acccumulation rate for current time step
      currently assumes that dt = dt_forcing (this checked in the constructor) 
      */
   long nt = _dm.getCurrentTimeStep();
   int idx = (int) nt;
   idx = idx % _acc_all.size();
   return _acc_all[idx] / _dt_forcing;
}

} // namespace

