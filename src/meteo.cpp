#include <iostream>
#include <sstream>
#include <netcdf>
#include <iomanip>
#include <numeric>
#include <string>
#include <cmath>

#include "config.h"
#include "logging.h"
#include "constants.h"
#include "meteo.h"
#include "dynamicmodel.h"

namespace DSM{ 

std::unique_ptr<Meteo> instantiate_meteo(DynamicModel& dm){
   const char * option_name = "forcing:which_forcing";
   int which_forcing = config.getInt(option_name, true, 0, 2, 0);

   switch (which_forcing) {
      case 0 :
         return { std::make_unique<MeteoIdealized>(dm) };
      case 1 :
         {
            std::unique_ptr<MeteoNetcdf> meteo = std::make_unique<MeteoNetcdfRacmoPoint>(dm);
            meteo->readForcing();
            return meteo;
         }
      case 2 :
         {
            std::unique_ptr<MeteoNetcdf> meteo = std::make_unique<MeteoNetcdfRacmoGridded>(dm);
            meteo->readForcing();
            return meteo;
         }
      default:
         logger << "ERROR: unknown value: " << which_forcing << " for config option " << option_name << std::endl;
         std::abort();
   }
}

/* class Meteo*/
Meteo::Meteo(DynamicModel& dm) : _dm(dm) { }

/* class MeteoIdealized*/
MeteoIdealized::MeteoIdealized(DynamicModel& dm) : Meteo(dm){
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

double MeteoIdealized::surfaceTemperature() {
   const long nt = _dm.getCurrentTimeInSeconds();
   const long t =  nt % sec_in_year;
   return _ideal_T_amp*cos(2*M_PI*(double)t/(double)sec_in_year) + _ideal_T_mean;
}

/* class MeteoNetcdf */
MeteoNetcdf::MeteoNetcdf(DynamicModel& dm) : Meteo(dm){
   logger << "INFO: using Netcdf forcing" << std::endl;

   const char * option_name = "forcing:netcdf_dt";
   _dt_forcing = config.getInt(option_name, true, 0, (int)1e9, 0);
   //if (_dm.getDt() != _dt_forcing) {
   //   logger << "ERROR: " << option_name << " [ " << _dt_forcing << " ] must equal model dt [ " << _dm.getDt() << " ]" << std::endl;
   //   std::abort();
   //}
   logger << "INFO: " << option_name << " = " << _dt_forcing << std::endl;
}

double MeteoNetcdf::surfaceTemperature() {
   /* return surface temperature at current time step
      currently assumes that dt = dt_forcing (this checked in the constructor) 
      */
   long nt = _dm.getCurrentTimeStep();
   long dt_model = _dm.getDt();
   int idx = (int) (nt * dt_model / _dt_forcing); // index in forcing array
   idx = idx % _tskin_all.size();
   return _tskin_all[idx];
}

double MeteoNetcdf::accumulationRate() {
   /* return acccumulation rate at current time step
      currently assumes that dt = dt_forcing (this checked in the constructor) 
      */
   long nt = _dm.getCurrentTimeStep();
   int idx = (int) nt;
   idx = idx % _acc_all.size();
   return _acc_all[idx] / _dt_forcing;
}

double MeteoNetcdf::surfaceWind(){
   /* return wind speed at current time step
      currently assumes that dt = dt_forcing (this checked in the constructor) 
      */
   long nt = _dm.getCurrentTimeStep();
   int idx = (int) nt;
   idx = idx % _w10m_all.size();
   return _w10m_all[idx];
}

/* class MeteoNetcdfRacmoPoint */
MeteoNetcdfRacmoPoint::MeteoNetcdfRacmoPoint(DynamicModel& dm) : MeteoNetcdf(dm){
   logger << "INFO: using RACMO point forcing" << std::endl;
}

void MeteoNetcdfRacmoPoint::readForcing(){
   /* Open NetCDF forcing files and read them into memory entirely 
      currently it is assumed that the files describe a single point, 
      thus they are dependent on time only. */

   std::string f_acc = config.getFilename("forcing:netcdf_acc", true, true, "");
   std::string f_w10m = config.getFilename("forcing:netcdf_w10m", true, true, "");
   std::string f_tskin = config.getFilename("forcing:netcdf_tskin", true, true, "");

   logger << "INFO: opening RACMO point file " << f_acc << std::endl;
   logger << "INFO: opening RACMO point file " << f_w10m << std::endl;
   logger << "INFO: opening RACMO point file " << f_tskin << std::endl;

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

   int ndims = var1.getDimCount();
   std::vector<netCDF::NcDim> dims = var1.getDims(); 
   if (ndims != 4 || dims[1].getSize() != 1 || dims[2].getSize() != 1 || dims[3].getSize() != 1) {
      logger << "ndims = " << ndims << " (required: 4)" << std::endl;
      for (int i=0; i<ndims; i++) {
         logger << "dim = " << i << ", len = " << dims[i].getSize() << std::endl;
      }
      logger << "ERROR: forcing file does not contain valid RACMO point data" << std::endl;
      std::abort();
   }

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

/* class MeteoNetcdfRacmoGridded */
MeteoNetcdfRacmoGridded::MeteoNetcdfRacmoGridded(DynamicModel& dm) : MeteoNetcdf(dm){
   logger << "INFO: using RACMO gridded forcing" << std::endl;

}

void MeteoNetcdfRacmoGridded::readForcing(){
   /* Open NetCDF forcing files and read them into memory entirely 
      RACMO variable layout : (time, height, x, y) which is assumed to simplify to
      (time, 1, x, y), i.e. only a single height is present */

   std::string f_acc = config.getFilename("forcing:netcdf_acc", true, true, "");
   std::string f_w10m = config.getFilename("forcing:netcdf_w10m", true, true, "");
   std::string f_tskin = config.getFilename("forcing:netcdf_tskin", true, true, "");

   logger << "INFO: opening RACMO gridded data file " << f_acc << std::endl;
   logger << "INFO: opening RACMO gridded data file " << f_w10m << std::endl;
   logger << "INFO: opening RACMO gridded data file " << f_tskin << std::endl;

   const double my_lat = config.getDouble("forcing:netcdf_lat", true, -90., 90., -999.);
   const double my_lon = config.getDouble("forcing:netcdf_lon", true, -180., 180., -999.);

   netCDF::NcFile nc_file_acc(f_acc, netCDF::NcFile::read);
   netCDF::NcFile nc_file_w10m(f_w10m, netCDF::NcFile::read);
   netCDF::NcFile nc_file_tskin(f_tskin, netCDF::NcFile::read);

   netCDF::NcDim time1 = nc_file_acc.getDim("time");
   netCDF::NcDim time2 = nc_file_w10m.getDim("time");
   netCDF::NcDim time3 = nc_file_tskin.getDim("time");

   const unsigned nt = time1.getSize();
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

   netCDF::NcVar nclon = nc_file_acc.getVar("lon");
   netCDF::NcVar nclat = nc_file_acc.getVar("lat");

   int ndims = nclon.getDimCount();
   //logger << "INFO: variable lon has " << ndims << " dimensions " << std::endl;
   if (ndims != 2) {
      logger << "ndims = " << ndims << std::endl;
      logger << "ERROR: forcing file is not in RACMO format" << std::endl;
      std::abort();
   }
   std::vector<netCDF::NcDim> dims = nclon.getDims(); 
   const unsigned nrlat = dims[0].getSize();
   const unsigned nrlon = dims[1].getSize();
   logger << "INFO: rlon/rlat first dimension : " << nrlat << std::endl;
   logger << "INFO: rlon/rlat second dimension : " << nrlon << std::endl;

   netCDF::NcVar var1 = nc_file_acc.getVar("acc");
   netCDF::NcVar var2 = nc_file_w10m.getVar("w10m");
   netCDF::NcVar var3 = nc_file_tskin.getVar("tskin");

   ndims = var1.getDimCount();
   dims = var1.getDims(); 
   if (ndims != 4) { 
      logger << "ndims = " << ndims << " (required: 4)" << std::endl;
      logger << "ERROR: forcing file does not contain valid RACMO point data" << std::endl;
      std::abort();
   }
   if (dims[0].getSize() != nt || dims[1].getSize() != 1 || dims[2].getSize() != nrlat || dims[3].getSize() != nrlon) {
      logger << "dim 0, len = " << dims[0].getSize() << " (required: " << nt << ")" << std::endl;
      logger << "dim 1, len = " << dims[1].getSize() << " (required: " << 1 << ")" << std::endl;
      logger << "dim 2, len = " << dims[2].getSize() << " (required: " << nrlat << ")" << std::endl;
      logger << "dim 3, len = " << dims[3].getSize() << " (required: " << nrlon << ")" << std::endl;
      logger << "ERROR: forcing file does not contain valid RACMO point data" << std::endl;
      std::abort();
   }

   std::vector<double> rlats, rlons;
   rlats.resize(nrlat*nrlon);
   rlons.resize(nrlat*nrlon);

   nclat.getVar(&rlats.front());
   nclon.getVar(&rlons.front());

   double mindist, dist;
   int idx1, idx2;
   mindist = 1e5; 
   for (int ilat = 0; ilat < nrlat; ++ilat) {
      for (int ilon = 0; ilon < nrlon; ++ilon) {
         dist = sqrt(pow(rlats[ilat*nrlon+ilon]-my_lat,2) + pow(rlons[ilat*nrlon+ilon]-my_lon,2));
         if (dist < mindist) {
            mindist = dist;
            idx1 = ilat;
            idx2 = ilon;
         }
      }
   }
   logger << "INFO: Best rlat/ilat indices: (" << idx1 << ", " << idx2 << ")"<< std::endl;
   logger << "INFO: requested lat = " << my_lat << ", closest lat =" << rlats[idx1*nrlon+idx2] << std::endl;
   logger << "INFO: requested lon = " << my_lon << ", closest lon =" << rlons[idx1*nrlon+idx2] << std::endl;
   logger << "INFO: dist = " << mindist << std::endl;

   // Naming follows : https://www.unidata.ucar.edu/software/netcdf/docs/cxx4/classnetCDF_1_1NcVar.html#a74d273a0d5f572d04b78c92cd0f5d41b
   std::vector<size_t> start = {0, 0, idx1, idx2};
   std::vector<size_t> count = {nt, 1, 1, 1};

   var1.getVar(start, count, &_acc_all.front());
   var2.getVar(start, count, &_w10m_all.front());
   var3.getVar(start, count, &_tskin_all.front());
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

} // namespace

