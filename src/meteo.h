#ifndef METEO_H
#define METEO_H

#include <memory>
#include <vector>

#include "constants.h"

namespace DSM{ 

class Meteo;
class DynamicModel;
std::unique_ptr<Meteo> instantiate_meteo(DynamicModel& dm); // based on the INI file, instantiate the appropriate Meteo-subtype

class Meteo{
   /* abstract class
      represents meteorological forcing at the surface */
 public: 
   Meteo(DynamicModel& dm);
   ~Meteo() {}; 

   void updateAnnualStatistics();

   double annualAcc() { /* annually integrated accumulation */ return _acc_ann_mean; }
   double annualTskin() { /* annual average surface temperature */ return _Ts_ann_mean; }
   double annualW10m() { /* annual average wind speed */ return _w10m_ann_mean; }

   virtual double surfaceTemperature() =0;
   virtual double accumulationRate() =0;
   virtual double surfaceWind() =0;

 protected:
   DynamicModel& _dm; // reference to the host model

   // some annual meteorological statistics (used by Helsen2008)
   double _acc_ann_mean;
   double _w10m_ann_mean;
   double _Ts_ann_mean;
};

class MeteoIdealized : public Meteo{
   /* idealized forcing */
 public:
   MeteoIdealized(DynamicModel& dm);
   ~MeteoIdealized() {}; 
   double surfaceTemperature() { return _ideal_T_mean; }
   double accumulationRate() { return _ideal_acc / sec_in_year; }
   double surfaceWind() { return _ideal_w10m; }
 private:
   double _ideal_T_mean;
   double _ideal_T_amp;
   double _ideal_acc;
   double _ideal_w10m;
};

class MeteoNetcdf : public Meteo{
   /* idealized forcing */
 public:
   MeteoNetcdf(DynamicModel& dm);
   ~MeteoNetcdf() {}; 
   double surfaceTemperature();
   double accumulationRate();
   double surfaceWind();

 private:
   // raw data from NetCDF files
   int _dt_forcing; // time step for forcing files [s]
   std::vector<double> _acc_all;
   std::vector<double> _tskin_all;
   std::vector<double> _w10m_all;

   void readForcing();
};

} // namespace

#endif
