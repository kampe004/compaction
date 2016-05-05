#ifndef SETTINGS_H
#define SETTINGS_H

namespace Densification{

enum class DensificationMethod {Ligtenberg2011, Anderson1976, Barnola1991, Spencer2001, BarnolaSpencer, HerronLangway};
std::ostream& operator<<(std::ostream& os, DensificationMethod dm); // overload the '<<' operator for function toString()

enum class ForcingMethod {Idealized, Netcdf};
std::ostream& operator<<(std::ostream& os, ForcingMethod dm); // overload the '<<' operator for function toString()

typedef struct Settings Settings; // data container to hold user settings
struct Settings {

   /* general stuff (physics used etc) */
   bool have_diffusion;
   DensificationMethod dm;
   double max_depth; // maximum depth of 1D firn model at which to consider the simulation failed (prevents excessive runtimes)
   int max_year; // maximum number of years at which to consider the simulation failed (prevents excessive runtimes)
   ForcingMethod fm; 

   /* fresh snow density */
   double rho_s;

   /* Anderson 1976 parameters */
   double eta0;
   double c5;
   double c6;

   /* meteorological forcing */
   int forcing_dt; 
   std::string f_acc; 
   std::string f_wind10m;
   std::string f_tskin;

   /* idealized forcing */
   double T_mean; // annual mean temperature [K]
   double Tamp; // amplitude of annual cycle  [K]
   double acc_mean; // annual accumulation [mm]
   double v10m_mean; // 10m wind speed [m/s]
};
}

#endif
