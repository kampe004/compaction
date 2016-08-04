
#include <netcdf>
#include "history.h"
#include "modelstate.h"
#include "logging.h"
#include "config.h"


namespace DSM{ 
History::History(ModelState& mstate) : _mstate(mstate) {
   _ro1m.sum      = 0.0;
   _ro1m.N        = 0;
   _ro1mavg.sum   = 0.0;
   _ro1mavg.N     = 0;
   _nrec          = 0;

   const char * option_name;
   option_name = "history:netcdf_output";
   _have_netcdf_output = config.getInt(option_name, false, 0, 1, 0) == 0 ? false : true;

   if (_have_netcdf_output)
      initNetcdfOutput();
}

void History::initNetcdfOutput(){
   const char * option_name;
   option_name = "history:freq_days";
   _hist_freq = config.getInt(option_name, false, 1, 1, 1); // output frequency in days

   const char * filepath = "daily.nc";

   try{
      // this is a trick for older version of Netcdf4 library that does not have the open() function yet...
      // could be replaced by making daily_output_nc a data member of History directly, and use the open(..) function
      datafile_ref = std::make_unique<netCDF::NcFile>(filepath, netCDF::NcFile::replace);
      netCDF::NcFile& daily_output_nc = *datafile_ref;

//      logger << "netCDF id = " << daily_output_nc.getId() << std::endl;

      netCDF::NcDim d1 = daily_output_nc.addDim("time", 0);
      netCDF::NcDim d2 = daily_output_nc.addDim("lev", NLEV);

      std::vector<netCDF::NcDim> dims;
      dims.push_back(d1);
      dims.push_back(d2);

      netCDF::NcVar rho_nc = daily_output_nc.addVar("rho", netCDF::ncFloat, dims);
      netCDF::NcVar T_nc = daily_output_nc.addVar("T", netCDF::ncFloat, dims);

      rho_nc.putAtt("_FillValue", netCDF::ncFloat, FILL_VALUE);
      T_nc.putAtt("_FillValue", netCDF::ncFloat, FILL_VALUE);

      dims.clear();
      dims.push_back(d2);
      netCDF::NcVar zlev_nc = daily_output_nc.addVar("zlev", netCDF::ncFloat, dims);

      _zlev[0] = 0.0;
      for (unsigned i = 1; i<NLEV; i++){
         _zlev[i]=_zlev[i-1] + 0.5;   // 1/2 meter spacing for now
      }

      zlev_nc.putVar(_zlev);
   }
   catch(netCDF::exceptions::NcException& e) 
   {
      logger << "ERROR: netcdf exception in History::initNetcdfOutput()" << std::endl;
      e.what();
      throw;
   }
}

History::~History(){
//   logger << "History::~History()" << std::endl;
}

void History::update(bool start_of_day) {
   /* with every call to this function, the running means are updated */
   Grid& grid = _mstate.getGrid();

   // Density at 1 meter depth (if found)
   static const double d1m = 1.0;
   double ro1m, ro1m_avg;
   double dm, dp; // depth+ and depth-
   // TODO: use nodal depths, not interface depths? 
   dp = 0.0;
   ro1m_avg = 0.0;
   for (int i = grid.size()-1; i >= 0; i--) {
      dp += grid[i].dz;
      if (dp > d1m ) { 
         if (i == grid.size()-1) {
            // unlikely: top layer exceeds one meter thickness: 
            ro1m = grid[i].dens;
         } else {
            // linear interpolation between depth dm and dp
            dm = dp - grid[i].dz;
            ro1m = grid[i].dens * (d1m-dm)/(dp-dm) + grid[i+1].dens * (dp-d1m)/(dp-dm);
         }
         ro1m_avg    += (d1m-dm) * grid[i].dens;

         _ro1m.sum   += ro1m;
         _ro1m.N     += 1;

         _ro1mavg.sum   += ro1m_avg / d1m;
         _ro1mavg.N     += 1;
         break;
      }
      ro1m_avg += grid[i].dz * grid[i].dens;
   }


   if (_have_netcdf_output && start_of_day) {
      if (_hist_freq != 1) {
         logger << "ERROR: programmer error, hist_freq != 1 not supported" << std::endl;
         std::abort();
      }

      // compute nodal depths
      const int Np = grid.size();
      double zn[Np-1];
      zn[Np-1] = grid[Np-1].dz * 0.5;
      for (int i = Np-2; i >= 0; i--) {
         zn[i] = zn[i+1] + 0.5 * (grid[i+1].dz + grid[i].dz);
      }

      /* 
         Append one record to the netCDF output.
         The record written is a linear interpolation of the current model grid
         to the netCDF output levels (fixed depths).
         The first output level is the surface (0.0 m depth) which is not interpolated
         but assumed to be equal to the top level. 
   
         Note that the levels are number 0...NLEV top to bottom, whereas the model 
         layering is reversed (0 is the bottom layer).
      */ 
      float rho_out[NLEV]; 
      float T_out[NLEV];
      rho_out[0]  = grid[grid.size()-1].dens;
      T_out[0]    = grid[grid.size()-1].T;
      dp = 0.0;
      int i = Np-1;
      for (int k = 1; k < NLEV; k++) {
         while(i >= 0 && zn[i] < _zlev[k]) {
            i--; // traverse downwards the model state snow pack
         }
         if (i < 0) {
            // output level depth has not been reached by model
            rho_out[k] = FILL_VALUE;
            T_out[k] = FILL_VALUE;
         } else {
            // linear interpolation between depth dm and dp
            dp = zn[i];
            dm = zn[i+1];
            rho_out[k]  = grid[i].dens * (_zlev[k]-dm)/(dp-dm) + grid[i+1].dens * (dp-_zlev[k])/(dp-dm);
            T_out[k]    = grid[i].T * (_zlev[k]-dm)/(dp-dm) + grid[i+1].T * (dp-_zlev[k])/(dp-dm);
         }
      }
   
      std::vector<size_t> startp,countp;
      startp.push_back(_nrec);
      startp.push_back(0);
   
      countp.push_back(1);
      countp.push_back(NLEV);
   
      netCDF::NcFile& daily_output_nc = *datafile_ref;
      //logger << "netCDF id = " << daily_output_nc.getId() << std::endl;
      netCDF::NcVar rho_nc = daily_output_nc.getVar("rho");
      netCDF::NcVar T_nc = daily_output_nc.getVar("T");
   
      rho_nc.putVar(startp,countp,rho_out);
      T_nc.putVar(startp,countp,T_out);
   
      // increase record counter
      _nrec++; 
   } // _have_netcdf_output
}

void History::writeHistory() {
   /* writes running means to file*/
   std::ofstream f_ro1m;
   f_ro1m.open("ro1m.txt");
   if (_ro1m.N > 0) {
      f_ro1m << (_ro1m.sum / _ro1m.N) << std::endl ;
   } else {
      f_ro1m << -1 << std::endl;
   }
   f_ro1m.close();

   std::ofstream f_ro1mavg;
   f_ro1mavg.open("ro1mavg.txt");
   if (_ro1mavg.N > 0) {
      f_ro1mavg << (_ro1mavg.sum / _ro1mavg.N) << std::endl ;
   } else {
      f_ro1mavg << -1 << std::endl;
   }
   f_ro1mavg.close();
}



} // namespace
