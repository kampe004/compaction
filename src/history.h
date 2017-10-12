#ifndef HISTORY_H
#define HISTORY_H

#include <memory>

namespace netCDF{
class NcFile;
}


namespace DSM{ 
class ModelState;

struct HistVar {
   double sum; 
   long N;
   bool missing;
};

static const int NLEV = 101; // number of output levels
static const float FILL_VALUE = 1.e30; 

class History{
 public:
   History(ModelState& mstate);
   ~History();
   void update(bool start_of_day);
   void writeHistory();

 private:
   void initNetcdfOutput();
   void addNetcdfRecord();

   bool _have_netcdf_output;
   int _hist_freq;
 
   ModelState& _mstate;
   HistVar _ro1cm; // density at 1 centimetre depth
   HistVar _ro1m; // density at 1 metre depth
   HistVar _ro1mavg; // average interval density (0,1) metre depth
   HistVar _z550; // depth of 550 kg/m3 density
   HistVar _z830; // depth of 830 kg/m3 density
   HistVar _Tice; // deep ice temperature

//   netCDF::NcFile daily_output_nc;
   std::unique_ptr<netCDF::NcFile> datafile_ref;

   float _zlev[NLEV]; // depth of netCDF output levels
   int _nrec; // current number of records in netCDF
};

} // namespace

#endif
