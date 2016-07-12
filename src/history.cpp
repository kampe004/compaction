
#include "history.h"
#include "modelstate.h"
#include "logging.h"


namespace DSM{ 
History::History(ModelState& mstate) : _mstate(mstate) {
   _ro1m.sum      = 0.0;
   _ro1m.N        = 0;
   _ro1mavg.sum   = 0.0;
   _ro1mavg.N     = 0;
}

void History::update() {
   /* with every call to this function, the running means are updated */
   Grid& grid = _mstate.getGrid();

   // Density at 1 meter depth (if found)
   static const double d1m = 1.0;
   double ro1m, ro1m_avg;
   double dm, dp; // depth+ and depth-
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
