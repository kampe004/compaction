
#include <iostream>

#include "../config.h"
#include "../constants.h"
#include "../logging.h"
#include "../meteo.h"
#include "../surfacedensity.h"

#include "test.h"

namespace DSM{
/* globally unique instances */
std::ostream logger(NULL);
ConfigParser config;

void set_logger(const std::ostream& os){
   logger.rdbuf(os.rdbuf());
}


void Test::testSlater2016(){
   logger << "Test::testSlater2016()" << std::endl;
   MeteoSimple meteo(*this);
   SurfaceDensitySlater2016 surf(meteo);

   double rho;

   meteo.setTemperature(T0-5.);
   meteo.setWind(1.0);
   rho = surf.density();
   logger << "T = " << meteo.surfaceTemperature() << ", U = " << meteo.surfaceWind() << ", rho = " << rho << std::endl;

   meteo.setTemperature(T0-10.);
   meteo.setWind(1.0);
   rho = surf.density();
   logger << "T = " << meteo.surfaceTemperature() << ", U = " << meteo.surfaceWind() << ", rho = " << rho << std::endl;

   meteo.setTemperature(T0-10.);
   meteo.setWind(3.0);
   rho = surf.density();
   logger << "T = " << meteo.surfaceTemperature() << ", U = " << meteo.surfaceWind() << ", rho = " << rho << std::endl;

   meteo.setTemperature(T0-10.);
   meteo.setWind(6.0);
   rho = surf.density();
   logger << "T = " << meteo.surfaceTemperature() << ", U = " << meteo.surfaceWind() << ", rho = " << rho << std::endl;

   meteo.setTemperature(T0-10.);
   meteo.setWind(9.0);
   rho = surf.density();
   logger << "T = " << meteo.surfaceTemperature() << ", U = " << meteo.surfaceWind() << ", rho = " << rho << std::endl;
}
} // namespace
 
using namespace DSM; 

int main(int argc, char** argv){
   std::ostream consolelogger(NULL);
   consolelogger.rdbuf( std::cout.rdbuf() );
   std::ofstream filelogger("model.log", std::ofstream::out);

   /* Choose one of the following: */
   set_logger(consolelogger);
   //set_logger(filelogger);
   logger.precision(16);

   Test t;
   t.testSlater2016();

   return 0;
}

