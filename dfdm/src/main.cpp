/* --------------------------------------------------------------------------------
   
   DESCRIPTION
      simple dry snow model (DSM) with overburden compaction and heat diffusion
   
   DATE
      15-JAN-2016
   
   AUTHOR
      L.vankampenhout@uu.nl
   -------------------------------------------------------------------------------- */

#include <iostream>

#include "logging.h"
#include "config.h"
#include "dynamicmodel.h"

namespace DSM{
/* globally unique instances */
std::ostream logger(NULL); 
ConfigParser config;
} // namespace

using namespace DSM; 

int main(){
   //logger.open("densification.log");

   logger.rdbuf( std::cout.rdbuf() );
   logger.precision(16);

   DynamicModel dm;
   dm.run();

   //logger.close();
   return 0;
}
