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

int main(int argc, char** argv){
   //logger.open("densification.log");
   logger.rdbuf( std::cout.rdbuf() );
   logger.precision(16);

   if (argc < 2) {
      logger << "Usage: " << argv[0] << " <config_file_name>" << std::endl;   
      return -1;
   }
   config.init(argv[1]);

   DynamicModel dm;
   dm.run();

   //logger.close();
   return 0;
}
