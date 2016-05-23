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

void set_logger(const std::ostream& os){
   logger.rdbuf(os.rdbuf());
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
