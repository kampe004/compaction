#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include "iniparser.h"

namespace DSM{ 

static const char config_name[] = "settings.ini";
class Config; 
extern Config config; // globally unique config file

class Config{
   /* Wrapper class for our INI parser, to hide the actual library that is used
      and make future migration easier */
 public:
   Config();
   ~Config();

   int getInt(const char* name, bool required, int minval, int maxval, int defval);
   double getDouble(const char* name, bool required, double minval, double maxval, double defval);
   std::string getFilename(const char* name, bool required, bool check_existence, const char* defval);
   //std::string getFilename(const char* name, bool required, bool check_existence); // version without defval

 private:
   dictionary * _dict;
};

} // namespace

#endif
