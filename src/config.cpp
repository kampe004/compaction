
#include <cstring>

#include "logging.h"
#include "utils.h"
#include "config.h"

namespace DSM{ 

ConfigParser::ConfigParser(){ }

ConfigParser::~ConfigParser(){
   iniparser_freedict(_dict);
}

void ConfigParser::init(const char* config_name){
   _dict = iniparser_load(config_name);
}

int ConfigParser::getInt(const char* name, bool required, int minval, int maxval, int defval){
   const char * rtnnam = "ConfigParser::getInt()";
   int tmp;
   if (required){
      int missval = -9999;
      tmp = iniparser_getint(_dict, name, missval);
      if (tmp == missval){
         logger << "ERROR: " << rtnnam << " required value not found: " << name << std::endl;
         std::abort();
      }
   } else {
      tmp = iniparser_getint(_dict, name, defval);
   }
   if (tmp < minval || tmp > maxval) {
      logger << "ERROR: " << rtnnam << " value out of bounds " << name << std::endl;
      logger << "ERROR: " << rtnnam << " val = " << tmp << ", minval = " << minval << ", maxval = " << maxval << std::endl;
      std::abort();
   }
   logger << "INFO: " << name << " = " << tmp << std::endl;
   return tmp;
}

double ConfigParser::getDouble(const char* name, bool required, double minval, double maxval, double defval){
   const char * rtnnam = "ConfigParser::getDouble()";
   double tmp;
   if (required){
      double missval = -9999.;
      tmp = iniparser_getdouble(_dict, name, missval);
      if (doubles_equal(tmp, missval)){
         logger << "ERROR: " << rtnnam << " required value not found: " << name << std::endl;
         std::abort();
      }
   } else {
      tmp = iniparser_getdouble(_dict, name, defval);
   }
   if (tmp < minval || tmp > maxval) {
      logger << "ERROR: " << rtnnam << " value out of bounds " << name << std::endl;
      logger << "ERROR: " << rtnnam << " val = " << tmp << ", minval = " << minval << ", maxval = " << maxval << std::endl;
      std::abort();
   }
   logger << "INFO: " << name << " = " << tmp << std::endl;
   return tmp;
}

std::string ConfigParser::getFilename(const char* name, bool required, bool check_existence, const char* defval ){
   const char * rtnnam = "ConfigParser::getFileName()";
   const char * tmp;
   if (required){
      tmp = iniparser_getstring(_dict, name, "NONE");
      if (strcmp(tmp, "NONE") == 0){
         logger << "ERROR: " << rtnnam << " required value not found: " << name << std::endl;
         std::abort();
      }
   } else {
      tmp = iniparser_getstring(_dict, name, defval);
   }
   logger << "INFO: " << name << " = " << tmp << std::endl;
   if (!file_exists(tmp)){
      logger << "ERROR: " << rtnnam << " file does not exist: " << tmp << std::endl;
      std::abort();
   }
   return std::string(tmp);
}

} // namespace
