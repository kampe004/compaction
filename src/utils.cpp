#include <fstream>
#include <cmath>

#include "utils.h"

namespace DSM{ 

bool file_exists(const char *fileName) {
   std::ifstream infile(fileName);
   return infile.good();
}

bool doubles_equal(double d1, double d2){
      return (std::abs(d1-d2) < 2.*eps ? true : false);
}


} // namespace
