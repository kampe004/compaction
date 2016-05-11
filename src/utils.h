#ifndef UTILS_H
#define UTILS_H

#include <limits>

namespace DSM{ 

/* Contains some general purpose routines */

static const double eps = std::numeric_limits<double>::epsilon(); // machine precision

bool file_exists(const char *fileName);
bool doubles_equal(double d1, double d2); 

} // namespace
#endif
