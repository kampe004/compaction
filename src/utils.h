#ifndef UTILS_H
#define UTILS_H

#include <limits>

namespace DSM{ 

/* Contains some general purpose routines */

static const double eps = std::numeric_limits<double>::epsilon(); // machine precision

bool file_exists(const char *fileName);
bool doubles_equal(const double d1, const double d2); 

} // namespace
#endif
