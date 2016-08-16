#ifndef MODELSTATE_H
#define MODELSTATE_H

#include <string>
#include <vector>
#include <algorithm>

namespace DSM{ 

class Meteo;
class SurfaceDensity;
class DynamicModel;

struct Layer {
   double T;    // temperature [K]
   double dens; // density [kg/m3]
   double dz;   // thickness [m]
   double dnd;    // dendricity [-]
   double sph;    // sphericity [-]
   double gs;   // grain size [m]
};

/* We chose the data structure to be std::vector.
   Layers are added when they are accumulating;
   this means the newest layer has the highest index and the deepest layer has index 0 */
typedef std::vector<Layer> Grid;

const double ref_height = 0.1;            // reference height of layer in meters
const double dzmax = 2*ref_height;        // maximum thickness of layer
const double dzmin = 0.8*ref_height;      // minimum thickness of layer
const double dzmin_top = dzmax / 10; 
const unsigned int initial_layers = 1;


/* Constants that determine the 'ideal' layer thickness profile that is used in the 
   layer checking scheme. For both the minimum and the maximum thickness an exponential
   curve is defined by the constants a and b, given as:
      dz = a*exp(1/zk * log(b/a)) 
   where zk is the kink-depth, below which dz assumes a constant value (equal to that 
   at the kink depth).  
   */
static const double ref_height_min_a = 0.005; // max height of first layer
static const double ref_height_min_b = 1.00; // max height of first layer
static const double ref_height_min_kink = 50.; // depth of kink

static const double ref_height_max_a = 0.01; // max height of first layer
static const double ref_height_max_b = 3.00; // max height of first layer
static const double ref_height_max_kink = 50.; // depth of kink


class ModelState {
   /* Holds physical state parameters per snow layer and methods to update layers */
 public:
   ModelState(Meteo& meteo, SurfaceDensity& surf, DynamicModel& dm);
   ~ModelState() {}; 

   void writeModelState();
   void printSummary();
   void combineLayers(Layer& lay1, Layer& lay2);

   bool hasReachedDensity(const double dens);

   int gridsize() { return _grid.size(); };

   double getDepthOfDensity(const double dens);
   double totalDepth();
   double maxTemp();
   double minTemp();
   double maxDens();

   Grid& getGrid() { return _grid; };
   Meteo& getMeteo() { return _meteo; };
   SurfaceDensity& getSurf() { return _surf; };

 private:
   Meteo& _meteo;
   DynamicModel& _dm;
   SurfaceDensity& _surf;
   Grid _grid; 
   double getZ550();
   double getZ830();
};

} // namespace
#endif
