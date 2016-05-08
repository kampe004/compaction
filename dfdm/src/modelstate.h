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
   //double d;    // dendricity
   //double s;    // sphericity
   //double gs;   // grain size [mm]
};

/* We chose the data structure to be std::vector.
   Layers are added when they are accumulating;
   this means the newest layer has the highest index and the deepest layer has index 0 */
typedef std::vector<Layer> Grid;

const double ref_height = 0.1;            // reference height of layer in meters
const double dzmax = 2*ref_height;        // maximum thickness of layer
const double dzmin = 0.8*ref_height;      // minimum thickness of layer
const unsigned int initial_layers = 1;

class ModelState {
   /* Holds physical state parameters per snow layer and methods to update layers */
 public:
   ModelState(Meteo& meteo, SurfaceDensity& surf, DynamicModel& dm);
   ~ModelState() {}; 

   void writeModelState();
   void printSummary();

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
