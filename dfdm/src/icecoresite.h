#ifndef ICECORESITE_H
#define ICECORESITE_H

#include <string>
#include <vector>
#include <algorithm>
#include "constants.h"
#include "settings.h"

namespace Densification { 

struct Layer {
    double T;
    double dens;
    double mass;
    double dz;
};

class IceCoreSite {
    /// abstract base class for all ice core classes
private:
    bool is_initialized = false;
    const double ref_height = 0.4; //  		# reference height of layer in meters
    const double dzmax = 2*ref_height; //  # maximum thickness of layer
    const double dzmin = 0.5*ref_height; // # minimum thickness of layer
    const unsigned int initial_layers = 1;

    // layers are ordered bottom [grid.front()] to top [grid.back()]
    std::vector<Layer> grid;

    // initialize constants in case the default constructor is used
    bool have_diffusion = false;
    DensificationMethod densification_method = DensificationMethod::Ar10T;
    // constants for overburden compaction 
    double eta0 = 9e5;
    double c5 = 0.08;
    double c6 = 0.023;

    double getDepthOfDensity(double dens);

    long current_time = 0; // time in seconds since start of simulation

protected:

    void accumulate(long dt);
    void compact(long dt);
    void compactAr10T(long dt);
    void compactOverburden(long dt);
    void heatDiffusion(long dt);
    double getZ550();
    double getZ830();
 
    // pure virtual functions (must be implemented by derived class)
    virtual double surfaceDensity(long time) =0;
    virtual double surfaceTemperature(long time) =0;
    virtual double accumulationRate(long time) =0;
    virtual double annualIntegratedAccumulation() =0;
    virtual double annualMeanSurfaceTemperature() =0;
//    virtual double annualSurfaceDensity() =0;

public:
    IceCoreSite() {};
    IceCoreSite(Settings& settings);
    ~IceCoreSite() {}; 

    void init();
    void runTimeStep(long dt);
    void writeFiles();

    long getSecondsSinceStart();

    // member functions with implementation
    int gridsize() { return grid.size(); };
    double maxTemp() {
       double retval = grid[0].T;
       for (int i = 1; i < gridsize(); i++){
          retval = std::max(retval,grid[i].T);
       }
       return retval;
    }
    double minTemp() {
       double retval = grid[0].T;
       for (int i = 1; i < gridsize(); i++){
          retval = std::min(retval,grid[i].T);
       }
       return retval;
    }
    double maxDens() {
       double retval = grid[0].dens;
       for (int i = 1; i < gridsize(); i++){
          retval = std::max(retval,grid[i].dens);
       }
       return retval;
    }

    // polymorphic functions (should be overridden by derived class)
    virtual std::string toString();
    virtual void printIceCoreSummary();
};


} 
#endif
