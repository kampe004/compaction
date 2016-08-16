#ifndef SURFACEDENSITY_H
#define SURFACEDENSITY_H

#include <memory>

namespace DSM{ 

class SurfaceDensity;
class Meteo;
std::unique_ptr<SurfaceDensity> instantiate_surfacedensity(Meteo& meteo); // based on the INI file, instantiate the appropriate Meteo-subtype

class SurfaceDensity{
 public: 
   SurfaceDensity(Meteo& meteo);
   ~SurfaceDensity() {};
   virtual double density()=0;
 protected:
   Meteo& _meteo;
};

class SurfaceDensityConstant : public SurfaceDensity { 
 public:
   SurfaceDensityConstant(Meteo& meteo);
   double density();
 private:
   double _val;
};

class SurfaceDensityHelsen2008 : public SurfaceDensity { 
 public:
   SurfaceDensityHelsen2008(Meteo& meteo);
   double density();
};

class SurfaceDensityLenaerts2012 : public SurfaceDensity { 
 public:
   SurfaceDensityLenaerts2012(Meteo& meteo);
   double density();
};

class SurfaceDensityCROCUS : public SurfaceDensity { 
 public:
   SurfaceDensityCROCUS(Meteo& meteo);
   double density();
};

class SurfaceDensityAnderson : public SurfaceDensity { 
 public:
   SurfaceDensityAnderson(Meteo& meteo);
   double density();
};

class SurfaceDensityAndersonListon : public SurfaceDensity { 
 public:
   SurfaceDensityAndersonListon(Meteo& meteo);
   double density();
 private:
   std::unique_ptr<SurfaceDensityAnderson> _sda; // needed to use the Anderson density() function to avoid code duplication
};

class SurfaceDensitySlater2016 : public SurfaceDensity { 
 public:
   SurfaceDensitySlater2016(Meteo& meteo);
   double density();
};

} // namespace

#endif
