#ifndef METAMORPHISM_H
#define METAMORPHISM_H

#include <memory>

namespace DSM{ 

class Metamorphism;
class ModelState;
class DynamicModel;
std::unique_ptr<Metamorphism> instantiate_metamorphism(ModelState& mstate, DynamicModel& dm); // based on the INI file, instantiate the appropriate Meteo-subtype

class Metamorphism{
 public:
   Metamorphism(ModelState& mstate, DynamicModel& dm);
   ~Metamorphism() {};
   virtual void metamorphism() =0;
 protected:
   ModelState& _mstate;
   DynamicModel& _dm;
};

class MetamorphismNone : public Metamorphism{
 public:
   MetamorphismNone(ModelState& mstate, DynamicModel& dm);
   void metamorphism();
};

class MetamorphismAnderson1976 : public Metamorphism{
 public:
   MetamorphismAnderson1976(ModelState& mstate, DynamicModel& dm);
   void metamorphism();
};

class MetamorphismCROCUS : public Metamorphism{
 public:
   MetamorphismCROCUS(ModelState& mstate, DynamicModel& dm);
   void metamorphism();
};

} // namespace

#endif
