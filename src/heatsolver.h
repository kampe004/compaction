#ifndef HEATSOLVER_H
#define HEATSOLVER_H

#include <memory>

namespace DSM{ 

class HeatSolver;
class ModelState;
class DynamicModel;
std::unique_ptr<HeatSolver> instantiate_heatsolver(ModelState& mstate, DynamicModel& dm); // based on the INI file, instantiate the appropriate Meteo-subtype

extern "C" void dgtsv_(int* dim, int* nrhs, double* dl, double* d, double* du, double* b, int* ldb, int* info);

class HeatSolver{
 public:
   HeatSolver(ModelState& mstate, DynamicModel& dm);
   ~HeatSolver() {};
   virtual void heatdiffusion() =0;
 protected:
   ModelState& _mstate;
   DynamicModel& _dm;
};

class HeatSolverNone : public HeatSolver{
 public:
   HeatSolverNone(ModelState& mstate, DynamicModel& dm);
   void heatdiffusion();
};

class HeatSolverExplicit : public HeatSolver{
 public:
   HeatSolverExplicit(ModelState& mstate, DynamicModel& dm);
   void heatdiffusion();
};

class HeatSolverImplicit : public HeatSolver{
 public:
   HeatSolverImplicit(ModelState& mstate, DynamicModel& dm);
   void heatdiffusion();
};

} // namespace

#endif
