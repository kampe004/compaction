#ifndef DYNAMICMODEL_H
#define DYNAMICMODEL_H

namespace DSM{ 

class ModelState;

class DynamicModel{
 public:
   DynamicModel();
   ~DynamicModel() {};

   // getters and setters
   long getCurrentTimeInSeconds() { return {_nt * _dt}; };
   long getCurrentTimeStep() { return _nt; };
   long getDt() { return _dt; };

   // methods
   void run();
   void runTimeStep(ModelState& mstate);
   void accumulate(ModelState& mstate);
   void doGridChecks(ModelState& mstate);
   void heatDiffusion(ModelState& mstate);

 private:
   long _nt;   // number of timesteps since start of simulation
   long _dt;   // simulation timestep [s]
   bool _has_heat; // heat diffusion
};

} // namespace

#endif
