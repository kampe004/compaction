#ifndef DYNAMICMODEL_H
#define DYNAMICMODEL_H

namespace DSM{ 

class ModelState;
class Metamorphism;
class Compaction;

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

 private:
   long _nt;   // number of timesteps since start of simulation
   long _dt;   // simulation timestep [s]
   bool _has_heat; // heat diffusion

   void runTimeStep(ModelState& mstate, Metamorphism& mm, Compaction& comp);
   void accumulate(ModelState& mstate);
   void doGridChecks(ModelState& mstate);
   void heatDiffusion(ModelState& mstate);
   void heatDiffusionShallow(ModelState& mstate);
};

} // namespace

#endif
