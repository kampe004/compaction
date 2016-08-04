#ifndef DYNAMICMODEL_H
#define DYNAMICMODEL_H

namespace DSM{ 

class ModelState;
class Metamorphism;
class Compaction;
class HeatSolver;

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

   void runTimeStep(ModelState& mstate, Metamorphism& mm, Compaction& comp, HeatSolver& hs);
   void doGridChecks(ModelState& mstate);
   void accumulate(ModelState& mstate);
   void accumulatePositive(ModelState& mstate);
   void accumulateNegative(ModelState& mstate);

   bool startOfDay(long curr_sec);
};

} // namespace

#endif
