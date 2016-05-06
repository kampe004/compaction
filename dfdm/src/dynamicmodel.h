#ifndef DYNAMICMODEL_H
#define DYNAMICMODEL_H

namespace DSM{ 

class DynamicModel{
 public:
   DynamicModel();
   ~DynamicModel() {};
   void run();
   long getCurrentTimeInSeconds() { return {_nt * _dt}; };
   long getCurrentTimeStep() { return _nt; };
   long getDt() { return _dt; };

 private:
   long _nt;   // number of timesteps since start of simulation
   long _dt;   // simulation timestep [s]
};

} // namespace

#endif
