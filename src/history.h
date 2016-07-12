#ifndef HISTORY_H
#define HISTORY_H

namespace DSM{ 

class ModelState;

struct HistVar {
   double sum; 
   long N;
};

class History{
 public:
   History(ModelState& mstate);
   ~History() {};
   void update();
   void writeHistory();
 private:
   ModelState& _mstate;
   HistVar _ro1m; // density at 1 metre depth
   HistVar _ro1mavg; // average interval density (0,1) metre depth
};

} // namespace

#endif
