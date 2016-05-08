#ifndef COMPACTION_H
#define COMPACTION_H

#include <memory>

namespace DSM{ 

class Compaction;
class ModelState;
class DynamicModel;
std::unique_ptr<Compaction> instantiate_compaction(ModelState& mstate, DynamicModel& dm); // based on the INI file, instantiate the appropriate Meteo-subtype

class Compaction{
 public:
   Compaction(ModelState& mstate, DynamicModel& dm);
   ~Compaction() {}
   virtual void compaction() =0;
 protected:
   ModelState& _mstate;
   DynamicModel& _dm;
};

class CompactionNone : public Compaction{
 public:
   CompactionNone(ModelState& mstate, DynamicModel& dm);
   void compaction() {};
};

class CompactionHerronLangway : public Compaction{
 public:
   CompactionHerronLangway(ModelState& mstate, DynamicModel& dm);
   void compaction();
};

class CompactionAnderson : public Compaction{
 public:
   CompactionAnderson(ModelState& mstate, DynamicModel& dm);
   void compaction();
};

class CompactionBarnolaPimienta : public Compaction{
 public:
   CompactionBarnolaPimienta(ModelState& mstate, DynamicModel& dm);
   void compaction();
};

class CompactionLigtenberg : public Compaction{
 public:
   CompactionLigtenberg(ModelState& mstate, DynamicModel& dm);
   void compaction();
};

class CompactionCROCUS : public Compaction{
 public:
   CompactionCROCUS(ModelState& mstate, DynamicModel& dm);
   void compaction();
};

} // namespace

#endif
