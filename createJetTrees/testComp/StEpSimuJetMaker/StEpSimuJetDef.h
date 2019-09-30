// Jet Algo Container

#ifndef StEpSimuJetDef_def
#define StEpSimuJetDef_def

#include "TObject.h"
#include "TRefArray.h"
#include "TVector3.h"

class StEpSimuJetEvent;
class StEpSimuJet;

class StEpSimuJetDef : public TObject {
 public:

  StEpSimuJetDef();

  enum ALGO { AKT, KT };
  enum FRAME { BEAM, GAMMA, BREIT};

  // Declare Setters
  void setRadius(float rad);
  void setMinPt(float pt);
  void setAlgo(int algo);
  void setFrame(int frame);

  // Declare and Implement Getters
  float radius() const { return mRadius; }
  float minPt() const { return mMinPt; }
  int algo() const { return mALGO; }
  int frame() const { return mFRAME; }


  // Get Number of Jets
  int numberOfJets() const;

  // Get Jet
  StEpSimuJet* jet(int i) const; //{ return (StForwardJet*)mJets.At(i); }

  // Add Jet to Dijet
  StEpSimuJet* addJet(StEpSimuJet* jet); //{ mJets.Add((TObject*)jet); return (StForwardJet*)mJets.Last(); }
  
 private:

  float mRadius;
  float mMinPt;
  int mALGO;
  int mFRAME;

  TRefArray mJets;
  
  ClassDef(StEpSimuJetDef,3);
};

#endif
