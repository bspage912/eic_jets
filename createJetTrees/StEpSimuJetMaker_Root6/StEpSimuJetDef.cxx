// Test of Dijet

#include "StEpSimuJetDef.h"
#include "StEpSimuJet.h"
#include "TVector3.h"

ClassImp(StEpSimuJetDef);

StEpSimuJetDef::StEpSimuJetDef()
{
}

void StEpSimuJetDef::setRadius(float rad) { mRadius = rad; }
void StEpSimuJetDef::setMinPt(float pt) { mMinPt = pt; }
void StEpSimuJetDef::setAlgo(int algo) { mALGO = algo; }
void StEpSimuJetDef::setFrame(int frame) { mFRAME = frame; }


// Get Number of Dijets
int StEpSimuJetDef::numberOfJets() const { return mJets.GetEntriesFast(); }

// Get Jet
StEpSimuJet* StEpSimuJetDef::jet(int i) const { return (StEpSimuJet*)mJets.At(i); }

// Add Jet to Dijet
StEpSimuJet* StEpSimuJetDef::addJet(StEpSimuJet* jet) { mJets.Add((TObject*)jet); return (StEpSimuJet*)mJets.Last(); }

