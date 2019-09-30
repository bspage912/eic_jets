// Test of Jet Container

#include "StEpSimuSubJet.h"

ClassImp(StEpSimuSubJet);

StEpSimuSubJet::StEpSimuSubJet(Double_t pt, Double_t eta, Double_t phi, Double_t E)
  : mPt(pt)
  , mEta(eta)
  , mPhi(phi)
  , mE(E)
{
}

// Setters
void StEpSimuSubJet::setRadius(Double_t rad)
{
  mRadius = rad;
}

void StEpSimuSubJet::setNumParticles(Double_t particles)
{
  mParticles = particles;
}

// Getters
Double_t StEpSimuSubJet::rap()
{
  return StEpSimuSubJet::fourMomentum().Rapidity();
}


// Set Parent Jet Definition
void StEpSimuSubJet::setParentJet(const StEpSimuJet* jet) { mParentJets = (TObject*)jet; }


