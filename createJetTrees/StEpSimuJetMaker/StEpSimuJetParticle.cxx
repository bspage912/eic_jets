// Test of Track Container

#include "StEpSimuJetParticle.h"

ClassImp(StEpSimuJetParticle);

StEpSimuJetParticle::StEpSimuJetParticle(Double_t pt, Double_t eta, Double_t phi, Double_t e)
  : mPt(pt)
  , mEta(eta)
  , mPhi(phi)
  , mE(e)
{
}

// Implement Setters
void StEpSimuJetParticle::setIndex(int id) { mIndex = id; }

Double_t StEpSimuJetParticle::rap()
{
  return StEpSimuJetParticle::fourMomentum().Rapidity();
}

