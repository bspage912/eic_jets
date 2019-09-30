// Test of Dijet Tree

#include "TClonesArray.h"
#include "TVector3.h"
#include "StEpSimuJetEvent.h"
#include "StEpSimuJetDef.h"
#include "StEpSimuJet.h"
#include "StEpSimuParticle.h"
#include "StEpSimuJetParticle.h"
#include "StEpSimuSubJet.h"


ClassImp(StEpSimuJetEvent);

  StEpSimuJetEvent::StEpSimuJetEvent()
    : mJetDefs(new TClonesArray("StEpSimuJetDef",50))
    , mJets(new TClonesArray("StEpSimuJet",50))
    , mJetParticles(new TClonesArray("StEpSimuJetParticle",50))
    , mParticles(new TClonesArray("StEpSimuParticle",50))
    , mSubJets(new TClonesArray("StEpSimuSubJet",50))
{
}

void StEpSimuJetEvent::Clear()
{
  mJetDefs->Clear();
  mJets->Clear();
  mParticles->Clear();
  mJetParticles->Clear();
  mSubJets->Clear();
}

// Implement Event Setters
void StEpSimuJetEvent::setGenEvent(Int_t genEvn) { mGenEvent = genEvn; }
void StEpSimuJetEvent::setTgtPartonX(Double32_t tgtX) { mTgtPartonX = tgtX; }
void StEpSimuJetEvent::setBeamPartonX(Double32_t beamX) { mBeamPartonX = beamX; }
void StEpSimuJetEvent::setBeamPartonTheta(Double32_t beamTheta) { mBeamPartonTheta = beamTheta; }
void StEpSimuJetEvent::setLeptonPhi(Double32_t lepPhi) { mLeptonPhi = lepPhi; }
void StEpSimuJetEvent::setF1(Double32_t f1) { mF1 = f1; }
void StEpSimuJetEvent::setSigmaRad(Double32_t sigRad) { mSigmaRad = sigRad; }
void StEpSimuJetEvent::setTHat(Double32_t tHat) { mTHat = tHat; }
void StEpSimuJetEvent::setUHat(Double32_t uHat) { mUHat = uHat; }
void StEpSimuJetEvent::setQ2Hat(Double32_t q2Hat) { mQ2Hat = q2Hat; }
void StEpSimuJetEvent::setSigRadCor(Double32_t sigRadCor) { mSigRadCor = sigRadCor; }
void StEpSimuJetEvent::setEBrems(Double32_t eBrems) { mEBrems = eBrems; }
void StEpSimuJetEvent::setPhotonFlux(Double32_t photonFlux) { mPhotonFlux = photonFlux; }
void StEpSimuJetEvent::setTrueY(Double32_t truey) { mTrueY = truey; }
void StEpSimuJetEvent::setTrueQ2(Double32_t trueq2) { mTrueQ2 = trueq2; }
void StEpSimuJetEvent::setTrueX(Double32_t truex) { mTrueX = truex; }
void StEpSimuJetEvent::setTrueW2(Double32_t truew) { mTrueW2 = truew; }
void StEpSimuJetEvent::setTrueNu(Double32_t truenu) { mTrueNu = truenu; }
void StEpSimuJetEvent::setF2(Double32_t f2) { mF2 = f2; }
void StEpSimuJetEvent::setR(Double32_t r) { mR = r; }
void StEpSimuJetEvent::setPt2Hat(Double32_t pt2Hat) { mPt2Hat = pt2Hat; }
void StEpSimuJetEvent::setSHat(Double32_t sHat) { mSHat = sHat; }

// Implement Event MC Setters
void StEpSimuJetEvent::setN(Int_t number) { mN = number; }
void StEpSimuJetEvent::setProcessID(Int_t process) { mProcessID = process; }
void StEpSimuJetEvent::setNTracks(Int_t numTracks) { mNTracks = numTracks; }

// Implement Event DIS Setters
void StEpSimuJetEvent::setX(Double32_t x) { mX = x; }
void StEpSimuJetEvent::setQ2(Double32_t q2) { mQ2 = q2; }
void StEpSimuJetEvent::setY(Double32_t y) { mY = y; }
void StEpSimuJetEvent::setW2(Double32_t w2) { mW2 = w2; }
void StEpSimuJetEvent::setNu(Double32_t nu) { mNu = nu; }
void StEpSimuJetEvent::setXJB(Double32_t xJB) { mXJB = xJB; }
void StEpSimuJetEvent::setQ2JB(Double32_t q2JB) { mQ2JB = q2JB; }
void StEpSimuJetEvent::setYJB(Double32_t yJB) { mYJB = yJB; }
void StEpSimuJetEvent::setW2JB(Double32_t w2JB) { mW2JB = w2JB; }
void StEpSimuJetEvent::setXDA(Double32_t xDA) { mXDA = xDA; }
void StEpSimuJetEvent::setQ2DA(Double32_t q2DA) { mQ2DA = q2DA; }
void StEpSimuJetEvent::setYDA(Double32_t yDA) { mYDA = yDA; }
void StEpSimuJetEvent::setW2DA(Double32_t w2DA) { mW2DA = w2DA; }

// Implement Boost Setters
void StEpSimuJetEvent::setBoost(Double_t px, Double_t py, Double_t pz, Double_t e)
{
  mBoostPx = px;
  mBoostPy = py;
  mBoostPz = pz;
  mBoostE = e;
}
void StEpSimuJetEvent::setBoostAngles(Double_t phi, Double_t theta)
{
  mBoostPhi = phi;
  mBoostTheta = theta;
}
void StEpSimuJetEvent::setGoodBoost(bool isGoodBoost)
{
  mGoodBoost = isGoodBoost;
}


// Get Number of Total Particles
int StEpSimuJetEvent::numberOfTotalParticles() const { return mParticles->GetEntriesFast(); }

// Get Number of Jet Definitions
int StEpSimuJetEvent::numberOfJetDefs() const { return mJetDefs->GetEntriesFast(); }

// Get Jet Def
StEpSimuJetDef* StEpSimuJetEvent::jetDef(int i) const { return (StEpSimuJetDef*)mJetDefs->At(i); }

// Get Full Particle
StEpSimuParticle* StEpSimuJetEvent::particle(int i) const { return (StEpSimuParticle*)mParticles->At(i); }


// Declare New Objects
StEpSimuJet* StEpSimuJetEvent::newJet(Double_t pt, Double_t eta, Double_t phi, Double_t energy) 
{ 
  return new ((*mJets)[mJets->GetEntriesFast()]) StEpSimuJet(pt, eta, phi, energy); 
}
StEpSimuJetDef* StEpSimuJetEvent::newJetDef()
{
  return new ((*mJetDefs)[mJetDefs->GetEntriesFast()]) StEpSimuJetDef;
}
StEpSimuParticle* StEpSimuJetEvent::newParticle()
{
  return new ((*mParticles)[mParticles->GetEntriesFast()]) StEpSimuParticle;
}
StEpSimuJetParticle* StEpSimuJetEvent::newJetParticle(Double_t pt, Double_t eta, Double_t phi, Double_t energy)
{
  return new ((*mJetParticles)[mJetParticles->GetEntriesFast()]) StEpSimuJetParticle(pt, eta, phi, energy);
}
StEpSimuSubJet* StEpSimuJetEvent::newSubJet(Double_t pt, Double_t eta, Double_t phi, Double_t energy)
{
  return new ((*mSubJets)[mSubJets->GetEntriesFast()]) StEpSimuSubJet(pt, eta, phi, energy);
}



