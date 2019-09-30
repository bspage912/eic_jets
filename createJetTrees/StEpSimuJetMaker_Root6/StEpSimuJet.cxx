// Test of Jet Container

#include "StEpSimuJet.h"

ClassImp(StEpSimuJet);

StEpSimuJet::StEpSimuJet(Double_t pt, Double_t eta, Double_t phi, Double_t E)
  : mPt(pt)
  , mEta(eta)
  , mPhi(phi)
  , mE(E)
{
}

// Setters
void StEpSimuJet::setArea(Double_t area)
{
  mArea = area;
}

void StEpSimuJet::setAreaError(Double_t error)
{
  mAreaError = error;
}

// Getters
Double_t StEpSimuJet::rap()
{
  return StEpSimuJet::fourMomentum().Rapidity();
}


// Implement Utility Functions
Double_t StEpSimuJet::sumParticlePt() const
{
  //float s = 0.;
  TVector3 mom;
  for(int i=0; i<numberOfJetParticles(); ++i)
    {
      const StEpSimuJetParticle* t = jetParticle(i);
      //s += t->pt();
      mom += t->momentum();
    }
  //return s;
  return mom.Pt();
}

Double_t StEpSimuJet::sumParticlePt(Double_t radius) const
{
  //float s = 0.;
  TVector3 mom;
  for(int i=0; i<numberOfJetParticles(); ++i)
    {
      const StEpSimuJetParticle* t = jetParticle(i);
      //if(deltaR(t) < radius) s += t->pt();
      if(deltaR(t) < radius) mom += t->momentum();
    }
  //return s;
  return mom.Pt();
}

Double_t StEpSimuJet::deltaR(const StEpSimuJetParticle* part) const
{
  Double_t jetRap = fourMomentum().Rapidity();
  Double_t partRap = part->fourMomentum().Rapidity();

  Double_t jetPhi = StEpSimuJet::phi();
  Double_t partPhi = part->phi();

  Double_t dRap = jetRap-partRap;
  Double_t dPhi = TVector2::Phi_mpi_pi(jetPhi-partPhi);

  return TMath::Sqrt( dRap*dRap + dPhi*dPhi );
}

// Set Parent Jet Definition
void StEpSimuJet::setJetDef(const StEpSimuJetDef* jetDef) { mJetDefs = (TObject*)jetDef; }

// Get Number of Particles
int StEpSimuJet::numberOfJetParticles() const { return mJetParticles.GetEntriesFast(); }
//int StForwardJet::numberOfTowers() const { return mTowers.GetEntriesFast(); }

// Get Number of Subjets
int StEpSimuJet::numberOfSubJets() const { return mSubJets.GetEntriesFast(); }

// Get Particle
StEpSimuJetParticle* StEpSimuJet::jetParticle(int i) const { return (StEpSimuJetParticle*)mJetParticles.At(i); }
//StForwardTower* StForwardJet::tower(int i) const { return (StForwardTower*)mTowers.At(i); }

// Get Subjet
StEpSimuSubJet* StEpSimuJet::subJet(int i) const { return (StEpSimuSubJet*)mSubJets.At(i); }

// Add Particle to Jet
StEpSimuJetParticle* StEpSimuJet::addJetParticle(StEpSimuJetParticle* particle) { mJetParticles.Add((TObject*)particle); return (StEpSimuJetParticle*)mJetParticles.Last(); }

// Add Subjet to Jet
StEpSimuSubJet* StEpSimuJet::addSubJet(StEpSimuSubJet* subjet) { mSubJets.Add((TObject*)subjet); return (StEpSimuSubJet*)mSubJets.Last(); }
//StForwardTower* StForwardJet::addTower(StForwardTower* tower) { mTowers.Add((TObject*)tower); return (StForwardTower*)mTowers.Last(); }
