// Test of Track Used with StForwardJet

#ifndef StEpSimuJetParticle_def
#define StEpSimuJetParticle_def

#include "TObject.h"
#include "TLorentzVector.h"
#include "TRef.h"
#include "TVector3.h"

class StEpSimuJetEvent;
class StEpSimuJetDef;
class StEpSimuJet;

class StEpSimuJetParticle : public TObject { 

 public:
  StEpSimuJetParticle()
    : mPt(0)
    , mEta(0)
    , mPhi(0)
    , mE(0)
    {
    }

  StEpSimuJetParticle(Double_t pt, Double_t eta, Double_t phi, Double_t energy);

  TVector3 momentum() const;
  TLorentzVector fourMomentum() const;

  // Declare Setters
  void setIndex(int id);
  

  // Declare and Implement Getters
  Double_t pt()                      const { return mPt; }
  Double_t eta()                     const { return mEta; }
  Double_t phi()                     const { return mPhi; }
  Double_t energy()                  const { return mE; }
  Double_t index()                   const { return mIndex; }

  Double_t rap();

  // Set Jet
  void setJet(const StEpSimuJet* jet) {mJets = (TObject*)jet; }

 private:

  TRef mJets;

  Double_t    mPt;
  Double_t    mEta;
  Double_t    mPhi;
  Double_t    mE;

  int mIndex;

  ClassDef(StEpSimuJetParticle,3);
};

inline TVector3 StEpSimuJetParticle::momentum() const
{
  TVector3 mom;
  mom.SetPtEtaPhi(mPt,mEta,mPhi);
  return mom;
}

inline TLorentzVector StEpSimuJetParticle::fourMomentum() const
{
  TLorentzVector fourMom;
  fourMom.SetPtEtaPhiE(mPt,mEta,mPhi,mE);
  return fourMom;
}

#endif
