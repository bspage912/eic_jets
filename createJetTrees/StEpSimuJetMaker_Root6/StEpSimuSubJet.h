// Test of Jet Used with StForwardDijetEvent

#ifndef StEpSimuSubJet_def
#define StEpSimuSubJet_def

#include "TObject.h"
#include "TLorentzVector.h"
#include "TRef.h"
#include "TRefArray.h"

#include "StEpSimuJetParticle.h"

class StEpSimuJetEvent;
class StEpSimuJetDef;
class StEpSimuJet;

class StEpSimuSubJet : public TObject {
 public:
  StEpSimuSubJet()
    : mPt(0)
    , mEta(0)
    , mPhi(0)
    , mE(0)
    {
    }

  StEpSimuSubJet(Double_t pt, Double_t eta, Double_t phi, Double_t energy);

  //StEpSimuJet(const TLorentzVector& fourMomentum);

  TVector3 momentum() const;
  TLorentzVector fourMomentum() const;
  
  // Declare Setters
  void setRadius(Double_t rad);
  void setNumParticles(Double_t part);
  
  // Declare and Implement Getters
  Double_t pt()                const { return mPt; }
  Double_t eta()               const { return mEta; }
  Double_t phi()               const { return mPhi; }
  Double_t E()                 const { return mE; }
  Double_t radius()            const { return mRadius; }
  Double_t numberOfParticles() const { return mParticles; }
  //Double_t area()       const { return mArea; }
  
  Double_t rap();

  
  // Set Parent Jet
  void setParentJet(const StEpSimuJet* jet); //{ mDijet = (TObject*)dijet; }
  
  
 private:
  
  TRef mParentJets;
  
  Double_t mPt;
  Double_t mEta;
  Double_t mPhi;
  Double_t mE;
  //Double_t mArea;
  Double_t mRadius;
  Double_t mParticles;
  
  
  ClassDef(StEpSimuSubJet,3);
};


inline TVector3 StEpSimuSubJet::momentum() const
{
  TVector3 mom;
  mom.SetPtEtaPhi(mPt,mEta,mPhi);
  return mom;
}

inline TLorentzVector StEpSimuSubJet::fourMomentum() const
{
  TLorentzVector fourMom;
  fourMom.SetPtEtaPhiE(mPt,mEta,mPhi,mE);
  return fourMom;
}

#endif
