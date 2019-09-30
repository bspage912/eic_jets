// Container for all Particles

#ifndef StEpSimuParticle_def
#define StEpSimuParticle_def

#include "TObject.h"
#include "TRefArray.h"
#include "TVector3.h"
#include "TLorentzVector.h"
//#include "TParticle.h"
#include "TClonesArray.h"

class StEpSimuJetEvent;

class StEpSimuParticle : public TObject {
 public:

  StEpSimuParticle();

  TVector3 momentum() const;
  TLorentzVector fourMomentum() const;

  TVector3 momentum_hb() const;
  TLorentzVector fourMomentum_hb() const;

  TVector3 momentum_breit() const;
  TLorentzVector fourMomentum_breit() const;

  // Declare Setters
  void setFourMom(Double_t px, Double_t py, Double_t pz, Double_t e);
  void setFourMom_hb(Double_t px, Double_t py, Double_t pz, Double_t e); // Hadron Boson Frame
  void setFourMom_breit(Double_t px, Double_t py, Double_t pz, Double_t e); // Breit Frame
  void setMass(Double_t mass);
  void setIndex(int id);
  void setStatus(int status);
  void setPdgCode(int code);
  void setParentIndex(int id);
  void setNChildren(int n);
  void setChild1Index(int id);
  void setChildNIndex(int id);
  void setZ(Double_t z);
  void setXF(Double_t xf);
  void setThetaVsGamma(Double_t theta);
  void setPtVsGamma(Double_t pt);
  void setVertex(Double_t x, Double_t y, Double_t z);
  

  // Declare and Implement Getters
  Double_t pt();
  Double_t eta();
  Double_t phi();
  Double_t rap();

  // Return Kinematics in Hadron Boson Frame
  Double_t pt_hb();
  Double_t eta_hb();
  Double_t phi_hb();
  Double_t rap_hb();

  // Return Kinematics in Breit Frame
  Double_t pt_breit();
  Double_t eta_breit();
  Double_t phi_breit();
  Double_t rap_breit();

  Double_t mass()         { return mMass; }
  int index()             { return mIndex; }
  int status()            { return mStatus; }
  int pdgCode()           { return mPDG; }
  int parentIndex()       { return mParentIndex; }
  int numChildren()       { return mNChildren; }
  int firstChild()        { return mFirstChild; }
  int lastChild()         { return mLastChild; }
  Double_t z()            { return mZ; }
  Double_t xF()           { return mXF; }
  Double_t thetaVsGamma() { return mThetaVsGamma; }
  Double_t ptVsGamma()    { return mPtVsGamma; }
  Double_t vx()           { return mVx; }
  Double_t vy()           { return mVy; }
  Double_t vz()           { return mVz; }


 private:

  Double_t mPx, mPy, mPz, mE;
  Double_t mPx_hb, mPy_hb, mPz_hb, mE_hb;
  Double_t mPx_breit, mPy_breit, mPz_breit, mE_breit;
  Double_t mMass;
  int mIndex;
  int mStatus;
  int mPDG;
  int mParentIndex;
  int mNChildren;
  int mFirstChild;
  int mLastChild;
  Double_t mZ;
  Double_t mXF;
  Double_t mThetaVsGamma;
  Double_t mPtVsGamma;
  Double_t mVx, mVy, mVz;

  ClassDef(StEpSimuParticle,3);
};

inline TVector3 StEpSimuParticle::momentum() const
{
  TVector3 mom;
  mom.SetXYZ(mPx,mPy,mPz);
  return mom;
}

inline TLorentzVector StEpSimuParticle::fourMomentum() const
{
  TLorentzVector fourMom;
  fourMom.SetPxPyPzE(mPx,mPy,mPz,mE);
  return fourMom;
}

inline TVector3 StEpSimuParticle::momentum_hb() const
{
  TVector3 mom;
  mom.SetXYZ(mPx_hb,mPy_hb,mPz_hb);
  return mom;
}

inline TLorentzVector StEpSimuParticle::fourMomentum_hb() const
{
  TLorentzVector fourMom;
  fourMom.SetPxPyPzE(mPx_hb,mPy_hb,mPz_hb,mE_hb);
  return fourMom;
}

inline TVector3 StEpSimuParticle::momentum_breit() const
{
  TVector3 mom;
  mom.SetXYZ(mPx_breit,mPy_breit,mPz_breit);
  return mom;
}

inline TLorentzVector StEpSimuParticle::fourMomentum_breit() const
{
  TLorentzVector fourMom;
  fourMom.SetPxPyPzE(mPx_breit,mPy_breit,mPz_breit,mE_breit);
  return fourMom;
}

#endif
