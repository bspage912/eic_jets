// EP Simu Jet Treese

#ifndef StEpSimuJetEvent_def
#define StEpSimuJetEvent_def

class TLorentzVector;
class TClonesArray;
class StEpSimuJetDef;
class StEpSimuJet;
class StEpSimuJetParticle;
class StEpSimuParticle;
class StEpSimuSubJet;

#include "TObject.h"
#include "TVector3.h"
#include "TDatime.h"


class StEpSimuJetEvent : public TObject {

 public:

  StEpSimuJetEvent();

  void Clear();

  // Declare Event Setters
  void setGenEvent(Int_t genEvn);
  void setTgtPartonX(Double32_t tgtX);
  void setBeamPartonX(Double32_t beamX);
  void setBeamPartonTheta(Double32_t beamTheta);
  void setLeptonPhi(Double32_t lepPhi);
  void setF1(Double32_t f1);
  void setSigmaRad(Double32_t sigRad);
  void setTHat(Double32_t tHat);
  void setUHat(Double32_t uHat);
  void setQ2Hat(Double32_t q2Hat);
  void setSigRadCor(Double32_t sigRadCor);
  void setEBrems(Double32_t eBrems);
  void setPhotonFlux(Double32_t photonFlux);
  void setTrueY(Double32_t truey);
  void setTrueQ2(Double32_t trueq2);
  void setTrueX(Double32_t truex);
  void setTrueW2(Double32_t truew);
  void setTrueNu(Double32_t truenu);
  void setF2(Double32_t f2);
  void setR(Double32_t r);
  void setPt2Hat(Double32_t pt2Hat);
  void setSHat(Double32_t sHat);

  // Declare Event MC Setters
  void setN(Int_t number);
  void setProcessID(Int_t process);
  void setNTracks(Int_t numTracks);

  // Declare Event DIS Setters
  void setX(Double32_t x);
  void setQ2(Double32_t q2);
  void setY(Double32_t y);
  void setW2(Double32_t w2);
  void setNu(Double32_t nu);
  void setXJB(Double32_t xJB);
  void setQ2JB(Double32_t q2JB);
  void setYJB(Double32_t yJB);
  void setW2JB(Double32_t w2JB);
  void setXDA(Double32_t xDA);
  void setQ2DA(Double32_t q2DA);
  void setYDA(Double32_t yDA);
  void setW2DA(Double32_t w2DA);

  // Declare Boost Parameter Setters
  void setBoost(Double_t px, Double_t py, Double_t pz, Double_t e);
  void setBoostAngles(Double_t phi, Double_t theta);
  void setGoodBoost(bool isGoodBoost);
  
  
  // Declare and Implement Event Getters
  Int_t genEvent()              const { return mGenEvent; }
  Double32_t targetPartonX()    const { return mTgtPartonX; }
  Double32_t beamPartonX()      const { return mBeamPartonX; }
  Double32_t beamPartonTheta()  const { return mBeamPartonTheta; }
  Double32_t leptonPhi()        const { return mLeptonPhi; }
  Double32_t F1()               const { return mF1; }
  Double32_t sigmaRad()         const { return mSigmaRad; }
  Double32_t tHat()             const { return mTHat; }
  Double32_t uHat()             const { return mUHat; }
  Double32_t q2Hat()            const { return mQ2Hat; }
  Double32_t sigRadCor()        const { return mSigRadCor; }
  Double32_t EBrems()           const { return mEBrems; }
  Double32_t photonFlux()       const { return mPhotonFlux; }
  Double32_t trueY()            const { return mTrueY; }
  Double32_t trueQ2()           const { return mTrueQ2; }
  Double32_t trueX()            const { return mTrueX; }
  Double32_t trueW2()           const { return mTrueW2; }
  Double32_t trueNu()           const { return mTrueNu; }
  Double32_t F2()               const { return mF2; }
  Double32_t R()                const { return mR; }
  Double32_t pt2Hat()           const { return mPt2Hat; }
  Double32_t sHat()             const { return mSHat; }

  // Declare and Implement Event MC Getters
  Int_t N()          const { return mN; }
  Int_t processID()  const { return mProcessID; }
  Int_t nTracks()    const { return mNTracks; }

  // Declare and Implement Event DIS Getters
  Double32_t x()  const { return mX; }
  Double32_t q2()  const { return mQ2; }
  Double32_t y()  const { return mY; }
  Double32_t w2()  const { return mW2; }
  Double32_t nu()  const { return mNu; }
  Double32_t xJB()  const { return mXJB; }
  Double32_t q2JB()  const { return mQ2JB; }
  Double32_t yJB()  const { return mYJB; }
  Double32_t w2JB()  const { return mW2JB; }
  Double32_t xDA()  const { return mXDA; }
  Double32_t q2DA()  const { return mQ2DA; }
  Double32_t yDA()  const { return mYDA; }
  Double32_t w2DA()  const { return mW2DA; }

  // Declare and Implement Boost Getters
  Double_t boostPx()    const { return mBoostPx; }
  Double_t boostPy()    const { return mBoostPy; }
  Double_t boostPz()    const { return mBoostPz; }
  Double_t boostE()     const { return mBoostE; }
  Double_t boostPhi()   const { return mBoostPhi; }
  Double_t boostTheta() const { return mBoostTheta; }
  bool goodBoost()      const { return mGoodBoost; }


  // Get Total Number of Particles
  int numberOfTotalParticles() const;

  // Get Number of Jet Definitions
  int numberOfJetDefs() const;

  // Get Jet Finder Configuration
  StEpSimuJetDef* jetDef(int i) const;

  // Get Full Particle
  StEpSimuParticle* particle(int i) const;

  // Declare New Objects
  StEpSimuJetDef* newJetDef();
  StEpSimuJet* newJet(Double_t pt, Double_t eta, Double_t phi, Double_t energy);
  StEpSimuParticle* newParticle();
  StEpSimuJetParticle* newJetParticle(Double_t pt, Double_t eta, Double_t phi, Double_t energy);
  StEpSimuSubJet* newSubJet(Double_t pt, Double_t eta, Double_t phi, Double_t energy);

  

 private:

  // Event Members
  Int_t mGenEvent;
  Double32_t mTgtPartonX;
  Double32_t mBeamPartonX;
  Double32_t mBeamPartonTheta;
  Double32_t mLeptonPhi;
  Double32_t mF1;
  Double32_t mSigmaRad;
  Double32_t mTHat;
  Double32_t mUHat;
  Double32_t mQ2Hat;
  Double32_t mSigRadCor;
  Double32_t mEBrems;
  Double32_t mPhotonFlux;
  Double32_t mTrueY;
  Double32_t mTrueQ2;
  Double32_t mTrueX;
  Double32_t mTrueW2;
  Double32_t mTrueNu;
  Double32_t mF2;
  Double32_t mR;
  Double32_t mPt2Hat;
  Double32_t mSHat;

  // Event MC Members
  Int_t mN;
  Int_t mProcessID;
  Int_t mNTracks;  

  // Event DIS Members
  Double32_t mX;
  Double32_t mQ2;
  Double32_t mY;
  Double32_t mW2;
  Double32_t mNu;
  Double32_t mXJB;
  Double32_t mQ2JB;
  Double32_t mYJB;
  Double32_t mW2JB;
  Double32_t mXDA;
  Double32_t mQ2DA;
  Double32_t mYDA;
  Double32_t mW2DA;

  // Boost Members
  Double_t mBoostPx, mBoostPy, mBoostPz, mBoostE;
  Double_t mBoostPhi, mBoostTheta;
  bool mGoodBoost;

  // Branches
  TClonesArray* mJetDefs;
  TClonesArray* mJets;
  TClonesArray* mJetParticles;
  TClonesArray* mParticles;
  TClonesArray* mSubJets;

  ClassDef(StEpSimuJetEvent,3);
};

#endif

