
#include "fastjet/ClusterSequence.hh"
//#include "fastjet/ClusterSequenceArea.hh"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include "/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker/StEpSimuJetEvent.h"
#include "/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker/StEpSimuJetDef.h"
#include "/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker/StEpSimuJet.h"
#include "/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker/StEpSimuParticle.h"
#include "/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker/StEpSimuJetParticle.h"
#include "/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker/StEpSimuSubJet.h"

#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TRandom3.h"

#include "/afs/rhic.bnl.gov/eic/PACKAGES/EicRoot/eic-smear/include/eicsmear/erhic/EventBase.h"
#include "/afs/rhic.bnl.gov/eic/PACKAGES/EicRoot/eic-smear/include/eicsmear/erhic/EventPythia.h"
#include "/afs/rhic.bnl.gov/eic/PACKAGES/EicRoot/eic-smear/include/eicsmear/erhic/Particle.h"
#include "/afs/rhic.bnl.gov/eic/PACKAGES/EicRoot/eic-smear/include/eicsmear/smear/EventSmear.h"
#include "/afs/rhic.bnl.gov/eic/PACKAGES/EicRoot/eic-smear/include/eicsmear/smear/EventS.h"
#include "/afs/rhic.bnl.gov/eic/PACKAGES/EicRoot/eic-smear/include/eicsmear/smear/Smear.h"
#include "/afs/rhic.bnl.gov/eic/PACKAGES/EicRoot/eic-smear/include/eicsmear/smear/ParticleMCS.h"

using namespace fastjet;
using namespace std;

PseudoJet operator-(const PseudoJet &p)
{
  return PseudoJet(-p.px(), -p.py(), -p.pz(), p.E());
}

//-------------------------------------------------------------------
// Lorentz transformations (Use For Breit Frame)
//-------------------------------------------------------------------
// boost the vector p with the boost given by (bx,by,bz)
PseudoJet boost(const PseudoJet p, 
		const double &bx, const double &by, const double &bz){
  double b2 = bx*bx + by*by + bz*bz;
  //assert(b2 < 1.0);
  if(b2 >= 1.0)
    {
      b2 = 0.999;
      //cout << "Bad in boost" << endl;
      //cout << bx << " " << by << " " << bz << endl;
      //cout << " " << endl;
    }
  //cout << "Matt's b2 = " << b2 << endl;
  double gamma = 1.0/sqrt(1.0 - b2);
  double bp = bx*p.px() + by*p.py() + bz*p.pz();
  double gamma2 = (b2 > 0.0 ? (gamma - 1.0)/b2 : 0.0);

  return PseudoJet(p.px() + gamma2*bp*bx + gamma*bx*p.E(),
		   p.py() + gamma2*bp*by + gamma*by*p.E(),
		   p.pz() + gamma2*bp*bz + gamma*bz*p.E(),
		   gamma*(p.E() + bp));
}

// boost the vector p with the boost given by b (i.e. (bx/bE,by/bE,bz/bE))
PseudoJet boost(const PseudoJet p, const PseudoJet b){
  return boost(p, b.px()/b.E(), b.py()/b.E(), b.pz()/b.E());
}

// rotation around the x axis
PseudoJet rotateX(const PseudoJet p, const double& psi){
  double sp = sin(psi);
  double cp = cos(psi);
  return PseudoJet(p.px(),
		   cp*p.py()-sp*p.pz(),
		   sp*p.py()+cp*p.pz(),
		   p.E());
}

// rotation around the y axis
PseudoJet rotateY(const PseudoJet p, const double& psi){
  double sp = sin(psi);
  double cp = cos(psi);
  return PseudoJet(cp*p.px()+sp*p.pz(),
		   p.py(),
		   cp*p.pz()-sp*p.px(),
		   p.E());
}

// rotation around the z axis
PseudoJet rotateZ(const PseudoJet p, const double& psi){
  double sp = sin(psi);
  double cp = cos(psi);
  return PseudoJet(cp*p.px()-sp*p.py(),
		   sp*p.px()+cp*p.py(),
		   p.pz(),
		   p.E());
}


int main(int argc, char* argv[]) {

  const int nevents = atoi(argv[1]);
  const int beginEvent = atoi(argv[2]);
  const int endEvent = atoi(argv[3]);
  const int keepNeutrals = atoi(argv[4]);
  const double trackEff = atof(argv[5]);
  const char* inFileName1 = argv[6];
  const char* inFileName2 = argv[7];
  const char* outFileName = argv[8];

  cout << "nevents = " << nevents << endl;
  cout << "Process Events " << beginEvent << " - " << endEvent << endl;
  cout << "inFileName1 (Unsmeared) = " << inFileName1 << endl;
  cout << "inFileName2 (Smeared) = " << inFileName2 << endl;
  cout << "outFileName = " << outFileName << endl;

  cout << endl;
  if(keepNeutrals == 0) cout << "####### Reject Neutral Hadrons for eta < 1.0 #######" << endl;
  if(keepNeutrals == 1) cout << "####### Accept Neutral Hadrons for eta < 1.0 #######" << endl;

  cout << endl;
  cout << "####### Reject " << (1.0-trackEff)*100.0 << "% of Accepted Tracks #######" << endl;

  // Chain for Simu Tree
  TChain* inTree = new TChain("EICTree");
  inTree->Add(inFileName1);

  inTree->AddFriend("Smeared",inFileName2);

  // Setup Input Event Buffer
  erhic::EventPythia* inEvent(NULL);
  Smear::Event* inEventS(NULL);
  inTree->SetBranchAddress("event",&inEvent);
  inTree->SetBranchAddress("eventS",&inEventS);

  // Open Output File
  TFile *ofile = TFile::Open(outFileName,"recreate");
  assert(ofile);

  // Create Jet Tree
  TTree *mTree = new TTree("simuJetsSmear","SimuJetTrees");

  // Set Jet Tree Structure
  StEpSimuJetEvent *event = 0;
  mTree->Branch("myJetsSmear","StEpSimuJetEvent",&event);

  TH1D *efficiency = new TH1D("efficiency","",5,0.,5.);

  //Counters
  Int_t badInputTrack = 0;
  Int_t unassignedTrack = 0;

  Int_t numRejectTracks = 0;
  Int_t numAcceptTracks = 0;
  Int_t numNeutrals = 0;

  // Set up random number for track efficiency
  TRandom3 r;

  // Loop Over Events in Simu Trees
  for(int iEvent=0; iEvent<nevents; iEvent++)
    {
      //Read Next Event
      //inTree->GetEntry(iEvent);
      if(inTree->GetEntry(iEvent) <=0) break;

      if(iEvent%10000 == 0) cout << "Event " << iEvent << endl;

      if((iEvent < beginEvent) || (iEvent >= endEvent)) continue;

      //cout << "Event " << iEvent << " Particles = " << inEvent->GetNTracks() << endl;

      event->Clear();

      // Event Level Quantities
      event->setGenEvent(inEvent->GetGenEvent());
      event->setTgtPartonX(inEvent->GetTargetPartonX());
      event->setBeamPartonX(inEvent->GetBeamPartonX());
      event->setBeamPartonTheta(inEvent->GetBeamPartonTheta());
      event->setLeptonPhi(inEvent->GetLeptonPhi());
      event->setF1(inEvent->GetF1());
      event->setSigmaRad(inEvent->GetSigmaRad());
      event->setTHat(inEvent->GetHardT());
      event->setUHat(inEvent->GetHardU());
      event->setQ2Hat(inEvent->GetHardQ2());
      event->setSigRadCor(inEvent->GetSigRadCor());
      event->setEBrems(inEvent->GetEBrems());
      event->setPhotonFlux(inEvent->GetPhotonFlux());
      event->setTrueY(inEvent->GetTrueY());
      event->setTrueQ2(inEvent->GetTrueQ2());
      event->setTrueX(inEvent->GetTrueX());
      event->setTrueW2(inEvent->GetTrueW2());
      event->setTrueNu(inEvent->GetTrueNu());
      event->setF2(inEvent->GetF2());
      event->setR(inEvent->GetR());
      event->setPt2Hat(inEvent->GetHardPt2());
      event->setSHat(inEvent->GetHardS());

      // Event MC Quantities
      event->setN(inEvent->GetN());
      event->setProcessID(inEvent->GetProcess());
      event->setNTracks(inEvent->GetNTracks());
      
      // Event DIS Quantities
      event->setX(inEventS->GetX());
      event->setQ2(inEventS->GetQ2());
      event->setY(inEventS->GetY());
      event->setW2(inEventS->GetW2());
      event->setNu(inEventS->GetNu());
      event->setXJB(inEventS->GetXJacquetBlondel());
      event->setQ2JB(inEventS->GetQ2JacquetBlondel());
      event->setYJB(inEventS->GetYJacquetBlondel());
      event->setW2JB(inEventS->GetW2JacquetBlondel());
      event->setXDA(inEventS->GetXDoubleAngle());
      event->setQ2DA(inEventS->GetQ2DoubleAngle());
      event->setYDA(inEventS->GetYDoubleAngle());
      event->setW2DA(inEventS->GetW2DoubleAngle());


      // Create Containers for Different Jet Defs
      //StEpSimuJetDef *akt10Pt = event->newJetDef();
      //StEpSimuJetDef *akt10Pt_hbv = event->newJetDef();
      //StEpSimuJetDef *akt10Pt_breitE = event->newJetDef();
      //StEpSimuJetDef *akt10Pt_breitJB = event->newJetDef();
      //StEpSimuJetDef *akt10Pt_breitDA = event->newJetDef();
      StEpSimuJetDef *akt10_breit = event->newJetDef();
      StEpSimuJetDef *akt10Pt_breit = event->newJetDef();

      // AKT BEAM
      //akt10Pt->setRadius(1.0);
      //akt10Pt->setMinPt(0.250);
      //akt10Pt->setAlgo(StEpSimuJetDef::AKT);
      //akt10Pt->setFrame(StEpSimuJetDef::BEAM);

      // AKT HBV
      //akt10Pt_hbv->setRadius(1.0);
      //akt10Pt_hbv->setMinPt(0.250);
      //akt10Pt_hbv->setAlgo(StEpSimuJetDef::AKT);
      //akt10Pt_hbv->setFrame(StEpSimuJetDef::GAMMA);

      // AKT Breit Electron Method X
      //akt10Pt_breitE->setRadius(1.0);
      //akt10Pt_breitE->setMinPt(0.250);
      //akt10Pt_breitE->setAlgo(StEpSimuJetDef::AKT);
      //akt10Pt_breitE->setFrame(StEpSimuJetDef::BREIT);

      // AKT Breit JB Method X
      //akt10Pt_breitJB->setRadius(1.0);
      //akt10Pt_breitJB->setMinPt(0.250);
      //akt10Pt_breitJB->setAlgo(StEpSimuJetDef::AKT);
      //akt10Pt_breitJB->setFrame(StEpSimuJetDef::BREIT);

      // AKT Breit DA Method X
      //akt10Pt_breitDA->setRadius(1.0);
      //akt10Pt_breitDA->setMinPt(0.250);
      //akt10Pt_breitDA->setAlgo(StEpSimuJetDef::AKT);
      //akt10Pt_breitDA->setFrame(StEpSimuJetDef::BREIT);

      // AKT Breit True X
      akt10_breit->setRadius(1.0);
      akt10_breit->setMinPt(0.250);
      akt10_breit->setAlgo(StEpSimuJetDef::AKT);
      akt10_breit->setFrame(StEpSimuJetDef::BREIT);

      akt10Pt_breit->setRadius(1.0);
      akt10Pt_breit->setMinPt(0.500);
      akt10Pt_breit->setAlgo(StEpSimuJetDef::AKT);
      akt10Pt_breit->setFrame(StEpSimuJetDef::BREIT);


      // Create Lab and Hadron Boson Jet Vectors
      //vector<PseudoJet> particlesPt;
      //vector<PseudoJet> particlesPt_hbv;
      //vector<PseudoJet> particlesPt_breitE;
      //vector<PseudoJet> particlesPt_breitJB;
      //vector<PseudoJet> particlesPt_breitDA;
      vector<PseudoJet> particles_breit;
      vector<PseudoJet> particlesPt_breit;


      // Set Up Boost to Breit Frame
      const Particle* part2 = inEvent->GetTrack(1);
      const Particle* part4 = inEvent->GetTrack(3);

      PseudoJet proton(part2->GetPx(),part2->GetPy(),part2->GetPz(),part2->GetE());
      PseudoJet gamma(part4->GetPx(),part4->GetPy(),part4->GetPz(),part4->GetE());

      //PseudoJet boost_vectorE = -(gamma + 2.0*inEvent->GetX()*proton);
      //PseudoJet boost_vectorJB = -(gamma + 2.0*inEvent->GetXJacquetBlondel()*proton);
      //PseudoJet boost_vectorDA = -(gamma + 2.0*inEvent->GetXDoubleAngle()*proton);
      PseudoJet boost_vector = -(gamma + 2.0*inEvent->GetTrueX()*proton);

      double bxL = boost_vector.px()/boost_vector.E(); // Check for valid boost
      double byL = boost_vector.py()/boost_vector.E(); 
      double bzL = boost_vector.pz()/boost_vector.E();

      if(bxL*bxL + byL*byL + bzL*bzL < 1.0) event->setGoodBoost(true);
      if(bxL*bxL + byL*byL + bzL*bzL >= 1.0) 
	{
	  event->setGoodBoost(false);
	  //cout << "Bad in Check: Event = " << iEvent << endl;
	  //cout << bxL << " " << byL << " " << bzL << endl;
	}

      //PseudoJet boosted_protonE = boost(proton,boost_vectorE);
      //Double_t phi_pE = boosted_protonE.phi();
      //Double_t theta_pE = TMath::ATan2(boosted_protonE.perp(),boosted_protonE.pz());

      //event->setBoost(boost_vectorE.px(),boost_vectorE.py(),boost_vectorE.pz(),boost_vectorE.e());
      //event->setBoostAngles(phi_pE,theta_pE);

      //PseudoJet boosted_protonJB = boost(proton,boost_vectorJB);
      //Double_t phi_pJB = boosted_protonJB.phi();
      //Double_t theta_pJB = TMath::ATan2(boosted_protonJB.perp(),boosted_protonJB.pz());

      //PseudoJet boosted_protonDA = boost(proton,boost_vectorDA);
      //Double_t phi_pDA = boosted_protonDA.phi();
      //Double_t theta_pDA = TMath::ATan2(boosted_protonDA.perp(),boosted_protonDA.pz());

      PseudoJet boosted_proton = boost(proton,boost_vector);
      Double_t phi_p = boosted_proton.phi();
      Double_t theta_p = TMath::ATan2(boosted_proton.perp(),boosted_proton.pz());

      event->setBoost(boost_vector.px(),boost_vector.py(),boost_vector.pz(),boost_vector.e());
      event->setBoostAngles(phi_p,theta_p);


      // Loop over Particles
      for(unsigned int j=0; j<inEventS->GetNTracks(); j++)
	{
	  const Smear::ParticleMCS* inParticleS = inEventS->GetTrack(j); // Smeared Particle

	  const Particle* inParticle = inEvent->GetTrack(j); // Unsmeared Particle

	  Int_t goodSmear = 1;

	  Double_t px = -999.;
	  Double_t py = -999.;
	  Double_t pz = -999.;
	  Double_t E = -999.;

	  if(inParticleS == NULL) // Particle was not smeared, use original
	    {	  
	      px = inParticle->GetPx();
	      py = inParticle->GetPy();
	      pz = inParticle->GetPz();
	      E = inParticle->GetE();
	    }
	  
	  if(inParticleS != NULL) // Particle was smeared
	    {
	      Double_t tmpX = inParticleS->GetPx();
	      Double_t tmpY = inParticleS->GetPy();
	      Double_t tmpZ = inParticleS->GetPz();
	      Double_t tmpE = inParticleS->GetE();
	      //Double_t tmpOrigE = inParticle->GetE();

	      if(TMath::Abs(tmpX) < 0.0001 || TMath::Abs(tmpY) < 0.0001 || TMath::Abs(tmpE) < 0.0001)
		{
		  goodSmear = 0;
		  badInputTrack++;
		  //cout << "Event " << iEvent << " Particle " << j << " Smeared Bad" << " " << tmpX << " " << tmpY << " " << tmpZ << " " << tmpE << endl;
		}

	      // EM particles have only their energy smeared, momentum is original
	      // Currently will treat clusters (pi0s and electrons) as massless
	      // Need to alter P so that P^2 = E^2
	      // Do for Neutral Hadrons too
	      Int_t localCode = TMath::Abs(inParticle->GetPdgCode());
	      if(localCode == 11 || localCode == 22 || localCode == 2112 || localCode == 130)
		{
		  Double_t factor2 = (tmpE*tmpE)/(tmpX*tmpX + tmpY*tmpY + tmpZ*tmpZ);

		  px = tmpX*TMath::Sqrt(factor2);
		  py = tmpY*TMath::Sqrt(factor2);
		  pz = tmpZ*TMath::Sqrt(factor2);
		  E = tmpE;
		}

	      // Charged Hadrons with 3.5 < |eta| < 4 will not be smeared by the tracking system
	      // Will be detected by the hadron calorimeter
	      // Need to alter P for these events same as for neutral hadrons above
	      if(localCode > 110 && localCode != 2112 && localCode != 130 && TMath::Abs(inParticle->GetEta()) > 3.5 && TMath::Abs(inParticle->GetEta()) <= 4.0)
		{
		  Double_t factor2 = (tmpE*tmpE)/(tmpX*tmpX + tmpY*tmpY + tmpZ*tmpZ);
		  
		  px = tmpX*TMath::Sqrt(factor2);
		  py = tmpY*TMath::Sqrt(factor2);
		  pz = tmpZ*TMath::Sqrt(factor2);
		  E = tmpE;
		}

	      // Charged Hadrons have only their momentum smeared, energy is original
	      // Currently will treat tracks as having pion mass
	      // Need to alter E so that E^2 = P^2 + mp^2
	      if(localCode > 110 && localCode != 2112 && localCode != 130 && TMath::Abs(inParticle->GetEta()) < 3.5)
		{
		  Double_t P2 = tmpX*tmpX + tmpY*tmpY + tmpZ*tmpZ;

		  px = tmpX;
		  py = tmpY;
		  pz = tmpZ;
		  E = TMath::Sqrt(P2 + 0.139570*0.139570);
		}
	    }
	  
	  if(px == -999. || py == -999. || pz == -999. || E == -999.)
	    {
	      goodSmear = 0;
	      unassignedTrack++;
	      //cout << "Event " << iEvent << " Particle " << j << " PDG " << inParticle->GetPdgCode() << " Unassigned" << " " << inParticle->GetPx() << " " << inParticle->GetPy() << " " << inParticle->GetPz() << endl;
	    }

	  // Define Hadron Boson Frame Kinematics
	  TLorentzVector hadBosVec = inParticle->Get4VectorInHadronBosonFrame();

	  Double_t px_hbv = hadBosVec.Px();
	  Double_t py_hbv = hadBosVec.Py();
	  Double_t pz_hbv = hadBosVec.Pz();
	  Double_t E_hbv = hadBosVec.E();

	  // Define Breit Frame Kinematics Using Electron Method
	  //fastjet::PseudoJet p_labE(px,py,pz,E);
	  //fastjet::PseudoJet boosted_particleE = boost(p_labE,boost_vectorE);
	  //boosted_particleE = rotateZ(boosted_particleE, -phi_pE);
	  //boosted_particleE = rotateY(boosted_particleE, -theta_pE);

	  // Define Breit Frame Kinematics Using JB Method
	  //fastjet::PseudoJet p_labJB(px,py,pz,E);
	  //fastjet::PseudoJet boosted_particleJB = boost(p_labJB,boost_vectorJB);
	  //boosted_particleJB = rotateZ(boosted_particleJB, -phi_pJB);
	  //boosted_particleJB = rotateY(boosted_particleJB, -theta_pJB);

	  // Define Breit Frame Kinematics Using DA Method
	  //fastjet::PseudoJet p_labDA(px,py,pz,E);
	  //fastjet::PseudoJet boosted_particleDA = boost(p_labDA,boost_vectorDA);
	  //boosted_particleDA = rotateZ(boosted_particleDA, -phi_pDA);
	  //boosted_particleDA = rotateY(boosted_particleDA, -theta_pDA);

	  // Define Breit Frame Kinematics Using True Kinematics
	  fastjet::PseudoJet p_lab(px,py,pz,E);
	  fastjet::PseudoJet boosted_particle = boost(p_lab,boost_vector);
	  boosted_particle = rotateZ(boosted_particle, -phi_p);
	  boosted_particle = rotateY(boosted_particle, -theta_p);

	  // Not Sure if set_user_index can be called twice on same object, so make a copy
	  //fastjet::PseudoJet boosted_particlePt = boost(p_lab,boost_vector);
	  //boosted_particlePt = rotateZ(boosted_particlePt, -phi_p);
	  //boosted_particlePt = rotateY(boosted_particlePt, -theta_p);

	  
	  // Count number of rejected particles
	  if(j>10 && inParticle->GetStatus() == 1 && inParticle->GetParentIndex() != 3 && goodSmear == 0)
	    {
	      if(inParticleS->GetEta() <= 4.0 && inParticleS->GetPt() >= 0.250)
		{
		  numRejectTracks++;
		  //cout << inParticle->GetPdgCode() << endl;
		}
	    }

	  // Select Particles for Jets
	  if(j>10 && inParticle->GetStatus() == 1 && inParticle->GetParentIndex() != 3 && goodSmear == 1)
	    {
	      // Create random number for track rejection
	      Double_t ranFactor = r.Uniform(0.0,1.0);

	      if(inParticleS == NULL) // Not Smeared
		{
		  if(inParticle->GetEta() <= 4.0 && inParticle->GetPt() >= 0.250)
		    {
		      Int_t localCode = TMath::Abs(inParticle->GetPdgCode());

		      numAcceptTracks++;
		      if(localCode == 2112 || localCode == 130) numNeutrals++;

		      int neutralHadronAccept = 1;
		      if((localCode == 2112 || localCode == 130) && TMath::Abs(inParticle->GetEta()) < 1.0)
			{
			  neutralHadronAccept = 0;
			}
		      if(keepNeutrals == 1) neutralHadronAccept = 1;

		      int trackAccept = 1;
		      if(localCode > 110 && localCode != 2112 && localCode != 130 && TMath::Abs(inParticle->GetEta()) <= 3.5)
			{
			  if(ranFactor > trackEff) trackAccept = 0;
			}

		      if(neutralHadronAccept == 1 && trackAccept == 1)
			{
			  //fastjet::PseudoJet pPt(px,py,pz,E);
			  //pPt.set_user_index(inParticle->GetIndex());
			  //particlesPt.push_back(pPt);
			  
			  //fastjet::PseudoJet pPt_hbv(px_hbv,py_hbv,pz_hbv,E_hbv);
			  //pPt_hbv.set_user_index(inParticle->GetIndex());
			  //particlesPt_hbv.push_back(pPt_hbv);
			  
			  //boosted_particlePt.set_user_index(inParticle->GetIndex());
			  //particlesPt_breit.push_back(boosted_particlePt);
			  
			  //boosted_particleE.set_user_index(inParticle->GetIndex());
			  //particlesPt_breitE.push_back(boosted_particleE);
			  
			  //boosted_particleJB.set_user_index(inParticle->GetIndex());
			  //particlesPt_breitJB.push_back(boosted_particleJB);
			  
			  //boosted_particleDA.set_user_index(inParticle->GetIndex());
			  //particlesPt_breitDA.push_back(boosted_particleDA);
			  
			  boosted_particle.set_user_index(inParticle->GetIndex());
			  particles_breit.push_back(boosted_particle);
			  if(inParticle->GetPt() >= 0.500) particlesPt_breit.push_back(boosted_particle);

			  //cout << "!!! Not Smeared " << inParticle->GetIndex() << " " << localCode << " " << inParticle->GetEta() << " !!!" << endl;
			}
		    }
		}

	      if(inParticleS != NULL) // Smeared
		{
		  if(inParticleS->GetEta() <= 4.0 && inParticleS->GetPt() >= 0.250)
		    {
		      Int_t localCode = TMath::Abs(inParticle->GetPdgCode());

		      numAcceptTracks++;
		      if(localCode == 2112 || localCode == 130) numNeutrals++;

		      int neutralHadronAccept = 1;
		      if((localCode == 2112 || localCode == 130) && TMath::Abs(inParticle->GetEta()) < 1.0)
			{
			  neutralHadronAccept = 0;
			}
		      if(keepNeutrals == 1) neutralHadronAccept = 1;

		      int trackAccept = 1;
		      if(localCode > 110 && localCode != 2112 && localCode != 130 && TMath::Abs(inParticle->GetEta()) <= 3.5)
			{
			  if(ranFactor > trackEff) trackAccept = 0;
			  efficiency->Fill(trackAccept);
			}

		      if(neutralHadronAccept == 1 && trackAccept == 1)
			{
			  //fastjet::PseudoJet pPt(px,py,pz,E);
			  //pPt.set_user_index(inParticle->GetIndex());
			  //particlesPt.push_back(pPt);
			  
			  //fastjet::PseudoJet pPt_hbv(px_hbv,py_hbv,pz_hbv,E_hbv);
			  //pPt_hbv.set_user_index(inParticle->GetIndex());
			  //particlesPt_hbv.push_back(pPt_hbv);
			  
			  //boosted_particlePt.set_user_index(inParticle->GetIndex());
			  //particlesPt_breit.push_back(boosted_particlePt);
			  
			  //boosted_particleE.set_user_index(inParticle->GetIndex());
			  //particlesPt_breitE.push_back(boosted_particleE);
			  
			  //boosted_particleJB.set_user_index(inParticle->GetIndex());
			  //particlesPt_breitJB.push_back(boosted_particleJB);
			  
			  //boosted_particleDA.set_user_index(inParticle->GetIndex());
			  //particlesPt_breitDA.push_back(boosted_particleDA);
			  
			  boosted_particle.set_user_index(inParticle->GetIndex());
			  particles_breit.push_back(boosted_particle);
			  if(inParticleS->GetPt() >= 0.500) particlesPt_breit.push_back(boosted_particle);
			}
		    }
		}
	    }
	  
	  // Store All Particles
	  StEpSimuParticle *myParticle = event->newParticle();
	  myParticle->setFourMom(px,py,pz,E);
	  myParticle->setFourMom_hb(px_hbv,py_hbv,pz_hbv,E_hbv);
	  myParticle->setFourMom_breit(boosted_particle.px(),boosted_particle.py(),boosted_particle.pz(),boosted_particle.e());
	  myParticle->setMass(inParticle->GetM());
	  myParticle->setIndex(inParticle->GetIndex());
	  myParticle->setStatus(inParticle->GetStatus());
	  myParticle->setPdgCode(inParticle->GetPdgCode());
	  myParticle->setParentIndex(inParticle->GetParentIndex());
	  myParticle->setNChildren(inParticle->GetNChildren());
	  myParticle->setChild1Index(inParticle->GetChild1Index());
	  myParticle->setChildNIndex(inParticle->GetChildNIndex());
	  myParticle->setZ(inParticle->GetZ());
	  myParticle->setXF(inParticle->GetXFeynman());
	  myParticle->setThetaVsGamma(inParticle->GetThetaVsGamma());
	  myParticle->setPtVsGamma(inParticle->GetPtVsGamma());
	  myParticle->setVertex(inParticle->GetVertex().X(),inParticle->GetVertex().Y(),inParticle->GetVertex().Z());
	}

      // Set Jet Definitions
      double R_10 = 1.0;
      JetDefinition jet_def_akt_10(antikt_algorithm,R_10);

      // Select Area Parameters
      //double maxrap = 5.5;
      //unsigned int n_repeat = 1; // was 3
      //double ghost_area = 0.01;

      //GhostedAreaSpec area_spec(maxrap,n_repeat,ghost_area);
      //AreaDefinition area_def(fastjet::active_area,area_spec);
      //AreaDefinition area_def(fastjet::passive_area,area_spec);

      // Run Clustering and Extract the Jets
      double ptmin = 1.0;

      // Cluster in Lab Frame
      //ClusterSequence csPt_akt_10(particlesPt, jet_def_akt_10);
      //ClusterSequenceArea csPt_akt_10(particlesPt, jet_def_akt_10, area_def);

      // Cluster in Hadron Boson Frame
      //ClusterSequence csPt_akt_10_hbv(particlesPt_hbv, jet_def_akt_10);
      //ClusterSequenceArea csPt_akt_10_hbv(particlesPt_hbv, jet_def_akt_10, area_def);

      // Cluster in Breit Frame Electron Method
      //ClusterSequence csPt_akt_10_breitE(particlesPt_breitE, jet_def_akt_10);
      //ClusterSequenceArea csPt_akt_10_breit(particlesPt_breit, jet_def_akt_10, area_def);

      // Cluster in Breit Frame JB Method
      //ClusterSequence csPt_akt_10_breitJB(particlesPt_breitJB, jet_def_akt_10);
      //ClusterSequenceArea csPt_akt_10_breit(particlesPt_breit, jet_def_akt_10, area_def);

      // Cluster in Breit Frame DA Method
      //ClusterSequence csPt_akt_10_breitDA(particlesPt_breitDA, jet_def_akt_10);
      //ClusterSequenceArea csPt_akt_10_breit(particlesPt_breit, jet_def_akt_10, area_def);

      // Cluster in Breit Frame True Kinematics
      ClusterSequence cs_akt_10_breit(particles_breit, jet_def_akt_10);
      ClusterSequence csPt_akt_10_breit(particlesPt_breit, jet_def_akt_10);
      //ClusterSequenceArea csPt_akt_10_breit(particlesPt_breit, jet_def_akt_10, area_def);


      // Lab Frame Jets
      //vector<PseudoJet> jetsPt_akt_10 = sorted_by_pt(csPt_akt_10.inclusive_jets(ptmin));

      // Hadron Boson Frame Jets
      //vector<PseudoJet> jetsPt_akt_10_hbv = sorted_by_pt(csPt_akt_10_hbv.inclusive_jets(ptmin));
      
      // Breit Frame Jets Electron Method
      //vector<PseudoJet> jetsPt_akt_10_breitE = sorted_by_pt(csPt_akt_10_breitE.inclusive_jets(ptmin));

      // Breit Frame Jets JB Method
      //vector<PseudoJet> jetsPt_akt_10_breitJB = sorted_by_pt(csPt_akt_10_breitJB.inclusive_jets(ptmin));

      // Breit Frame Jets DA Method
      //vector<PseudoJet> jetsPt_akt_10_breitDA = sorted_by_pt(csPt_akt_10_breitDA.inclusive_jets(ptmin));

      // Breit Frame Jets True Kinematics
      vector<PseudoJet> jets_akt_10_breit = sorted_by_pt(cs_akt_10_breit.inclusive_jets(ptmin));
      vector<PseudoJet> jetsPt_akt_10_breit = sorted_by_pt(csPt_akt_10_breit.inclusive_jets(ptmin));
      

      // Set SubJet Definitions
      //double R_Sub = 0.25;
      double R_Sub = 0.2;
      JetDefinition jet_def_sub(antikt_algorithm,R_Sub);

      /*
      // AKT 10 Pt Cut
      for(int i=0; i<jetsPt_akt_10.size(); i++)
	{
	  StEpSimuJet *jet = event->newJet(jetsPt_akt_10[i].pt(),jetsPt_akt_10[i].eta(),jetsPt_akt_10[i].phi(),jetsPt_akt_10[i].e());
	  akt10Pt->addJet(jet);
	  jet->setJetDef(akt10Pt);
	  //jet->setArea(jetsPt_akt_10[i].area());
	  //jet->setAreaError(jetsPt_akt_10[i].area_error());
	  jet->setArea(-1.0);
	  jet->setAreaError(-1.0);

	  vector<PseudoJet> constituents = jetsPt_akt_10[i].constituents();

	  for(int j=0; j<constituents.size(); j++)
	    {
	      StEpSimuJetParticle *localPart = event->newJetParticle(constituents[j].pt(),constituents[j].eta(),constituents[j].phi(),constituents[j].E());
	      jet->addJetParticle(localPart);
	      localPart->setJet(jet);

	      localPart->setIndex(constituents[j].user_index());
	    }

	  double subPtMin = 0.5;
	  ClusterSequence cs_Sub(constituents, jet_def_sub);
	  vector<PseudoJet> subjets = sorted_by_pt(cs_Sub.inclusive_jets(subPtMin));
	  for(int k=0; k<subjets.size(); k++)
	    {
	      StEpSimuSubJet *localSub = event->newSubJet(subjets[k].pt(),subjets[k].eta(),subjets[k].phi(),subjets[k].e());
	      jet->addSubJet(localSub);
	      localSub->setParentJet(jet);
	      localSub->setRadius(R_Sub);
	      localSub->setNumParticles(subjets[k].constituents().size());
	    }
	}
      */

      // AKT 10 HB Pt Cut
      /*
      for(int i=0; i<jetsPt_akt_10_hbv.size(); i++)
	{
	  StEpSimuJet *jet = event->newJet(jetsPt_akt_10_hbv[i].pt(),jetsPt_akt_10_hbv[i].eta(),jetsPt_akt_10_hbv[i].phi(),jetsPt_akt_10_hbv[i].e());
	  akt10Pt_hbv->addJet(jet);
	  jet->setJetDef(akt10Pt_hbv);
	  //jet->setArea(jetsPt_akt_10_hbv[i].area());
	  //jet->setAreaError(jetsPt_akt_10_hbv[i].area_error());
	  jet->setArea(-1.0);
	  jet->setAreaError(-1.0);

	  vector<PseudoJet> constituents = jetsPt_akt_10_hbv[i].constituents();

	  for(int j=0; j<constituents.size(); j++)
	    {
	      StEpSimuJetParticle *localPart = event->newJetParticle(constituents[j].pt(),constituents[j].eta(),constituents[j].phi(),constituents[j].E());
	      jet->addJetParticle(localPart);
	      localPart->setJet(jet);

	      localPart->setIndex(constituents[j].user_index());
	    }

	  double subPtMin = 0.5;
	  ClusterSequence cs_Sub(constituents, jet_def_sub);
	  vector<PseudoJet> subjets = sorted_by_pt(cs_Sub.inclusive_jets(subPtMin));
	  for(int k=0; k<subjets.size(); k++)
	    {
	      StEpSimuSubJet *localSub = event->newSubJet(subjets[k].pt(),subjets[k].eta(),subjets[k].phi(),subjets[k].e());
	      jet->addSubJet(localSub);
	      localSub->setParentJet(jet);
	      localSub->setRadius(R_Sub);
	      localSub->setNumParticles(subjets[k].constituents().size());
	    }
	}
      */
      /*
      // AKT 10 Breit Pt Cut Electron Method
      for(int i=0; i<jetsPt_akt_10_breitE.size(); i++)
	{
	  StEpSimuJet *jet = event->newJet(jetsPt_akt_10_breitE[i].pt(),jetsPt_akt_10_breitE[i].eta(),jetsPt_akt_10_breitE[i].phi(),jetsPt_akt_10_breitE[i].e());
	  akt10Pt_breitE->addJet(jet);
	  jet->setJetDef(akt10Pt_breitE);
	  //jet->setArea(jetsPt_akt_10_breit[i].area());
	  //jet->setAreaError(jetsPt_akt_10_breit[i].area_error());
	  jet->setArea(-1.0);
	  jet->setAreaError(-1.0);

	  vector<PseudoJet> constituents = jetsPt_akt_10_breitE[i].constituents();

	  for(int j=0; j<constituents.size(); j++)
	    {
	      StEpSimuJetParticle *localPart = event->newJetParticle(constituents[j].pt(),constituents[j].eta(),constituents[j].phi(),constituents[j].E());
	      jet->addJetParticle(localPart);
	      localPart->setJet(jet);

	      localPart->setIndex(constituents[j].user_index());
	    }

	  double subPtMin = 0.5;
	  ClusterSequence cs_Sub(constituents, jet_def_sub);
	  vector<PseudoJet> subjets = sorted_by_pt(cs_Sub.inclusive_jets(subPtMin));
	  for(int k=0; k<subjets.size(); k++)
	    {
	      StEpSimuSubJet *localSub = event->newSubJet(subjets[k].pt(),subjets[k].eta(),subjets[k].phi(),subjets[k].e());
	      jet->addSubJet(localSub);
	      localSub->setParentJet(jet);
	      localSub->setRadius(R_Sub);
	      localSub->setNumParticles(subjets[k].constituents().size());
	    }
	}
      */
      /*
      // AKT 10 Breit Pt Cut JB Method
      for(int i=0; i<jetsPt_akt_10_breitJB.size(); i++)
	{
	  StEpSimuJet *jet = event->newJet(jetsPt_akt_10_breitJB[i].pt(),jetsPt_akt_10_breitJB[i].eta(),jetsPt_akt_10_breitJB[i].phi(),jetsPt_akt_10_breitJB[i].e());
	  akt10Pt_breitJB->addJet(jet);
	  jet->setJetDef(akt10Pt_breitJB);
	  //jet->setArea(jetsPt_akt_10_breit[i].area());
	  //jet->setAreaError(jetsPt_akt_10_breit[i].area_error());
	  jet->setArea(-1.0);
	  jet->setAreaError(-1.0);

	  vector<PseudoJet> constituents = jetsPt_akt_10_breitJB[i].constituents();

	  for(int j=0; j<constituents.size(); j++)
	    {
	      StEpSimuJetParticle *localPart = event->newJetParticle(constituents[j].pt(),constituents[j].eta(),constituents[j].phi(),constituents[j].E());
	      jet->addJetParticle(localPart);
	      localPart->setJet(jet);

	      localPart->setIndex(constituents[j].user_index());
	    }

	  double subPtMin = 0.5;
	  ClusterSequence cs_Sub(constituents, jet_def_sub);
	  vector<PseudoJet> subjets = sorted_by_pt(cs_Sub.inclusive_jets(subPtMin));
	  for(int k=0; k<subjets.size(); k++)
	    {
	      StEpSimuSubJet *localSub = event->newSubJet(subjets[k].pt(),subjets[k].eta(),subjets[k].phi(),subjets[k].e());
	      jet->addSubJet(localSub);
	      localSub->setParentJet(jet);
	      localSub->setRadius(R_Sub);
	      localSub->setNumParticles(subjets[k].constituents().size());
	    }
	}

      // AKT 10 Breit Pt Cut DA Method
      for(int i=0; i<jetsPt_akt_10_breitDA.size(); i++)
	{
	  StEpSimuJet *jet = event->newJet(jetsPt_akt_10_breitDA[i].pt(),jetsPt_akt_10_breitDA[i].eta(),jetsPt_akt_10_breitDA[i].phi(),jetsPt_akt_10_breitDA[i].e());
	  akt10Pt_breitDA->addJet(jet);
	  jet->setJetDef(akt10Pt_breitDA);
	  //jet->setArea(jetsPt_akt_10_breit[i].area());
	  //jet->setAreaError(jetsPt_akt_10_breit[i].area_error());
	  jet->setArea(-1.0);
	  jet->setAreaError(-1.0);

	  vector<PseudoJet> constituents = jetsPt_akt_10_breitDA[i].constituents();

	  for(int j=0; j<constituents.size(); j++)
	    {
	      StEpSimuJetParticle *localPart = event->newJetParticle(constituents[j].pt(),constituents[j].eta(),constituents[j].phi(),constituents[j].E());
	      jet->addJetParticle(localPart);
	      localPart->setJet(jet);

	      localPart->setIndex(constituents[j].user_index());
	    }

	  double subPtMin = 0.5;
	  ClusterSequence cs_Sub(constituents, jet_def_sub);
	  vector<PseudoJet> subjets = sorted_by_pt(cs_Sub.inclusive_jets(subPtMin));
	  for(int k=0; k<subjets.size(); k++)
	    {
	      StEpSimuSubJet *localSub = event->newSubJet(subjets[k].pt(),subjets[k].eta(),subjets[k].phi(),subjets[k].e());
	      jet->addSubJet(localSub);
	      localSub->setParentJet(jet);
	      localSub->setRadius(R_Sub);
	      localSub->setNumParticles(subjets[k].constituents().size());
	    }
	}
      */

      // AKT 10 Breit True Kinematics
      for(unsigned int i=0; i<jets_akt_10_breit.size(); i++)
	{
	  StEpSimuJet *jet = event->newJet(jets_akt_10_breit[i].pt(),jets_akt_10_breit[i].eta(),jets_akt_10_breit[i].phi(),jets_akt_10_breit[i].e());
	  akt10_breit->addJet(jet);
	  jet->setJetDef(akt10_breit);
	  //jet->setArea(jetsPt_akt_10_breit[i].area());
	  //jet->setAreaError(jetsPt_akt_10_breit[i].area_error());
	  jet->setArea(-1.0);
	  jet->setAreaError(-1.0);

	  vector<PseudoJet> constituents = jets_akt_10_breit[i].constituents();

	  for(unsigned int j=0; j<constituents.size(); j++)
	    {
	      StEpSimuJetParticle *localPart = event->newJetParticle(constituents[j].pt(),constituents[j].eta(),constituents[j].phi(),constituents[j].E());
	      jet->addJetParticle(localPart);
	      localPart->setJet(jet);

	      localPart->setIndex(constituents[j].user_index());
	    }

	  double subPtMin = 0.5;
	  ClusterSequence cs_Sub(constituents, jet_def_sub);
	  vector<PseudoJet> subjets = sorted_by_pt(cs_Sub.inclusive_jets(subPtMin));
	  for(unsigned int k=0; k<subjets.size(); k++)
	    {
	      StEpSimuSubJet *localSub = event->newSubJet(subjets[k].pt(),subjets[k].eta(),subjets[k].phi(),subjets[k].e());
	      jet->addSubJet(localSub);
	      localSub->setParentJet(jet);
	      localSub->setRadius(R_Sub);
	      localSub->setNumParticles(subjets[k].constituents().size());
	    }
	}
      
      // AKT 10 Breit Pt Cut True Kinematics
      for(unsigned int i=0; i<jetsPt_akt_10_breit.size(); i++)
	{
	  StEpSimuJet *jet = event->newJet(jetsPt_akt_10_breit[i].pt(),jetsPt_akt_10_breit[i].eta(),jetsPt_akt_10_breit[i].phi(),jetsPt_akt_10_breit[i].e());
	  akt10Pt_breit->addJet(jet);
	  jet->setJetDef(akt10Pt_breit);
	  //jet->setArea(jetsPt_akt_10_breit[i].area());
	  //jet->setAreaError(jetsPt_akt_10_breit[i].area_error());
	  jet->setArea(-1.0);
	  jet->setAreaError(-1.0);

	  vector<PseudoJet> constituents = jetsPt_akt_10_breit[i].constituents();

	  for(unsigned int j=0; j<constituents.size(); j++)
	    {
	      StEpSimuJetParticle *localPart = event->newJetParticle(constituents[j].pt(),constituents[j].eta(),constituents[j].phi(),constituents[j].E());
	      jet->addJetParticle(localPart);
	      localPart->setJet(jet);

	      localPart->setIndex(constituents[j].user_index());
	    }

	  double subPtMin = 0.5;
	  ClusterSequence cs_Sub(constituents, jet_def_sub);
	  vector<PseudoJet> subjets = sorted_by_pt(cs_Sub.inclusive_jets(subPtMin));
	  for(unsigned int k=0; k<subjets.size(); k++)
	    {
	      StEpSimuSubJet *localSub = event->newSubJet(subjets[k].pt(),subjets[k].eta(),subjets[k].phi(),subjets[k].e());
	      jet->addSubJet(localSub);
	      localSub->setParentJet(jet);
	      localSub->setRadius(R_Sub);
	      localSub->setNumParticles(subjets[k].constituents().size());
	    }
	}

      
      mTree->Fill();
    }

  cout << "Number of Tracks with too small momentum components = " << badInputTrack << endl;
  cout << "Number of unassigned tracks = " << unassignedTrack << endl;

  cout << "Number of Tracks Smeared Out of Bounds = " << numRejectTracks << endl;
  cout << "Number of Accepted Tracks = " << numAcceptTracks << endl;
  cout << "Number of Neutrons and K0_L = " << numNeutrals << endl;
      
  ofile->Write();
  ofile->Close();
  
}  
