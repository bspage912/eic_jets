// Analyze the My Simu Jet Trees

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <map>
#include <vector>

using namespace std;

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TStopwatch.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"
#include "TVectorD.h"

//#include "LHAPDF/LHAPDF.h"

#include "/eic/u/bpage/epJets/createJetTrees/testComp/StEpSimuJetMaker/StEpSimuJetEvent.h"
#include "/eic/u/bpage/epJets/createJetTrees/testComp/StEpSimuJetMaker/StEpSimuJetDef.h"
#include "/eic/u/bpage/epJets/createJetTrees/testComp/StEpSimuJetMaker/StEpSimuJet.h"
#include "/eic/u/bpage/epJets/createJetTrees/testComp/StEpSimuJetMaker/StEpSimuParticle.h"
#include "/eic/u/bpage/epJets/createJetTrees/testComp/StEpSimuJetMaker/StEpSimuJetParticle.h"
#include "/eic/u/bpage/epJets/createJetTrees/testComp/StEpSimuJetMaker/StEpSimuSubJet.h"


//-------------------------------------------------------------------
// Lorentz transformations (Boost to CM Frame)
//-------------------------------------------------------------------
TLorentzRotation computeBoost(const TLorentzVector& rest, const TLorentzVector* z) {
  TLorentzRotation toRest(-(rest.BoostVector()));
  if(z) 
    {
      TRotation rotate;
      TLorentzVector boostedZ(*z);
      boostedZ *= toRest;
      rotate.SetZAxis(boostedZ.Vect());
      // We need the rotation of the frame, so take the inverse.
      // See the TRotation documentation for more explanation.
      rotate = rotate.Inverse();
      toRest.Transform(rotate);
    } // if
  return toRest;
}


//-------------------------------------------------------------------
// Lorentz transformations (Use For Breit Frame)
//-------------------------------------------------------------------
// boost the vector p with the boost given by (bx,by,bz)
TLorentzVector boost(const TLorentzVector p, 
		const double &bx, const double &by, const double &bz){
  double b2 = bx*bx + by*by + bz*bz;
  assert(b2 < 1.0);
  //cout << "Matt's b2 = " << b2 << endl;
  double gamma = 1.0/sqrt(1.0 - b2);
  double bp = bx*p.Px() + by*p.Py() + bz*p.Pz();
  double gamma2 = (b2 > 0.0 ? (gamma - 1.0)/b2 : 0.0);

  //return PseudoJet(p.px() + gamma2*bp*bx + gamma*bx*p.E(),
  //p.py() + gamma2*bp*by + gamma*by*p.E(),
  //p.pz() + gamma2*bp*bz + gamma*bz*p.E(),
  //gamma*(p.E() + bp));

  TLorentzVector out;
  out.SetXYZT(p.Px()+gamma2*bp*bx+gamma*bx*p.E(),p.Py()+gamma2*bp*by+gamma*by*p.E(),p.Pz()+gamma2*bp*bz+gamma*bz*p.E(),gamma*(p.E()+bp));
  return out;
}

// boost the vector p with the boost given by b (i.e. (bx/bE,by/bE,bz/bE))
TLorentzVector boost(const TLorentzVector p, const TLorentzVector b){
  return boost(p, b.Px()/b.E(), b.Py()/b.E(), b.Pz()/b.E());
}

// rotation around the y axis
TLorentzVector rotateY(const TLorentzVector p, const double& psi){
  double sp = sin(psi);
  double cp = cos(psi);
  //return PseudoJet(cp*p.px()+sp*p.pz(),
  //p.py(),
  //cp*p.pz()-sp*p.px(),
  //p.E());
  TLorentzVector out;
  out.SetXYZT(cp*p.Px()+sp*p.Pz(),p.Py(),cp*p.Pz()-sp*p.Px(),p.E());
  return out;
}

// rotation around the z axis
TLorentzVector rotateZ(const TLorentzVector p, const double& psi){
  double sp = sin(psi);
  double cp = cos(psi);
  //return PseudoJet(cp*p.px()-sp*p.py(),
  //sp*p.px()+cp*p.py(),
  //p.pz(),
  //p.E());
  TLorentzVector out;
  out.SetXYZT(cp*p.Px()-sp*p.Py(),sp*p.Px()+cp*p.Py(),p.Pz(),p.E());
  return out;
}


Double_t calcDijetMass(StEpSimuJet* hiJet, StEpSimuJet* loJet)
{
  // Calculate Invariant Mass of Dijet System
  Double_t hiRapidity = hiJet->rap();
  Double_t loRapidity = loJet->rap();
  
  Double_t highPtDiJetMass = hiJet->fourMomentum().M();
  Double_t highPtDiJetMass2 = highPtDiJetMass*highPtDiJetMass;
  
  Double_t lowPtDiJetMass = loJet->fourMomentum().M();
  Double_t lowPtDiJetMass2 = lowPtDiJetMass*lowPtDiJetMass;
  
  Double_t highPtDiJetTMass = TMath::Sqrt(highPtDiJetMass2 + (hiJet->pt()*hiJet->pt()));
  Double_t lowPtDiJetTMass = TMath::Sqrt(lowPtDiJetMass2 + (loJet->pt()*loJet->pt()));
  
  Double_t term1 = 2*highPtDiJetTMass*lowPtDiJetTMass*(TMath::CosH(hiRapidity-loRapidity));
  Double_t term2 = 2*(hiJet->pt())*(loJet->pt())*(TMath::Cos(hiJet->phi()-loJet->phi()));
  
  Double_t invMass = TMath::Sqrt(highPtDiJetMass2 + lowPtDiJetMass2 + term1 - term2);
  
  return invMass;
}


Double_t calcDeltaR(StEpSimuJet *jet1, StEpSimuJet *jet2)
{
  //Double_t deltaEta = part->eta_breit() - jet->eta();
  //Double_t deltaPhi = TVector2::Phi_mpi_pi(part->phi_breit() - TVector2::Phi_mpi_pi(jet->phi()));

  Double_t deltaEta = jet1->eta() - jet2->eta();
  Double_t deltaPhi = TVector2::Phi_mpi_pi(jet1->phi() - TVector2::Phi_mpi_pi(jet2->phi()));

  return TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
}


Double_t calcDeltaR(StEpSimuParticle *part, StEpSimuJet *jet)
{
  Double_t deltaEta = part->eta_breit() - jet->eta();
  Double_t deltaPhi = TVector2::Phi_mpi_pi(part->phi_breit() - TVector2::Phi_mpi_pi(jet->phi()));

  //Double_t deltaEta = part->eta() - jet->eta();
  //Double_t deltaPhi = TVector2::Phi_mpi_pi(part->phi() - TVector2::Phi_mpi_pi(jet->phi()));

  return TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
}


void jetPartonAssoc(StEpSimuJetEvent *event, StEpSimuJet *jet, Double_t *deltaR, Int_t *index)
{
  // Generating Particles
  StEpSimuParticle *partOne = event->particle(9);
  StEpSimuParticle *partTwo = event->particle(10);

  // Find Particle Matched to Jet
  Double_t partOneDeltaR = calcDeltaR(partOne,jet);
  Double_t partTwoDeltaR = calcDeltaR(partTwo,jet);

  Double_t dR = 999.0;
  int code = -999;

  if(partOneDeltaR < partTwoDeltaR)
    {
      dR = partOneDeltaR;
      //code = partOne->pdgCode();
      code = partOne->index();
    }

  if(partOneDeltaR > partTwoDeltaR)
    {
      dR = partTwoDeltaR;
      //code = partTwo->pdgCode();
      code = partTwo->index();
    }

  *deltaR = dR;
  *index = code;
}




// Main Program
int main(int argc, char* argv[])
{
  if(argc != 4) // 4
    {
      cerr << "WRONG" << endl;
      exit(EXIT_FAILURE);
    }

  const int q2Set = atoi(argv[1]);
  const char* infile = argv[2];
  const char* outfile = argv[3];

  //const int q2Set = 0;
  //const char* infile = "/gpfs02/eic/bpage/jetTrees_partonRad/jetTree.ep.20x250.5Mevents.RadCor=0_1.T.root";
  //const char* outfile = "t.hist.root";
  //const char* infile = "/eicdata/eic0009/bpage/jetTrees_akt10_v3/Q2_10-100/jetTree.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_4.C.root";
  //const char* infile="/eicdata/eic0009/bpage/jetTrees_akt10_v3/lowQ2/jetTree.ep.20x250.1Mevents.RadCor=0.Q2=0.01-0.1.kT=1.0_4.A.root";
  //const char* outfile = "test.hist.root";

  // Open Tree
  TChain* jetChain = new TChain("simuJets");
  jetChain->Add(infile);

  // Set Buffer
  StEpSimuJetEvent* simuJetEvent = 0;
  jetChain->SetBranchAddress("myJets",&simuJetEvent);

  // Open Output ROOT File
  TFile* ofile = TFile::Open(outfile,"recreate");
  assert(ofile);

  // List Q2 Range
  if(q2Set == 0) cout << "Q2 Range = All" << endl;
  if(q2Set == 1) cout << "Q2 Range: 10^-9 <= Q2 < 10^-8" << endl;
  if(q2Set == 2) cout << "Q2 Range: 10^-8 <= Q2 < 10^-7" << endl;
  if(q2Set == 3) cout << "Q2 Range: 10^-7 <= Q2 < 10^-6" << endl;
  if(q2Set == 4) cout << "Q2 Range: 10^-6 <= Q2 < 10^-5" << endl;
  if(q2Set == 5) cout << "Q2 Range: 10^-5 <= Q2 < 10^-4" << endl;
  if(q2Set == 6) cout << "Q2 Range: 10^-4 <= Q2 < 10^-3" << endl;
  if(q2Set == 7) cout << "Q2 Range: 10^-3 <= Q2 < 10^-2" << endl;
  if(q2Set == 8) cout << "Q2 Range: 10^-2 <= Q2 < 10^-1" << endl;
  if(q2Set == 9) cout << "Q2 Range: 10^-1 <= Q2 < 10^0" << endl;
  if(q2Set == 10) cout << "Q2 Range: 10^0 <= Q2 < 10^1" << endl;
  if(q2Set == 11) cout << "Q2 Range: 10^1 <= Q2 < 10^2" << endl;
  

  //TH2D *testDeltaR = new TH2D("testDeltaR","",500,0.,5.,500,0.,5.);

  TH2D *jetPartonDR_strict = new TH2D("jetPartonDR_strict","Jet 2 Vs Jet 1 Parton DR",500,0.,5.,500,0.,5.);
  TH1D *massOverRootS_strict = new TH1D("massOverRootS_strict","Dijet Mass / Root S",500,0.,5.);

  TH2D *jetPartonDR_loose = new TH2D("jetPartonDR_loose","Jet 2 Vs Jet 1 Parton DR",500,0.,5.,500,0.,5.);
  TH1D *massOverRootS_loose = new TH1D("massOverRootS_loose","Dijet Mass / Root S",500,0.,5.);


  TH2D *jetPartonDR = new TH2D("jetPartonDR","Jet 2 Vs Jet 1 Parton DR",500,0.,5.,500,0.,5.);
  TH2D *jetWTAPartonDR = new TH2D("jetWTAPartonDR","Jet 2 Vs Jet 1 Parton DR",500,0.,5.,500,0.,5.);

  TH1D *jet1PartonDRRatio = new TH1D("jet1PartonDRRatio","WTA/Thrust Parton-Jet DR",500,0.,5.);
  TH1D *jet2PartonDRRatio = new TH1D("jet2PartonDRRatio","WTA/Thrust Parton-Jet DR",500,0.,5.);

  TH1D *massOverRootS = new TH1D("massOverRootS","Dijet Mass / Root S",500,0.,5.);
  TH1D *massWTAOverRootS = new TH1D("massWTAOverRootS","Dijet Mass / Root S",500,0.,5.);
  
  
  TH1D *subProcessHist = new TH1D("subProcessHist","Subprocess Classes",200,0.,200.);

  // Kinematic Cut Counter
  TH1D *subProcessClass = new TH1D("subProcessClass","Subprocess classes before kin cuts",5,0.,5.);
  //TH1D *subProcessClass_q2 = new TH1D("subProcessClass_q2","Subprocess classes failing Q2 cut",5,0.,5.);
  TH1D *subProcessClass_y = new TH1D("subProcessClass_y","Subprocess classes failing Y cut",5,0.,5.);
  TH1D *subProcessClass_cut = new TH1D("subProcessClass_cut","Subprocess classes passing kin cuts",5,0.,5.);


  //int totalCounts = 0;
  //int abnormalRecord = 0;
  int badBoosts = 0;


  // Event Loop
  for(int iEvent=0; iEvent<10000000000; iEvent++) //100000000
    {
      if(jetChain->GetEvent(iEvent) <= 0) break;

      //if(iEvent > 5) break;

      // Progress Indicator
      if(iEvent % 100000 == 0) cout << "Event " << iEvent << endl;
      //if(iEvent % 1000 == 0) cout << "Event " << iEvent << endl;

      // Should Not be Null
      assert(simuJetEvent);

      // Check Number of Jet Definitions
      if(simuJetEvent->numberOfJetDefs() != 4) cout << "BAD NUM JET DEFS" << endl;

      // Check for Valid Boost
      double bxL = simuJetEvent->boostPx()/simuJetEvent->boostE();
      double byL = simuJetEvent->boostPy()/simuJetEvent->boostE();
      double bzL = simuJetEvent->boostPz()/simuJetEvent->boostE();
      if(bxL*bxL + byL*byL + bzL*bzL >= 1.0)
	{
	  //cout << iEvent << " " << simuJetEvent->trueQ2() << " " << simuJetEvent->trueX() << " " << simuJetEvent->trueY() << " " << bxL*bxL + byL*byL + bzL*bzL << endl;
	  //printf("",bxL*bxL + byL*byL + bzL*bzL);
	  badBoosts++;
	  continue;
	}
      /*
      if(!simuJetEvent->goodBoost()) 
	{
	  badBoosts++;
	  continue;
	}
      */

      // Set Q2 Analysis Range
      if(q2Set != 0)
	{
	  Double_t loQ2Limit = -1.;
	  Double_t hiQ2Limit = -1.;

	  switch ( q2Set ) {
	  case 1:
	    loQ2Limit = 0.000000001;
	    hiQ2Limit = 0.00000001;
	    break;
	  case 2:
	    loQ2Limit = 0.00000001;
	    hiQ2Limit = 0.0000001;
	    break;
	  case 3:
	    loQ2Limit = 0.0000001;
	    hiQ2Limit = 0.000001;
	    break;
	  case 4:
	    loQ2Limit = 0.000001;
	    hiQ2Limit = 0.00001;
	    break;
	  case 5:
	    loQ2Limit = 0.00001;
	    hiQ2Limit = 0.0001;
	    break;
	  case 6:
	    loQ2Limit = 0.0001;
	    hiQ2Limit = 0.001;
	    break;
	  case 7:
	    loQ2Limit = 0.001;
	    hiQ2Limit = 0.01;
	    break;
	  case 8:
	    loQ2Limit = 0.01;
	    hiQ2Limit = 0.1;
	    break;
	  case 9:
	    loQ2Limit = 0.1;
	    hiQ2Limit = 1.0;
	    break;
	  case 10:
	    loQ2Limit = 1.0;
	    hiQ2Limit = 10.0;
	    break;
	  case 11:
	    loQ2Limit = 10.0;
	    hiQ2Limit = 100.0;
	    break;
	  default:
	    loQ2Limit = 0.0000000001;
	    hiQ2Limit = 1000000000.0;
	    cout << "Invalid Q2 Range" << endl;
	  }

	  if(simuJetEvent->trueQ2() < loQ2Limit || simuJetEvent->trueQ2() >= hiQ2Limit) continue;
	}
      

      // Create a Jet Definitions Vector
      vector<StEpSimuJetDef*> jetDefs;

      for(int iJetDef=0; iJetDef<simuJetEvent->numberOfJetDefs(); iJetDef++)
	{
	  StEpSimuJetDef *def = simuJetEvent->jetDef(iJetDef);
	  jetDefs.push_back(def);

	  // Jet Definition Consistency Check
	  float rad10 = TMath::Abs(def->radius() - 1.0);
	  float rad04 = TMath::Abs(def->radius() - 0.4);
	  float minpt25 = TMath::Abs(def->minPt() - 0.25);

	  //cout << iJetDef << " " << def->radius() << " " << def->minPt() << endl;

	  switch ( iJetDef ) {
	  case 0:
	    if(rad10 > 0.001 || minpt25 > 0.001 || def->algo() != 0 || def->frame() != 2) cout << "BAD DEF 0" << endl;
	    break;
	  case 1:
	    if(rad04 > 0.001 || minpt25 > 0.001 || def->algo() != 0 || def->frame() != 2) cout << "BAD DEF 1" << endl;
	    break;
	  case 2:
	    if(rad10 > 0.001 || minpt25 > 0.001 || def->algo() != 0 || def->frame() != 2) cout << "BAD DEF 2" << endl;
	    //cout << rad08 << " " << def->radius() << " " << minpt25 << " " << def->algo() << " " << def->frame() << endl;
	    break;
	  case 3:
	    if(rad04 > 0.001 || minpt25 > 0.001 || def->algo() != 0 || def->frame() != 2) cout << "BAD DEF 3" << endl;
	    break;
	  default:
	    cout << "Something wrong with iJetDef" << endl;
	  }

	}
      
      
      // Look at Subprocesses
      subProcessHist->Fill(simuJetEvent->processID());

      int subIndex = -1;
      if(simuJetEvent->processID() >= 11 && simuJetEvent->processID() <= 68) subIndex = 0;
      if(simuJetEvent->processID() >= 91 && simuJetEvent->processID() <= 95) subIndex = 1;
      if(simuJetEvent->processID() == 131 || simuJetEvent->processID() == 132) subIndex = 2;
      if(simuJetEvent->processID() == 135 || simuJetEvent->processID() == 136) subIndex = 3;
      if(simuJetEvent->processID() == 99) subIndex = 4;

      //testQ2->Fill(simuJetEvent->trueQ2());
      
      // Place Event Kinematics Cuts
      int kinCut = 0;
      subProcessClass->Fill(subIndex);
      if(simuJetEvent->trueY() < 0.01 || simuJetEvent->trueY() > 0.95) 
	{
	  subProcessClass_y->Fill(subIndex);
	  kinCut = 1;
	}
      if(kinCut == 1) continue;
      subProcessClass_cut->Fill(subIndex);

      if(subIndex != 2 && subIndex != 3) continue;

      StEpSimuParticle *part10 = simuJetEvent->particle(9);
      StEpSimuParticle *part11 = simuJetEvent->particle(10);

      //cout << part10->parentIndex() << endl;
      //cout << part11->parentIndex() << endl;

      /*
      cout << "####################" << endl;
      cout << "Parton 1 pT eta rap phi: " << part10->pt_breit() << " " << part10->eta_breit() << " " << part10->rap_breit() << " " << part10->phi_breit() << endl;
      cout << "Parton 2 pT eta rap phi: " << part11->pt_breit() << " " << part11->eta_breit() << " " << part11->rap_breit() << " " << part11->phi_breit() << endl;
      for(int i=0; i<jetDefs[0]->numberOfJets(); i++)
	{
	  StEpSimuJet *jet = jetDefs[0]->jet(i);

	  if(jet->pt() < 5.0) continue;

	  cout << endl;
	  cout << "Jet " << i << " pT eta rap phi mass: " << jet->pt() << " " << jet->eta() << " " << jet->rap() << " " << TVector2::Phi_mpi_pi(jet->phi()) << " " << jet->fourMomentum().M() << endl;

	  TLorentzVector vec;
	  for(int j=0; j<jet->numberOfJetParticles(); j++)
	    {
	      StEpSimuJetParticle *part = jet->jetParticle(j);
	      
	      cout << "Particle " << j << " pT eta rap phi: " << part->pt() << " " << part->eta() << " " << part->rap() << " " << part->phi() << endl;

	      vec += part->fourMomentum();
	    }
	  cout << "Reco Thrust pT eta rap phi mass: " << vec.Perp() << " " << vec.PseudoRapidity() << " " << vec.Rapidity() << " " << vec.Phi() << " " << vec.M() << endl;
	}
      cout << endl;

      for(int i=0; i<jetDefs[2]->numberOfJets(); i++)
	{
	  StEpSimuJet *jet = jetDefs[2]->jet(i);

	  if(jet->pt() < 5.0) continue;

	  cout << endl;
	  cout << "WTA Jet " << i << " pT eta rap phi mass: " << jet->pt() << " " << jet->eta() << " " << jet->rap() << " " << TVector2::Phi_mpi_pi(jet->phi()) << " " << jet->fourMomentum().M() << endl;

	  TLorentzVector vec;
	  for(int j=0; j<jet->numberOfJetParticles(); j++)
	    {
	      StEpSimuJetParticle *part = jet->jetParticle(j);
	      
	      cout << "Particle " << j << " pT eta rap phi: " << part->pt() << " " << part->eta() << " " << part->rap() << " " << part->phi() << endl;

	      vec += part->fourMomentum();
	    }
	  cout << "Reco Thrust pT eta rap phi mass: " << vec.Perp() << " " << vec.PseudoRapidity() << " " << vec.Rapidity() << " " << vec.Phi() << " " << vec.M() << endl;
	}
      cout << endl;
      */

      // Dijet Kin
      if(jetDefs[0]->numberOfJets() == 2)
	{
	  StEpSimuJet *jetHiStrict = jetDefs[0]->jet(0);
	  StEpSimuJet *jetLoStrict = jetDefs[0]->jet(1);

	  if(TMath::Cos(jetHiStrict->phi() - jetLoStrict->phi()) > -0.5) continue;
	  if(jetHiStrict->pt() < 5.0 || jetLoStrict->pt() < 4.0) continue;

	  Double_t jetHiStrictDR, jetLoStrictDR;
	  Int_t jetHiStrictPartI, jetLoStrictPartI;

	  jetPartonAssoc(simuJetEvent,jetHiStrict,&jetHiStrictDR,&jetHiStrictPartI);
	  jetPartonAssoc(simuJetEvent,jetLoStrict,&jetLoStrictDR,&jetLoStrictPartI);

	  jetPartonDR_strict->Fill(jetHiStrictDR,jetLoStrictDR);

	  Double_t mass_strict = calcDijetMass(jetHiStrict,jetLoStrict);

	  massOverRootS_strict->Fill(mass_strict/TMath::Sqrt(simuJetEvent->sHat()));
	}

      if(jetDefs[0]->numberOfJets() > 1)
	{
	  StEpSimuJet *jetHiLoose = jetDefs[0]->jet(0);
	  StEpSimuJet *jetLoLoose = jetDefs[0]->jet(1);

	  if(TMath::Cos(jetHiLoose->phi() - jetLoLoose->phi()) > -0.5) continue;
	  if(jetHiLoose->pt() < 5.0 || jetLoLoose->pt() < 4.0) continue;

	  Double_t jetHiLooseDR, jetLoLooseDR;
	  Int_t jetHiLoosePartI, jetLoLoosePartI;

	  jetPartonAssoc(simuJetEvent,jetHiLoose,&jetHiLooseDR,&jetHiLoosePartI);
	  jetPartonAssoc(simuJetEvent,jetLoLoose,&jetLoLooseDR,&jetLoLoosePartI);

	  jetPartonDR_loose->Fill(jetHiLooseDR,jetLoLooseDR);

	  Double_t mass_loose = calcDijetMass(jetHiLoose,jetLoLoose);

	  massOverRootS_loose->Fill(mass_loose/TMath::Sqrt(simuJetEvent->sHat()));
	}


      // Compare Axes
      if(jetDefs[0]->numberOfJets() == 2 && jetDefs[2]->numberOfJets() == 2)
	{
	  StEpSimuJet *jet1 = jetDefs[0]->jet(0);
	  StEpSimuJet *jet2 = jetDefs[0]->jet(1);
	  StEpSimuJet *jet1WTA = jetDefs[2]->jet(0);
	  StEpSimuJet *jet2WTA = jetDefs[2]->jet(1);

	  if(TMath::Cos(jet1->phi() - jet2->phi()) > -0.5) continue;
	  if(jet1->pt() < 5.0 || jet2->pt() < 4.0) continue;

	  //Double_t jet1DeltaR = calcDeltaR(jet1,jet1WTA);
	  //Double_t jet2DeltaR = calcDeltaR(jet2,jet2WTA);

	  Double_t jet1PartDR, jet2PartDR, jet1WTAPartDR, jet2WTAPartDR;
	  Int_t jet1PartI, jet2PartI, jet1WTAPartI, jet2WTAPartI;

	  jetPartonAssoc(simuJetEvent,jet1,&jet1PartDR,&jet1PartI);
	  jetPartonAssoc(simuJetEvent,jet2,&jet2PartDR,&jet2PartI);
	  jetPartonAssoc(simuJetEvent,jet1WTA,&jet1WTAPartDR,&jet1WTAPartI);
	  jetPartonAssoc(simuJetEvent,jet2WTA,&jet2WTAPartDR,&jet2WTAPartI);

	  if(jet1PartI == jet2PartI) 
	    {
	      cout << "THRUST JETS POINT TO SAME PARTON " << jet1PartI << endl;
	      cout << "####################" << endl;
	      cout << "Parton 1 pT eta rap phi: " << part10->pt_breit() << " " << part10->eta_breit() << " " << part10->rap_breit() << " " << part10->phi_breit() << endl;
	      cout << "Parton 2 pT eta rap phi: " << part11->pt_breit() << " " << part11->eta_breit() << " " << part11->rap_breit() << " " << part11->phi_breit() << endl;
	      cout << "Jet 1 pT eta rap phi: " << jet1->pt() << " " << jet1->eta() << " " << jet1->rap() << " " << jet1->phi() << endl;
	      cout << "Jet 2 pT eta rap phi: " << jet2->pt() << " " << jet2->eta() << " " << jet2->rap() << " " << jet2->phi() << endl;
	    }
	  //if(jet1WTAPartI == jet2WTAPartI) cout << "WTA JETS POINT TO SAME PARTON" << endl;

	  jetPartonDR->Fill(jet1PartDR,jet2PartDR);
	  jetWTAPartonDR->Fill(jet1WTAPartDR,jet2WTAPartDR);

	  jet1PartonDRRatio->Fill(jet1WTAPartDR/jet1PartDR);
	  jet2PartonDRRatio->Fill(jet2WTAPartDR/jet2PartDR);

	  Double_t mass = calcDijetMass(jet1,jet2);
	  Double_t massWTA = calcDijetMass(jet1WTA,jet2WTA); 


	  massOverRootS->Fill(mass/TMath::Sqrt(simuJetEvent->sHat()));
	  massWTAOverRootS->Fill(massWTA/TMath::Sqrt(simuJetEvent->sHat()));

	  //testDeltaR->Fill(jet1DeltaR,jet2DeltaR);
	}

    } // End Event Loop

  //cout << "Total Events = " << totalCounts << " Events Where 10 11 Don't Have Parent Index 0 = " << abnormalRecord << endl;
  cout << "Events with bad Boosts = " << badBoosts << endl;

  // Write and Close Output ROOT File
  ofile->Write();
  ofile->Close();

}
