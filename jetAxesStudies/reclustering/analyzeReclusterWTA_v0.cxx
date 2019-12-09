// Analyze the My Simu Jet Trees

#include "fastjet/ClusterSequence.hh"

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

#include "/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker_Root6/StEpSimuJetEvent.h"
#include "/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker_Root6/StEpSimuJetDef.h"
#include "/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker_Root6/StEpSimuJet.h"
#include "/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker_Root6/StEpSimuParticle.h"
#include "/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker_Root6/StEpSimuJetParticle.h"
#include "/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker_Root6/StEpSimuSubJet.h"

using namespace fastjet;
using namespace std;


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



// Main Program
int main(int argc, char* argv[])
{
  if(argc != 4)
    {
      cerr << "WRONG" << endl;
      exit(EXIT_FAILURE);
    }

  const int q2Set = atoi(argv[1]);
  const char* infile = argv[2];
  const char* outfile = argv[3];

  //const int q2Set = 0;
  //const char* infile = "/eicdata/eic0009/bpage/jetTrees_akt10_v3/Q2_10-100/jetTree.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_4.C.root";
  //const char* infile="/eicdata/eic0009/bpage/jetTrees_akt10_v3/lowQ2/jetTree.ep.20x250.1Mevents.RadCor=0.Q2=0.01-0.1.kT=1.0_4.A.root";
  //const char * infile="/gpfs02/eic/bpage/jetTrees_akt10Pt_v3/Root6/jetTree.ep.20x250.1Mevents.RadCor=0.Q2=0.00001-1.0.kT=1.0_23.A.root";
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


  // Create Histograms
  char algoNames[6][100];
  char pTNames[2][100];
  char subNames[5][100];

  sprintf(algoNames[0],"lab");
  sprintf(algoNames[1],"labPt");
  sprintf(algoNames[2],"hbv");
  sprintf(algoNames[3],"hbvPt");
  sprintf(algoNames[4],"breit");
  sprintf(algoNames[5],"breitPt");

  sprintf(pTNames[0],"allPt");
  sprintf(pTNames[1],"hiPt");
  
  sprintf(subNames[0],"hQCD");
  sprintf(subNames[1],"sQCD");
  sprintf(subNames[2],"QCDC");
  sprintf(subNames[3],"PGF");
  sprintf(subNames[4],"DIS");


  // Jet Kinematics
  TH1D *jetPtHist[6][2][5];
  TH1D *jetEtaHist[6][2][5];
  TH1D *jetRapHist[6][2][5];

  // WTA Jet Kinematics
  TH1D *wtaPtHist[6][2][5];
  TH1D *wtaEtaHist[6][2][5];
  TH1D *wtaRapHist[6][2][5];

  // Comparisons
  TH1D *jetWTADeltaRHist[6][2][5];
  TH1D *jetWTAPDeltaRHist[6][2][5];
  TH1D *wtaWTAPDeltaRHist[6][2][5];

  TH2D *jetWTADeltaRVsPtHist[6][2][5];
  TH2D *jetWTADeltaRVsEtaHist[6][2][5];

  TH2D *jetWTAPDeltaRVsPtHist[6][2][5];
  TH2D *jetWTAPDeltaRVsEtaHist[6][2][5];

  TH2D *wtaWTAPDeltaRVsPtHist[6][2][5];
  TH2D *wtaWTAPDeltaRVsEtaHist[6][2][5];


  //Inclusive Jets
  for(int i=0; i<6; i++)
    {
      for(int j=0; j<2; j++)
	{
	  for(int k=0; k<5; k++)
	    {
	      jetPtHist[i][j][k] = new TH1D(Form("jetPt_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("Jet pT: %s %s %s",algoNames[i],pTNames[j],subNames[k]),100,0.,50.);
	      jetPtHist[i][j][k]->Sumw2();
	      jetEtaHist[i][j][k] = new TH1D(Form("jetEta_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("Jet Eta: %s %s %s",algoNames[i],pTNames[j],subNames[k]),200,-5.,5.);
	      jetEtaHist[i][j][k]->Sumw2();
	      jetRapHist[i][j][k] = new TH1D(Form("jetRap_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("Jet Rap: %s %s %s",algoNames[i],pTNames[j],subNames[k]),400,-5.,20.);
	      jetRapHist[i][j][k]->Sumw2();

	      wtaPtHist[i][j][k] = new TH1D(Form("wtaPt_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("WTA pT: %s %s %s",algoNames[i],pTNames[j],subNames[k]),100,0.,50.);
	      wtaPtHist[i][j][k]->Sumw2();
	      wtaEtaHist[i][j][k] = new TH1D(Form("wtaEta_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("WTA Eta: %s %s %s",algoNames[i],pTNames[j],subNames[k]),200,-5.,5.);
	      wtaEtaHist[i][j][k]->Sumw2();
	      wtaRapHist[i][j][k] = new TH1D(Form("wtaRap_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("WTA Rap: %s %s %s",algoNames[i],pTNames[j],subNames[k]),400,-5.,20.);
	      wtaRapHist[i][j][k]->Sumw2();

	      jetWTADeltaRHist[i][j][k] = new TH1D(Form("jetWTADeltaR_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("DeltaR Between Jet WTA: %s %s %s",algoNames[i],pTNames[j],subNames[k]),500,0.,1.);
	      jetWTADeltaRHist[i][j][k]->Sumw2();
	      jetWTAPDeltaRHist[i][j][k] = new TH1D(Form("jetWTAPDeltaR_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("DeltaR Between Jet WTAP: %s %s %s",algoNames[i],pTNames[j],subNames[k]),500,0.,1.);
	      jetWTAPDeltaRHist[i][j][k]->Sumw2();
	      wtaWTAPDeltaRHist[i][j][k] = new TH1D(Form("wtaWTAPDeltaR_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("DeltaR Between WTA WTAP: %s %s %s",algoNames[i],pTNames[j],subNames[k]),500,0.,1.);
	      wtaWTAPDeltaRHist[i][j][k]->Sumw2();

	      jetWTADeltaRVsPtHist[i][j][k] = new TH2D(Form("jetWTADeltaRVsPt_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("DeltaR Between Jet WTA Vs Pt: %s %s %s",algoNames[i],pTNames[j],subNames[k]),100,0.,100.,500,0.,1.);
	      jetWTADeltaRVsPtHist[i][j][k]->Sumw2();
	      jetWTADeltaRVsEtaHist[i][j][k] = new TH2D(Form("jetWTADeltaRVsEta_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("DeltaR Between Jet WTA Vs Eta: %s %s %s",algoNames[i],pTNames[j],subNames[k]),200,-5.,5.,500,0.,1.);
	      jetWTADeltaRVsEtaHist[i][j][k]->Sumw2();

	      jetWTAPDeltaRVsPtHist[i][j][k] = new TH2D(Form("jetWTAPDeltaRVsPt_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("DeltaR Between Jet WTAP Vs Pt: %s %s %s",algoNames[i],pTNames[j],subNames[k]),100,0.,100.,500,0.,1.);
	      jetWTAPDeltaRVsPtHist[i][j][k]->Sumw2();
	      jetWTAPDeltaRVsEtaHist[i][j][k] = new TH2D(Form("jetWTAPDeltaRVsEta_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("DeltaR Between Jet WTAP Vs Eta: %s %s %s",algoNames[i],pTNames[j],subNames[k]),200,-5.,5.,500,0.,1.);
	      jetWTAPDeltaRVsEtaHist[i][j][k]->Sumw2();

	      wtaWTAPDeltaRVsPtHist[i][j][k] = new TH2D(Form("wtaWTAPDeltaRVsPt_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("DeltaR Between WTA WTAP Vs Pt: %s %s %s",algoNames[i],pTNames[j],subNames[k]),100,0.,100.,500,0.,1.);
	      wtaWTAPDeltaRVsPtHist[i][j][k]->Sumw2();
	      wtaWTAPDeltaRVsEtaHist[i][j][k] = new TH2D(Form("wtaWTAPDeltaRVsEta_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("DeltaR Between WTA WTAP Vs Eta: %s %s %s",algoNames[i],pTNames[j],subNames[k]),200,-5.,5.,500,0.,1.);
	      wtaWTAPDeltaRVsEtaHist[i][j][k]->Sumw2();
	    }
	}
    }

  
  TH1D *subProcessHist = new TH1D("subProcessHist","Subprocess Classes",200,0.,200.);

  //TH1D *highXEta = new TH1D("highXEta","",200,-5.,5.);

  // Kinematic Cut Counter
  TH1D *subProcessClass = new TH1D("subProcessClass","Subprocess classes before kin cuts",5,0.,5.);
  //TH1D *subProcessClass_q2 = new TH1D("subProcessClass_q2","Subprocess classes failing Q2 cut",5,0.,5.);
  TH1D *subProcessClass_y = new TH1D("subProcessClass_y","Subprocess classes failing Y cut",5,0.,5.);
  TH1D *subProcessClass_cut = new TH1D("subProcessClass_cut","Subprocess classes passing kin cuts",5,0.,5.);


  //int totalCounts = 0;
  //int abnormalRecord = 0;
  int badBoosts = 0;


  // Event Loop
  for(int iEvent=0; iEvent<100000000; iEvent++) //100000000
    {
      if(jetChain->GetEvent(iEvent) <= 0) break;

      //if(iEvent > 5) break;

      // Progress Indicator
      if(iEvent % 100000 == 0) cout << "Event " << iEvent << endl;
      //if(iEvent % 1000 == 0) cout << "Event " << iEvent << endl;

      // Should Not be Null
      assert(simuJetEvent);

      // Check Number of Jet Definitions
      if(simuJetEvent->numberOfJetDefs() != 6) cout << "BAD NUM JET DEFS" << endl;

      // Check for Valid Boost
      if(!simuJetEvent->goodBoost()) 
	{
	  badBoosts++;
	  continue;
	}

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
	  float minpt25 = TMath::Abs(def->minPt() - 0.25);
	  float minpt50 = TMath::Abs(def->minPt() - 0.50);

	  switch ( iJetDef ) {
	  case 0:
	    if(rad10 > 0.001 || minpt25 > 0.001 || def->algo() != 0 || def->frame() != 0) cout << "BAD DEF 0" << endl;
	    break;
	  case 1:
	    if(rad10 > 0.001 || minpt50 > 0.001 || def->algo() != 0 || def->frame() != 0) cout << "BAD DEF 0" << endl;
	    break;
	  case 2:
	    if(rad10 > 0.001 || minpt25 > 0.001 || def->algo() != 0 || def->frame() != 1) cout << "BAD DEF 1" << endl;
	    break;
	  case 3:
	    if(rad10 > 0.001 || minpt50 > 0.001 || def->algo() != 0 || def->frame() != 1) cout << "BAD DEF 1" << endl;
	    break;
	  case 4:
	    if(rad10 > 0.001 || minpt25 > 0.001 || def->algo() != 0 || def->frame() != 2) cout << "BAD DEF 2" << endl;
	    break;
	  case 5:
	    if(rad10 > 0.001 || minpt50 > 0.001 || def->algo() != 0 || def->frame() != 2) cout << "BAD DEF 2" << endl;
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


      // Set Up Boost From Lab to Breit and From Breit to Lab
      Double_t boostX = simuJetEvent->boostPx();
      Double_t boostY = simuJetEvent->boostPy();
      Double_t boostZ = simuJetEvent->boostPz();
      Double_t boostE = simuJetEvent->boostE();

      Double_t boostTheta = simuJetEvent->boostTheta();
      Double_t boostPhi = simuJetEvent->boostPhi();

      TLorentzVector labToBreitBoost;
      labToBreitBoost.SetXYZT(boostX,boostY,boostZ,boostE);

      TLorentzVector breitToLabBoost;
      breitToLabBoost.SetXYZT(-boostX,-boostY,-boostZ,boostE);


      // Look at Inclusive Jets
      for(int branch=0; branch<6; branch++)
	{
	  if(branch == 1 || branch == 2 || branch == 3 || branch == 5) continue; // Only look at Lab and Breit frame jets with min pt > 250 MeV
	  for(int i=0; i<jetDefs[branch]->numberOfJets(); i++)
	    {
	      StEpSimuJet *jet = jetDefs[branch]->jet(i);

	      // Transform from Breit to Lab
	      TLorentzVector breitToLab;
	      breitToLab = rotateY(jet->fourMomentum(),boostTheta);
	      breitToLab = rotateZ(breitToLab,boostPhi);
	      breitToLab = boost(breitToLab,breitToLabBoost);

	      if(jet->pt() >= 5.0)
		{
		  // Recluster Event Using WTA
		  vector<PseudoJet> particlesJet;

		  for(int j=0; j<jet->numberOfJetParticles(); j++)
		    {
		      StEpSimuJetParticle *p = jet->jetParticle(j);

		      fastjet::PseudoJet pPt(p->fourMomentum().Px(),p->fourMomentum().Py(),p->fourMomentum().Pz(),p->fourMomentum().E());
		      pPt.set_user_index(p->index());
		      particlesJet.push_back(pPt);
		    }

		  // Set Jet Definitions
		  double R_10 = 2.5;
		  JetDefinition jet_def_akt_wtaPt_10(antikt_algorithm,R_10,fastjet::WTA_pt_scheme,Best);
		  //JetDefinition jet_def_akt_wtaPt_10(kt_algorithm,R_10);
		  JetDefinition jet_def_akt_wtaMP_10(antikt_algorithm,R_10,fastjet::WTA_modp_scheme,Best);
		  
		  // Run Clustering and Extract the Jets
		  double ptmin = 1.0;

		  // Cluster
		  ClusterSequence csPt_akt_10_wtaPt(particlesJet, jet_def_akt_wtaPt_10);
		  ClusterSequence csPt_akt_10_wtaMP(particlesJet, jet_def_akt_wtaMP_10);

		  // WTA Jets
		  vector<PseudoJet> jetsPt_akt_10_wtaPt = sorted_by_pt(csPt_akt_10_wtaPt.inclusive_jets(ptmin));
		  vector<PseudoJet> jetsPt_akt_10_wtaMP = sorted_by_pt(csPt_akt_10_wtaMP.inclusive_jets(ptmin));
      
		  if(jetsPt_akt_10_wtaPt.size() > 1 || jetsPt_akt_10_wtaMP.size() > 1)
		    {
		      cout << "MORE JETS" << endl;
		    } 

		  // Construct WTA Jets
		  TLorentzVector wtaJet;
		  wtaJet.SetPtEtaPhiE(jetsPt_akt_10_wtaPt[0].pt(),jetsPt_akt_10_wtaPt[0].eta(),jetsPt_akt_10_wtaPt[0].phi(),jetsPt_akt_10_wtaPt[0].e());

		  TLorentzVector wtaPJet;
		  wtaPJet.SetPtEtaPhiE(jetsPt_akt_10_wtaMP[0].pt(),jetsPt_akt_10_wtaMP[0].eta(),jetsPt_akt_10_wtaMP[0].phi(),jetsPt_akt_10_wtaMP[0].e());

		  // Transform from Breit to Lab
		  TLorentzVector breitToLabWTA;
		  breitToLabWTA = rotateY(wtaJet,boostTheta);
		  breitToLabWTA = rotateZ(breitToLabWTA,boostPhi);
		  breitToLabWTA = boost(breitToLabWTA,breitToLabBoost);


		  // Plot Jet Quantities
		  jetPtHist[branch][0][subIndex]->Fill(jet->pt());
		  if(branch == 0) jetEtaHist[branch][0][subIndex]->Fill(jet->eta());
		  if(branch == 4) jetEtaHist[branch][0][subIndex]->Fill(breitToLab.PseudoRapidity());
		  jetRapHist[branch][0][subIndex]->Fill(jet->rap());

		  // Plot WTA Jet Quantities
		  wtaPtHist[branch][0][subIndex]->Fill(wtaJet.Perp());
		  if(branch == 0) wtaEtaHist[branch][0][subIndex]->Fill(wtaJet.PseudoRapidity());
		  if(branch == 4) wtaEtaHist[branch][0][subIndex]->Fill(breitToLabWTA.PseudoRapidity());
		  wtaRapHist[branch][0][subIndex]->Fill(wtaJet.Rapidity());
		  

		  // Calc Delta R
		  Double_t dRap1 = wtaJet.Rapidity() - jet->rap();
		  Double_t dPhi1 = TVector2::Phi_mpi_pi(wtaJet.Phi() - jet->phi());

		  Double_t dRap2 = wtaPJet.Rapidity() - jet->rap();
		  Double_t dPhi2 = TVector2::Phi_mpi_pi(wtaPJet.Phi() - jet->phi());

		  Double_t dRap3 = wtaJet.Rapidity() - wtaPJet.Rapidity();
		  Double_t dPhi3 = TVector2::Phi_mpi_pi(wtaJet.Phi() - wtaPJet.Phi());


		  // Plot Comparisons
		  jetWTADeltaRHist[branch][0][subIndex]->Fill(TMath::Sqrt(dRap1*dRap1 + dPhi1*dPhi1));
		  jetWTAPDeltaRHist[branch][0][subIndex]->Fill(TMath::Sqrt(dRap2*dRap2 + dPhi2*dPhi2));
		  wtaWTAPDeltaRHist[branch][0][subIndex]->Fill(TMath::Sqrt(dRap3*dRap3 + dPhi3*dPhi3));

		  jetWTADeltaRVsPtHist[branch][0][subIndex]->Fill(jet->pt(),TMath::Sqrt(dRap1*dRap1 + dPhi1*dPhi1));
		  if(branch == 0) jetWTADeltaRVsEtaHist[branch][0][subIndex]->Fill(jet->eta(),TMath::Sqrt(dRap1*dRap1 + dPhi1*dPhi1));
		  if(branch == 4) jetWTADeltaRVsEtaHist[branch][0][subIndex]->Fill(breitToLab.PseudoRapidity(),TMath::Sqrt(dRap1*dRap1 + dPhi1*dPhi1));

		  jetWTAPDeltaRVsPtHist[branch][0][subIndex]->Fill(jet->pt(),TMath::Sqrt(dRap2*dRap2 + dPhi2*dPhi2));
		  if(branch == 0) jetWTAPDeltaRVsEtaHist[branch][0][subIndex]->Fill(jet->eta(),TMath::Sqrt(dRap2*dRap2 + dPhi2*dPhi2));
		  if(branch == 4) jetWTAPDeltaRVsEtaHist[branch][0][subIndex]->Fill(breitToLab.PseudoRapidity(),TMath::Sqrt(dRap2*dRap2 + dPhi2*dPhi2));

		  wtaWTAPDeltaRVsPtHist[branch][0][subIndex]->Fill(jet->pt(),TMath::Sqrt(dRap3*dRap3 + dPhi3*dPhi3));
		  if(branch == 0) wtaWTAPDeltaRVsEtaHist[branch][0][subIndex]->Fill(jet->eta(),TMath::Sqrt(dRap3*dRap3 + dPhi3*dPhi3));
		  if(branch == 4) wtaWTAPDeltaRVsEtaHist[branch][0][subIndex]->Fill(breitToLab.PseudoRapidity(),TMath::Sqrt(dRap3*dRap3 + dPhi3*dPhi3));
		}
	    }
	}
    } // End Event Loop

  //cout << "Total Events = " << totalCounts << " Events Where 10 11 Don't Have Parent Index 0 = " << abnormalRecord << endl;
  cout << "Events with bad Boosts = " << badBoosts << endl;

  // Write and Close Output ROOT File
  ofile->Write();
  ofile->Close();

}
