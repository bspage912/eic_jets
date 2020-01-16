// Analyze the My Simu Jet Trees

#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/Recluster.hh"
#include "fastjet/contrib/SoftDrop.hh"

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


Int_t disPartIndexAfterRad(StEpSimuJetEvent *event)
{
  StEpSimuParticle *part10 = event->particle(9);

  Int_t part10PDG = part10->pdgCode();

  Int_t workingIndex = -1;
  //Double_t workingEta = -1.;

  for(int i=0; i<event->numberOfTotalParticles(); i++)
    {
      StEpSimuParticle *p = event->particle(i);

      if(p->parentIndex() == 10 && (p->pdgCode() == part10PDG))
	{
	  workingIndex = p->index();
	}
    }
  return workingIndex;
}


Double_t calcDeltaR(StEpSimuParticle *part, StEpSimuJet *jet)
{
  Double_t deltaEta = part->eta_breit() - jet->eta();
  Double_t deltaPhi = TVector2::Phi_mpi_pi(part->phi_breit() - TVector2::Phi_mpi_pi(jet->phi()));

  return TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
}


Int_t jetPartonAssoc(StEpSimuJetEvent *event, int defNum, StEpSimuJet *jet)
{
  // Jet Objects
  //StEpSimuJetDef *def = event->jetDef(defNum);

  // Assign Each Particle to a Jet
  StEpSimuParticle *partOne = event->particle(9);
  StEpSimuParticle *partTwo = event->particle(10);

  Double_t partOneDeltaR = calcDeltaR(partOne,jet);
  Double_t partTwoDeltaR = calcDeltaR(partTwo,jet);

  //Double_t jetDeltaR = 999.0;
  Int_t index = -1;

  if(partOneDeltaR < partTwoDeltaR)
    {
      //jetDeltaR = partOneDeltaR;
      index = partOne->index();
    }
  if(partTwoDeltaR < partOneDeltaR)
    {
      //jetDeltaR = partTwoDeltaR;
      index = partTwo->index();
    }

  return index;
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

  // Groming parameters
  //double beta = 0.0;
  //double zcut = 0.1;
  // Paper variations: 
  // fix zcut = 0.1 and vary beta = 0, 1, 2
  // fix beta = 1   and vary zcut = 0.05, 0.1, 0.2, 0.3
  // Not that's for LHC kinematics

  // Create Histograms
  char algoNames[2][100];
  char pTNames[2][100];
  char subNames[5][100];
  char grmNames[6][100];

  sprintf(algoNames[0],"lab");
  sprintf(algoNames[1],"breit");

  sprintf(pTNames[0],"allPt");
  sprintf(pTNames[1],"hiPt");
  
  sprintf(subNames[0],"hQCD");
  sprintf(subNames[1],"sQCD");
  sprintf(subNames[2],"QCDC");
  sprintf(subNames[3],"PGF");
  sprintf(subNames[4],"DIS");

  sprintf(grmNames[0],"z010b00");
  sprintf(grmNames[1],"z010b10");
  sprintf(grmNames[2],"z010b20");
  sprintf(grmNames[3],"z005b10");
  sprintf(grmNames[4],"z020b10");
  sprintf(grmNames[5],"z030b10");

  // Jet Kinematics
  TH1D *jetPtHist[2][2][5];
  TH1D *jetEtaLabHist[2][2][5];
  TH1D *jetRapHist[2][2][5];

  // WTA Jet Kinematics
  TH1D *wtaPtHist[2][2][5];
  TH1D *wtaEtaLabHist[2][2][5];
  TH1D *wtaRapHist[2][2][5];

  // Groomed Jet Kinematics
  TH1D *grmPtHist[2][2][5][6];
  TH1D *grmEtaLabHist[2][2][5][6];
  TH1D *grmRapHist[2][2][5][6];

  TH2D *jetVsGRMPtHist[2][2][5][6];
  TH2D *jetVsGRMNumPartHist[2][2][5][6];

  // Comparisons
  TH1D *jetWTADeltaRHist[2][2][5];
  TH1D *jetWTAPDeltaRHist[2][2][5];
  TH1D *wtaWTAPDeltaRHist[2][2][5];

  TH1D *jetWTADeltaRLogHist[2][2][5];
  TH1D *jetWTAPDeltaRLogHist[2][2][5];
  TH1D *wtaWTAPDeltaRLogHist[2][2][5];

  TH2D *jetWTADeltaRVsPtHist[2][2][5];
  TH2D *jetWTADeltaRVsEtaHist[2][2][5];

  TH2D *jetWTAPDeltaRVsPtHist[2][2][5];
  TH2D *jetWTAPDeltaRVsEtaHist[2][2][5];

  TH2D *wtaWTAPDeltaRVsPtHist[2][2][5];
  TH2D *wtaWTAPDeltaRVsEtaHist[2][2][5];

  // Groomed axes comparisons
  TH1D *jetGRMDeltaRHist[2][2][5][6];
  TH1D *wtaGRMDeltaRHist[2][2][5][6];

  TH1D *jetGRMDeltaRLogHist[2][2][5][6];
  TH1D *wtaGRMDeltaRLogHist[2][2][5][6];

  TH2D *jetGRMDeltaRVsPtHist[2][2][5][6];
  TH2D *jetGRMDeltaRVsEtaHist[2][2][5][6];

  TH2D *wtaGRMDeltaRVsPtHist[2][2][5][6];
  TH2D *wtaGRMDeltaRVsEtaHist[2][2][5][6];

  // Jet Parton DeltaR
  TH1D *jetPartonDeltaRHist[2][2][5];
  TH1D *wtaPartonDeltaRHist[2][2][5];

  TH2D *wtaVsJetPartonDeltaRHist[2][2][5];
  TH1D *wtaJetPartonDeltaRDiffHist[2][2][5];

  // Groomed Parton DeltaR
  TH1D *grmPartonDeltaRHist[2][2][5][6];

  TH2D *grmVsJetPartonDeltaRHist[2][2][5][6];
  TH1D *grmJetPartonDeltaRDiffHist[2][2][5][6];

  TH2D *grmVsWTAPartonDeltaRHist[2][2][5][6];
  TH1D *grmWTAPartonDeltaRDiffHist[2][2][5][6];

  TH2::SetDefaultSumw2(true);

  //Inclusive Jets
  for(int i=0; i<2; i++)
    {
      for(int j=0; j<2; j++)
	{
	  for(int k=0; k<5; k++)
	    {
	      jetPtHist[i][j][k] = new TH1D(Form("jetPt_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("Jet pT: %s %s %s",algoNames[i],pTNames[j],subNames[k]),100,0.,50.);
	      jetEtaLabHist[i][j][k] = new TH1D(Form("jetEtaLab_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("Jet Lab Eta: %s %s %s",algoNames[i],pTNames[j],subNames[k]),200,-5.,5.);
	      jetRapHist[i][j][k] = new TH1D(Form("jetRap_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("Jet Rap: %s %s %s",algoNames[i],pTNames[j],subNames[k]),400,-5.,20.);
	      
	      wtaPtHist[i][j][k] = new TH1D(Form("wtaPt_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("WTA pT: %s %s %s",algoNames[i],pTNames[j],subNames[k]),100,0.,50.);
	      wtaEtaLabHist[i][j][k] = new TH1D(Form("wtaEtaLab_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("WTA Lab Eta: %s %s %s",algoNames[i],pTNames[j],subNames[k]),200,-5.,5.);
	      wtaRapHist[i][j][k] = new TH1D(Form("wtaRap_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("WTA Rap: %s %s %s",algoNames[i],pTNames[j],subNames[k]),400,-5.,20.);
	      
	      for(int l=0; l<6; l++)
		{
		  grmPtHist[i][j][k][l] = new TH1D(Form("grmPt_%s_%s_%s_%s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),Form("GRM pT: %s %s %s %s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),100,0.,50.);
		  grmEtaLabHist[i][j][k][l] = new TH1D(Form("grmEtaLab_%s_%s_%s_%s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),Form("GRM Lab Eta: %s %s %s %s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),200,-5.,5.);
		  grmRapHist[i][j][k][l] = new TH1D(Form("grmRap_%s_%s_%s_%s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),Form("GRM Rap: %s %s %s %s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),400,-5.,20.);

		  jetVsGRMPtHist[i][j][k][l] = new TH2D(Form("jetVsGRMPt_%s_%s_%s_%s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),Form("Groomed Vs Ungroomed Jet pT: %s %s %s %s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),200,0.,50.,200,0.,50.);
		  jetVsGRMNumPartHist[i][j][k][l] = new TH2D(Form("jetVsGRMNumPart_%s_%s_%s_%s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),Form("Groomed Vs Ungroomed Jet Particles: %s %s %s %s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),20,0.,20.,20,0.,20.);
		}

	      jetWTADeltaRHist[i][j][k] = new TH1D(Form("jetWTADeltaR_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("DeltaR Between Jet WTA: %s %s %s",algoNames[i],pTNames[j],subNames[k]),500,0.,1.);
	      jetWTAPDeltaRHist[i][j][k] = new TH1D(Form("jetWTAPDeltaR_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("DeltaR Between Jet WTAP: %s %s %s",algoNames[i],pTNames[j],subNames[k]),500,0.,1.);
	      wtaWTAPDeltaRHist[i][j][k] = new TH1D(Form("wtaWTAPDeltaR_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("DeltaR Between WTA WTAP: %s %s %s",algoNames[i],pTNames[j],subNames[k]),500,0.,1.);

	      jetWTADeltaRLogHist[i][j][k] = new TH1D(Form("jetWTADeltaRLog_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("Log DeltaR Between Jet WTA: %s %s %s",algoNames[i],pTNames[j],subNames[k]),100,-9.,1.);
	      jetWTAPDeltaRLogHist[i][j][k] = new TH1D(Form("jetWTAPDeltaRLog_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("Log DeltaR Between Jet WTAP: %s %s %s",algoNames[i],pTNames[j],subNames[k]),100,-9.,1.);
	      wtaWTAPDeltaRLogHist[i][j][k] = new TH1D(Form("wtaWTAPDeltaRLog_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("Log DeltaR Between WTA WTAP: %s %s %s",algoNames[i],pTNames[j],subNames[k]),100,-9.,1.);

	      jetWTADeltaRVsPtHist[i][j][k] = new TH2D(Form("jetWTADeltaRVsPt_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("DeltaR Between Jet WTA Vs Pt: %s %s %s",algoNames[i],pTNames[j],subNames[k]),100,0.,100.,500,0.,1.);
	      jetWTADeltaRVsEtaHist[i][j][k] = new TH2D(Form("jetWTADeltaRVsEta_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("DeltaR Between Jet WTA Vs Eta: %s %s %s",algoNames[i],pTNames[j],subNames[k]),200,-5.,5.,500,0.,1.);

	      jetWTAPDeltaRVsPtHist[i][j][k] = new TH2D(Form("jetWTAPDeltaRVsPt_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("DeltaR Between Jet WTAP Vs Pt: %s %s %s",algoNames[i],pTNames[j],subNames[k]),100,0.,100.,500,0.,1.);
	      jetWTAPDeltaRVsEtaHist[i][j][k] = new TH2D(Form("jetWTAPDeltaRVsEta_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("DeltaR Between Jet WTAP Vs Eta: %s %s %s",algoNames[i],pTNames[j],subNames[k]),200,-5.,5.,500,0.,1.);
	      
	      wtaWTAPDeltaRVsPtHist[i][j][k] = new TH2D(Form("wtaWTAPDeltaRVsPt_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("DeltaR Between WTA WTAP Vs Pt: %s %s %s",algoNames[i],pTNames[j],subNames[k]),100,0.,100.,500,0.,1.);
	      wtaWTAPDeltaRVsEtaHist[i][j][k] = new TH2D(Form("wtaWTAPDeltaRVsEta_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("DeltaR Between WTA WTAP Vs Eta: %s %s %s",algoNames[i],pTNames[j],subNames[k]),200,-5.,5.,500,0.,1.);

	      for(int l=0; l<6; l++)
		{
		  jetGRMDeltaRHist[i][j][k][l] = new TH1D(Form("jetGRMDeltaR_%s_%s_%s_%s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),Form("DeltaR Between Jet GRM: %s %s %s %s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),500,0.,1.);
		  wtaGRMDeltaRHist[i][j][k][l] = new TH1D(Form("wtaGRMDeltaR_%s_%s_%s_%s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),Form("DeltaR Between WTA GRM: %s %s %s %s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),500,0.,1.);

		  jetGRMDeltaRLogHist[i][j][k][l] = new TH1D(Form("jetGRMDeltaRLog_%s_%s_%s_%s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),Form("Log DeltaR Between Jet GRM: %s %s %s %s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),100,-9.,1.);
		  wtaGRMDeltaRLogHist[i][j][k][l] = new TH1D(Form("wtaGRMDeltaRLog_%s_%s_%s_%s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),Form("Log DeltaR Between WTA GRM: %s %s %s %s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),100,-9.,1.);

		  jetGRMDeltaRVsPtHist[i][j][k][l] = new TH2D(Form("jetGRMDeltaRVsPt_%s_%s_%s_%s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),Form("DeltaR Between Jet GRM Vs Pt: %s %s %s %s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),100,0.,100.,500,0.,1.);
		  jetGRMDeltaRVsEtaHist[i][j][k][l] = new TH2D(Form("jetGRMDeltaRVsEta_%s_%s_%s_%s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),Form("DeltaR Between Jet GRM Vs Eta: %s %s %s %s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),200,-5.,5.,500,0.,1.);

		  wtaGRMDeltaRVsPtHist[i][j][k][l] = new TH2D(Form("wtaGRMDeltaRVsPt_%s_%s_%s_%s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),Form("DeltaR Between WTA GRM Vs Pt: %s %s %s %s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),100,0.,100.,500,0.,1.);
		  wtaGRMDeltaRVsEtaHist[i][j][k][l] = new TH2D(Form("wtaGRMDeltaRVsEta_%s_%s_%s_%s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),Form("DeltaR Between WTA GRM Vs Eta: %s %s %s %s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),200,-5.,5.,500,0.,1.);
		}

	      jetPartonDeltaRHist[i][j][k] = new TH1D(Form("jetPartonDeltaR_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("Jet Parton DeltaR: %s %s %s",algoNames[i],pTNames[j],subNames[k]),2500,0.,5.);
	      wtaPartonDeltaRHist[i][j][k] = new TH1D(Form("wtaPartonDeltaR_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("WTA Jet Parton DeltaR: %s %s %s",algoNames[i],pTNames[j],subNames[k]),2500,0.,5.);
	      wtaVsJetPartonDeltaRHist[i][j][k] = new TH2D(Form("wtaVsJetPartonDeltaR_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("WTA Vs Jet Parton DeltaR: %s %s %s",algoNames[i],pTNames[j],subNames[k]),2500,0.,5.,2500,0.,5.);
	      wtaJetPartonDeltaRDiffHist[i][j][k] = new TH1D(Form("wtaJetPartonDeltaRDiff_%s_%s_%s",algoNames[i],pTNames[j],subNames[k]),Form("WTA - Jet Parton DeltaR: %s %s %s",algoNames[i],pTNames[j],subNames[k]),200,-1.,1.);

	      for(int l=0; l<6; l++)
		{
		  grmPartonDeltaRHist[i][j][k][l] = new TH1D(Form("grmPartonDeltaR_%s_%s_%s_%s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),Form("GRM Jet Parton DeltaR: %s %s %s %s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),2500,0.,5.);

		  grmVsJetPartonDeltaRHist[i][j][k][l] = new TH2D(Form("grmVsJetPartonDeltaR_%s_%s_%s_%s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),Form("GRM Vs Jet Parton DeltaR: %s %s %s %s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),2500,0.,5.,2500,0.,5.);
		  grmJetPartonDeltaRDiffHist[i][j][k][l] = new TH1D(Form("grmJetPartonDeltaRDiff_%s_%s_%s_%s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),Form("GRM - Jet Parton DeltaR: %s %s %s %s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),200,-1.,1.);

		  grmVsWTAPartonDeltaRHist[i][j][k][l] = new TH2D(Form("grmVsWTAPartonDeltaR_%s_%s_%s_%s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),Form("GRM Vs WTA Parton DeltaR: %s %s %s %s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),2500,0.,5.,2500,0.,5.);
		  grmWTAPartonDeltaRDiffHist[i][j][k][l] = new TH1D(Form("grmWTAPartonDeltaRDiff_%s_%s_%s_%s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),Form("GRM - WTA Parton DeltaR: %s %s %s %s",algoNames[i],pTNames[j],subNames[k],grmNames[l]),200,-1.,1.);
		}
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


      // Find Init Parton for DIS Jet
      //StEpSimuParticle *disParton = simuJetEvent->particle(disPartIndexAfterRad(simuJetEvent)-1);

      // Look at Inclusive Jets
      for(int branchN=0; branchN<6; branchN++)
	{
	  if(branchN == 1 || branchN == 2 || branchN == 3 || branchN == 5) continue; // Only look at Lab and Breit frame jets with min pt > 250 MeV
	  for(int i=0; i<jetDefs[branchN]->numberOfJets(); i++)
	    {
	      StEpSimuJet *jet = jetDefs[branchN]->jet(i);

	      int branch = -1;
	      if(branchN == 0) branch = 0;
	      if(branchN == 4) branch = 1;

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
		  JetDefinition jet_def_akt_wtaMP_10(antikt_algorithm,R_10,fastjet::WTA_modp_scheme,Best);
		  
		  // Cambridge Aachen for grooming
		  JetDefinition jet_def_CA_10( fastjet::cambridge_algorithm, R_10 );		  

		  // Run Clustering and Extract the Jets
		  double ptmin = 1.0;

		  // Cluster
		  ClusterSequence csPt_akt_10_wtaPt(particlesJet, jet_def_akt_wtaPt_10);
		  ClusterSequence csPt_akt_10_wtaMP(particlesJet, jet_def_akt_wtaMP_10);
		  ClusterSequence csPt_CA_10(particlesJet, jet_def_CA_10);

		  // WTA Jets
		  vector<PseudoJet> jetsPt_akt_10_wtaPt = sorted_by_pt(csPt_akt_10_wtaPt.inclusive_jets(ptmin));
		  vector<PseudoJet> jetsPt_akt_10_wtaMP = sorted_by_pt(csPt_akt_10_wtaMP.inclusive_jets(ptmin));
		  
		  // CA jets for grooming 
		  vector<PseudoJet> jetsPt_CA_10 = sorted_by_pt(csPt_CA_10.inclusive_jets(ptmin));
     
		  if(jetsPt_akt_10_wtaPt.size() > 1 || jetsPt_akt_10_wtaMP.size() > 1 || jetsPt_CA_10.size() > 1 )
		    {
		      cout << "MORE JETS" << endl;
		      throw std::runtime_error("Shouldn't find more jets");
		    } 

		  // Construct WTA Jets
		  TLorentzVector wtaJet;
		  wtaJet.SetPtEtaPhiE(jetsPt_akt_10_wtaPt[0].pt(),jetsPt_akt_10_wtaPt[0].eta(),jetsPt_akt_10_wtaPt[0].phi(),jetsPt_akt_10_wtaPt[0].e());

		  TLorentzVector wtaPJet;
		  wtaPJet.SetPtEtaPhiE(jetsPt_akt_10_wtaMP[0].pt(),jetsPt_akt_10_wtaMP[0].eta(),jetsPt_akt_10_wtaMP[0].phi(),jetsPt_akt_10_wtaMP[0].e());


		  // Generate Groomed Pseudojets
		  vector<PseudoJet> groomedPJets;

		  double bVals[3] = {0.0, 1.0, 2.0};
		  for(int beta=0; beta<3; beta++)
		    {
		      double zcut = 0.1;
		      contrib::SoftDrop sd( bVals[beta], zcut);
		      // Add custom recluster or background subtractor here
		      PseudoJet sd_jet = sd( jetsPt_CA_10[0]);
		      if ( sd_jet == 0){
			std::cerr << "Something caused SoftDrop to return 0 ---" << endl;
			throw(-1);
			continue;
		      }
		      groomedPJets.push_back(sd_jet);
		    }

		  double zVals[3] = {0.05, 0.2, 0.3};
		  for(int zcut=0; zcut<3; zcut++)
		    {
		      double beta = 1.0;
		      contrib::SoftDrop sd( beta, zVals[zcut]);
		      // Add custom recluster or background subtractor here
		      PseudoJet sd_jet = sd( jetsPt_CA_10[0]);
		      if ( sd_jet == 0){
			std::cerr << "Something caused SoftDrop to return 0 ---" << endl;
			throw(-1);
			continue;
		      }
		      groomedPJets.push_back(sd_jet);
		    }

		  // could do double zg = sd_jet.structure_of<contrib::SoftDrop>().symmetry();

		  // Create Vector of Groomed Jet 4-vectors
		  // Vector 0: zcut = 0.10 beta = 0.0
		  // Vector 1: zcut = 0.10 beta = 1.0
		  // Vector 2: zcut = 0.10 beta = 2.0
		  // Vector 3: zcut = 0.05 beta = 1.0
		  // Vector 4: zcut = 0.20 beta = 1.0
		  // Vector 5: zcut = 0.30 beta = 1.0
		  vector<TLorentzVector> grmJet;
		  for(int j=0; j<6; j++)
		    {
		      TLorentzVector tmpJet;
		      tmpJet.SetPtEtaPhiE(groomedPJets[j].pt(),groomedPJets[j].eta(),groomedPJets[j].phi(),groomedPJets[j].e());
		      grmJet.push_back(tmpJet);
		    }

		  // Create Vector of Groomed Constituents
		  vector<vector<PseudoJet>> constituents;
		  for(int j=0; j<6; j++)
		    {
		      vector<PseudoJet> c = groomedPJets[j].constituents();
		      constituents.push_back(c);
		    }


		  // Transform from Breit to Lab
		  TLorentzVector breitToLabWTA;
		  breitToLabWTA = rotateY(wtaJet,boostTheta);
		  breitToLabWTA = rotateZ(breitToLabWTA,boostPhi);
		  breitToLabWTA = boost(breitToLabWTA,breitToLabBoost);

		  vector<TLorentzVector> breitToLabGRM;
		  for(int j=0; j<6; j++)
		    {
		      TLorentzVector tmp;
		      tmp = rotateY(grmJet[j],boostTheta);
		      tmp = rotateZ(tmp,boostPhi);
		      tmp = boost(tmp,breitToLabBoost);
		      breitToLabGRM.push_back(tmp);
		    }


		  // Plot Jet Quantities
		  jetPtHist[branch][0][subIndex]->Fill(jet->pt());
		  if(branchN == 0) jetEtaLabHist[branch][0][subIndex]->Fill(jet->eta());
		  if(branchN == 4) jetEtaLabHist[branch][0][subIndex]->Fill(breitToLab.PseudoRapidity());
		  jetRapHist[branch][0][subIndex]->Fill(jet->rap());

		  // Plot WTA Jet Quantities
		  wtaPtHist[branch][0][subIndex]->Fill(wtaJet.Perp());
		  if(branchN == 0) wtaEtaLabHist[branch][0][subIndex]->Fill(wtaJet.PseudoRapidity());
		  if(branchN == 4) wtaEtaLabHist[branch][0][subIndex]->Fill(breitToLabWTA.PseudoRapidity());
		  wtaRapHist[branch][0][subIndex]->Fill(wtaJet.Rapidity());

		  // Plot GRM Jet Quantities
		  for(int j=0; j<6; j++)
		    {
		      grmPtHist[branch][0][subIndex][j]->Fill(grmJet[j].Perp());
		      if(branchN == 0) grmEtaLabHist[branch][0][subIndex][j]->Fill(grmJet[j].PseudoRapidity());
		      if(branchN == 4) grmEtaLabHist[branch][0][subIndex][j]->Fill(breitToLabGRM[j].PseudoRapidity());
		      grmRapHist[branch][0][subIndex][j]->Fill(grmJet[j].Rapidity());
		      
		      jetVsGRMPtHist[branch][0][subIndex][j]->Fill(jet->pt(),grmJet[j].Perp());
		      jetVsGRMNumPartHist[branch][0][subIndex][j]->Fill(jet->numberOfJetParticles(),constituents[j].size());
		    }
		  

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

		  jetWTADeltaRLogHist[branch][0][subIndex]->Fill(TMath::Log10(TMath::Sqrt(dRap1*dRap1 + dPhi1*dPhi1)));
		  jetWTAPDeltaRLogHist[branch][0][subIndex]->Fill(TMath::Log10(TMath::Sqrt(dRap2*dRap2 + dPhi2*dPhi2)));
		  wtaWTAPDeltaRLogHist[branch][0][subIndex]->Fill(TMath::Log10(TMath::Sqrt(dRap3*dRap3 + dPhi3*dPhi3)));
		  
		  jetWTADeltaRVsPtHist[branch][0][subIndex]->Fill(jet->pt(),TMath::Sqrt(dRap1*dRap1 + dPhi1*dPhi1));
		  if(branchN == 0) jetWTADeltaRVsEtaHist[branch][0][subIndex]->Fill(jet->eta(),TMath::Sqrt(dRap1*dRap1 + dPhi1*dPhi1));
		  if(branchN == 4) jetWTADeltaRVsEtaHist[branch][0][subIndex]->Fill(breitToLab.PseudoRapidity(),TMath::Sqrt(dRap1*dRap1 + dPhi1*dPhi1));

		  jetWTAPDeltaRVsPtHist[branch][0][subIndex]->Fill(jet->pt(),TMath::Sqrt(dRap2*dRap2 + dPhi2*dPhi2));
		  if(branchN == 0) jetWTAPDeltaRVsEtaHist[branch][0][subIndex]->Fill(jet->eta(),TMath::Sqrt(dRap2*dRap2 + dPhi2*dPhi2));
		  if(branchN == 4) jetWTAPDeltaRVsEtaHist[branch][0][subIndex]->Fill(breitToLab.PseudoRapidity(),TMath::Sqrt(dRap2*dRap2 + dPhi2*dPhi2));

		  wtaWTAPDeltaRVsPtHist[branch][0][subIndex]->Fill(jet->pt(),TMath::Sqrt(dRap3*dRap3 + dPhi3*dPhi3));
		  if(branchN == 0) wtaWTAPDeltaRVsEtaHist[branch][0][subIndex]->Fill(jet->eta(),TMath::Sqrt(dRap3*dRap3 + dPhi3*dPhi3));
		  if(branchN == 4) wtaWTAPDeltaRVsEtaHist[branch][0][subIndex]->Fill(breitToLab.PseudoRapidity(),TMath::Sqrt(dRap3*dRap3 + dPhi3*dPhi3));

		  for(int j=0; j<6; j++)
		    {
		      Double_t dRap4 = grmJet[j].Rapidity() - jet->rap();
		      Double_t dPhi4 = TVector2::Phi_mpi_pi(grmJet[j].Phi() - jet->phi());

		      Double_t dRap5 = grmJet[j].Rapidity() - wtaJet.Rapidity();
		      Double_t dPhi5 = TVector2::Phi_mpi_pi(grmJet[j].Phi() - wtaJet.Phi());

		      jetGRMDeltaRHist[branch][0][subIndex][j]->Fill(TMath::Sqrt(dRap4*dRap4 + dPhi4*dPhi4));
		      wtaGRMDeltaRHist[branch][0][subIndex][j]->Fill(TMath::Sqrt(dRap5*dRap5 + dPhi5*dPhi5));

		      jetGRMDeltaRLogHist[branch][0][subIndex][j]->Fill(TMath::Log10(TMath::Sqrt(dRap4*dRap4 + dPhi4*dPhi4)));
		      wtaGRMDeltaRLogHist[branch][0][subIndex][j]->Fill(TMath::Log10(TMath::Sqrt(dRap5*dRap5 + dPhi5*dPhi5)));

		      jetGRMDeltaRVsPtHist[branch][0][subIndex][j]->Fill(jet->pt(),TMath::Sqrt(dRap4*dRap4 + dPhi4*dPhi4));
		      if(branchN == 0) jetGRMDeltaRVsEtaHist[branch][0][subIndex][j]->Fill(jet->eta(),TMath::Sqrt(dRap4*dRap4 + dPhi4*dPhi4));
		      if(branchN == 4) jetGRMDeltaRVsEtaHist[branch][0][subIndex][j]->Fill(breitToLab.PseudoRapidity(),TMath::Sqrt(dRap4*dRap4 + dPhi4*dPhi4));

		      wtaGRMDeltaRVsPtHist[branch][0][subIndex][j]->Fill(jet->pt(),TMath::Sqrt(dRap5*dRap5 + dPhi5*dPhi5));
		      if(branchN == 0) wtaGRMDeltaRVsEtaHist[branch][0][subIndex][j]->Fill(jet->eta(),TMath::Sqrt(dRap5*dRap5 + dPhi5*dPhi5));
		      if(branchN == 4) wtaGRMDeltaRVsEtaHist[branch][0][subIndex][j]->Fill(breitToLab.PseudoRapidity(),TMath::Sqrt(dRap5*dRap5 + dPhi5*dPhi5));
		    }

		  
		  // Compare to Parton - DIS
		  if(branchN == 0 && subIndex == 4 && i == 0) // Look at highest jet only
		    {
		      // Find Init Parton for DIS Jet
		      StEpSimuParticle *disParton = simuJetEvent->particle(disPartIndexAfterRad(simuJetEvent)-1);

		      Double_t dRap6 = disParton->fourMomentum().Rapidity() - jet->rap();
		      Double_t dPhi6 = TVector2::Phi_mpi_pi(disParton->fourMomentum().Phi() - jet->phi());

		      Double_t dRap7 = disParton->fourMomentum().Rapidity() - wtaJet.Rapidity();
		      Double_t dPhi7 = TVector2::Phi_mpi_pi(disParton->fourMomentum().Phi() - wtaJet.Phi());

		      jetPartonDeltaRHist[branch][0][subIndex]->Fill(TMath::Sqrt(dRap6*dRap6 + dPhi6*dPhi6));
		      wtaPartonDeltaRHist[branch][0][subIndex]->Fill(TMath::Sqrt(dRap7*dRap7 + dPhi7*dPhi7));

		      wtaVsJetPartonDeltaRHist[branch][0][subIndex]->Fill(TMath::Sqrt(dRap6*dRap6 + dPhi6*dPhi6),TMath::Sqrt(dRap7*dRap7 + dPhi7*dPhi7));
		      wtaJetPartonDeltaRDiffHist[branch][0][subIndex]->Fill(TMath::Sqrt(dRap7*dRap7 + dPhi7*dPhi7)-TMath::Sqrt(dRap6*dRap6 + dPhi6*dPhi6));

		      for(int j=0; j<6; j++)
			{
			  Double_t dRap8 = disParton->fourMomentum().Rapidity() - grmJet[j].Rapidity();
			  Double_t dPhi8 = TVector2::Phi_mpi_pi(disParton->fourMomentum().Phi() - grmJet[j].Phi());

			  grmPartonDeltaRHist[branch][0][subIndex][j]->Fill(TMath::Sqrt(dRap8*dRap8 + dPhi8*dPhi8));

			  grmVsJetPartonDeltaRHist[branch][0][subIndex][j]->Fill(TMath::Sqrt(dRap6*dRap6 + dPhi6*dPhi6),TMath::Sqrt(dRap8*dRap8 + dPhi8*dPhi8));
			  grmJetPartonDeltaRDiffHist[branch][0][subIndex][j]->Fill(TMath::Sqrt(dRap8*dRap8 + dPhi8*dPhi8)-TMath::Sqrt(dRap6*dRap6 + dPhi6*dPhi6));

			  grmVsWTAPartonDeltaRHist[branch][0][subIndex][j]->Fill(TMath::Sqrt(dRap7*dRap7 + dPhi7*dPhi7),TMath::Sqrt(dRap8*dRap8 + dPhi8*dPhi8));
			  grmWTAPartonDeltaRDiffHist[branch][0][subIndex][j]->Fill(TMath::Sqrt(dRap8*dRap8 + dPhi8*dPhi8)-TMath::Sqrt(dRap7*dRap7 + dPhi7*dPhi7));
			}
		    }

		  // Compare to Parton - Higher Order
		  if(branchN == 4 && subIndex != 1 && subIndex != 4 && i < 2)
		    {
		      Double_t jetIndex = jetPartonAssoc(simuJetEvent,4,jet);
		      //cout << jetINdex << endl;
		      StEpSimuParticle *hoParton = simuJetEvent->particle(jetIndex-1);

		      Double_t dRap9 = hoParton->fourMomentum().Rapidity() - jet->rap();
		      Double_t dPhi9 = TVector2::Phi_mpi_pi(hoParton->fourMomentum().Phi() - jet->phi());

		      Double_t dRap10 = hoParton->fourMomentum().Rapidity() - wtaJet.Rapidity();
		      Double_t dPhi10 = TVector2::Phi_mpi_pi(hoParton->fourMomentum().Phi() - wtaJet.Phi());

		      jetPartonDeltaRHist[branch][0][subIndex]->Fill(TMath::Sqrt(dRap9*dRap9 + dPhi9*dPhi9));
		      wtaPartonDeltaRHist[branch][0][subIndex]->Fill(TMath::Sqrt(dRap10*dRap10 + dPhi10*dPhi10));

		      wtaVsJetPartonDeltaRHist[branch][0][subIndex]->Fill(TMath::Sqrt(dRap9*dRap9 + dPhi9*dPhi9),TMath::Sqrt(dRap10*dRap10 + dPhi10*dPhi10));
		      wtaJetPartonDeltaRDiffHist[branch][0][subIndex]->Fill(TMath::Sqrt(dRap10*dRap10 + dPhi10*dPhi10)-TMath::Sqrt(dRap9*dRap9 + dPhi9*dPhi9));

		      for(int j=0; j<6; j++)
			{
			  Double_t dRap11 = hoParton->fourMomentum().Rapidity() - grmJet[j].Rapidity();
			  Double_t dPhi11 = TVector2::Phi_mpi_pi(hoParton->fourMomentum().Phi() - grmJet[j].Phi());

			  grmPartonDeltaRHist[branch][0][subIndex][j]->Fill(TMath::Sqrt(dRap11*dRap11 + dPhi11*dPhi11));

			  grmVsJetPartonDeltaRHist[branch][0][subIndex][j]->Fill(TMath::Sqrt(dRap9*dRap9 + dPhi9*dPhi9),TMath::Sqrt(dRap11*dRap11 + dPhi11*dPhi11));
			  grmJetPartonDeltaRDiffHist[branch][0][subIndex][j]->Fill(TMath::Sqrt(dRap11*dRap11 + dPhi11*dPhi11)-TMath::Sqrt(dRap9*dRap9 + dPhi9*dPhi9));

			  grmVsWTAPartonDeltaRHist[branch][0][subIndex][j]->Fill(TMath::Sqrt(dRap10*dRap10 + dPhi10*dPhi10),TMath::Sqrt(dRap11*dRap11 + dPhi11*dPhi11));
			  grmWTAPartonDeltaRDiffHist[branch][0][subIndex][j]->Fill(TMath::Sqrt(dRap11*dRap11 + dPhi11*dPhi11)-TMath::Sqrt(dRap10*dRap10 + dPhi10*dPhi10));
			}
		    }

		  // Look at High Pt Jets
		  if(jet->pt() >= 10.0)
		    {
		      // Plot Jet Quantities
		      jetPtHist[branch][1][subIndex]->Fill(jet->pt());
		      if(branchN == 0) jetEtaLabHist[branch][1][subIndex]->Fill(jet->eta());
		      if(branchN == 4) jetEtaLabHist[branch][1][subIndex]->Fill(breitToLab.PseudoRapidity());
		      jetRapHist[branch][1][subIndex]->Fill(jet->rap());
		      
		      // Plot WTA Jet Quantities
		      wtaPtHist[branch][1][subIndex]->Fill(wtaJet.Perp());
		      if(branchN == 0) wtaEtaLabHist[branch][1][subIndex]->Fill(wtaJet.PseudoRapidity());
		      if(branchN == 4) wtaEtaLabHist[branch][1][subIndex]->Fill(breitToLabWTA.PseudoRapidity());
		      wtaRapHist[branch][1][subIndex]->Fill(wtaJet.Rapidity());
		      
		      // Plot GRM Jet Quantities
		      for(int j=0; j<6; j++)
			{
			  grmPtHist[branch][1][subIndex][j]->Fill(grmJet[j].Perp());
			  if(branchN == 0) grmEtaLabHist[branch][1][subIndex][j]->Fill(grmJet[j].PseudoRapidity());
			  if(branchN == 4) grmEtaLabHist[branch][1][subIndex][j]->Fill(breitToLabGRM[j].PseudoRapidity());
			  grmRapHist[branch][1][subIndex][j]->Fill(grmJet[j].Rapidity());
			  
			  jetVsGRMPtHist[branch][1][subIndex][j]->Fill(jet->pt(),grmJet[j].Perp());
			  jetVsGRMNumPartHist[branch][1][subIndex][j]->Fill(jet->numberOfJetParticles(),constituents[j].size());
			}

		      // Plot Comparisons
		      jetWTADeltaRHist[branch][1][subIndex]->Fill(TMath::Sqrt(dRap1*dRap1 + dPhi1*dPhi1));
		      jetWTAPDeltaRHist[branch][1][subIndex]->Fill(TMath::Sqrt(dRap2*dRap2 + dPhi2*dPhi2));
		      wtaWTAPDeltaRHist[branch][1][subIndex]->Fill(TMath::Sqrt(dRap3*dRap3 + dPhi3*dPhi3));

		      jetWTADeltaRLogHist[branch][1][subIndex]->Fill(TMath::Log10(TMath::Sqrt(dRap1*dRap1 + dPhi1*dPhi1)));
		      jetWTAPDeltaRLogHist[branch][1][subIndex]->Fill(TMath::Log10(TMath::Sqrt(dRap2*dRap2 + dPhi2*dPhi2)));
		      wtaWTAPDeltaRLogHist[branch][1][subIndex]->Fill(TMath::Log10(TMath::Sqrt(dRap3*dRap3 + dPhi3*dPhi3)));
		      
		      jetWTADeltaRVsPtHist[branch][1][subIndex]->Fill(jet->pt(),TMath::Sqrt(dRap1*dRap1 + dPhi1*dPhi1));
		      if(branchN == 0) jetWTADeltaRVsEtaHist[branch][1][subIndex]->Fill(jet->eta(),TMath::Sqrt(dRap1*dRap1 + dPhi1*dPhi1));
		      if(branchN == 4) jetWTADeltaRVsEtaHist[branch][1][subIndex]->Fill(breitToLab.PseudoRapidity(),TMath::Sqrt(dRap1*dRap1 + dPhi1*dPhi1));
		      
		      jetWTAPDeltaRVsPtHist[branch][1][subIndex]->Fill(jet->pt(),TMath::Sqrt(dRap2*dRap2 + dPhi2*dPhi2));
		      if(branchN == 0) jetWTAPDeltaRVsEtaHist[branch][1][subIndex]->Fill(jet->eta(),TMath::Sqrt(dRap2*dRap2 + dPhi2*dPhi2));
		      if(branchN == 4) jetWTAPDeltaRVsEtaHist[branch][1][subIndex]->Fill(breitToLab.PseudoRapidity(),TMath::Sqrt(dRap2*dRap2 + dPhi2*dPhi2));
		      
		      wtaWTAPDeltaRVsPtHist[branch][1][subIndex]->Fill(jet->pt(),TMath::Sqrt(dRap3*dRap3 + dPhi3*dPhi3));
		      if(branchN == 0) wtaWTAPDeltaRVsEtaHist[branch][1][subIndex]->Fill(jet->eta(),TMath::Sqrt(dRap3*dRap3 + dPhi3*dPhi3));
		      if(branchN == 4) wtaWTAPDeltaRVsEtaHist[branch][1][subIndex]->Fill(breitToLab.PseudoRapidity(),TMath::Sqrt(dRap3*dRap3 + dPhi3*dPhi3));
		      
		      for(int j=0; j<6; j++)
			{
			  Double_t dRap4 = grmJet[j].Rapidity() - jet->rap();
			  Double_t dPhi4 = TVector2::Phi_mpi_pi(grmJet[j].Phi() - jet->phi());
			  
			  Double_t dRap5 = grmJet[j].Rapidity() - wtaJet.Rapidity();
			  Double_t dPhi5 = TVector2::Phi_mpi_pi(grmJet[j].Phi() - wtaJet.Phi());
			  
			  jetGRMDeltaRHist[branch][1][subIndex][j]->Fill(TMath::Sqrt(dRap4*dRap4 + dPhi4*dPhi4));
			  wtaGRMDeltaRHist[branch][1][subIndex][j]->Fill(TMath::Sqrt(dRap5*dRap5 + dPhi5*dPhi5));

			  jetGRMDeltaRLogHist[branch][1][subIndex][j]->Fill(TMath::Log10(TMath::Sqrt(dRap4*dRap4 + dPhi4*dPhi4)));
			  wtaGRMDeltaRLogHist[branch][1][subIndex][j]->Fill(TMath::Log10(TMath::Sqrt(dRap5*dRap5 + dPhi5*dPhi5)));
			  
			  jetGRMDeltaRVsPtHist[branch][1][subIndex][j]->Fill(jet->pt(),TMath::Sqrt(dRap4*dRap4 + dPhi4*dPhi4));
			  if(branchN == 0) jetGRMDeltaRVsEtaHist[branch][1][subIndex][j]->Fill(jet->eta(),TMath::Sqrt(dRap4*dRap4 + dPhi4*dPhi4));
			  if(branchN == 4) jetGRMDeltaRVsEtaHist[branch][1][subIndex][j]->Fill(breitToLab.PseudoRapidity(),TMath::Sqrt(dRap4*dRap4 + dPhi4*dPhi4));
			  
			  wtaGRMDeltaRVsPtHist[branch][1][subIndex][j]->Fill(jet->pt(),TMath::Sqrt(dRap5*dRap5 + dPhi5*dPhi5));
			  if(branchN == 0) wtaGRMDeltaRVsEtaHist[branch][1][subIndex][j]->Fill(jet->eta(),TMath::Sqrt(dRap5*dRap5 + dPhi5*dPhi5));
			  if(branchN == 4) wtaGRMDeltaRVsEtaHist[branch][1][subIndex][j]->Fill(breitToLab.PseudoRapidity(),TMath::Sqrt(dRap5*dRap5 + dPhi5*dPhi5));
			}

		      // Compare to Parton - DIS
		      if(branchN == 0 && subIndex == 4 && i == 0) // Look at highest jet only
			{
			  // Find Init Parton for DIS Jet
			  StEpSimuParticle *disParton = simuJetEvent->particle(disPartIndexAfterRad(simuJetEvent)-1);
			  
			  Double_t dRap6 = disParton->fourMomentum().Rapidity() - jet->rap();
			  Double_t dPhi6 = TVector2::Phi_mpi_pi(disParton->fourMomentum().Phi() - jet->phi());
			  
			  Double_t dRap7 = disParton->fourMomentum().Rapidity() - wtaJet.Rapidity();
			  Double_t dPhi7 = TVector2::Phi_mpi_pi(disParton->fourMomentum().Phi() - wtaJet.Phi());
			  
			  jetPartonDeltaRHist[branch][1][subIndex]->Fill(TMath::Sqrt(dRap6*dRap6 + dPhi6*dPhi6));
			  wtaPartonDeltaRHist[branch][1][subIndex]->Fill(TMath::Sqrt(dRap7*dRap7 + dPhi7*dPhi7));
			  
			  wtaVsJetPartonDeltaRHist[branch][1][subIndex]->Fill(TMath::Sqrt(dRap6*dRap6 + dPhi6*dPhi6),TMath::Sqrt(dRap7*dRap7 + dPhi7*dPhi7));
			  wtaJetPartonDeltaRDiffHist[branch][1][subIndex]->Fill(TMath::Sqrt(dRap7*dRap7 + dPhi7*dPhi7)-TMath::Sqrt(dRap6*dRap6 + dPhi6*dPhi6));
			  
			  for(int j=0; j<6; j++)
			    {
			      Double_t dRap8 = disParton->fourMomentum().Rapidity() - grmJet[j].Rapidity();
			      Double_t dPhi8 = TVector2::Phi_mpi_pi(disParton->fourMomentum().Phi() - grmJet[j].Phi());
			      
			      grmPartonDeltaRHist[branch][1][subIndex][j]->Fill(TMath::Sqrt(dRap8*dRap8 + dPhi8*dPhi8));
			      
			      grmVsJetPartonDeltaRHist[branch][1][subIndex][j]->Fill(TMath::Sqrt(dRap6*dRap6 + dPhi6*dPhi6),TMath::Sqrt(dRap8*dRap8 + dPhi8*dPhi8));
			      grmJetPartonDeltaRDiffHist[branch][1][subIndex][j]->Fill(TMath::Sqrt(dRap8*dRap8 + dPhi8*dPhi8)-TMath::Sqrt(dRap6*dRap6 + dPhi6*dPhi6));
			      
			      grmVsWTAPartonDeltaRHist[branch][1][subIndex][j]->Fill(TMath::Sqrt(dRap7*dRap7 + dPhi7*dPhi7),TMath::Sqrt(dRap8*dRap8 + dPhi8*dPhi8));
			      grmWTAPartonDeltaRDiffHist[branch][1][subIndex][j]->Fill(TMath::Sqrt(dRap8*dRap8 + dPhi8*dPhi8)-TMath::Sqrt(dRap7*dRap7 + dPhi7*dPhi7));
			    }
			}

		      // Compare to Parton - Higher Order
		      if(branchN == 4 && subIndex != 1 && subIndex != 4 && i < 2)
			{
			  Double_t jetIndex = jetPartonAssoc(simuJetEvent,4,jet);
			  //cout << jetINdex << endl;
			  StEpSimuParticle *hoParton = simuJetEvent->particle(jetIndex-1);
			  
			  Double_t dRap9 = hoParton->fourMomentum().Rapidity() - jet->rap();
			  Double_t dPhi9 = TVector2::Phi_mpi_pi(hoParton->fourMomentum().Phi() - jet->phi());
			  
			  Double_t dRap10 = hoParton->fourMomentum().Rapidity() - wtaJet.Rapidity();
			  Double_t dPhi10 = TVector2::Phi_mpi_pi(hoParton->fourMomentum().Phi() - wtaJet.Phi());
			  
			  jetPartonDeltaRHist[branch][1][subIndex]->Fill(TMath::Sqrt(dRap9*dRap9 + dPhi9*dPhi9));
			  wtaPartonDeltaRHist[branch][1][subIndex]->Fill(TMath::Sqrt(dRap10*dRap10 + dPhi10*dPhi10));
			  
			  wtaVsJetPartonDeltaRHist[branch][1][subIndex]->Fill(TMath::Sqrt(dRap9*dRap9 + dPhi9*dPhi9),TMath::Sqrt(dRap10*dRap10 + dPhi10*dPhi10));
			  wtaJetPartonDeltaRDiffHist[branch][1][subIndex]->Fill(TMath::Sqrt(dRap10*dRap10 + dPhi10*dPhi10)-TMath::Sqrt(dRap9*dRap9 + dPhi9*dPhi9));
			  
			  for(int j=0; j<6; j++)
			    {
			      Double_t dRap11 = hoParton->fourMomentum().Rapidity() - grmJet[j].Rapidity();
			      Double_t dPhi11 = TVector2::Phi_mpi_pi(hoParton->fourMomentum().Phi() - grmJet[j].Phi());
			      
			      grmPartonDeltaRHist[branch][1][subIndex][j]->Fill(TMath::Sqrt(dRap11*dRap11 + dPhi11*dPhi11));
			      
			      grmVsJetPartonDeltaRHist[branch][1][subIndex][j]->Fill(TMath::Sqrt(dRap9*dRap9 + dPhi9*dPhi9),TMath::Sqrt(dRap11*dRap11 + dPhi11*dPhi11));
			      grmJetPartonDeltaRDiffHist[branch][1][subIndex][j]->Fill(TMath::Sqrt(dRap11*dRap11 + dPhi11*dPhi11)-TMath::Sqrt(dRap9*dRap9 + dPhi9*dPhi9));
			      
			      grmVsWTAPartonDeltaRHist[branch][1][subIndex][j]->Fill(TMath::Sqrt(dRap10*dRap10 + dPhi10*dPhi10),TMath::Sqrt(dRap11*dRap11 + dPhi11*dPhi11));
			      grmWTAPartonDeltaRDiffHist[branch][1][subIndex][j]->Fill(TMath::Sqrt(dRap11*dRap11 + dPhi11*dPhi11)-TMath::Sqrt(dRap10*dRap10 + dPhi10*dPhi10));
			    } 
			}
		    } // Pt > 10
		  
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
