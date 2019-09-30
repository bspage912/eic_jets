
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

#include "/afs/rhic.bnl.gov/eic/PACKAGES/EicRoot/eic-smear/include/eicsmear/erhic/EventBase.h"
#include "/afs/rhic.bnl.gov/eic/PACKAGES/EicRoot/eic-smear/include/eicsmear/erhic/EventPythia.h"
#include "/afs/rhic.bnl.gov/eic/PACKAGES/EicRoot/eic-smear/include/eicsmear/erhic/Particle.h"

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
  const char* inFileName = argv[4];
  const char* outFileName = argv[5];

  cout << "nevents = " << nevents << endl;
  cout << "Process Events " << beginEvent << " - " << endEvent << endl;
  cout << "inFileName = " << inFileName << endl;
  cout << "outFileName = " << outFileName << endl;

  // Chain for Simu Tree
  TChain* inTree = new TChain("EICTree");
  inTree->Add(inFileName);

  // Setup Input Event Buffer
  erhic::EventPythia* inEvent(NULL);
  inTree->SetBranchAddress("event",&inEvent);

  // Open Output File
  TFile *ofile = TFile::Open(outFileName,"recreate");
  assert(ofile);

  // Create Jet Tree
  TTree *mTree = new TTree("simuJets","SimuJetTrees");

  // Set Jet Tree Structure
  StEpSimuJetEvent *event = 0;
  mTree->Branch("myJets","StEpSimuJetEvent",&event);


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
      event->setX(inEvent->GetX());
      event->setQ2(inEvent->GetQ2());
      event->setY(inEvent->GetY());
      event->setW2(inEvent->GetW2());
      event->setNu(inEvent->GetNu());
      event->setXJB(inEvent->GetXJacquetBlondel());
      event->setQ2JB(inEvent->GetQ2JacquetBlondel());
      event->setYJB(inEvent->GetYJacquetBlondel());
      event->setW2JB(inEvent->GetW2JacquetBlondel());
      event->setXDA(inEvent->GetXDoubleAngle());
      event->setQ2DA(inEvent->GetQ2DoubleAngle());
      event->setYDA(inEvent->GetYDoubleAngle());
      event->setW2DA(inEvent->GetW2DoubleAngle());


      // Create Containers for Different Jet Defs
      //StEpSimuJetDef *akt10Pt = event->newJetDef();
      //StEpSimuJetDef *akt10Pt_hbv = event->newJetDef();
      StEpSimuJetDef *akt10Pt_breit = event->newJetDef();
      StEpSimuJetDef *akt04Pt_breit = event->newJetDef();
      StEpSimuJetDef *akt10Pt_wta_breit = event->newJetDef();
      StEpSimuJetDef *akt04Pt_wta_breit = event->newJetDef();

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

      // AKT Breit
      akt10Pt_breit->setRadius(1.0);
      akt10Pt_breit->setMinPt(0.250);
      akt10Pt_breit->setAlgo(StEpSimuJetDef::AKT);
      akt10Pt_breit->setFrame(StEpSimuJetDef::BREIT);

      akt04Pt_breit->setRadius(0.4);
      akt04Pt_breit->setMinPt(0.250);
      akt04Pt_breit->setAlgo(StEpSimuJetDef::AKT);
      akt04Pt_breit->setFrame(StEpSimuJetDef::BREIT);

      akt10Pt_wta_breit->setRadius(1.0);
      akt10Pt_wta_breit->setMinPt(0.250);
      akt10Pt_wta_breit->setAlgo(StEpSimuJetDef::AKT);
      akt10Pt_wta_breit->setFrame(StEpSimuJetDef::BREIT);

      akt04Pt_wta_breit->setRadius(0.4);
      akt04Pt_wta_breit->setMinPt(0.250);
      akt04Pt_wta_breit->setAlgo(StEpSimuJetDef::AKT);
      akt04Pt_wta_breit->setFrame(StEpSimuJetDef::BREIT);


      // Create Lab and Hadron Boson Jet Vectors
      vector<PseudoJet> particlesPt;
      vector<PseudoJet> particlesPt_hbv;
      vector<PseudoJet> particlesPt_breit;


      // Set Up Boost to Breit Frame
      const Particle* part2 = inEvent->GetTrack(1);
      const Particle* part4 = inEvent->GetTrack(3);

      PseudoJet proton(part2->GetPx(),part2->GetPy(),part2->GetPz(),part2->GetE());
      PseudoJet gamma(part4->GetPx(),part4->GetPy(),part4->GetPz(),part4->GetE());

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

      PseudoJet boosted_proton = boost(proton,boost_vector);
      Double_t phi_p = boosted_proton.phi();
      Double_t theta_p = TMath::ATan2(boosted_proton.perp(),boosted_proton.pz());

      event->setBoost(boost_vector.px(),boost_vector.py(),boost_vector.pz(),boost_vector.e());
      event->setBoostAngles(phi_p,theta_p);


      // Loop over Particles
      for(int j=0; j<inEvent->GetNTracks(); j++)
	{
	  const Particle* inParticle = inEvent->GetTrack(j);
	  
	  Double_t px = inParticle->GetPx();
	  Double_t py = inParticle->GetPy();
	  Double_t pz = inParticle->GetPz();
	  Double_t E = inParticle->GetE();

	  // Define Hadron Boson Frame Kinematics
	  TLorentzVector hadBosVec = inParticle->Get4VectorInHadronBosonFrame();

	  Double_t px_hbv = hadBosVec.Px();
	  Double_t py_hbv = hadBosVec.Py();
	  Double_t pz_hbv = hadBosVec.Pz();
	  Double_t E_hbv = hadBosVec.E();

	  // Define Breit Fram Kinematics
	  fastjet::PseudoJet p_lab(px,py,pz,E);
	  fastjet::PseudoJet boosted_particle = boost(p_lab,boost_vector);
	  boosted_particle = rotateZ(boosted_particle, -phi_p);
	  boosted_particle = rotateY(boosted_particle, -theta_p);

	  // Not Sure if set_user_index can be called twice on same object, so make a copy
	  fastjet::PseudoJet boosted_particlePt = boost(p_lab,boost_vector);
	  boosted_particlePt = rotateZ(boosted_particlePt, -phi_p);
	  boosted_particlePt = rotateY(boosted_particlePt, -theta_p);


	  // Select Particles for Jets
	  if(j>10 && inParticle->GetStatus() == 1 && inParticle->GetParentIndex() != 3)
	    {
	      if(inParticle->GetEta() <= 4.5 && inParticle->GetPt() >= 0.250)
		{
		  fastjet::PseudoJet pPt(px,py,pz,E);
		  pPt.set_user_index(inParticle->GetIndex());
		  particlesPt.push_back(pPt);
		  
		  fastjet::PseudoJet pPt_hbv(px_hbv,py_hbv,pz_hbv,E_hbv);
		  pPt_hbv.set_user_index(inParticle->GetIndex());
		  particlesPt_hbv.push_back(pPt_hbv);

		  boosted_particlePt.set_user_index(inParticle->GetIndex());
		  particlesPt_breit.push_back(boosted_particlePt);
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
      double R_04 = 0.4;
      JetDefinition jet_def_akt_10(antikt_algorithm,R_10);
      JetDefinition jet_def_akt_04(antikt_algorithm,R_04);
      JetDefinition jet_def_akt_wta_10(antikt_algorithm,R_10,fastjet::WTA_pt_scheme,Best);
      JetDefinition jet_def_akt_wta_04(antikt_algorithm,R_04,fastjet::WTA_pt_scheme,Best);

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

      // Cluster in Breit Frame
      ClusterSequence csPt_akt_10_breit(particlesPt_breit, jet_def_akt_10);
      ClusterSequence csPt_akt_04_breit(particlesPt_breit, jet_def_akt_04);
      ClusterSequence csPt_akt_10_wta_breit(particlesPt_breit, jet_def_akt_wta_10);
      ClusterSequence csPt_akt_04_wta_breit(particlesPt_breit, jet_def_akt_wta_04);
      //ClusterSequenceArea csPt_akt_10_breit(particlesPt_breit, jet_def_akt_10, area_def);


      // Lab Frame Jets
      //vector<PseudoJet> jetsPt_akt_10 = sorted_by_pt(csPt_akt_10.inclusive_jets(ptmin));

      // Hadron Boson Frame Jets
      //vector<PseudoJet> jetsPt_akt_10_hbv = sorted_by_pt(csPt_akt_10_hbv.inclusive_jets(ptmin));
      
      // Breit Frame Jets
      vector<PseudoJet> jetsPt_akt_10_breit = sorted_by_pt(csPt_akt_10_breit.inclusive_jets(ptmin));
      vector<PseudoJet> jetsPt_akt_04_breit = sorted_by_pt(csPt_akt_04_breit.inclusive_jets(ptmin));
      vector<PseudoJet> jetsPt_akt_10_wta_breit = sorted_by_pt(csPt_akt_10_wta_breit.inclusive_jets(ptmin));
      vector<PseudoJet> jetsPt_akt_04_wta_breit = sorted_by_pt(csPt_akt_04_wta_breit.inclusive_jets(ptmin));
      

      // Set SubJet Definitions
      //double R_Sub = 0.25;
      double R_Sub = 0.2;
      JetDefinition jet_def_sub(antikt_algorithm,R_Sub);

      // AKT 10 Breit Pt Cut
      for(int i=0; i<jetsPt_akt_10_breit.size(); i++)
	{
	  StEpSimuJet *jet = event->newJet(jetsPt_akt_10_breit[i].pt(),jetsPt_akt_10_breit[i].eta(),jetsPt_akt_10_breit[i].phi(),jetsPt_akt_10_breit[i].e());
	  akt10Pt_breit->addJet(jet);
	  jet->setJetDef(akt10Pt_breit);
	  //jet->setArea(jetsPt_akt_10_breit[i].area());
	  //jet->setAreaError(jetsPt_akt_10_breit[i].area_error());
	  jet->setArea(-1.0);
	  jet->setAreaError(-1.0);

	  vector<PseudoJet> constituents = jetsPt_akt_10_breit[i].constituents();

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

      // AKT 04 Breit Pt Cut
      for(int i=0; i<jetsPt_akt_04_breit.size(); i++)
	{
	  StEpSimuJet *jet = event->newJet(jetsPt_akt_04_breit[i].pt(),jetsPt_akt_04_breit[i].eta(),jetsPt_akt_04_breit[i].phi(),jetsPt_akt_04_breit[i].e());
	  akt04Pt_breit->addJet(jet);
	  jet->setJetDef(akt04Pt_breit);
	  //jet->setArea(jetsPt_akt_10_breit[i].area());
	  //jet->setAreaError(jetsPt_akt_10_breit[i].area_error());
	  jet->setArea(-1.0);
	  jet->setAreaError(-1.0);

	  vector<PseudoJet> constituents = jetsPt_akt_04_breit[i].constituents();

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

      // AKT 10 Breit Pt Cut WTA
      for(int i=0; i<jetsPt_akt_10_wta_breit.size(); i++)
	{
	  StEpSimuJet *jet = event->newJet(jetsPt_akt_10_wta_breit[i].pt(),jetsPt_akt_10_wta_breit[i].eta(),jetsPt_akt_10_wta_breit[i].phi(),jetsPt_akt_10_wta_breit[i].e());
	  akt10Pt_wta_breit->addJet(jet);
	  jet->setJetDef(akt10Pt_wta_breit);
	  //jet->setArea(jetsPt_akt_10_breit[i].area());
	  //jet->setAreaError(jetsPt_akt_10_breit[i].area_error());
	  jet->setArea(-1.0);
	  jet->setAreaError(-1.0);

	  vector<PseudoJet> constituents = jetsPt_akt_10_wta_breit[i].constituents();

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

      // AKT 04 Breit Pt Cut WTA
      for(int i=0; i<jetsPt_akt_04_wta_breit.size(); i++)
	{
	  StEpSimuJet *jet = event->newJet(jetsPt_akt_04_wta_breit[i].pt(),jetsPt_akt_04_wta_breit[i].eta(),jetsPt_akt_04_wta_breit[i].phi(),jetsPt_akt_04_wta_breit[i].e());
	  akt04Pt_wta_breit->addJet(jet);
	  jet->setJetDef(akt04Pt_wta_breit);
	  //jet->setArea(jetsPt_akt_10_breit[i].area());
	  //jet->setAreaError(jetsPt_akt_10_breit[i].area_error());
	  jet->setArea(-1.0);
	  jet->setAreaError(-1.0);

	  vector<PseudoJet> constituents = jetsPt_akt_04_wta_breit[i].constituents();

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

      
      mTree->Fill();
    }
      
  ofile->Write();
  ofile->Close();
  
}  
