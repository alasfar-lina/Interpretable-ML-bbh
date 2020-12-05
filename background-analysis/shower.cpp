#include "Pythia8/Pythia.h"
#include "TSystem.h"
#include "TString.h"
#include "TMath.h"
#include"TROOT.h"
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm> // for finding elements in a vector
///// Fastjet includes /////
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/tools/Filter.hh"
///
#include"./analysis_tools.hpp"
///
using namespace fastjet;
using namespace Pythia8;


/////////////////////  class add pgd-id to the jets //////////////////////////
class MyUserInfo : public PseudoJet::UserInfoBase{
public:
   //  - pdg_id        the PDG id of the particle
   MyUserInfo(const Int_t & pdg_id_in ):
   _pdg_id(pdg_id_in){}
   /// access to the PDG id
   Int_t pdg_id() const { return _pdg_id;}
protected:
   Int_t _pdg_id;         // the associated pdg id
};

///////////////////// class to  b-tag the jets //////////////////////////

class SW_IsBeauty : public SelectorWorker{
public:
   // default ctor
   SW_IsBeauty(){}
   // the selector's description
   std::string description() const{
      return "beauty quarks ";
   }
   // keeps the ones that have the pdg id of the b
   bool pass(const PseudoJet &p) const{
      const int & btag = p.user_index();
      return (btag == 1); // beauty or anti-beauty quarks
   }
};
// the function that allows to write simply the selector
Selector SelectorIsBeauty(){
   return Selector(new SW_IsBeauty());
}
//----------------------------------------------------------------------
// set up a class to give standard (by default E-scheme)
// recombination, with additional tracking of flavour information in
// the user_index.
//
// b-tagged particles are assumed to have their user_index set to 1,
// and other particles should have user_index to 0.
//
//----------------------------------------------------------------------

typedef JetDefinition::DefaultRecombiner DefRecomb;
class FlavourRecombiner : public  DefRecomb {
public:
   FlavourRecombiner(RecombinationScheme recomb_scheme = E_scheme) :
   DefRecomb(recomb_scheme) {};

   virtual std::string description() const {
      return DefRecomb::description()+" (with user index addition)";}

      /// recombine pa and pb and put result into pab
      virtual void recombine(const PseudoJet & pa, const PseudoJet & pb,
         PseudoJet & pab) const {
            DefRecomb::recombine(pa,pb,pab);
            // Note: see the above discussion for the fact that we consider
            // negative user indices as "0"
            pab.set_user_index(max(pa.user_index(),0) + max(pb.user_index(),0));
         }
      };



      ///////////////////// class to  select gammas //////////////////////////

      class SW_Isgamma : public SelectorWorker{
      public:
         // default ctor
         SW_Isgamma(){}
         // the selector's description
         std::string description() const{
            return "photons";
         }
         // keeps the ones that have the pdg id of the b
         bool pass(const PseudoJet &p) const{
            const int & pdgid = p.user_info<MyUserInfo>().pdg_id();
            return (pdgid == 22); // beauty or anti-beauty quarks
         }
      };
      // the function that allows to write simply the selector
      Selector SelectorIsGamma(){
         return Selector(new SW_Isgamma());
      }


      ///////////////////// class to  select Higgses //////////////////////////

      class SW_IsHiggs : public SelectorWorker{
      public:
         // default ctor
         SW_IsHiggs(){}
         // the selector's description
         std::string description() const{
            return "Higgs";
         }
         // keeps the ones that have the pdg id of the b
         bool pass(const PseudoJet &p) const{
            const int & pdgid = p.user_info<MyUserInfo>().pdg_id();
            return (pdgid == 25|| (pdgid == 35)); // beauty or anti-beauty quarks
         }
      };
      // the function that allows to write simply the selector
      Selector SelectorIsHiggs(){
         return Selector(new SW_IsHiggs());
      }



int run_analysis() {
   TFile* rf = TFile::Open("./result_bbaa-14TeV.root","recreate");
   TTree *OutTree = new TTree("HHSM14","hh SM 14 TeV events");
   Pythia pythia;
   pythia.readString("PhaseSpace:pTHatMin = 20.");
   pythia.readString("Beams:frameType = 4");
   pythia.readString("Beams:LHEF = ./events.lhe");
///////////
pythia.readString("310:mayDecay = off");  //  K0s
pythia.readString("3112:mayDecay = off"); //  Sigma-
pythia.readString("3122:mayDecay = off"); //  Lambda
pythia.readString("3222:mayDecay = off"); //  Sigma+
pythia.readString("3312:mayDecay = off"); //  Xi-
pythia.readString("3322:mayDecay = off"); //  Xi0
pythia.readString("3334:mayDecay = off"); //  Omega-

 pythia.init();

 std::vector<PseudoJet> particles;
 std::vector<PseudoJet> Final_state_particles;
 Double_t weight;


 Int_t nbjet;
Double_t ptb1;
Double_t ptb2;
Double_t pta1 ;
Double_t pta2 ;
Double_t ptaa;
Double_t etab1;
Double_t etab2;
Double_t etaa1;
Double_t etaa2;
Double_t etaaa;
Double_t mbb;
Double_t maa;
Double_t mb1h;
Double_t mbbh;
Double_t met;
Double_t ht;
Double_t drbamin;
Double_t drba1;
Double_t dphiba1;
Double_t dphibb;

//////
OutTree->Branch("nbjet", &nbjet, "nbjet/I");
OutTree->Branch("ptb1", &ptb1, "ptb1/D");
OutTree->Branch("ptb2", &ptb2, "ptb2/D");
OutTree->Branch("pta1", &pta1, "pta1/D");
OutTree->Branch("pta2", &pta2, "pta2/D");
OutTree->Branch("ptaa", &ptaa, "ptaa/D");
OutTree->Branch("etab1", &etab1, "etab1/D");
OutTree->Branch("etab2", &etab2,"etab2/D");
OutTree->Branch("etaa1", &etaa1, "etaa1/D");
OutTree->Branch("etaa2", &etaa2,"etaa2/D");
OutTree->Branch("etaaa", &etaaa, "etaaa/D");
OutTree->Branch("mbb", &mbb,"mbb/D");
// OutTree->Branch("met", &met, "met/D”);
OutTree->Branch("ht", &ht, "ht/D");
OutTree->Branch("drbamin", &drbamin,"drbamin/D");
OutTree->Branch("drba1", &drba1,"drba1/D");
OutTree->Branch("dphiba1", &dphiba1,"dphiba1/D");
OutTree->Branch("dphibb", &dphibb,"dphibb/D");
// OutTree->Branch("weight", &weight, "weight/D");

int nevt=20000;
int Ninel =0;
for (int ievt=0; ievt<nevt; ++ievt) {
  while(!pythia.next());  // iterate until the next succesful event
  if (pythia.info.atEndOfFile()) break;

  // if((ievt+1)%100==0) fprintf(stderr," PYTHIA Status: %8d/%d events generated\r",ievt+1,nevt);

  int proc = pythia.info.code();
  int ntrk = pythia.event.size();
  Int_t flag_event = 0;
  particles.clear();
  // Int_t  flag[pid ->size()];
  Final_state_particles.clear();
  // skip elastic scattering events

  // if(proc==102) continue;
  if (ntrk != 0) {
  Ninel++ ;
  }

  // Find number of all final charged particles and fill histogram.

  for(int i=0; i <ntrk ; ++i) {



     Int_t id = pythia.event[i].id();
     Int_t mi = pythia.event[i].mother1();
     Double_t px =pythia.event[i].px();
     Double_t py =pythia.event[i].py() ;
     Double_t pz =pythia.event[i].pz();
     Double_t e =pythia.event[i].e();


    // accept only charged final state particles


    // check that final state particles are what they should be

    PseudoJet p(px,py,pz,e);
    p.set_user_info(new MyUserInfo(id));

     // pre cuts
    // double absp = pythia.event[i].pAbs();
    // if(absp  <2.) continue;
    double eta = pythia.event[i].eta();
    if(eta>5.) continue;
 if(!pythia.event[i].isFinal()  ) continue;
 Int_t mother_id = pythia.event[mi].id();
Int_t btag = Btag(mother_id) ;
p.set_user_index(btag);
 Final_state_particles.push_back(p);
 if(!checkFinal(id)) {
printf("unexpected final state particle id=%d\n",pythia.event[i].id());
}
} // end particle loop

Selector   sel_gamma = SelectorIsGamma();
// ATLAS CUTS
Selector select_gammapt = SelectorPtMin(20.0);
Selector atlas_select_gammaeta = SelectorAbsEtaMax(2.37);
Selector select_gammaeta = SelectorAbsEtaMax(3.0);
Selector exclude_crackreigon = SelectorAbsEtaRange(1.37,1.52);
//
Selector select_gamma =      sel_gamma && select_gammapt && select_gammaeta;
// && !exclude_crackreigon; for ATLAS
//

std::vector<PseudoJet>  gammas =  sel_gamma(Final_state_particles);
// The sorting so one can choose the hardest 2 gammas
std::vector<PseudoJet> selected_gammas =  sorted_by_pt(gammas);
if (gammas.size() <2 ) continue;







double R = 0.4;  // select wide jet radius
//for the mass drop algorithm use 1.2
FlavourRecombiner flav_recombiner; // for tracking flavour

JetDefinition jet_def(antikt_algorithm, R, &flav_recombiner);
// run the jet finding; find the hardest jet
ClusterSequence cs(Final_state_particles, jet_def);
vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
//
// Selector btagger =     SelectorIsBeauty();
Selector sel_hardest=   SelectorAbsRapMax(3.0)&& SelectorPtMin(20.0);
vector<PseudoJet> sel_jets = sel_hardest(jets);
vector<PseudoJet> sel_bjets;

// add b-tagged jets (TRUTH)
sel_bjets.clear();
for (size_t k = 0; k < sel_jets.size(); k++) {
   if (sel_jets[k].user_index() >0) {
      sel_bjets.push_back(sel_jets[k]);
   }
}
// if(sel_jets.size() > 5) continue;
// only use events with 2 b-tagged jets
nbjet = sel_bjets.size();
if(sel_bjets.size() < 2) continue;

// sort b-jets and select the 2 hardest
vector<PseudoJet> btagged =    sorted_by_pt(sel_bjets);



PseudoJet bjet1 = btagged[0];
PseudoJet bjet2 = btagged[1];
Double_t   Rbb = bjet1.delta_R(bjet2);



// next we "filter" it, to remove UE & pileup contamination
//----------------------------------------------------------



   if(bjet1.pt()<30.0) continue ;



PseudoJet gamma1 = selected_gammas[0];
PseudoJet gamma2 = selected_gammas[1];
Double_t delta_rgg = GetDeltaR(gamma1.eta(),gamma1.phi(),gamma2.eta(),gamma2.phi());
Double_t delta_rbg11 = GetDeltaR(bjet1.eta(),bjet1.phi(),gamma1.eta(),gamma1.phi());
Double_t delta_rbg22 = GetDeltaR(bjet2.eta(),bjet2.phi(),gamma2.eta(),gamma2.phi());
Double_t delta_rbg12 = GetDeltaR(bjet1.eta(),bjet1.phi(),gamma2.eta(),gamma2.phi());
Double_t delta_rbg21 = GetDeltaR(bjet2.eta(),bjet2.phi(),gamma1.eta(),gamma1.phi());
// if( Rbb > 2.0) continue;
// if( delta_rgg > 2.0) continue;
if(     (delta_rbg11 < 0.2) ||
(delta_rbg12 < 0.2) ||
(delta_rbg21 < 0.2) ||
(delta_rbg22 < 0.2))continue ;
   Double_t min_delta_rg = TMath::Min(TMath::Min(delta_rbg11,delta_rbg12),TMath::Min(delta_rbg21,delta_rbg22));
PseudoJet gammagamma = join(gamma1,gamma2);
PseudoJet bb = join(bjet1,bjet2);
PseudoJet HH  = join(bb,gammagamma) ;
PseudoJet Hb1  = join(bjet1,gammagamma) ;
   if(gamma1.pt()<20.0 ) continue ; // for ATLAS let the upper bound be 50


if ((gammagamma.m() >110.0) && (gammagamma.m() < 140.0)){
      Double_t phi11= bjet1.delta_phi_to(gamma1);
      Double_t phi12 =bjet1.delta_phi_to(gamma2);
      Double_t phi21= bjet2.delta_phi_to(gamma1);
      Double_t phi22 =bjet2.delta_phi_to(gamma2);
      Double_t dphibg= TMath::Min(TMath::Min(phi11,phi12) ,TMath::Min(phi21,phi22) );
      /////
      ptb1=bjet1.pt();
      ptb2=bjet2.pt();;
      pta1=gamma1.pt() ;
      pta2 =gamma2.pt();
      ptaa= gammagamma.pt();
      etab1=bjet1.eta();
      etab2=bjet2.eta();
      etaa1=gamma1.eta();
      etaa2=gamma2.eta();
      etaaa=gammagamma.eta();
      mbb= bb.m();
      maa=gammagamma.m();
      mb1h=Hb1.m();
      mbbh=HH.m();
      ht= bjet1.mperp()+bjet2.mperp()+gamma1.mperp()+gamma2.mperp();
      drbamin=min_delta_rg;
      drba1=TMath::Min(delta_rbg11,delta_rbg21);
      dphiba1 =dphibg;
      dphibb=bjet2.delta_phi_to(bjet1);

      OutTree->Fill();

}


Float_t  prog=  (Float_t)ievt/(Float_t) nevt *100;
printf("\r[Events Analysed  ::  %.1f   %c ]", prog, 37);
  }// end loop over events
///////////////// Normalise Histograms



OutTree->Print();
OutTree->SetDirectory(rf);
OutTree->Write();
rf->Close();
return 0;
}
