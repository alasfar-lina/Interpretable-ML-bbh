// @(#)root/net / Author: Lina alasfar 10.04.2019
// use rootlogon.cc first in order to load the Fastjet Library
//////////////////////////////////////////////////////////////////////////////
//               Double Higgs production analysis program                   //
//             This program uses Fastjet for jet clustering                 //
//             Analysis of hh -> gamma gmamma  2 b-tagged jets              //
//////////////////////////////////////////////////////////////////////////////
////ROOT includes////
#include "TString.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include"TROOT.h"
#include"TFile.h"
#include"TTree.h"
//////// C ++ STL  includes
#include <cstdio>
#include<iostream>
#include <cmath>
#include<algorithm>
//// some useful functions ///
#include"./analysis_tools.hpp"
///// Fastjet includes /////
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/tools/Filter.hh"
///
using namespace fastjet;
///

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




      //
      int run_analysis() {


         std::string x = "./root_files/SM-ggF-0TeV.root ";
         std::string y= "./result_SM-ggF-14TeV.root";
         TString infile = x;
         TString outfile = y;


         TFile *f = new TFile(infile);
         TFile* rf = TFile::Open(outfile,"recreate");
         TTree *t1 = (TTree*)f->Get("RootTuple");
         TTree *OutTree = new TTree("HHSM14","hh SM 14 TeV events");
         // t1->Print();
         std::vector<PseudoJet> particles;
         std::vector<PseudoJet> Final_state_particles;
         Double_t weight;
         // pointer to a vector
         std::vector<double>* vpx =0 ;
         std::vector<double>* vpy =0 ;
         std::vector<double>* vpz =0 ;
         std::vector<double>* ve =0 ;
         std::vector<int>* pid =0;
         std::vector<int>* mother_index =0;
         // acccessing branches
         t1->SetBranchAddress("Px",&vpx);
         t1->SetBranchAddress("Py",&vpy);
         t1->SetBranchAddress("Pz",&vpz);
         t1->SetBranchAddress("E",&ve);
         t1->SetBranchAddress("PID",&pid);
         t1->SetBranchAddress("mother_index",&mother_index);
         t1->SetBranchAddress("weight",&weight); // event weight in fb  ( d \sigma)
         /////////////////////////////////////////////
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
                OutTree->Branch("weight", &weight, "weight/D");


         // all entries and fill the histograms
         Long64_t nentries =t1->GetEntries();
         std::cout << nentries << std::endl;
       Double_t n_tot = 0.0; Double_t  n_trig = 0.0;
         for (Long64_t i=0;i<nentries;i++) {  // loop over events
            t1->GetEntry(i);
            Int_t flag_event = 0;
            particles.clear();
            // Int_t  flag[pid ->size()];
            Final_state_particles.clear();

            for (size_t j = 12; j < pid ->size(); j++) { // loop over particle
               Int_t id = pid->at(j);
               Int_t mi = mother_index->at(j);
               Double_t px =vpx->at(j);
               Double_t py =vpy->at(j);
               Double_t pz =vpz->at(j);
               Double_t e =ve->at(j);


               // an event with  particles: px py pz  E and PDG ID
               PseudoJet p(px,py,pz,e);
               p.set_user_info(new MyUserInfo(id));
               // veto events with leptons having low eta nd high pt
               if( (id==11 ) || (id==13 )|| (id==15 ) ) {
                  if( (p.pt() > 20.0) && std::abs(p.eta())>3.0) flag_event = 1 ;
               }
         //
               // remove neutrinos
               if( (id==12 ) || (id==14 )|| (id==16 ) ) continue;


               // remove particles outside the acceptence reigon
               if( std::abs(p.eta())>5.0 ) continue;


   // B-tag particles (3 generations)
               if(checkFinal(id))
               {
                  Int_t mother_id = pid->at(mi);
                  Int_t gmother_id = pid->at( mother_index->at(mi));
                  Int_t ggmother_id = pid->at(mother_index->at(mother_index->at(mi)));
                  Int_t gggmother_id = pid->at(mother_index->at(mother_index->at(mother_index->at(mi))));
                  Int_t btag = Btag(mother_id) ;
                  Int_t btag_mother = Btag(gmother_id) ;
                  Int_t btag_gmother = Btag(ggmother_id) ;
                  Int_t btag_ggmother = Btag(gggmother_id) ;

                  p.set_user_index(((btag | btag_mother )| btag_gmother) |btag_ggmother);
                  // Int_t tag =  (btag | btag_mother )| btag_gmother;
                  Final_state_particles.push_back(p);
                  // this vector will contain the  event final state  particles,
               }

               particles.push_back(p);
               // this vector will contain the whole event particles,
            }
            // end loop over particle
            n_tot = n_tot +1 ;
            if(flag_event)  continue;  // to preform the leptonic veto

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
            n_trig = n_trig +1 ;

            // obtain the mhh distribution
            Selector   sel_higgs = SelectorIsHiggs();
            std::vector<PseudoJet>  Higgses =  sel_higgs(particles);
            if(Higgses.size() < 2) continue;
            // for (size_t i = 0; i < Higgses.size(); i++) {
            //    for (size_t j = i+1; j < Higgses.size(); j++) {
            PseudoJet H1 = Higgses[0];
            PseudoJet H2 = Higgses[1];
            PseudoJet HHorg = join(H1,H2);
            //    }
            // }




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
            if(sel_bjets.size() != 2) continue;

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


            Float_t  prog=  (Float_t)i/(Float_t) nentries *100;
            printf("\r[Events Analysed  ::  %.1f   %c ]", prog, 37);
         } // end loop over events
         ///////////////// Normalise Histograms



         OutTree->Print();
         OutTree->SetDirectory(rf);
         OutTree->Write();
         rf->Close();
         return 0;
      }
