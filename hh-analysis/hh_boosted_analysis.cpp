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
// Watch out however that, by default, the user_index of a particle is
// set to -1 and you may not have control over that (e.g. if you
// compute the jet area using explicit ghosts, the ghosts will have a
// default user_index of -1). For that reason, if one of the particle
// being combined has a user index of -1, we assume it is not b-tagged
// (i.e. we count it as 0 in the recombination)
//
// This will work for native algorithms, but not for all plugins
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
            return "Higgs ";
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
         // progress bar
         // for (size_t ii = 1; ii < 4; ii++) {
         //    for (size_t jj = 0; jj < 7; jj++) {
         std::vector<Int_t> W {
            13,
            14,
            27,
            100
         };
         std::vector<Int_t> vals{
            25,50,100,150,175,250,200
         };
         // (/2d0, 4d0,6d0,8d0,10d0,14d0,20d0 /)
         //(/5d0, 15d0,25d0,50d0,100d0,150d0,200d0 /)
         // 250d0, 500d0,1d3,1.5d3,1.75d3,150d0,2d3
         // Form(%d-%,,)
          std::string dirc1 = "./root_files/strange/cs";
          std::string dirc2 = "./result_root_files/strange/results_cs";
          std::string exte1 ="TeV.root";

         char *x = new char[dirc1.length()+ exte1.length()+50];
         char *y = new char[dirc2.length()+ exte1.length()+50];

         sprintf(x, "%s%d-%d%s", dirc1.c_str(), vals.at(jj),W.at(ii), exte1.c_str() );
         sprintf(y, "%s%d-%d%s", dirc2.c_str(), vals.at(jj),W.at(ii), exte1.c_str() );

         TString infile = x;
         TString outfile = y;

         TFile *f = new TFile(infile);
         TFile* rf = TFile::Open(outfile,"recreate");
         TTree *t1 = (TTree*)f->Get("RootTuple");
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
         //  Histograms
         TH1::SetDefaultSumw2(kTRUE);
         TH1::AddDirectory(kFALSE);
         // //
         TH2::SetDefaultSumw2(kTRUE);
         TH2::AddDirectory(kFALSE);
         // //
         TH3::SetDefaultSumw2(kTRUE);
         TH3::AddDirectory(kFALSE);
         //
         // WARNING  ::  You need to weight histograms if you wan to obtain correct kin distributions
         //
         // b-jet b-jet
         auto hmbb   = new TH1D("hmbb"," B-tagged di-jets invariant mass  ",60,110,140);
         hmbb->GetXaxis()->SetTitle("m_{jj} (b-tag) [GeV]");
         hmbb->GetYaxis()->SetTitle("Entries / 0.5 GeV");
         // gamma gamma
         auto hmgg   = new TH1D("hmgg"," #gamma #gamma invariant mass  ",60,110,140);
         hmgg->GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
         hmgg->GetYaxis()->SetTitle("Entries / 0.5 GeV");
         // jet radius seaperation between the 2 b-jets
         auto hrbb   = new TH1D("hrbb"," #Delta R(bb)  ",60,.1,6.1);
         hrbb->GetXaxis()->SetTitle("#Delta R(bb)");
         hrbb->GetYaxis()->SetTitle("Entries");

         auto hrgaga   = new TH1D("hrgg"," #Delta R(#gamma #gamma)  ",60,.1,6.1);
         hrbb->GetXaxis()->SetTitle("#Delta R(#gamma #gamma)");
         hrbb->GetYaxis()->SetTitle("Entries");

         auto hrbga   = new TH1D("hrbg"," #Delta R(b #gamma)  ",60,.1,6.1);
         hrbb->GetXaxis()->SetTitle("#Delta R(b #gamma)");
         hrbb->GetYaxis()->SetTitle("Entries");

         auto hmbjet = new TH1D("hmbjet","number of tagged b jets",10,0,10);


         //////////// Important Histograms////////////
         auto hmhh   = new TH1D("hmhh","  Truth Di-Higgs invariant mass  ",20,250,1000);
         hmhh->GetXaxis()->SetTitle("m_{h h} [GeV]");
         hmhh->GetYaxis()->SetTitle("Entries");
         auto hmbbgg   = new TH1D("hmbbgg"," m_{bb #gamma #gamma}  ",20,250,1000);
         auto hpthh = new TH1D("hpthh"," truth ",10,10,700);
         auto hptbbgg = new TH1D("hptbbgg","  ",10,10,700);
         //////////////////////////////////////////////

         // all entries and fill the histograms
         Long64_t nentries =t1->GetEntries();
         std::cout << nentries << std::endl;

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
               if( (id==11 ) || (id==13 )|| (id==15 ) ) {
                  if( (p.pt() > 20.0)) flag_event = 1 ;
               }
               if( (id==12 ) || (id==14 )|| (id==16 ) ) continue;
               if(j== 14) {
                  Int_t id2 = pid->at(j+1);
                  Double_t px2 =vpx->at(j+1);
                  Double_t py2 =vpy->at(j+1);
                  Double_t pz2 =vpz->at(j+1);
                  Double_t e2 =ve->at(j+1);
                  PseudoJet p2(px2,py2,pz2,e2);
                  PseudoJet pgg = join(p,p2);
                  // hmgg ->Fill(pgg.m(),weight) ;
               }
               if(j== 20) {
                  Int_t id2 = pid->at(j+1);
                  Double_t px2 =vpx->at(j+1);
                  Double_t py2 =vpy->at(j+1);
                  Double_t pz2 =vpz->at(j+1);
                  Double_t e2 =ve->at(j+1);
                  PseudoJet p2(px2,py2,pz2,e2);
                  PseudoJet pgg = join(p,p2);
                  // hmgg ->Fill(pgg.m(),weight) ;
               }
               if( std::abs(p.eta())>5.0 ) continue ;
               // B-tag particles


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
                  Final_state_particles.push_back(p);
               }

               particles.push_back(p);
               // this vector will contain the whole event particles,
            }
            // end loop over particle
            if(flag_event)  continue;
            // obtain the mhh distribution
            Selector   sel_higgs = SelectorIsHiggs();
            std::vector<PseudoJet>  Higgses =  sel_higgs(particles);
            if(Higgses.size() < 2) continue;
            PseudoJet H1 = Higgses[0];
            PseudoJet H2 = Higgses[1];
            PseudoJet HHorg = join(H1,H2);
            hmhh-> Fill(HHorg.m(),weight);
            hpthh->Fill(H1.pt(),weight);
            //    }
            // }
            Selector   sel_gamma = SelectorIsGamma();
            // ATLAS CUTS
            Selector select_gammapt = SelectorPtMin(20.0);
            Selector atlas_select_gammaeta = SelectorAbsEtaMax(2.37);
            Selector select_gammaeta = SelectorAbsEtaMax(3.0);
            Selector exclude_crackreigon = SelectorAbsEtaRange(1.37,1.52);
            //
            Selector select_gamma = sel_gamma && select_gammapt && select_gammaeta;
            std::vector<PseudoJet>  gammas =  select_gamma(Final_state_particles);
            std::vector<PseudoJet> selected_gammas =  sorted_by_pt(gammas);
///////////////////////////  Mass drop tagger /////////////////////////////////////
            //
            // double R = 2;  // select wide jet radius for the mass drop algorithm
            // FlavourRecombiner flav_recombiner; // for tracking flavour
            // JetDefinition jet_def(cambridge_algorithm, R, &flav_recombiner);
            // // run the jet finding; find the hardest jet
            // ClusterSequence cs(Final_state_particles, jet_def);
            // vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
            // //
            // // Selector btagger =     SelectorIsBeauty();
            // Selector sel_2hardest=   SelectorAbsRapMax(2.5)&& SelectorPtMin(25.0);
            // vector<PseudoJet> sel_jets = sel_2hardest(jets);
            // vector<PseudoJet> sel_bjets;
            // sel_bjets.clear();
            // for (size_t k = 0; k < sel_jets.size(); k++) {
            //    if (sel_jets[k].user_index() >0) {
            //       sel_bjets.push_back(sel_jets[k]);
            //    }
            // }
            //
            // MassDropTagger md_tagger(0.667, 0.09);
            //  PseudoJet btagged = md_tagger(jets[0]);
            // // //  CASubJetTagger ca_tagger;
            //
            //
            // if (btagged == 0) continue;
            //  PseudoJet parent1 = btagged.pieces()[0];
            //  PseudoJet parent2 = btagged.pieces()[1];
            // double   Rbb = parent1.delta_R(parent2);
            //
            //
            // // next we "filter" it, to remove UE & pileup contamination
            // //----------------------------------------------------------
            //   double   Rfilt =  min(Rbb/2,0.3); // somewhat arbitrary choice
            //   unsigned nfilt = 2;               // number of pieces we'll take
            // //
            // //
            //   Filter filter(JetDefinition(cambridge_algorithm, Rfilt, &flav_recombiner),
            //                 SelectorNHardest(nfilt));
            //   PseudoJet filtered = filter(btagged);
            //    PseudoJet bjet1=      filtered.pieces()[0];
            //    PseudoJet bjet2=      filtered.pieces()[1];
///////////////////////////////////////////ORDINARY tagger ////////////////

               double R = 0.4;  // select wide jet radius
               //for the mass drop algorithm use 1.2
               FlavourRecombiner flav_recombiner; // for tracking flavour

               JetDefinition jet_def(antikt_algorithm, R, &flav_recombiner);
               // run the jet finding; find the hardest jet
               ClusterSequence cs(Final_state_particles, jet_def);
               vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
               //SelectorAbsRapMax(2.5)&&
               Selector btagger =     SelectorIsBeauty();
               Selector sel_hardest=    SelectorPtMin(20.0);
               vector<PseudoJet> sel_jets = sel_hardest(jets);
               vector<PseudoJet> sel_bjets;

               // add b-tagged jets (TRUTH)
               sel_bjets.clear();
               for (size_t k = 0; k < sel_jets.size(); k++) {
                  if (sel_jets[k].user_index() >0) {
                     sel_bjets.push_back(sel_jets[k]);
                  }
               }
                hmbjet->Fill(sel_bjets.size(),weight);
               // if(sel_jets.size() > 5) continue;
               // only use events with 2 b-tagged jets
               // if(sel_bjets.size() != 2) continue;

               // sort b-jets and select the 2 hardest
               vector<PseudoJet> btagged =    sorted_by_pt(sel_bjets);



               PseudoJet bjet1 = btagged[0];
               PseudoJet bjet2 = btagged[1];
///////////////////////////////////////////////////////////
            PseudoJet gamma1 = selected_gammas[0];
            PseudoJet gamma2 = selected_gammas[1];
            Double_t delta_rgg = GetDeltaR(gamma1.eta(),gamma1.phi(),gamma2.eta(),gamma2.phi());
            Double_t delta_rbg11 = GetDeltaR(bjet1.eta(),bjet1.phi(),gamma1.eta(),gamma1.phi());
            Double_t delta_rbg22 = GetDeltaR(bjet2.eta(),bjet2.phi(),gamma2.eta(),gamma2.phi());
            Double_t delta_rbg12 = GetDeltaR(bjet1.eta(),bjet1.phi(),gamma2.eta(),gamma2.phi());
            Double_t delta_rbg21 = GetDeltaR(bjet2.eta(),bjet2.phi(),gamma1.eta(),gamma1.phi());
            // if( Rbb > 2.0) continue;
            // if( delta_rgg > 2.0) continue;
            if((delta_rbg11 < 0.2) ||
            (delta_rbg12 < 0.2) ||
            (delta_rbg21 < 0.2) ||
            (delta_rbg22 < 0.2))continue ;
            Double_t min_delta_rg = TMath::Min(TMath::Min(delta_rbg11,delta_rbg12),TMath::Min(delta_rbg21,delta_rbg22));

            PseudoJet gammagamma = join(gamma1,gamma2);
            PseudoJet bb = join(bjet1,bjet2);
            PseudoJet HH  = join(bb,gammagamma) ;
            if (gamma1.pt() > gamma2.pt()) {
               if(gamma1.pt()<20.0)  continue ;
            }
            if (gamma1.pt() < gamma2.pt()) {
               if(gamma2.pt()<20.0) continue ;
            }
         ////////////////////////////
         if (bjet1.pt() > bjet2.pt()) {
            if(bjet1.pt()<30.0)  continue ;
         }
         if (bjet1.pt() < bjet2.pt()) {
            if(bjet2.pt()<30.0) continue ;
         }

            if ((gammagamma.m() >110.0) && (gammagamma.m() < 140.0)){
                  if((bjet1.user_index() != 0) || (bjet2.user_index()) != 0)
                  {
                     hptbbgg->Fill(HH.pt(),weight);
                     hmbbgg->Fill(HH.m(),weight);
                     hmbb->Fill(bb.m(),weight);
                     hmgg->Fill(gammagamma.m(),weight);
                     hetagg ->Fill(gammagamma.eta(),weight);
                     hrbb->Fill(Rbb,weight);
                     hrgaga->Fill(delta_rgg,weight);
                     hrbga->Fill(min_delta_rg,weight);
                     Double_t max_pt_b = TMath:: Max(bjet1.pt(),bjet2.pt());
                     Double_t min_pt_b = TMath:: Min(bjet1.pt(),bjet2.pt());
                     Double_t max_pt_g = TMath:: Max(gamma1.pt(),gamma2.pt());
                     Double_t min_phi_b = TMath:: Min(bjet1.phi(),bjet2.phi());



                  }

         }




            Float_t  prog=  (Float_t)i/(Float_t) nentries *100;
            // std::cout << weight << std::endl;
            printf("\r[Events Analysed  ::  %.1f   %c ]", prog, 37);
         } // end loop over events


         // Double_t drbb_scale =1.0/hrbb->GetBinWidth(1)/(Float_t)nentries/21.337;
         // hmhh->Scale(scale/hmhh->GetBinWidth(1));
         // hetahh->Scale(scale/hetahh->GetBinWidth(1)) ;
         // hpthh->Scale(scale/hpthh->GetBinWidth(1)) ;
         // hmbbgg->Scale(scale/hmbbgg->GetBinWidth(1));
         // hmbbgg_mistag->Scale(scale/hmbbgg->GetBinWidth(1));
         // hptbbgg->Scale(scale/hptbbgg->GetBinWidth(1));
         // hptbbgg_mistag->Scale(scale/hptbbgg->GetBinWidth(1));
         // Double_t mxs =hmhh->Integral(-1,21,"width" );

         ///////////////////////////////////////////
         // Write Histograms to file
         // hrbb ->SetDirectory(rf);
         // hrbb ->Write();
         hmhh ->SetDirectory(rf);
         hmhh ->Write();
         // hdrbbgg->SetDirectory(rf);
         // hdrbbgg ->Write();
         hpthh->SetDirectory(rf);
         hpthh ->Write();
         hptbbgg->SetDirectory(rf);
         hptbbgg ->Write();
         hetahh->SetDirectory(rf);
         hetahh ->Write();
         hmbbgg->SetDirectory(rf);
         hmbbgg ->Write();
         // hmh->SetDirectory(rf);
         // hmh ->Write();
         rf->Close();
         printf( " XS  = %.6f\n",mxs);
   //    }
   // }
         return 0;
      }
