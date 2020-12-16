#include "classes/DelphesClasses.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootUtilities.h"
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


int run_analysis()
{
  gSystem->Load("libDelphes");
	const char *inputFile ="/users/pep/alasfarl/HEP_tools/Delphes-3.4.2/hh/HL-LHC14HH-kd750.root";

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");
	TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
	TClonesArray *branchHT = treeReader->UseBranch("ScalarHT");
  // Loop over all events


     TFile* rf = TFile::Open("./kd750-1btag-HLLHC14.root","recreate");
     TTree *OutTree = new TTree("HHSM14","hh   14 TeV events");
     // t1->Print();
     std::vector<PseudoJet> jjets;
     std::vector<PseudoJet> bjets;
     std::vector<PseudoJet> photons;
     Double_t weight;
     // pointer to a vector
     /////////////////////////////////////////////
            Double_t njjet;
             Double_t nbjet;
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
             OutTree->Branch("njjet", &njjet, "njjet/D");
             OutTree->Branch("nbjet", &nbjet, "nbjet/D");
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
             OutTree->Branch("maa", &maa,"maa/D");
             OutTree->Branch("mb1h", &mb1h,"mb1h/D");
             OutTree->Branch("mbbh", &mbbh,"mbbh/D");
             OutTree->Branch("met", &met, "met/D");
             OutTree->Branch("ht", &ht, "ht/D");
             OutTree->Branch("drbamin", &drbamin,"drbamin/D");
             OutTree->Branch("drba1", &drba1,"drba1/D");
             OutTree->Branch("dphiba1", &dphiba1,"dphiba1/D");
             OutTree->Branch("dphibb", &dphibb,"dphibb/D");
            OutTree->Branch("weight", &weight, "weight/D");



            Double_t PartonPTMin = 30.0;
            Double_t PartonEtaMax = 2.5;
            // Other Cuts
            Double_t EtaObsMax = 2.5; // jet selection
            Double_t PTObsMin = 30.; // jet selection
            Double_t PTAObsMin = 20.; // photon jet selection

  for(Long64_t entry = 0; entry < numberOfEntries; ++entry) // loop over events
  {

    treeReader->ReadEntry(entry);
		if (entry!=0) {
			jjets.clear();
			bjets.clear();
			photons.clear();
		}




    // If event contains at least 1 jet

    for (size_t i = 0; i < branchJet->GetEntries(); i++) {
       Jet *jet = (Jet*) branchJet->At(i);
       TLorentzVector dumjet;
       dumjet.SetPtEtaPhiM(jet->PT,jet->Eta,jet->Phi,jet->Mass); // dummy 4 vector to convert from pt eta to p_i and e to feed into the Psedojet
       PseudoJet pjet(dumjet.Px(),dumjet.Py(),dumjet.Pz(),dumjet.E());
       Int_t bflag = jet->BTag;
       pjet.set_user_index(bflag);
       if(jet->PT>PTObsMin && abs(jet->Eta)<EtaObsMax){
				  if(bflag==0)  jjets.push_back(pjet);
	    if(bflag==1) bjets.push_back(pjet);
    }
		// delete jet;
	}
		for (size_t k = 0; k < branchPhoton->GetEntries(); k++) {
			 Photon *ph = (Photon*) branchPhoton->At(k);
			 TLorentzVector dumphoton;
			 dumphoton.SetPtEtaPhiM(ph->PT,ph->Eta,ph->Phi,0.0); // dummy 4 vector to convert from pt eta to p_i and e to feed into the Psedojet
			 PseudoJet p(dumphoton.Px(),dumphoton.Py(),dumphoton.Pz(),dumphoton.E());
			 if(ph->PT>PTAObsMin && abs(ph->Eta)<EtaObsMax){
	        photons.push_back(p);
			 }
			 // delete ph;
		 }
           std::vector<PseudoJet> selected_gammas =  sorted_by_pt(photons);
		             if (selected_gammas.size() <2 ) continue;
								 PseudoJet gamma1 = selected_gammas[0];
				         PseudoJet gamma2 = selected_gammas[1];
								  PseudoJet gammagamma = join(gamma1,gamma2);
								 if ((gammagamma.m() >110.0) && (gammagamma.m() < 140.0)){
					std::vector<PseudoJet> sel_jets  = sorted_by_pt(jjets);
						std::vector<PseudoJet> sel_bjets  = sorted_by_pt(bjets);
										njjet = jjets.size();
										if (bjets.size() <1) continue;
										nbjet = bjets.size();

										PseudoJet bjet1 = sel_bjets[0];
										PseudoJet bjet2(0.0,0.0,0.0,0.0);
										if (nbjet>1) {
											bjet2=sel_bjets[1];
										}

									Double_t delta_rgg = GetDeltaR(gamma1.eta(),gamma1.phi(),gamma2.eta(),gamma2.phi());
					            Double_t delta_rbg11 = GetDeltaR(bjet1.eta(),bjet1.phi(),gamma1.eta(),gamma1.phi());
					            Double_t delta_rbg22 = GetDeltaR(bjet2.eta(),bjet2.phi(),gamma2.eta(),gamma2.phi());
					            Double_t delta_rbg12 = GetDeltaR(bjet1.eta(),bjet1.phi(),gamma2.eta(),gamma2.phi());
					            Double_t delta_rbg21 = GetDeltaR(bjet2.eta(),bjet2.phi(),gamma1.eta(),gamma1.phi());
					            // if( Rbb > 2.0) continue;
					            // if( delta_rgg > 2.0) continue;
					            // if(     (delta_rbg11 < 0.2) ||
					            // (delta_rbg12 < 0.2) ||
					            // (delta_rbg21 < 0.2) ||
					            // (delta_rbg22 < 0.2))continue ;
					              Double_t min_delta_rg = TMath::Min(TMath::Min(delta_rbg11,delta_rbg12),TMath::Min(delta_rbg21,delta_rbg22));
					              PseudoJet bb = join(bjet1,bjet2);
					              PseudoJet HH  = join(bb,gammagamma) ;
												PseudoJet Hb1  = join(bjet1,gammagamma) ;

					                  // Double_t phi11= bjet1.delta_phi_to(gamma1);
					                  // Double_t phi12 =bjet1.delta_phi_to(gamma2);
					                  // Double_t phi21= bjet2.delta_phi_to(gamma1);
					                  // Double_t phi22 =bjet2.delta_phi_to(gamma2);
					                  // Double_t dphibg= TMath::Min(TMath::Min(phi11,phi12) ,TMath::Min(phi21,phi22) );
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
					                  maa=gammagamma.m();
														mbb=bb.m();
														mbbh=HH.m();
														mb1h=Hb1.m();
                            ScalarHT * scht= (ScalarHT*) branchHT->At(0);
					                  ht= scht->HT;
					                  drbamin=min_delta_rg;
					                  drba1=TMath::Min(delta_rbg11,delta_rbg21);
					                  dphiba1 = TMath::Abs(bjet1.delta_phi_to(gamma1));
					                  dphibb=TMath::Abs(bjet2.delta_phi_to(bjet1));
														MissingET *Misset = (MissingET*) branchMissingET->At(0);
														met = Misset->MET;
														Float_t  prog=  (Float_t)entry/(Float_t) numberOfEntries *100;
														printf("\r[Events Analysed  ::  %.1f   %c ]", prog, 37);
					                  OutTree->Fill();
												}
}
OutTree->Print();
OutTree->SetDirectory(rf);
OutTree->Write();
rf->Close();
return 0;
}
