void rootlogon()
//=======================================================================
// setup root to use fastjet
//
// usage> root -l -n rootlogon.cc
//=======================================================================
{
  // specify the path to the fastjet version to be used

  const char *fastjet = "/users/pep/alasfarl/fastjet-install";

  // make pythia available in ROOT

  gSystem->Setenv("FASTJET",fastjet);
  // gSystem->Setenv("PYTHIA8DATA",Form("%s/share/Pythia8/xmldoc",pythia));
  gInterpreter->AddIncludePath("/users/pep/alasfarl/fastjet-install/include");
  int flag = gSystem->Load("$FASTJET/lib64/libfastjet.so");
  gSystem->Load("$FASTJET/lib64/libfastjettools.so");
  if(flag <0) printf("failed to load %s\n",fastjet);
  if(flag==0) printf("loaded %s\n",fastjet);
  if(flag==1) printf("already loaded %s\n",fastjet);



const char *delphes = "/users/pep/alasfarl/HEP_tools/Delphes-3.4.2";
  gSystem->Setenv("DELPHES",delphes);
  gInterpreter->AddIncludePath("/users/pep/alasfarl/HEP_tools/Delphes-3.4.2");
  gInterpreter->AddIncludePath("/users/pep/alasfarl/HEP_tools/Delphes-3.4.2/external");
  int flag2 = gSystem->Load("$DELPHES/libDelphes.so");
  gSystem->Load("$DELPHES/libDelphes.so");
  if(flag2 <0) printf("failed to load %s\n",delphes);
  if(flag2==0) printf("loaded %s\n",delphes);
  if(flag2==1) printf("already loaded %s\n",delphes);

// style things

Int_t TimeNewRoman        = 132;  // Old LHCb style: 62;
 // Line thickness
 Double_t lhcbWidth    = 2.00; // Old LHCb style: 3.00;
 // Text size
 Double_t lhcbTSize    = 0.04;

 // gStyle->SetOptTitle(kFALSE);  // do  not show title
 //  gStyle->SetOptStat(0);       // do not show stats

  gStyle->SetTextFont(TimeNewRoman);
  gStyle->SetPalette(kCandy);
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatColor(10);
  gStyle->SetStatX(0.9);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.2);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetFrameFillColor(10);
  gStyle->SetCanvasColor(10);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetLegendBorderSize(1);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetTextSize(lhcbTSize);
  gStyle->SetLabelFont(TimeNewRoman,"x");
  gStyle->SetLabelFont(TimeNewRoman,"y");
  gStyle->SetLabelFont(TimeNewRoman,"z");
  gStyle->SetLabelSize(lhcbTSize,"x");
  gStyle->SetLabelSize(lhcbTSize,"y");
  gStyle->SetLabelSize(lhcbTSize,"z");
  gStyle->SetTitleFont(TimeNewRoman);
  gStyle->SetTitleFont(TimeNewRoman,"x");
  gStyle->SetTitleFont(TimeNewRoman,"y");
  gStyle->SetTitleFont(TimeNewRoman,"z");
  gStyle->SetTitleSize(1.2*lhcbTSize,"x");
  gStyle->SetTitleSize(1.2*lhcbTSize,"y");
  gStyle->SetTitleSize(1.2*lhcbTSize,"z");
   gStyle->SetLineWidth(1) ;
  return;
}
