Int_t Btag(Int_t pid){
  //==============================================================================
  // select beauty hadrons
  //==============================================================================
  const  int N = 51;
  static int list[N] = {
    511,     // B0
    513,     // B*0
    521,     // B+-
    523,     // B*+-
    531,     // B_s0
    533,     // B_s*0
    541,    // B_c+-
    543,    // B_c*+-
    551,    // eta_b
    553,
    557,
    555,    // Upsilon
    515,   // B_2*0
    525,   // B_2+-
    535,   //  B_s2*0
    545,   // B_c2*+-
    5122,   // Lambda_b
    5112,   // Sigma_b-
    5212,   // Sigma_b0
    5222, // Sigma_b+
    5114,  // Sigma_b- *
    5214,  // Sigma_b0 *
    5224,  // Sigma_b+ *
    5132, // Xi_b particles
    5232,
    5312,
    5322,
    5314,
    5324,
    5332, // Omega_b particles
    5334,
    5142,
    5242,
    5412,
    5422,
    5414,
    5424,
    5324,
    5432,
    5434,
    5442,
    5444,
    5512,
    5522,
    5514,
    5524,
    5532,
    5534,
    5542,
    5544,
    5554};
  for(int i=0; i<N; ++i) {
    if(abs(pid)==list[i]) return 1;
  }
  // if(abs(pid) == 5) return 1;
  return 0;


}

Double_t GetDeltaR(Double_t eta1, Double_t phi1 , Double_t eta2, Double_t phi2) {
  //==============================================================================
  // Obtain  the jet seaperation
  //==============================================================================
  Double_t r = std::sqrt(  std::pow((eta1-eta2),2)
                           +std::pow((phi1-phi2),2));
  return r;
}


Double_t GetpT( Double_t px, Double_t py  ) {
  //==============================================================================
  // obtain invariant mass of particles given their 4 momenta
  //==============================================================================

Double_t pT = std::sqrt(  std::pow(px,2)
                         +std::pow(py,2));
return pT;
}




Double_t GetInvMass( Double_t px1, Double_t py1, Double_t pz1 , Double_t e1 ,
                    Double_t px2, Double_t py2, Double_t pz2 , Double_t e2  ) {
  //==============================================================================
  // obtain invariant mass of particles given their 4 momenta
  //==============================================================================

Double_t m = std::sqrt(  -std::pow(px1+px2,2)
                         -std::pow(py1+py2,2)
                         -std::pow(pz1+pz2,2)
                         +std::pow(e1+e2,2) );
return m;
}


int checkFinal(int pid)
//==============================================================================
// check whether a particle is counted as final state particle
//==============================================================================
{
  // list of particles with tau>0.03 ns (PIDcode, name, lifetime/ns)

  const  int N = 18;
  static int list[N] = {
    11,     // e-        inf
    12,     // nu_e      inf
    13,     // mu        2170
    14,     // nu_mu     inf
    16,     // nu_tau    inf
    22,     // gamma     inf
    130,    // K0L       51.2
    211,    // pi+       26
    310,    // K0S       0.089
    321,    // K+        12.4
    2112,   // n         880e9
    2212,   // p         inf
    3112,   // Sigma-    0.15
    3122,   // Lambda    0.26
    3222,   // Sigma+    0.08
    3312,   // Xi-       0.16
    3322,   // Xi0       0.29
    3334};  // Omega-    0.08

  for(int i=0; i<N; ++i) {
    if(abs(pid)==list[i]) return 1;
  }
  return 0;
}
