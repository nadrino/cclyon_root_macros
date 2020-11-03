


enum ATMPDEventType{
  SubGeV_elike_0dcy = 1,
  SubGeV_elike_1dcy,
  SubGeV_SingleRing_pi0like,
  SubGeV_mulike_0dcy,
  SubGeV_mulike_1dcy,
  SubGeV_mulike_2dcy,
  SubGeV_pi0like,
  MultiGeV_elike_nue,
  MultiGeV_elike_nuebar,
  MultiGeV_mulike,
  MultiRing_elike_nue,
  MultiRing_elike_nuebar,
  MultiRing_mulike,
  MultiRingOther_1,
  PCStop, // not yet in data
  PCThru,
  UpStop_mu,
  UpThruNonShower_mu,
  UpThruShower_mu
};

struct fqEvent{
  float*** fqmrdir;

  float** fq1rnll;
  float** fqmrmom;

  float* fqmrnll;
}

std::vector<double> AtmSubGeVRecoEdges;
std::vector<double> CosZenithRecoEdges;

std::vector< std::vector<std::string> > cutObservableListForSample;
std::vector< std::string > cutsList;

void init();


void detSystReweight(){

  init();



}


void init(){

  // CUTS DEFINITION
  // https://www.t2k.org/t2ksk/code/ntuplevariables
  cutsList.emplace_back("evis"); // Reconstructed Visible Energy
  cutsList.emplace_back("fq1rnll"); // LL for 1 ring
  cutsList.emplace_back("fqmrnll"); // LL for multi rings
  cutsList.emplace_back("fqpi0nll"); // LL for Pi0 id

  // fqnse is the number of sub events, with the first being the lepton, and subsequent decay electrons
  cutsList.emplace_back("fqnse - 1"); // Nb of decay electrons

  cutsList.emplace_back("PID");
  cutsList.emplace_back("Ndecay");
  cutsList.emplace_back("PID of MER");
  cutsList.emplace_back("MME LL");
  cutsList.emplace_back("NueNuebarSeparationLL");
  cutsList.emplace_back("Both rings are e-like and Minv in π0 range");

  // SAMPLES DEFINITION
  int nbSamples = UpThruShower_mu; // last one
  cutObservableListForSample.resize(nbSamples);

  // // SubGeV Samples
  // cutObservableListForSample[ SubGeV_elike_0dcy ].emplace_back("Evis");
  // cutObservableListForSample[ SubGeV_elike_0dcy ].emplace_back("Nrings");
  // cutObservableListForSample[ SubGeV_elike_0dcy ].emplace_back("PID");
  // cutObservableListForSample[ SubGeV_elike_0dcy ].emplace_back("Ndecay");
  //
  // cutObservableListForSample[ SubGeV_elike_1dcy ] = cutObservableListForSample[ SubGeV_elike_0dcy ];
  //
  // cutObservableListForSample[ SubGeV_mulike_0dcy ] = cutObservableListForSample[ SubGeV_elike_0dcy ];
  // cutObservableListForSample[ SubGeV_mulike_1dcy ] = cutObservableListForSample[ SubGeV_elike_0dcy ];
  // cutObservableListForSample[ SubGeV_mulike_2dcy ] = cutObservableListForSample[ SubGeV_elike_0dcy ];
  //
  // cutObservableListForSample[ MultiGeV_elike_nue ] = cutObservableListForSample[ SubGeV_elike_0dcy ];
  // cutObservableListForSample[ MultiGeV_elike_nuebar ] = cutObservableListForSample[ SubGeV_elike_0dcy ];
  //
  // cutObservableListForSample[ SubGeV_pi0like ] = cutObservableListForSample[ SubGeV_elike_0dcy ];
  // cutObservableListForSample[ SubGeV_pi0like ].pop_back();
  // cutObservableListForSample[ SubGeV_pi0like ].emplace_back("Both rings are e-like and Minv in π0 range");
  //
  // // MultiGeV / SR samples
  // cutObservableListForSample[ MultiGeV_mulike ].emplace_back("Evis");
  // cutObservableListForSample[ MultiGeV_mulike ].emplace_back("Nrings");
  // cutObservableListForSample[ MultiGeV_mulike ].emplace_back("PID");
  //
  // // MultiGeV / MR samples
  // cutObservableListForSample[ MultiRing_mulike ].emplace_back("Evis");
  // cutObservableListForSample[ MultiRing_mulike ].emplace_back("Nrings");
  // cutObservableListForSample[ MultiRing_mulike ].emplace_back("PID of MER");
  //
  // cutObservableListForSample[ MultiRingOther ] = cutObservableListForSample[ MultiRing_mulike ];
  // cutObservableListForSample[ MultiRingOther ].emplace_back("MME LL");
  //
  // cutObservableListForSample[ MultiRing_elike_nue ] = cutObservableListForSample[ MultiRingOther ];
  // cutObservableListForSample[ MultiRing_elike_nue ].emplace_back("NueNuebarSeparationLL");
  //
  // cutObservableListForSample[ MultiRing_elike_nuebar ] = cutObservableListForSample[ MultiRing_elike_nue ];


  // MultiGeV > 1.33 GeV
  cutObservableListForSample[ MultiGeV_elike_nue ].emplace_back("evis");
  cutObservableListForSample[ MultiGeV_elike_nue ].emplace_back("nrings"); // Single ring
  // elike
  cutObservableListForSample[ MultiGeV_elike_nue ].emplace_back("pid");
  // 1 or more decay electron
  // LLrat_nDcy_nueneb: number of decay electrons for nue/neubar separation
  cutObservableListForSample[ MultiGeV_elike_nue ].emplace_back("LLrat_nDcy_nueneb");


  // OSCILLATION FIT BINNING
  // binning edges from Fig. 7.2 in m.jiang's thesis
  const double tenTo2p0 =    100.0;
  const double tenTo2p4 =    251.2;
  const double tenTo2p6 =    398.1;
  const double tenTo2p8 =    631.0;
  const double tenTo3p0 =   1000.;
  const double tenTo3p2 =   1585.;
  const double tenTo3p4 =   2512.;
  const double tenTo3p7 =   5012.;
  const double tenTo4p0 =  10000.;
  const double tenTo5p0 = 100000.;

  AtmSubGeVRecoEdges = { // in MeV
    // 5 bins
    tenTo2p0, tenTo2p4, tenTo2p6, tenTo2p8, tenTo3p0, tenTo3p2
  };

  // 10 bins
  CosZenithRecoEdges = {
    // 0.20 steps all the way
    -1.00, -0.80, -0.60, -0.40, -0.20,
    0.00,  0.20,  0.40,  0.60,  0.80,  1.00
  };


}
