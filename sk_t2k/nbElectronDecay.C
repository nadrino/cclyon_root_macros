
std::string mcFilePath = "$RESOURCES_DIR/SK-T2K-Joint/sk4_fcmc_18a_fQv6r0_minituple_100yr_05.root";

float fq1rmom[10][7];
float fq1rt0[10][7];
float fq1rnll[10][7];
float fqmrmom[28][6];
int fqmrpid[28][6];
int fqmrnring[28];
int fqnse;
int fqnmrfit;
float genmom;
float MReIncLVal;
float MRnuenuebarLVal;
int atmpdEventType;

TFile* mcFile;
TTree* mcTree;

TH2D* hDecay = new TH2D("hDecay", ";#DeltaT (ns);ATMPDEventType;Counts", 100, 0., 200., 15, 0.5, 15.5);
TH2D* hEvis = new TH2D("hEvis", ";Reconstructed Momentum (MeV);;Counts", 200, 0, 2500, 15, 0.5, 15.5);
TH2D* hNbDecay = new TH2D("hNbDecay", ";# of Decay Electrons;;Counts", 5, -0.5, 4.5, 15, 0.5, 15.5);
TH2D* hmRingPid = new TH2D("hmRingPid", ";MR PID (MER);;Counts", 7, -0.5, 6.5, 15, 0.5, 15.5);
TH2D* h1RingPid = new TH2D("h1RingPid", ";1R PID (e or mu);;Counts", 200, -10000+10, 10010, 15, 0.5, 15.5);
TH2D* hmmeLl = new TH2D("hmmeLl", ";MME Likelihood;;Counts", 200, -12-0.25, 12-0.25, 15, 0.5, 15.5);
TH2D* hnueNuebarSeparation = new TH2D("hnueNuebarSeparation", ";#nu_{e} #bar{#nu}_{e} Separation;;Counts", 160, -5, 15, 15, 0.5, 15.5);

vector<string> atmpdLabels = {
  "SubGeV_elike_0dcy", // 1
  "SubGeV_elike_1dcy",
  "SubGeV_SingleRing_pi0like",
  "SubGeV_mulike_0dcy",
  "SubGeV_mulike_1dcy",
  "SubGeV_mulike_2dcy",
  "SubGeV_pi0like",
  "MultiGeV_elike_nue",
  "MultiGeV_elike_nuebar",
  "MultiGeV_mulike",
  "MultiRing_elike_nue",
  "MultiRing_elike_nuebar",
  "MultiRing_mulike",
  "MultiRingOther_1",
  "PCStop", // not yet in data
  "PCThru",
  "UpStop_mu",
  "UpThruNonShower_mu",
  "UpThruShower_mu"
};

vector<string> pidLabels = {
  "GAMMA", // 0
  "ELECTRON",
  "MUON",
  "PION",
  "KAON",
  "PROTON",
  "CONE GENERATOR"
};

map<string, TCanvas*> canvasMap;

void fixTH2display(TH2 *histogram_);

void nbElectronDecay(){

  mcFile = TFile::Open(mcFilePath.c_str());
  mcTree = (TTree*) mcFile->Get("atm_minituple");

  mcTree->SetBranchAddress("fq1rmom", fq1rmom);
  mcTree->SetBranchAddress("fq1rt0", fq1rt0);
  mcTree->SetBranchAddress("fq1rnll", fq1rnll);
  mcTree->SetBranchAddress("fqmrmom", fqmrmom);
  mcTree->SetBranchAddress("fqmrpid", fqmrpid);
  mcTree->SetBranchAddress("fqmrnring", fqmrnring);

  mcTree->SetBranchAddress("fqnse",    &fqnse);
  mcTree->SetBranchAddress("fqnmrfit", &fqnmrfit);
  mcTree->SetBranchAddress("genmom",   &genmom);
  mcTree->SetBranchAddress("MReIncLVal",   &MReIncLVal);
  mcTree->SetBranchAddress("MRnuenuebarLVal",   &MRnuenuebarLVal);
  mcTree->SetBranchAddress("ATMPDEventType", &atmpdEventType);

  // ATMPDEventTypes
  for(int iBinY = 1 ; iBinY <= atmpdLabels.size() ; iBinY++){
    hNbDecay ->GetYaxis()->SetBinLabel(iBinY, atmpdLabels[iBinY-1].c_str());
    hEvis    ->GetYaxis()->SetBinLabel(iBinY, atmpdLabels[iBinY-1].c_str());
    hmmeLl   ->GetYaxis()->SetBinLabel(iBinY, atmpdLabels[iBinY-1].c_str());
    hmRingPid->GetYaxis()->SetBinLabel(iBinY, atmpdLabels[iBinY-1].c_str());
    h1RingPid->GetYaxis()->SetBinLabel(iBinY, atmpdLabels[iBinY-1].c_str());
    hnueNuebarSeparation->GetYaxis()->SetBinLabel(iBinY, atmpdLabels[iBinY-1].c_str());
  }

  // PIDs
  for(int iBinX = 1 ; iBinX <= pidLabels.size() ; iBinX++){
    hmRingPid->GetXaxis()->SetBinLabel(iBinX, pidLabels[iBinX-1].c_str());
  }

  int nEvents = mcTree->GetEntries();
  float deltaTBuffer;
  int pidIndex;
  int nbDecayElec;
  float mom = 0;
  for(int iEvent = 0 ; iEvent < nEvents ; iEvent++){

//    GenericToolbox::displayProgressBar(iEvent, nEvents, "Reading SK Tree...");
    mcTree->GetEntry(iEvent);

    nbDecayElec = 0;
    pidIndex = 1; // electron
    bool isElectron = fq1rnll[0][1] - fq1rnll[0][3] < 10;
    if(not isElectron) pidIndex = 2; // muon

    for( int iSubPart = 1 ; iSubPart < fqnse ; iSubPart++ ){
      deltaTBuffer = fq1rt0[iSubPart][1] - fq1rt0[0][pidIndex];
      if(deltaTBuffer > 100 and deltaTBuffer < 10000){
        nbDecayElec++;
      }
    }
    hNbDecay->Fill(nbDecayElec, atmpdEventType);

    // hEvis->Fill(fq1rmom[0][pidIndex], atmpdEventType);
    // hEvis->Fill(fq1rmom[0][1], atmpdEventType);
    // hEvis->Fill(fq1rmom[0][3], atmpdEventType);
    // hEvis->Fill(fq1rmom[0][3], atmpdEventType);
    // hEvis->Fill(genmom, atmpdEventType);
    mom = 0;
    for(int iRing = 0 ; iRing < 6 ; iRing++){
      // mom += fqmrmom[0][iRing] - 155.9*(fqmrpid[0][iRing] == 2);
      mom += fqmrmom[0][iRing];
    }
    hEvis->Fill(mom, atmpdEventType);


    hmRingPid->Fill(fqmrpid[0][TMath::LocMax(fqmrnring[0], fqmrmom[0])], atmpdEventType);
    h1RingPid->Fill(fq1rnll[0][1] - fq1rnll[0][3], atmpdEventType);
    hmmeLl->Fill(MReIncLVal, atmpdEventType);
    hnueNuebarSeparation->Fill(MRnuenuebarLVal, atmpdEventType);

  }

   canvasMap["hNbDecay"] = new TCanvas("hNbDecay", "hNbDecay", 800, 600);
   hNbDecay->Draw("COLZ TEXT");
   hNbDecay->GetXaxis()->SetNdivisions(5);
   fixTH2display(hNbDecay);
   gPad->SetLogz();
   gPad->SetLeftMargin(0.25);

  canvasMap["hEvis"] = new TCanvas("hEvis", "hEvis", 800, 600);
  hEvis->Draw("COLZ");
  fixTH2display(hEvis);
  gPad->SetLogz();
  gPad->SetLeftMargin(0.25);

   canvasMap["hmmeLl"] = new TCanvas("hmmeLl", "hmmeLl", 800, 600);
   hmmeLl->Draw("COLZ");
   fixTH2display(hmmeLl);
   gPad->SetLogz();
   gPad->SetLeftMargin(0.25);

   canvasMap["hnueNuebarSeparation"] = new TCanvas("hnueNuebarSeparation", "hnueNuebarSeparation", 800, 600);
   hnueNuebarSeparation->Draw("COLZ");
   fixTH2display(hnueNuebarSeparation);
   gPad->SetLogz();
   gPad->SetLeftMargin(0.25);

   canvasMap["h1RingPid"] = new TCanvas("h1RingPid", "h1RingPid", 800, 600);
   h1RingPid->Draw("COLZ");
   h1RingPid->GetXaxis()->SetMaxDigits(3);
   fixTH2display(h1RingPid);
   gPad->SetLogz();
   gPad->SetLeftMargin(0.25);

   canvasMap["hmRingPid"] = new TCanvas("hmRingPid", "hmRingPid", 800, 600);
   hmRingPid->Draw("COLZ TEXT");
   hmRingPid->GetXaxis()->SetNdivisions(5);
   fixTH2display(hmRingPid);
   gPad->SetLogz();
   gPad->SetLeftMargin(0.25);

}

void fixTH2display(TH2 *histogram_){

  gPad->SetRightMargin(0.15);
  histogram_->GetZaxis()->SetTitleOffset(0.8);
  auto* pal = (TPaletteAxis*) histogram_->GetListOfFunctions()->FindObject("palette");
  // TPaletteAxis* pal = (TPaletteAxis*) histogram_->GetListOfFunctions()->At(0);
  if(pal != nullptr){
    pal->SetX1NDC(1 - 0.15 + 0.01);
    pal->SetX2NDC(1 - 0.15 + 0.05);
    pal->GetAxis()->SetMaxDigits(2);
    pal->Draw();
  }

}
