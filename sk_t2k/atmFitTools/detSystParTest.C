
#include "fqEvent.h"
#include "preProcess.cxx"

std::string mcFilePath = "$RESOURCES_DIR/SK-T2K-Joint/sk4_fcmc_18a_fQv6r0_minituple_100yr_05.root";


int atmpdEventType;

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

map<string, TCanvas*> canvasMap;

TH1D* hist = new TH1D("getRCParameter", "getRCParameter", 100, -200, 100);
TH2D* hist2D = new TH2D("hist2D", ";RC Parameter;;Counts", 100, -200, 100, 15, 0.5, 15.5);
int selectedEventType = -1;

preProcess* preprocess;
fqEvent* fqevent;
TFile* mcFile;
TTree* mcTree;


void fillHist();

void detSystParTest(){

  mcFile = TFile::Open(mcFilePath.c_str());
  mcTree = (TTree*) mcFile->Get("atm_minituple");
  fqevent = new fqEvent(mcTree);
  preprocess = new preProcess();

  mcTree->SetBranchAddress("ATMPDEventType", &atmpdEventType);

  fillHist();

}

void fillHist(){

  // ATMPDEventTypes
  for(int iBinY = 1 ; iBinY <= atmpdLabels.size() ; iBinY++){
    hist2D ->GetYaxis()->SetBinLabel(iBinY, atmpdLabels[iBinY-1].c_str());
  }

  hist->Reset();
  hist2D->Reset();

  cout << "Generating Hist..." << endl;
  int nEvents = mcTree->GetEntries();
  for(int iEvent = 0 ; iEvent < nEvents ; iEvent++){
    GenericToolbox::displayProgressBar(iEvent, nEvents, "READING TREE...");
    fqevent->GetEntry(iEvent);
    if(atmpdEventType == selectedEventType or selectedEventType == -1){
      // cout << GET_VAR_NAME_VALUE(preprocess->getRCParameter(fqevent)) << endl;
      hist->Fill(preprocess->getRCParameter(fqevent));
      hist2D->Fill(preprocess->getRCParameter(fqevent), atmpdEventType);
    }
  }

  hist->Draw();

  canvasMap["hist2D"] = new TCanvas("hist2D", "hist2D", 800, 600);
  hist2D->Draw("COLZ");
  GenericToolbox::fixTH2display(hist2D);
  gPad->SetLogz();
  gPad->SetLeftMargin(0.25);


}
