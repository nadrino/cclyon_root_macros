
#include "fqEvent.h"
#include "preProcess.cxx"

std::string mcFilePath = "$RESOURCES_DIR/SK-T2K-Joint/sk4_fcmc_18a_fQv6r0_minituple_100yr_05.root";

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

TH1D* hist = new TH1D("getRCParameter", "getRCParameter", 100, -200, 100);
ATMPDEventType selectedEventType = SubGeV_elike_0dcy;

void detSystParTest(){

  TFile* mcFile = TFile::Open(mcFilePath.c_str());
  TTree* mcTree = (TTree*) mcFile->Get("atm_minituple");
  fqEvent* fqevent = new fqEvent(mcTree);
  preProcess* preprocess = new preProcess();

  ATMPDEventType eventType;
  mcTree->SetBranchAddress("ATMPDEventType", &eventType);

  fillHist();

}

void fillHist(){

  cout << "READING TREE..." << endl;
  int nEvents = mcTree->GetEntries();
  for(int iEvent = 0 ; iEvent < nEvents ; iEvent++){
    GenericToolbox::displayProgressBar(iEvent, nEvents, "READING TREE...");
    fqevent->GetEntry(iEvent);
    if(eventType == selectedEventType){
      // cout << GET_VAR_NAME_VALUE(preprocess->getRCParameter(fqevent)) << endl;
      hist->Fill(preprocess->getRCParameter(fqevent));
    }
  }

  hist->Draw();

}
