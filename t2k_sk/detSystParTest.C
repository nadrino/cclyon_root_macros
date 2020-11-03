
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

int maxFqnmrfit = 25;
int maxFqnse = 5;

float fqmrdir[maxFqnmrfit][6][3];
float fq1rnll[maxFqnse][7];
float fqmrmom[maxFqnmrfit][6];

float fqmrnll[maxFqnmrfit];
int fqmrifit[maxFqnmrfit];

int fqnmrfit;

//  struct fqEvent{
//
//   int maxFqnmrfit = 25;
//   int maxFqnse = 5;
//
//   std::vector<std::vector<std::vector<float>>> fqmrdir;
//
//   std::vector<std::vector<float>> fq1rnll;
//   std::vector<std::vector<float>> fqmrmom;
//
//   float* fqmrnll = nullptr;
//   int* fqmrifit = nullptr;
//
//   int fqnmrfit; // array of one
//
//
//   void init(){
//
//     this->fqmrifit = new int[this->maxFqnmrfit];
//     this->fqmrnll = new float[this->maxFqnmrfit];
//
//     fq1rnll.resize(maxFqnse, std::vector<float>(7));
//     fqmrmom.resize(maxFqnmrfit, std::vector<float>(6));
//     fqmrdir.resize( maxFqnmrfit );
//     for(auto& fqmrdirRow : fqmrdir){
//       fqmrdirRow.resize(6);
//       for(auto& fqmrdirCol : fqmrdirRow){
//         fqmrdirCol.resize(3);
//       }
//     }
//
//
//
//
//
//     // this->fqmrmom = ;
//
//     // this->fqmrmom = new float*[this->maxFqnmrfit];
//     // for(int iRow = 0 ; iRow < this->maxFqnmrfit ; iRow++){
//     //   this->fqmrmom[iRow] = new float[6];
//     // }
//     //
//     // this->fq1rnll = new float*[this->maxFqnse];
//     // for(int iRow = 0 ; iRow < this->maxFqnse ; iRow++){
//     //   this->fq1rnll[iRow] = new float[7];
//     // }
//     //
//     // this->fqmrdir = new float**[this->maxFqnmrfit];
//     // for(int iRow = 0 ; iRow < this->maxFqnmrfit ; iRow++){
//     //   this->fqmrdir[iRow] = new float*[6];
//     //   for(int iCol = 0 ; iCol < 6 ; iCol++){
//     //     this->fqmrdir[iRow][iCol] = new float[3];
//     //   }
//     // }
//
//   }
// };

int best2RID;
float fqmrdot;

int getBest2RFitID();
float getRCParameter();

void detSystParTest(){

  TFile* mcFile = TFile::Open(mcFilePath.c_str());
  TTree* mcTree = (TTree*) mcFile->Get("atm_minituple");

  ATMPDEventType eventType;

  // fqEvent fqevent;
  // fqevent.init();

  mcTree->SetBranchAddress("fqmrdir", fqmrdir);

  mcTree->SetBranchAddress("fq1rnll", fq1rnll);
  mcTree->SetBranchAddress("fqmrmom", fqmrmom);

  mcTree->SetBranchAddress("fqmrnll", fqmrnll);
  mcTree->SetBranchAddress("fqmrifit", fqmrifit);

  mcTree->SetBranchAddress("fqnmrfit", &fqevent.fqnmrfit);

  mcTree->SetBranchAddress("ATMPDEventType", &eventType);

  cout << "READING TREE..." << endl;
  int nEvents = mcTree->GetEntries();
  for(int iEvent = 0 ; iEvent < nEvents ; iEvent++){
    GenericToolbox::displayProgressBar(iEvent, nEvents, "READING TREE...");
    mcTree->GetEntry(iEvent);
    if(eventType == SubGeV_elike_0dcy){
      cout << GET_VAR_NAME_VALUE(getRCParameter(&fqevent)) << endl;
    }
  }

}


///////////////////////////////////////////
float getRCParameter(){

  // get best 2R ID
  int ibest = getBest2RFitID();

  // get best 1R Likelihood
  float best1Rnglnl = TMath::Min(fq1rnll[0][1],fq1rnll[0][2] );

  // get mom of 2nd ring
  float ringmom = (float)TMath::Min(fqmrmom[best2RID][0],fqmrmom[best2RID][1]);

  // likelihood difference between 1R and 2R
  float deltaLnL = best1Rnglnl - fqmrnll[best2RID];

  // cut values from fiTQun.cc v4r0
  float a0 = 150.;
  float a1 = -0.6;

  // these values determine the cut line
  float cthresh = a0 + a1*(ringmom);
  if (!(cthresh>0.)) cthresh=0.;


  float rcpar = deltaLnL - cthresh;


  // also calculate angle between 2R rings (should be moved to separate function)
  fqmrdot =   fqmrdir[best2RID][0][0]*fqmrdir[best2RID][1][0]
            + fqmrdir[best2RID][0][1]*fqmrdir[best2RID][1][1]
            + fqmrdir[best2RID][0][2]*fqmrdir[best2RID][1][2];


  // signed sq root ///////////
//  if (rcpar<0.){
//    rcpar =  -1.*TMath::Sqrt(-1.*rcpar);
//  }
//  else{
//    rcpar = TMath::Sqrt(rcpar);
//  }
  ////////////////////////////
  // tricks to make peak widths in distribution more comparable
//  return 500.*(TMath::Log(rcpar+40.)-TMath::Log(40));

  return rcpar;
}

///////////////////////////////////////
//returns the index of the best 2R fit
// ! temporary change to return fit 20000033 !
int getBest2RFitID(){

  // total number of MR fits
  int nfits = (int)fqnmrfit;

  // loop to find the best likelihood
  double ngLnLBest = 10000000.;
  int bestindex = 0;
  for (int ifit=0;ifit<nfits;ifit++){

    int fitID = TMath::Abs(fqmrifit[ifit]); //< fit fit ID code

    // pick out the fits we want to compare to
    if ( TMath::Abs((TMath::Abs(fitID)-20000000))<50){
      // check if it's the best
      if (fqmrnll[ifit] < ngLnLBest){
        ngLnLBest = fqmrnll[ifit];
        bestindex = ifit;
      }
    }
  }
  best2RID = bestindex;
  return bestindex;
}
