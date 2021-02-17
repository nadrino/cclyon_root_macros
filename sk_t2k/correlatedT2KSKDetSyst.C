//
// Created by Adrien BLANCHET on 14/12/2020.
//

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
  UpThruShower_mu // 19
};

std::string outputPath = "./throws.root";

std::string filePath = "/Users/ablanche/Documents/Work/Repositories/cclyon_root_macros/sk_t2k/skDetCovMatrix.root";
std::string matrixName = "covMatrix";

std::string skFilePath = "/Users/ablanche/Documents/Work/Resources/SK-T2K-Joint/sk4_fcmc_18a_fQv6r0_minituple_100yr_05.root";
std::vector<std::string> atmMatrixBinCutsList;
std::vector<TTreeFormula*> atmFormulaList;
TFile* skFile;
TTree* skTree;

std::string skToT2kSamplePath = "/Users/ablanche/Documents/Work/Resources/SK-T2K-Joint/T2KforATMcuts.root";
std::vector<int> t2kSampleIntList = { 11, 13, 14 };
std::vector<std::string> t2kSampleNameList = { "1Re", "1Rmu", "1ReDe" };
std::vector<std::pair<int, std::string>> t2kSubSamplesList; // pair is <t2kSample, CutString>
std::vector<TTreeFormula*> t2kSubSampleFormulaList;
std::map<int,int> entryT2kSampleMapping;
std::map<int,int> entryT2kSubSampleMapping;
std::vector<std::vector<int>> entrySkSampleMapping; // entrySkSampleMapping[iEntry] = {skSample0, skSample4, ...}
std::vector<std::vector<double>> t2kSampleSkEntries;    // Count the number of events for each [T2KSample][SkSample]
std::vector<std::vector<double>> t2kSubSampleSkEntries; // Count the number of events for each [T2KSubSample][SkSample]

// output
TFile* throwsFile;
TTree* throwsTree;
std::vector<std::string> binTitlesList;

int nbThrows = 1E4;

TFile* covFile;
TMatrixD* covMatrix;
TMatrixDSym* covMatrixReg;
TMatrixD* fCholeskyCovMatrix;

TRandom3* fGenerator;
int fNSystParameters;

void init();
void mapSkEventSamples();
void fillThrowTTree();
void generatePlots();
std::vector<double> throwParameters();
TMatrixD* getCovarianceMatrixOfTree(TTree* tree_);

// MAIN
void correlatedT2KSKDetSyst(){

  init();
  mapSkEventSamples();

  // Init output
  throwsFile = TFile::Open(outputPath.c_str(), "RECREATE");
  throwsTree = new TTree("throws", "throws");

  fillThrowTTree();
  generatePlots();

  throwsFile->Close();
  std::cout << "OUTPUT FILE: " << outputPath << std::endl;

  exit(0);

}


void init(){

  cout << __METHOD_NAME__ << endl;

  covFile = TFile::Open(filePath.c_str());
  covMatrix = (TMatrixD*) covFile->Get(matrixName.c_str());

  fNSystParameters = covMatrix->GetNrows();

  fGenerator = new TRandom3();
  fGenerator->SetSeed(0);

  covMatrixReg = GenericToolbox::convertToSymmetricMatrix(covMatrix);

  //Cholesky decomposition of the covariance matrix
  TDecompChol *choleskyDecomposer = new TDecompChol((*covMatrixReg));
  bool decompSucceeded = choleskyDecomposer->Decompose();

  cout << "Computing Cholesky decomposition of the covariance matrix" << endl;

  int trial = 0;
  if(not decompSucceeded){
    cout << "Adding to the diagonal terms to make matrix positive definite..." << endl;

    while(not decompSucceeded and trial < 50){
      for (int i = 0; i < covMatrixReg->GetNcols(); i++)
      {
        (*covMatrixReg)[i][i] += 1E-19; // original: 0.000000000001 // 1E-12
      }
      choleskyDecomposer->SetMatrix((*covMatrixReg));
      decompSucceeded = choleskyDecomposer->Decompose();
      trial++;
    }

    if(decompSucceeded){
      cout << "Cholesky decomposition succeeded after " << trial+1 << " trials." << endl;
    }
    else{
      cout << "Cholesky decomposition failed." << endl;
      exit(1);
    }

  }

  cout << "choleskyDecomposer->GetU() = " << (((TMatrixD *)(choleskyDecomposer->GetU()).Clone())->T()).GetNcols() << endl;
  fCholeskyCovMatrix = (TMatrixD *)(((TMatrixD *)(choleskyDecomposer->GetU()).Clone())->T()).Clone();
  cout << "fCholeskyCovMatrix->GetNcols() = " << fCholeskyCovMatrix->GetNcols() << endl;

  delete choleskyDecomposer;

  TList* cutStrings_TList = (TList*) covFile->Get("resources/cutStrings_TList");
  for(int iComp = 0 ; iComp < cutStrings_TList->GetEntries() ; iComp++){
    atmMatrixBinCutsList.emplace_back(cutStrings_TList->At(iComp)->GetTitle());
  }


  cout << "Listing cuts for atmospheric samples" << endl;

  // Formula
  skFile = TFile::Open(skFilePath.c_str());
  skTree = (TTree*) skFile->Get("atm_minituple");

  for(int iPar = 0 ; iPar < int(atmMatrixBinCutsList.size()) ; iPar++){
    atmFormulaList.emplace_back(
      new TTreeFormula(
        Form("skAtmSample_%i", int(iPar)),
        atmMatrixBinCutsList[iPar].c_str(),
        skTree
      )
    );
  }

  //////////////////////////////////
  //> T2K selection sub-samples


  ///////////////
  // FHC 1Re: Osc-nue
  ///////////////
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (0, "fabs(mode)<=30 && ipnu[0] == 12 && genmom >= 0   && genmom < 300"));
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (0, "fabs(mode)<=30 && ipnu[0] == 12 && genmom >= 300 && genmom < 700"));
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (0, "fabs(mode)<=30 && ipnu[0] == 12 && genmom >= 700"));

  // FHC 1Re: numu
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (0, "fabs(mode)<=30 && ipnu[0] == 14 && genmom >= 0   && genmom < 300"));
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (0, "fabs(mode)<=30 && ipnu[0] == 14 && genmom >= 300 && genmom < 700"));
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (0, "fabs(mode)<=30 && ipnu[0] == 14 && genmom >= 700"));

  // FHC 1Re: beam-nue
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (0, "fabs(mode)<=30 && ipnu[0] == 12 && genmom >= 0   && genmom < 300"));
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (0, "fabs(mode)<=30 && ipnu[0] == 12 && genmom >= 300 && genmom < 700"));
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (0, "fabs(mode)<=30 && ipnu[0] == 12 && genmom >= 700"));

  // FHC 1Re: NC 
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (0, "fabs(mode)>30 && ipnu[0] > 0 && genmom >= 0   && genmom < 300"));
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (0, "fabs(mode)>30 && ipnu[0] > 0 && genmom >= 300 && genmom < 700"));
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (0, "fabs(mode)>30 && ipnu[0] > 0 && genmom >= 700"));


  ///////////////
  // FHC 1Rmu: numu CCQE
  ///////////////
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (1, "fabs(mode)==1 && ipnu[0] == 14 && pnu[2]*1000. >= 0   && pnu[2]*1000. < 400"));
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (1, "fabs(mode)==1 && ipnu[0] == 14 && pnu[2]*1000. >= 400 && pnu[2]*1000. < 1100"));
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (1, "fabs(mode)==1 && ipnu[0] == 14 && pnu[2]*1000. >= 1100 && pnu[2]*1000. < 30000"));

  // FHC 1Rmu: numu non-CCQE
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (1, "fabs(mode)!=1 && fabs(mode)<=30 && ipnu[0] == 14 && pnu[2]*1000. >= 0 && pnu[2]*1000. < 30000"));

  // FHC 1Rmu: nue
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (1, "fabs(mode)<=30 && ipnu[0] == 12 && pnu[2]*1000. >= 0 && pnu[2]*1000. < 30000"));

  // FHC 1Rmu: NC 
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (1, "fabs(mode)>30 && ipnu[0] > 0 && pnu[2]*1000. >= 0 && pnu[2]*1000. < 30000"));


  ///////////////
  // RHC 1Re: Osc-nue
  ///////////////
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (0, "fabs(mode)<=30 && ipnu[0] == -12 && genmom >= 0   && genmom < 300"));
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (0, "fabs(mode)<=30 && ipnu[0] == -12 && genmom >= 300 && genmom < 700"));
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (0, "fabs(mode)<=30 && ipnu[0] == -12 && genmom >= 700"));

  // RHC 1Re: numu
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (0, "fabs(mode)<=30 && ipnu[0] == -14 && genmom >= 0   && genmom < 300"));
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (0, "fabs(mode)<=30 && ipnu[0] == -14 && genmom >= 300 && genmom < 700"));
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (0, "fabs(mode)<=30 && ipnu[0] == -14 && genmom >= 700"));

  // RHC 1Re: beam-nue
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (0, "fabs(mode)<=30 && ipnu[0] == -12 && genmom >= 0   && genmom < 300"));
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (0, "fabs(mode)<=30 && ipnu[0] == -12 && genmom >= 300 && genmom < 700"));
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (0, "fabs(mode)<=30 && ipnu[0] == -12 && genmom >= 700"));

  // RHC 1Re: NC 
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (0, "fabs(mode)>30 && ipnu[0] < 0 && genmom >= 0   && genmom < 300"));
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (0, "fabs(mode)>30 && ipnu[0] < 0 && genmom >= 300 && genmom < 700"));
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (0, "fabs(mode)>30 && ipnu[0] < 0 && genmom >= 700"));


  ///////////////
  // RHC 1Rmu: numu CCQE
  ///////////////
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (1, "fabs(mode)==1 && ipnu[0] == -14 && pnu[2]*1000. >= 0   && pnu[2]*1000. < 400"));
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (1, "fabs(mode)==1 && ipnu[0] == -14 && pnu[2]*1000. >= 400 && pnu[2]*1000. < 1100"));
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (1, "fabs(mode)==1 && ipnu[0] == -14 && pnu[2]*1000. >= 1100 && pnu[2]*1000. < 30000"));

  // RHC 1Rmu: numu non-CCQE
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (1, "fabs(mode)!=1 && fabs(mode)<=30 && ipnu[0] == -14 && pnu[2]*1000. >= 0 && pnu[2]*1000. < 30000"));

  // RHC 1Rmu: nue
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (1, "fabs(mode)<=30 && ipnu[0] == -12 && pnu[2]*1000. >= 0 && pnu[2]*1000. < 30000"));

  // RHC 1Rmu: NC 
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (1, "fabs(mode)>30 && ipnu[0] < 0 && pnu[2]*1000 >= 0 && pnu[2]*1000. < 30000"));


  ///////////////
  // FHC 1Re1De: Osc-nue
  ///////////////
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (2, "fabs(mode)<=30 && ipnu[0] == 12 && genmom >= 0   && genmom < 300"));
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (2, "fabs(mode)<=30 && ipnu[0] == 12 && genmom >= 300 && genmom < 700"));
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (2, "fabs(mode)<=30 && ipnu[0] == 12 && genmom >= 700"));

  // FHC 1Re1De: numu
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (2, "fabs(mode)<=30 && ipnu[0] == 14 && genmom >= 0   && genmom < 300"));
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (2, "fabs(mode)<=30 && ipnu[0] == 14 && genmom >= 300 && genmom < 700"));
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (2, "fabs(mode)<=30 && ipnu[0] == 14 && genmom >= 700"));

  // FHC 1Re1De: beam-nue
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (2, "fabs(mode)<=30 && ipnu[0] == 12 && genmom >= 0   && genmom < 300"));
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (2, "fabs(mode)<=30 && ipnu[0] == 12 && genmom >= 300 && genmom < 700"));
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (2, "fabs(mode)<=30 && ipnu[0] == 12 && genmom >= 700"));

  // FHC 1Re1De: NC 
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (2, "fabs(mode)>30 && ipnu[0] > 0 && genmom >= 0   && genmom < 300"));
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (2, "fabs(mode)>30 && ipnu[0] > 0 && genmom >= 300 && genmom < 700"));
  t2kSubSamplesList.emplace_back(std::pair<int,std::string>
    (2, "fabs(mode)>30 && ipnu[0] > 0 && genmom >= 700"));

  for(size_t iSubSample = 0 ; iSubSample < t2kSubSamplesList.size() ; iSubSample++){
    t2kSubSampleFormulaList.emplace_back(
      new TTreeFormula(
        Form("t2kSubSample_%i", int(iSubSample)),
        t2kSubSamplesList[iSubSample].second.c_str(),
        skTree
      )
    );
  }


  //////////////////////////
  // Lucile files
  cout << "Mapping T2K selection on Atm files..." << endl;
  TFile* atmToT2kFile = TFile::Open(skToT2kSamplePath.c_str());
  for(int iT2kSample = 0 ; iT2kSample < t2kSampleIntList.size() ; iT2kSample++){
    TTree* t2kSampleTree = (TTree*) atmToT2kFile->Get(Form("fOutputTree%i", t2kSampleIntList[iT2kSample]));
    for(int iEntry = 0 ; iEntry < t2kSampleTree->GetEntries() ; iEntry++){
      t2kSampleTree->GetEntry(iEntry);
      entryT2kSampleMapping[t2kSampleTree->GetLeaf("SKEntry_d")->GetValue(0)] = iT2kSample;
    }
  }

}
void mapSkEventSamples(){

  cout << "Mapping SK Atmospheric Event Samples (SK + T2K cuts)" << std::endl;

  entrySkSampleMapping.resize(skTree->GetEntries());

  t2kSampleSkEntries.resize(t2kSampleIntList.size());
  for(auto& t2kSample : t2kSampleSkEntries){
    t2kSample.resize(atmMatrixBinCutsList.size());
  }
  t2kSubSampleSkEntries.resize(t2kSubSamplesList.size());
  for(auto& t2kSubSample : t2kSubSampleSkEntries){
    t2kSubSample.resize(atmMatrixBinCutsList.size());
  }

  for(int iEntry = 0 ; iEntry < skTree->GetEntries() ; iEntry++){

    GenericToolbox::displayProgressBar(iEntry, skTree->GetEntries(), "reading skTree...");

    skTree->GetEntry(iEntry);

    // iEntry is in which Atm samples?
    entrySkSampleMapping[iEntry] = {};
    for(int iSkSample = 0 ; iSkSample < int(atmMatrixBinCutsList.size()) ; iSkSample++){

      // Does the event passes the cut of this sample
      bool doEventPassCut = false;
      skTree->SetNotify(atmFormulaList[iSkSample]);
      for(int jInstance = 0; jInstance < atmFormulaList[iSkSample]->GetNdata(); jInstance++) {
        if ( atmFormulaList[iSkSample]->EvalInstance(jInstance) ) {
          doEventPassCut = true;
          break;
        }
      }

      if(doEventPassCut) {
        // Attaching this sample to iEntry:
        entrySkSampleMapping[iEntry].emplace_back(iSkSample);
      }

    }

    // iEntry is in which T2K samples?
    if(not entrySkSampleMapping[iEntry].empty()){

      for(const auto& entrySkSample : entrySkSampleMapping[iEntry]){
        t2kSampleSkEntries[entryT2kSampleMapping[iEntry]][entrySkSample]++; // [iT2kSample][iSkSample]
      }

      // T2K Sub-Samples
      for(size_t iSubSample = 0 ; iSubSample < t2kSubSamplesList.size() ; iSubSample++){

        if(entryT2kSampleMapping[iEntry] != t2kSubSamplesList[iSubSample].first){
          continue;
        }

        // Does the event passes the cut of this sample
        bool doEventPassCut = false;
        skTree->SetNotify(t2kSubSampleFormulaList[iSubSample]);
        for(int jInstance = 0; jInstance < t2kSubSampleFormulaList[iSubSample]->GetNdata(); jInstance++) {
          if ( t2kSubSampleFormulaList[iSubSample]->EvalInstance(jInstance) ) {
            doEventPassCut = true;
            break;
          }
        }

        if(doEventPassCut){
          for(const auto& entrySkSample : entrySkSampleMapping[iEntry]){
            t2kSubSampleSkEntries[iSubSample][entrySkSample]++; // [iT2kSample][iSkSample]
          }
        }

      }
    }

  }



}
void fillThrowTTree(){

  cout << "Filling up TTree with throws (" << nbThrows << ")" << endl;

  std::vector<double> skThrownParameters(covMatrix->GetNrows()); // reweighted norm SK

  std::vector<double> t2kPropagatedNorm(t2kSampleIntList.size()); // reweighted norm T2K
  std::vector<int> t2kPropagatedNormCount(t2kSampleIntList.size()); // nominal norm T2K

  std::vector<double> t2kSubPropagatedNorm(t2kSubSamplesList.size()); // reweighted norm T2K Sub sample
  std::vector<int> t2kSubPropagatedNormCount(t2kSubSamplesList.size()); // nominal norm T2K Sub sample

  // Hooking tree
  for(int iAtmSample = 0 ; iAtmSample < int(skThrownParameters.size()) ; iAtmSample++){
    throwsTree->Branch(Form("Atm_%i", iAtmSample), &skThrownParameters[iAtmSample]);
    binTitlesList.emplace_back(atmMatrixBinCutsList[iAtmSample]);
  }
  for(int iT2kSample = 0 ; iT2kSample < int(t2kSampleNameList.size()) ; iT2kSample++){
    throwsTree->Branch(Form("T2K_%s", t2kSampleNameList[iT2kSample].c_str()), &t2kPropagatedNorm[iT2kSample]);
    binTitlesList.emplace_back(t2kSampleNameList[iT2kSample]);
  }
  for(int iT2kSubSample = 0 ; iT2kSubSample < int(t2kSubSamplesList.size()) ; iT2kSubSample++){
    throwsTree->Branch(Form("T2K_Sub_%i", iT2kSubSample), &t2kSubPropagatedNorm[iT2kSubSample]);
    binTitlesList.emplace_back(Form("%s && %s",
                                    t2kSampleNameList[t2kSubSamplesList[iT2kSubSample].first].c_str(),
                                    t2kSubSamplesList[iT2kSubSample].second.c_str()
                                    )
                              );
  }

  // Make throws
  cout << "Throwing toys..." << endl;
  for(int iThrow = 0 ; iThrow < nbThrows ; iThrow++){

    GenericToolbox::displayProgressBar(iThrow, nbThrows, "throwing...");

    auto throws = throwParameters();

    // Propagating on Atm Samples
    for(int iAtmSample = 0 ; iAtmSample < int(throws.size()) ; iAtmSample++){
      skThrownParameters[iAtmSample] = throws[iAtmSample];
    }

    // Propagating on T2k Samples
    for(int iT2kSample = 0 ; iT2kSample < int(t2kSampleSkEntries.size()) ; iT2kSample++){
      t2kPropagatedNorm[iT2kSample] = 0; // reset
      double norm = 0;
      for(int iAtmSample = 0 ; iAtmSample < int(throws.size()) ; iAtmSample++){
        t2kPropagatedNorm[iT2kSample] += t2kSampleSkEntries[iT2kSample][iAtmSample]*skThrownParameters[iAtmSample];
        norm += t2kSampleSkEntries[iT2kSample][iAtmSample];
      } // iAtmSample
      t2kPropagatedNorm[iT2kSample] /= norm;
    } // iT2kSample

    // Propagating on T2k Sub-Samples
    for(int iT2kSubSample = 0 ; iT2kSubSample < int(t2kSubSampleSkEntries.size()) ; iT2kSubSample++){
      t2kSubPropagatedNorm[iT2kSubSample] = 0; // reset
      double norm = 0;
      for(int iAtmSample = 0 ; iAtmSample < int(throws.size()) ; iAtmSample++){
        t2kSubPropagatedNorm[iT2kSubSample] += t2kSubSampleSkEntries[iT2kSubSample][iAtmSample]*skThrownParameters[iAtmSample];
        norm += t2kSubSampleSkEntries[iT2kSubSample][iAtmSample];
      } // iAtmSample
      t2kSubPropagatedNorm[iT2kSubSample] /= norm;
    } // iT2kSubSample

    throwsTree->Fill();

  }

  throwsTree->Write("throws");
  covMatrix->Write("skCovMatrix");

}
void generatePlots(){

  cout << "Generating plots" << endl;

  TMatrixD* matrixTemp;
  TH2D* histTemp;

  TMatrixD* covSkT2kMatrix_TMatrixD = getCovarianceMatrixOfTree(throwsTree);
  covSkT2kMatrix_TMatrixD->Write("covSkT2kMatrix_TMatrixD");
  histTemp = GenericToolbox::convertTMatrixDtoTH2D(covSkT2kMatrix_TMatrixD, "covSkT2kMatrix");
  for(int iLeaf = 0 ; iLeaf < histTemp->GetNbinsY() ; iLeaf++){
//    histTemp->GetYaxis()->SetBinLabel(iLeaf+1, throwsTree->GetListOfLeaves()->At(iLeaf)->GetName());
    histTemp->GetYaxis()->SetBinLabel(iLeaf+1, binTitlesList[iLeaf].c_str());
  }
  histTemp->GetYaxis()->SetLabelSize(histTemp->GetYaxis()->GetLabelSize()/3.);
  histTemp->GetYaxis()->SetTitle("");
  histTemp->Write("covSkT2kMatrix_TH2D");

  // Sigmas List
  TVectorD* sigmaSkT2kList_TVectorD = new TVectorD(covSkT2kMatrix_TMatrixD->GetNcols());
  for(int iComp = 0 ; iComp < covSkT2kMatrix_TMatrixD->GetNcols() ; iComp++ ){
    (*sigmaSkT2kList_TVectorD)[iComp] = TMath::Sqrt((*covSkT2kMatrix_TMatrixD)[iComp][iComp]);
  }
  sigmaSkT2kList_TVectorD->Write("sigmaSkT2kList_TVectorD");
  TH1D* sigmaSkT2kList_TH1D = GenericToolbox::convertTVectorDtoTH1D(sigmaSkT2kList_TVectorD, "sigmaSkT2kList_TH1D");
  for(int iLeaf = 0 ; iLeaf < sigmaSkT2kList_TH1D->GetNbinsX() ; iLeaf++){
//    sigmaSkT2kList_TH1D->GetXaxis()->SetBinLabel(iLeaf+1, throwsTree->GetListOfLeaves()->At(iLeaf)->GetName());
    sigmaSkT2kList_TH1D->GetXaxis()->SetBinLabel(iLeaf+1, binTitlesList[iLeaf].c_str());
  }
  sigmaSkT2kList_TH1D->GetXaxis()->SetLabelSize(sigmaSkT2kList_TH1D->GetXaxis()->GetLabelSize()/3.);
  sigmaSkT2kList_TH1D->Write("sigmaSkT2kList_TH1D");

  // Correlation Matrix
  matrixTemp = GenericToolbox::convertToCorrelationMatrix(covSkT2kMatrix_TMatrixD);
  matrixTemp->Write("corSkT2kMatrix_TMatrixD");

  histTemp = GenericToolbox::convertTMatrixDtoTH2D(matrixTemp, "corMatrix");
  for(int iLeaf = 0 ; iLeaf < histTemp->GetNbinsY() ; iLeaf++){
//    histTemp->GetYaxis()->SetBinLabel(iLeaf+1, throwsTree->GetListOfLeaves()->At(iLeaf)->GetName());
    histTemp->GetYaxis()->SetBinLabel(iLeaf+1, binTitlesList[iLeaf].c_str());
  }
  histTemp->GetYaxis()->SetLabelSize(histTemp->GetYaxis()->GetLabelSize()/3.);
  histTemp->GetYaxis()->SetTitle("");
  histTemp->Write("corSkT2kMatrix_TH2D");

  matrixTemp = new TMatrixD(atmMatrixBinCutsList.size(), atmMatrixBinCutsList.size()); // only atm part
  for(int iCol = 0 ; iCol < matrixTemp->GetNcols() ; iCol++){
    for(int iRow = 0 ; iRow < matrixTemp->GetNrows() ; iRow++){
      (*matrixTemp)[iCol][iRow] =
        ( (*covSkT2kMatrix_TMatrixD)[iCol][iRow] - (*covMatrix)[iCol][iRow] )
        / ( TMath::Sqrt((*covMatrix)[iCol][iCol])*TMath::Sqrt((*covMatrix)[iRow][iRow]) );
      (*matrixTemp)[iCol][iRow] = TMath::Abs((*matrixTemp)[iCol][iRow])*100;
    } // iRow
  } // iCol

  histTemp = GenericToolbox::convertTMatrixDtoTH2D(matrixTemp, "covDeltaAtm");
  for(int iLeaf = 0 ; iLeaf < histTemp->GetNbinsY() ; iLeaf++){
//    histTemp->GetYaxis()->SetBinLabel(iLeaf+1, throwsTree->GetListOfLeaves()->At(iLeaf)->GetName());
    histTemp->GetYaxis()->SetBinLabel(iLeaf+1, binTitlesList[iLeaf].c_str());
  }
  histTemp->GetYaxis()->SetLabelSize(histTemp->GetYaxis()->GetLabelSize()/3.);
  histTemp->GetYaxis()->SetTitle("");
  histTemp->GetZaxis()->SetTitle("Deviation (%)");
  histTemp->GetZaxis()->SetRangeUser(0,100);
  histTemp->Write("covDiffAtm_TH2D");

//  histTemp = GenericToolbox::convertTMatrixDtoTH2D(GenericToolbox::convertToCorrelationMatrix(matrixTemp), "corDeltaAtm");
//  for(int iLeaf = 0 ; iLeaf < histTemp->GetNbinsY() ; iLeaf++){
////    histTemp->GetYaxis()->SetBinLabel(iLeaf+1, throwsTree->GetListOfLeaves()->At(iLeaf)->GetName());
//    histTemp->GetYaxis()->SetBinLabel(iLeaf+1, binTitlesList[iLeaf].c_str());
//  }
//  histTemp->GetYaxis()->SetLabelSize(histTemp->GetYaxis()->GetLabelSize()/3.);
//  histTemp->GetYaxis()->SetTitle("");
//  histTemp->Write("corDiffAtm_TH2D");
}
std::vector<double> throwParameters(){

  std::vector<double> throws(fNSystParameters);

  TVectorD Rn(fNSystParameters);

  for (int i = 0; i < fNSystParameters; i++){
    Rn[i] = fGenerator->Gaus();
  }

  TVectorD Incre = Rn;
  Incre *= (*fCholeskyCovMatrix);
  for (int i = 0; i < fNSystParameters; i++){
    throws[i] = 1 + Incre[i];
  }

//  for(int iThrow = 0 ; iThrow < int(throws.size()) ; iThrow++){
//    throws[iThrow] /= TMath::Sqrt((*covMatrixReg)[iThrow][iThrow]);
//  }

  return throws;

}
TMatrixD* getCovarianceMatrixOfTree(TTree* tree_){

  cout << "Computing covariance matrix of TTree: " << tree_ << std::endl;

  // Hook to tree
  cout << " > Hook to tree..." << std::endl;
  std::vector<double> leafValueList(tree_->GetListOfLeaves()->GetEntries(), 0);
  for(int iLeaf = 0 ; iLeaf < tree_->GetListOfLeaves()->GetEntries() ; iLeaf++){
    tree_->SetBranchAddress(tree_->GetListOfLeaves()->At(iLeaf)->GetName(), &leafValueList[iLeaf]);
  }

  // Compute mean of every variable
  std::vector<double> meanValueLeafList(tree_->GetListOfLeaves()->GetEntries(),0);
  for(int iEntry = 0 ; iEntry < tree_->GetEntries() ; iEntry++){
    GenericToolbox::displayProgressBar(iEntry, tree_->GetEntries(), " > Compute mean of every variable...");
    tree_->GetEntry(iEntry);
    for(int iLeaf = 0 ; iLeaf < tree_->GetListOfLeaves()->GetEntries() ; iLeaf++){
      meanValueLeafList[iLeaf] += leafValueList[iLeaf];
    }
  }
  for(int iLeaf = 0 ; iLeaf < tree_->GetListOfLeaves()->GetEntries() ; iLeaf++){
    meanValueLeafList[iLeaf] /= tree_->GetEntries();
  }

  // Compute covariance
  TMatrixD* outCovMatrix = new TMatrixD(tree_->GetListOfLeaves()->GetEntries(), tree_->GetListOfLeaves()->GetEntries());
  for(int iCol = 0 ; iCol < tree_->GetListOfLeaves()->GetEntries() ; iCol++){
    for(int iRow = 0 ; iRow < tree_->GetListOfLeaves()->GetEntries() ; iRow++){
      (*outCovMatrix)[iCol][iRow] = 0;
    }
  }
  for(int iEntry = 0 ; iEntry < tree_->GetEntries() ; iEntry++){
    GenericToolbox::displayProgressBar(iEntry, tree_->GetEntries(), " > Compute covariance...");
    tree_->GetEntry(iEntry);
    for(int iCol = 0 ; iCol < tree_->GetListOfLeaves()->GetEntries() ; iCol++){
      for(int iRow = 0 ; iRow < tree_->GetListOfLeaves()->GetEntries() ; iRow++){
        (*outCovMatrix)[iCol][iRow] += (leafValueList[iCol] - meanValueLeafList[iCol])*(leafValueList[iRow] - meanValueLeafList[iRow]);
      } // iRow
    } // iCol
  }
  for(int iCol = 0 ; iCol < tree_->GetListOfLeaves()->GetEntries() ; iCol++){
    for(int iRow = 0 ; iRow < tree_->GetListOfLeaves()->GetEntries() ; iRow++){
      (*outCovMatrix)[iCol][iRow] /= tree_->GetEntries();
    }
  }

  return outCovMatrix;

}

