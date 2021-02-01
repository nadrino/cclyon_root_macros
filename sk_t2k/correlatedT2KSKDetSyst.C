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

std::string filePath = "/Users/ablanche/Documents/Work/Repositories/cclyon_root_macros/sk_t2k/skDetCovMatrix.root";
std::string matrixName = "covMatrix";

std::string skFilePath = "/Users/ablanche/Documents/Work/Resources/SK-T2K-Joint/sk4_fcmc_18a_fQv6r0_minituple_100yr_05.root";
std::vector<std::string> matrixBinCutsList;
std::vector<TTreeFormula*> formulaList;
TFile* skFile;
TTree* skTree;

std::string skToT2kSamplePath = "/Users/ablanche/Documents/Work/Resources/SK-T2K-Joint/T2KforATMcuts.root";
std::vector<int> skSampleIntList = {11, 13, 14}; // nue / nuebar / nue_edecay
std::map<int,int> entryMapping;
std::map<int,int> entrySkMapping;
std::vector<std::vector<double>> t2kSampleSkEntries;

int nbThrows = 1000;

TFile* covFile;
TMatrixD* covMatrix;
TMatrixDSym* covMatrixReg;
TMatrixD* fCholeskyCovMatrix;

TRandom3* fGenerator;
int fNSystParameters;

void init();
std::vector<double> throwParameters();
TMatrixD* getCovarianceMatrixOfTree(TTree* tree_);

void correlatedT2KSKDetSyst(){

  init();

  std::vector<double> skThrownParameters(covMatrix->GetNrows());
  std::vector<double> t2kPropagatedNorm(3);
  std::vector<int> t2kPropagatedNormCount(3);

  TFile* throwsFile = TFile::Open("throws.root", "RECREATE");
  TTree* throwsTree = new TTree("throws", "throws");

  t2kSampleSkEntries.resize(skSampleIntList.size());
  for(auto& t2kSample : t2kSampleSkEntries){
    t2kSample.resize(skThrownParameters.size());
  }

  cout << "Counting atmospheric events in T2K samples..." << endl;
  for(int iEntry = 0 ; iEntry < skTree->GetEntries() ; iEntry++){

    GenericToolbox::displayProgressBar(iEntry, skTree->GetEntries(), "reading skTree...");

    skTree->GetEntry(iEntry);
    // which SK weight?
    if(not GenericToolbox::doesKeyIsInMap(iEntry, entrySkMapping)){
      for(int iPar = 0 ; iPar < int(skThrownParameters.size()) ; iPar++){

        bool doEventPassCut = false;
        skTree->SetNotify(formulaList[iPar]);
        for(int jInstance = 0; jInstance < formulaList[iPar]->GetNdata(); jInstance++) {
          if ( formulaList[iPar]->EvalInstance(jInstance) ) {
            doEventPassCut = true;
            break;
          }
        }

        if(doEventPassCut) {
          entrySkMapping[iEntry] = iPar;
          t2kSampleSkEntries[entryMapping[iEntry]][iPar]++;
          break;
        }

      }
    }
  }



  for(int iPar = 0 ; iPar < int(skThrownParameters.size()) ; iPar++){
    throwsTree->Branch(Form("Atm_%i", iPar), &skThrownParameters[iPar]);
  }
  for(int iSample = 0 ; iSample < int(t2kPropagatedNorm.size()) ; iSample++){
    throwsTree->Branch(Form("T2K_%i", iSample), &t2kPropagatedNorm[iSample]);
  }

  cout << "Throwing toys..." << endl;
  for(int iThrow = 0 ; iThrow < nbThrows ; iThrow++){

    GenericToolbox::displayProgressBar(iThrow, nbThrows, "throwing...");

    auto throws = throwParameters();

    for(int iPar = 0 ; iPar < int(throws.size()) ; iPar++){
      skThrownParameters[iPar] = throws[iPar];
    }

    for(int iSample = 0 ; iSample < int(t2kSampleSkEntries.size()) ; iSample++){
      t2kPropagatedNorm[iSample] = 0; // reset
      double norm = 0;
      for(int iPar = 0 ; iPar < int(throws.size()) ; iPar++){
        t2kPropagatedNorm[iSample] += t2kSampleSkEntries[iSample][iPar]*skThrownParameters[iPar];
        norm += t2kSampleSkEntries[iSample][iPar];
      } // iPar
      t2kPropagatedNorm[iSample] /= norm;
    } // iSample

    throwsTree->Fill();

  }

  throwsTree->Write("throws");
  covMatrix->Write("skCovMatrix");

  auto* covSkT2kMatrix = getCovarianceMatrixOfTree(throwsTree);
  covSkT2kMatrix->Write("covSkT2kMatrix");

  auto* covHistTemp = GenericToolbox::convertTMatrixDtoTH2D(covSkT2kMatrix, "covSkT2kMatrix");
  for(int iLeaf = 0 ; iLeaf < throwsTree->GetListOfLeaves()->GetEntries() ; iLeaf++){
    covHistTemp->GetYaxis()->SetBinLabel(iLeaf+1, throwsTree->GetListOfLeaves()->At(iLeaf)->GetName());
  }
  covHistTemp->GetYaxis()->SetTitle("");
  covHistTemp->Write("covSkT2kMatrixTH2D");

  auto* corTemp = GenericToolbox::convertToCorrelationMatrix(covSkT2kMatrix);
  corTemp->Write("corSkT2kMatrix");

  auto* corHistTemp = GenericToolbox::convertTMatrixDtoTH2D(corTemp, "corMatrix");
  for(int iLeaf = 0 ; iLeaf < throwsTree->GetListOfLeaves()->GetEntries() ; iLeaf++){
    corHistTemp->GetYaxis()->SetBinLabel(iLeaf+1, throwsTree->GetListOfLeaves()->At(iLeaf)->GetName());
  }
  corHistTemp->GetYaxis()->SetTitle("");
  corHistTemp->Write("corSkT2kMatrixTH2D");

  TMatrixD* covSkT2kMatrixSub = (TMatrixD*) covSkT2kMatrix->Clone();
  for(int iCol = 0 ; iCol < covMatrix->GetNcols() ; iCol++){
    for(int iRow = 0 ; iRow < covMatrix->GetNrows() ; iRow++){
      (*covSkT2kMatrixSub)[iCol][iRow] -= (*covMatrix)[iCol][iRow];
    } // iRow
  } // iCol

  corHistTemp = GenericToolbox::convertTMatrixDtoTH2D(covSkT2kMatrixSub, "covSkT2kMatrixSub");
  for(int iLeaf = 0 ; iLeaf < throwsTree->GetListOfLeaves()->GetEntries() ; iLeaf++){
    corHistTemp->GetYaxis()->SetBinLabel(iLeaf+1, throwsTree->GetListOfLeaves()->At(iLeaf)->GetName());
  }
  corHistTemp->GetYaxis()->SetTitle("");
  corHistTemp->Write("covSkT2kMatrixSubTH2D");

  corHistTemp = GenericToolbox::convertTMatrixDtoTH2D(GenericToolbox::convertToCorrelationMatrix(covSkT2kMatrixSub), "corSkT2kMatrixSub");
  for(int iLeaf = 0 ; iLeaf < throwsTree->GetListOfLeaves()->GetEntries() ; iLeaf++){
    corHistTemp->GetYaxis()->SetBinLabel(iLeaf+1, throwsTree->GetListOfLeaves()->At(iLeaf)->GetName());
  }
  corHistTemp->GetYaxis()->SetTitle("");
  corHistTemp->Write("corSkT2kMatrixSubTH2D");

  throwsFile->Close();

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

  // Hook to tree
  std::vector<double> leafValueList(tree_->GetListOfLeaves()->GetEntries(), 0);
  for(int iLeaf = 0 ; iLeaf < tree_->GetListOfLeaves()->GetEntries() ; iLeaf++){
    tree_->SetBranchAddress(tree_->GetListOfLeaves()->At(iLeaf)->GetName(), &leafValueList[iLeaf]);
  }

  // Compute mean of every variable
  std::vector<double> meanValueLeafList(tree_->GetListOfLeaves()->GetEntries(),0);
  for(int iEntry = 0 ; iEntry < tree_->GetEntries() ; iEntry++){
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

  int trial = 0;
  if(not decompSucceeded){
    cout << "Adding to the diagonal terms to make matrix positive definite..." << endl;

    while(not decompSucceeded and trial < 50){
      for (int i = 0; i < covMatrixReg->GetNcols(); i++)
      {
        (*covMatrixReg)[i][i] += 0.000000000001;
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

  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_elike_0dcy) + " && genmom <  400 && 50 < fqwall < 200");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_elike_0dcy) + " && genmom <  400 && 200 < fqwall");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_elike_0dcy) + " && genmom >= 400 && 50 < fqwall < 200");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_elike_0dcy) + " && genmom >= 400 && 200 < fqwall");

  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_elike_1dcy) + " && genmom <  400 && 50 < fqwall < 200");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_elike_1dcy) + " && genmom <  400 && 200 < fqwall");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_elike_1dcy) + " && genmom >= 400 && 50 < fqwall < 200");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_elike_1dcy) + " && genmom >= 400 && 200 < fqwall");

  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_SingleRing_pi0like));

  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_mulike_0dcy) + " && genmom <  400 && 50 < fqwall < 200");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_mulike_0dcy) + " && genmom <  400 && 200 < fqwall");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_mulike_0dcy) + " && genmom >= 400 && 50 < fqwall < 200");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_mulike_0dcy) + " && genmom >= 400 && 200 < fqwall");

  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_mulike_1dcy) + " && genmom <  400 && 50 < fqwall < 200");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_mulike_1dcy) + " && genmom <  400 && 200 < fqwall");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_mulike_1dcy) + " && genmom >= 400 && 50 < fqwall < 200");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_mulike_1dcy) + " && genmom >= 400 && 200 < fqwall");

  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_mulike_2dcy) + " && genmom <  400 && 50 < fqwall < 200");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_mulike_2dcy) + " && genmom <  400 && 200 < fqwall");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_mulike_2dcy) + " && genmom >= 400 && 50 < fqwall < 200");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_mulike_2dcy) + " && genmom >= 400 && 200 < fqwall");

//  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_pi0like));

  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(MultiGeV_elike_nue) + " && 50 < fqwall < 200");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(MultiGeV_elike_nue) + " && 200 < fqwall");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(MultiGeV_elike_nuebar) + " && 50 < fqwall < 200");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(MultiGeV_elike_nuebar) + " && 200 < fqwall");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(MultiGeV_mulike) + " && 50 < fqwall < 200");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(MultiGeV_mulike) + " && 200 < fqwall");


  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(MultiRing_elike_nue) + " && 50 < fqwall < 200");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(MultiRing_elike_nue) + " && 200 < fqwall");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(MultiRing_elike_nuebar) + " && 50 < fqwall < 200");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(MultiRing_elike_nuebar) + " && 200 < fqwall");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(MultiRing_mulike) + " && 50 < fqwall < 200");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(MultiRing_mulike) + " && 200 < fqwall");
//  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(MultiRingOther_1));

//  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(PCStop));
//  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(PCThru));

  // matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(UpStop_mu));
  // matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(UpThruNonShower_mu));
  // matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(UpThruShower_mu));

  // Formula
  skFile = TFile::Open(skFilePath.c_str());
  skTree = (TTree*) skFile->Get("atm_minituple");

  for(int iPar = 0 ; iPar < int(matrixBinCutsList.size()) ; iPar++){
    formulaList.emplace_back(
      new TTreeFormula(
        Form("skSample_%i", int(iPar)),
        matrixBinCutsList[iPar].c_str(),
        skTree
        )
      );
  }

  //
  TFile* atmToT2kFile = TFile::Open(skToT2kSamplePath.c_str());

  for(int iT2kSample = 0 ; iT2kSample < skSampleIntList.size() ; iT2kSample++){
    TTree* t2kSampleTree = (TTree*) atmToT2kFile->Get(Form("fOutputTree%i", skSampleIntList[iT2kSample]));
//    t2kSampleTree->Print();
//    cout << t2kSampleTree->GetLeaf("SKEntry_d") << endl;
    for(int iEntry = 0 ; iEntry < t2kSampleTree->GetEntries() ; iEntry++){
      t2kSampleTree->GetEntry(iEntry);
      entryMapping[t2kSampleTree->GetLeaf("SKEntry_d")->GetValue(0)] = iT2kSample;
    }
  }



}
