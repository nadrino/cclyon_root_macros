//
// Created by Adrien BLANCHET on 14/12/2020.
//

std::string filePath = "/Users/ablanche/Documents/Work/Repositories/cclyon_root_macros/sk_t2k/skDetCovMatrix.root";
std::string matrixName = "covMatrix";

TFile* covFile;
TMatrixD* covMatrix;
TMatrixDSym* covMatrixReg;
TMatrixD* fCholeskyCovMatrix;

TRandom3* fGenerator;
int fNSystParameters;

void init();
std::vector<double> throwParameters();

void throwCholesky(){

  init();

  auto throws = throwParameters();

  for(auto& throwVal: throws){
    cout << "throw = " << throwVal << endl;
  }

  cout << "Printing throws..." << std::endl;
  TH1D* hThrows = GenericToolbox::convertTVectorDtoTH1D(throws, "throws");
  hThrows->GetXaxis()->SetTitle("Parameter #");
  hThrows->GetYaxis()->SetTitle("Throw in unit of sigmas");
  hThrows->Draw("");

}


std::vector<double> throwParameters(){

  std::vector<double> throws(fNSystParameters);

  TVectorD Rn(fNSystParameters);
  bool IsOK = false;
  while (!IsOK) //Throw until we get values within limits
  {
    IsOK = true;
    for (int i = 0; i < fNSystParameters; i++)
      Rn[i] = fGenerator->Gaus();
    TVectorD Incre = Rn;
    cout << "Incre = " << Incre.GetNrows() << " / fCholeskyCovMatrix = " << fCholeskyCovMatrix->GetNcols() << endl;
    Incre *= (*fCholeskyCovMatrix);
    for (int i = 0; i < fNSystParameters; i++)
    {
      throws[i] = 0 + Incre[i];
//      if (throws[i] < fMinVals[i] || throws[i] > fMaxVals[i]){
//        IsOK = false;
//      }

    }
  }

  for(int iThrow = 0 ; iThrow < int(throws.size()) ; iThrow++){
    throws[iThrow] /= TMath::Sqrt((*covMatrixReg)[iThrow][iThrow]);
  }

  return throws;

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

}