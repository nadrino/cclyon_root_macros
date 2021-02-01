

std::string filePath = "/Users/ablanche/Documents/Work/Output/results/xsLLhFitter/BANFF_Fit/xsllhGenWeightsFormater/xsec_covariance_2020a_v5b.root";
std::string matrixName = "xsec_cov";

TFile* covFile;
TMatrixD* covMatrix;

TF1* gausFunction;

TVectorD *eigenValueList;
TMatrixD *eigenVectorList;

void init();
std::vector<double> getThrownParameters();

TMatrixDSym *convertToSymmetricMatrix(TMatrixD *covMatrix);

void generateToyParamatersFromCov(){

  init();

  cout << "Get throws..." << endl;
  std::vector<double> throws = getThrownParameters();

  GenericToolbox::convertTVectorDtoTH1D(throws, "throws")->Draw();

}


void init(){

  cout << "Init..." << endl;

  covFile = TFile::Open(filePath.c_str());
  covMatrix = (TMatrixD*) covFile->Get(matrixName.c_str());

  // Covariance matrices are symetric :
  auto *symmetric_matrix = convertToSymmetricMatrix(covMatrix);
  auto *eigenMatrixDecomposer = new TMatrixDSymEigen(*symmetric_matrix);
  eigenValueList = (TVectorD*) (eigenMatrixDecomposer->GetEigenValues()).Clone();
  eigenVectorList = (TMatrixD*) (eigenMatrixDecomposer->GetEigenVectors()).Clone();

  if(gausFunction == nullptr){
    gausFunction = new TF1("gausFunction", "gaus", -5, 5);
    gausFunction->SetParameter("Constant", 1);
    gausFunction->SetParameter("Mean", 0);
    gausFunction->SetParameter("Sigma", 1);
  }

}

std::vector<double> getThrownParameters(){

  std::vector<double> eigenThrowList(covMatrix->GetNcols());
  std::vector<double> covThrowList(covMatrix->GetNcols());

  // throw parameters
  for(int iEigen = 0 ; iEigen < covMatrix->GetNcols() ; iEigen++){
    eigenThrowList[iEigen] = gausFunction->GetRandom()*TMath::Sqrt((*eigenValueList)[iEigen]);
    eigenThrowList[iEigen] *= eigenThrowList[iEigen]; // squared in the matrix
  }

  // switch from the eigen vector space to the original space of the cov matrix
  for(int iPar = 0 ; iPar < covMatrix->GetNcols() ; iPar++){
    covThrowList[iPar] = 0;
    for(int iEigen = 0 ; iEigen < covMatrix->GetNcols() ; iEigen++){
      covThrowList[iPar] += eigenThrowList[iEigen] * (*eigenVectorList)[iPar][iEigen];
    }
     covThrowList[iPar] /= (*covMatrix)[iPar][iPar];
  }

  return covThrowList;

}

TMatrixDSym *convertToSymmetricMatrix(TMatrixD *covMatrix) {

  auto *symmetric_matrix = (TMatrixD *) covMatrix->Clone();
  auto *transposed_symmetric_matrix = new TMatrixD(*covMatrix);

  transposed_symmetric_matrix->Transpose(*covMatrix);
  *symmetric_matrix += *transposed_symmetric_matrix;
  for (int i_col = 0; i_col < covMatrix->GetNcols(); i_col++) {
    for (int i_row = 0; i_row < covMatrix->GetNrows(); i_row++) {
      (*symmetric_matrix)[i_row][i_col] /= 2.;
    }
  }

  auto *result = (TMatrixDSym *) symmetric_matrix->Clone(); // Convert to TMatrixDSym

  delete transposed_symmetric_matrix;
  delete symmetric_matrix;

  return result;

}
