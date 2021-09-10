

//std::string filePath = "/Users/ablanche/Documents/Work/Output/results/xsLLhFitter/BANFF_Fit/xsllhGenWeightsFormater/xsec_covariance_2020a_v5b.root";
std::string filePath = "/Users/ablanche/Documents/Work/Repositories/cclyon_root_macros/sk_t2k/skDetCovMatrix.root";
//std::string matrixName = "xsec_cov";
std::string matrixName = "covMatrix";

TFile* covFile;
TMatrixD* covMatrix;
TMatrixD* covSqrtMatrix;

void init();

void fillSqrtCovMatrix( const TMatrixD& V, TMatrixD& C ); // V = covMat, C = sqrtCovMat
std::vector<double> throwParameters( const TMatrixD& sqrtCovMat_ );

// copy/paste from GenericToolbox
TH1D* convertTVectorDtoTH1D(TVectorD *Y_values_, std::string histTitle_ = "", std::string Y_title_ = "", std::string X_title_ = "Entry #", TVectorD *Y_errors_ = nullptr);
TH1D* convertTVectorDtoTH1D(std::vector<double> Y_values_, std::string histTitle_ = "", std::string Y_title_ = "", std::string X_title_ = "Entry #", TVectorD *Y_errors_ = nullptr);


////////////////////
// MAIN
void throwPars(){

  init();

  cout << "Throwing parameters..." << endl;
  vector< double > throws = throwParameters (*covSqrtMatrix);

  for(size_t iThrow = 0 ; iThrow < throws.size() ; iThrow++){
    cout << throws[iThrow] << endl;
  }

  cout << "Printing throws..." << std::endl;
  TH1D* hThrows = convertTVectorDtoTH1D(throws, "throws");
  hThrows->GetXaxis()->SetTitle("Parameter #");
  hThrows->GetYaxis()->SetTitle("Throw in unit of sigmas");
  hThrows->Draw("");

  TFile* testF = TFile::Open("test.root", "RECREATE");
  TTree* testT = new TTree("test", "test");
  vector< double > pars(throws.size(), 0);
  for( int iPar = 0 ; iPar < pars.size() ; iPar++ ){
    testT->Branch(Form("%i", iPar), &pars.at(iPar));
  }

  for( size_t iThrow = 0 ; iThrow < 1000. ; iThrow++ ){
    vector< double > newThrows = throwParameters (*covSqrtMatrix);
    for( int iPar = 0 ; iPar < pars.size() ; iPar++ ){
      pars.at(iPar) = newThrows.at(iPar);
    }
    testT->Fill();
  }
  testF->Close();

}

void init(){
  cout << "Init..." << endl;

  covFile = TFile::Open(filePath.c_str());
  covMatrix = (TMatrixD*) covFile->Get(matrixName.c_str());
  covSqrtMatrix = new TMatrixD( covMatrix->GetNrows(), covMatrix->GetNcols() );
//  covSqrtMatrix = (TMatrixD*) &covMatrix->Sqrt();
  fillSqrtCovMatrix( *covMatrix, *covSqrtMatrix );

}

void fillSqrtCovMatrix( const TMatrixD& V, TMatrixD& C ){
  // calculate sqrt(V) as lower diagonal matrix

  // reset
  for( int i = 0; i < covMatrix->GetNrows(); ++i ) {
    for( int j = 0; j < covMatrix->GetNrows(); ++j ) {
      C[i][j] = 0;
    }
  }

  for( int j = 0; j < covMatrix->GetNrows(); ++j ) {
    // diagonal terms first
    double Ck = 0;
    for( int k = 0; k < j; ++k ) {
      Ck += C[j][k] * C[j][k];
    } // k
    C[j][j] = sqrt( fabs( V[j][j] - Ck ) );

    // off-diagonal terms
    for( int i = j+1; i < covMatrix->GetNrows(); ++i ) {
      Ck = 0;
      for( int k = 0; k < j; ++k ) {
	      Ck += C[i][k] * C[j][k];
      } //k
      C[i][j] = ( V[i][j] - Ck ) / C[j][j];
    }// i
  } // j
}

std::vector<double> throwParameters( const TMatrixD& sqrtCovMat_ ){
  std::vector<double> outParThrows (covMatrix->GetNrows(), 0);
  std::vector<double> randomThrows (covMatrix->GetNrows(), 0);

  gRandom->SetSeed(0);

  // np random numbers from unit Gaussian
  for( int i = 0; i < covMatrix->GetNrows(); ++i ) {
    randomThrows[i] = gRandom->Gaus( 0.0, 1.0 );
  }

  for( int i = 0; i < covMatrix->GetNrows(); ++i ) {
    outParThrows[i] = 0;
    for( int j = 0; j <= i; ++j ) {
      outParThrows[i] += sqrtCovMat_[i][j] * randomThrows[j];
    } // j

    outParThrows[i] /= TMath::Sqrt((*covMatrix)[i][i]);

  } // i

  return outParThrows;

}

TH1D* convertTVectorDtoTH1D(TVectorD* Y_values_, std::string histTitle_, std::string Y_title_, std::string X_title_, TVectorD* Y_errors_){

  auto* th1_histogram = new TH1D(histTitle_.c_str(), histTitle_.c_str(),
                                 Y_values_->GetNrows(), -0.5, Y_values_->GetNrows() - 0.5);

  for(int i_row = 0; i_row < Y_values_->GetNrows(); i_row++)
  {
    th1_histogram->SetBinContent(i_row + 1, (*Y_values_)[i_row]);
    if(Y_errors_ != nullptr)
      th1_histogram->SetBinError(i_row + 1, (*Y_errors_)[i_row]);
  }

  th1_histogram->SetLineWidth(2);
  th1_histogram->SetLineColor(kBlue);
  th1_histogram->GetXaxis()->SetTitle(X_title_.c_str());
  th1_histogram->GetYaxis()->SetTitle(Y_title_.c_str());

  return th1_histogram;
}

TH1D* convertTVectorDtoTH1D(std::vector<double> Y_values_, std::string histTitle_, std::string Y_title_, std::string X_title_, TVectorD *Y_errors_){
  TH1D* out = nullptr;
  auto* tVectorHandler = new TVectorD(Y_values_.size(), &Y_values_[0]);
  out = convertTVectorDtoTH1D(tVectorHandler, histTitle_, Y_title_, X_title_, Y_errors_);
  delete tVectorHandler;
  return out;
}
