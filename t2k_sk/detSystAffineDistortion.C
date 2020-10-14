

TF1* gausFunction = nullptr;
TFile* outFile = nullptr;
TTree* outTree = nullptr;
bool varMapIsHooked = false;
std::map<std::string, Float_t> varMap;
std::vector<Float_t> observableValueList;

double pickToyParameter(double oneSigma_);
void init();
void hookToTree();
double evalGaus(double x);

void detSystAffineDistortion()
{

  init();

  TFile* atmFile = TFile::Open("/sps/t2k/ablanche/resources/SK-T2K-Joint/sk4_fcmc_18a_fQv6r0_minituple_100yr_05.root");
  TTree* atm_minituple = (TTree*) atmFile->Get("atm_minituple");
  // atm_minituple->SetBranchStatus("*", false);
  // atm_minituple->SetBranchStatus("MReIncLVal", true);
  //
  // Float_t MReIncLVal;
  // atm_minituple->SetBranchAddress("MReIncLVal", &MReIncLVal);
  // for(int iEntry = 0 ; iEntry < atm_minituple->GetEntries() ; iEntry++){
  //     GenericToolbox::displayProgressBar(iEntry, atm_minituple->GetEntries(), "Reading...");
  //   atm_minituple->GetEntry(iEntry);
  //   observableValueList.emplace_back(MReIncLVal);
  // }

  double nominalCounts = atm_minituple->Draw("", "MReIncLVal>-0.25", "goff");

  std::map<std::string, std::pair<double, double>> throwingRanges;

  throwingRanges["b"].first = -15;
  throwingRanges["b"].second = 15;

  throwingRanges["a"].first = 0.1;
  throwingRanges["a"].second = 10;

  TH1D* deviation_hist = new TH1D("deviation_hist", "deviation_hist", 100, -0.2, 0.2);

  // L -> a*L + b
  string rcFormulae;
  double gausNorm = gausFunction->Integral(-5,5);
  double sigma = 0.06;
  int nbThrows = 10000;
  std::vector<std::map<std::string, Float_t>> throwsData;
  for(int iThrow = 0 ; iThrow < nbThrows ; iThrow++){

    GenericToolbox::displayProgressBar(iThrow, nbThrows, "Throwing parameters...");
    // bool bSignIsPositive = 0.5 > gRandom->Rndm();
    for(auto &throwingRange : throwingRanges){
      varMap[throwingRange.first] = (throwingRange.second.second - throwingRange.second.first)*gRandom->Rndm() + throwingRange.second.first;
    }
    // varMap["a"] = TMath::Power(10, varMap["log_a"]);
    // varMap["a"] = 1;
    // varMap["b"] = TMath::Power(10, varMap["log_b"]);
    // if(!bSignIsPositive) varMap["b"] *= -1;

    rcFormulae = std::to_string(varMap["a"]);
    rcFormulae += "*MReIncLVal+";
    rcFormulae += std::to_string(varMap["b"]);
    rcFormulae += ">-0.25";
    varMap["counts"] = atm_minituple->Draw("", rcFormulae.c_str(), "goff");
    varMap["delta_counts"] = varMap["counts"] - nominalCounts;
    varMap["delta_counts_over_counts"] = varMap["delta_counts"]/nominalCounts;
    deviation_hist->Fill(varMap["delta_counts_over_counts"]);
    throwsData.emplace_back(varMap);

  }

  varMap["weight"] = 0;
  hookToTree();

  for(int iThrow = 0 ; iThrow < nbThrows ; iThrow++){
    GenericToolbox::displayProgressBar(iThrow, nbThrows, "Writting throws...");
    for(auto &var : varMap){
      var.second = throwsData[iThrow][var.first];
    }
    varMap["weight"] = evalGaus( varMap["delta_counts_over_counts"] / sigma ) /
      deviation_hist->GetBinContent(deviation_hist->FindBin(varMap["delta_counts_over_counts"]));
    outTree->Fill();
  }

  outFile->WriteTObject(outTree, "outTree");
  outFile->Close();

}

void hookToTree()
{
    outFile = TFile::Open("outFile.root", "RECREATE");
    outTree = new TTree("outTree", "outTree");
    for(auto& varMapPair : varMap){
      outTree->Branch(varMapPair.first.c_str(), &varMapPair.second);
    }
    varMapIsHooked = true;
}

void init()
{
  if(gausFunction == nullptr){
    gausFunction = new TF1("gausFunction", "gaus", -5, 5);
    gausFunction->SetParameter("Constant", 1);
    gausFunction->SetParameter("Mean", 0);
    gausFunction->SetParameter("Sigma", 1);
  }
}

double evalGaus(double x){
  return exp(-0.5*x*x)/(sqrt(2*TMath::Pi()));
}

double pickToyParameter(double oneSigma_)
{



  return gausFunction->GetRandom()*oneSigma_;

}
