

TF1* gausFunction = nullptr;
TFile* outFile = nullptr;
TTree* outTree = nullptr;
bool varMapIsHooked = false;
std::map<std::string, double> varMap;

double pickToyParameter(double oneSigma_);
void init();
void hookToTree();

void detCovMatrixGenerator()
{

  init();

  TFile* atmFile = TFile::Open("/sps/t2k/ablanche/resources/T2KSKjoint/sk4_fcmc_18a_fQv6r0_minituple_100yr_05.root");
  TTree* atm_minituple = (TTree*) atmFile->Get("atm_minituple");

  double nominalCounts = atm_minituple->Draw("", "MReIncLVal>0", "goff");

  std::map<std::string, std::pair<double, double>> throwingRanges;
  throwingRanges["b"].first = -10;
  throwingRanges["b"].second = 10;
  throwingRanges["log_a"].first = -1;
  throwingRanges["log_a"].second = 1;

  // L -> a*L + b
  stringstream rcFormulaeSS;
  double gausNorm = gausFunction->Integral(-5,5);
  double sigma = 0.06;
  int nbThrows = 1000;
  for(int iThrow = 0 ; iThrow < nbThrows ; iThrow++){
  GenericToolbox::getElapsedTimeSinceLastCallStr(0);
  GenericToolbox::getElapsedTimeSinceLastCallStr(1);
    GenericToolbox::displayProgressBar(iThrow, nbThrows);
    varMap["b"] = (throwingRanges["b"].second - throwingRanges["b"].first)*gRandom->Rndm() + throwingRanges["b"].first;
    cout << GenericToolbox::getElapsedTimeSinceLastCallStr(1) << std::endl;
    varMap["log_a"] = (throwingRanges["log_a"].second - throwingRanges["log_a"].first)*gRandom->Rndm() + throwingRanges["log_a"].first;
    cout << GenericToolbox::getElapsedTimeSinceLastCallStr(1) << std::endl;
    varMap["a"] = TMath::Power(10, varMap["log_a"]);
    cout << GenericToolbox::getElapsedTimeSinceLastCallStr(1) << std::endl;

    rcFormulaeSS.str("");
    cout << GenericToolbox::getElapsedTimeSinceLastCallStr(1) << std::endl;
    rcFormulaeSS.clear();
    cout << GenericToolbox::getElapsedTimeSinceLastCallStr(1) << std::endl;
    rcFormulaeSS << varMap["a"] << "*MReIncLVal+" << varMap["b"] << ">0";
    cout << GenericToolbox::getElapsedTimeSinceLastCallStr(1) << std::endl;
    varMap["counts"] = atm_minituple->Draw("", rcFormulaeSS.str().c_str(), "goff");
    cout << GenericToolbox::getElapsedTimeSinceLastCallStr(1) << std::endl;
    varMap["delta_counts"] = varMap["counts"] - nominalCounts;
    cout << GenericToolbox::getElapsedTimeSinceLastCallStr(1) << std::endl;
    varMap["delta_counts_over_counts"] = varMap["delta_counts"]/nominalCounts;
cout << GenericToolbox::getElapsedTimeSinceLastCallStr(1) << std::endl;
    varMap["weight"] = gausFunction->Eval(varMap["delta_counts_over_counts"]/sigma)/gausNorm;
cout << GenericToolbox::getElapsedTimeSinceLastCallStr(1) << std::endl;
    if(not varMapIsHooked){
      hookToTree();
    }
    cout << GenericToolbox::getElapsedTimeSinceLastCallStr(1) << std::endl;
    outTree->Fill();
    cout << GenericToolbox::getElapsedTimeSinceLastCallStr(1) << std::endl;
    cout << "it=" << GenericToolbox::getElapsedTimeSinceLastCallStr(0) << std::endl;
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

double pickToyParameter(double oneSigma_)
{



  return gausFunction->GetRandom()*oneSigma_;

}
