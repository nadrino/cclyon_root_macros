

// http://www-sk.icrr.u-tokyo.ac.jp/indico/event/5454/contributions/15556/attachments/15622/18206/wendell_atm_info_20201117.pdf

std::string outFilePath = "skDetCovMatrix.root";

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

std::map<std::string, int> eventTypesList;
std::map<std::string, TMatrixD*> covMatrixComponents;
std::map<std::string, TMatrixD*> corMatrixComponents;

std::vector<std::string> matrixBinCutsList;
TMatrixD* covMatrixPtr;

void init();
void fillBinCutsList();

void processFcReduction();
void processPcReduction();
void processFcPcSeparation();

void processCosmicRayBackground();
void processFlasherbackgrounds();
void processFiducialVolume();

void processRingSeparation();
void processRingSeparationMR();

void processSingleRingPid();
void processMultiRingPid();

void squareDiag();
void propagateOffDiagCov();
void mergeComponents();
void writeToFile();
void hookLabels(TH2D* hist_);

string getAtmpdEventTypeName(int atmpdEventTypeId_);
int getAtmpdEventType(std::string cutsStr_);
bool isCloseToWall(std::string cutsStr_);
bool isFarToWall(std::string cutsStr_);
bool isLowEnergy(std::string cutsStr_);
bool isHighEnergy(std::string cutsStr_);

bool isFC(int ATMPDEventTypeInt);
bool isPC(int ATMPDEventTypeInt);
bool isSubGeV(int ATMPDEventTypeInt);
bool isMultiGeV(int ATMPDEventTypeInt);
bool isMultiRing(int ATMPDEventTypeInt);
bool isMuLike(int ATMPDEventTypeInt);
bool isELike(int ATMPDEventTypeInt);


////////////////
// MAIN
///////////////
void detSystCovGenerator(){

  cout << "Init..." << endl;
  init();

  cout << "Creating each sub-matrix..." << endl;
  processFcReduction();
  processPcReduction();
  // processFcPcSeparation();

  processCosmicRayBackground();
  processFlasherbackgrounds();

  processRingSeparation();
  processSingleRingPid();
  processMultiRingPid();

  squareDiag();
  propagateOffDiagCov();
  mergeComponents();
  writeToFile();

  exit(0);

}

void processFcReduction(){

  cout << __METHOD_NAME__ << endl;

  // FC Reduction
  // Uncertainty from fully-contained event selection
  // 1.3% normalization error, fully correlated across all FC samples
  covMatrixComponents["FC_Reduction"] = (TMatrixD*) covMatrixPtr->Clone();
  corMatrixComponents["FC_Reduction"] = (TMatrixD*) covMatrixPtr->Clone();

  // Cov diagonal
  for(int iComp = 0 ; iComp < matrixBinCutsList.size() ; iComp++){
    if(isFC(getAtmpdEventType(matrixBinCutsList[iComp]))){
      (*covMatrixComponents["FC_Reduction"])[iComp][iComp] = 1.3;
    }
  }

  // Correlations
  for(int iComp = 0 ; iComp < matrixBinCutsList.size() ; iComp++){
    (*corMatrixComponents["FC_Reduction"])[iComp][iComp] = 1; // obvious
    for(int jComp = 0 ; jComp < matrixBinCutsList.size() ; jComp++){
      if(
            (*covMatrixComponents["FC_Reduction"])[iComp][iComp] != 0
        and (*covMatrixComponents["FC_Reduction"])[jComp][jComp] != 0
      ){
        double sign = (*covMatrixComponents["FC_Reduction"])[iComp][iComp]*(*covMatrixComponents["FC_Reduction"])[jComp][jComp];
        sign /= TMath::Abs(sign);
        (*corMatrixComponents["FC_Reduction"])[iComp][jComp] = 1*sign; // 100 % correlated
      }
    }
  }

}
void processPcReduction(){

  cout << __METHOD_NAME__ << endl;

  // FC Reduction
  // Similarly but for partially-contained reduction
  // 1.0% normalization error, fully correlated across all PC samples
  covMatrixComponents["PC_Reduction"] = (TMatrixD*) covMatrixPtr->Clone();
  corMatrixComponents["PC_Reduction"] = (TMatrixD*) covMatrixPtr->Clone();

  // Cov diagonal
  for(int iComp = 0 ; iComp < matrixBinCutsList.size() ; iComp++){
    if(isPC(getAtmpdEventType(matrixBinCutsList[iComp]))){
      (*covMatrixComponents["PC_Reduction"])[iComp][iComp] = 1.0;
    }
  }

  // Correlations
  for(int iComp = 0 ; iComp < matrixBinCutsList.size() ; iComp++){
    (*corMatrixComponents["PC_Reduction"])[iComp][iComp] = 1; // obvious
    for(int jComp = 0 ; jComp < matrixBinCutsList.size() ; jComp++){
      if(
            (*covMatrixComponents["PC_Reduction"])[iComp][iComp] != 0
        and (*covMatrixComponents["PC_Reduction"])[jComp][jComp] != 0
      ){
        (*corMatrixComponents["PC_Reduction"])[iComp][jComp] = 100*0.01; // 100 % correlated
      }
    }
  }

}
void processFcPcSeparation(){

  cout << __METHOD_NAME__ << endl;

  // FC/PC Separation
  // Uncertainty in relative reductions
  // 0.02% normalization, fully anti-correlated between
  // - FC Single-ring multi-GeV mu-like and PC sample
  // - Total number of FC+PC events
  double sigmaValue = 1.*0.01;
  covMatrixComponents["PC_PC_Separation"] = (TMatrixD*) covMatrixPtr->Clone();
  corMatrixComponents["PC_PC_Separation"] = (TMatrixD*) covMatrixPtr->Clone();

}

void processCosmicRayBackground(){

  // Cosmic Ray background – 1 error parameter
  // Fully correlated normalization change between
  // Sub-GeV mu-like – 0.02%
  // Multi-GeV mu-like – 0.02%
  // Multi-Ring mu-like – 0.07%
  // PC – 0.49%

  // Sub-GeV mu-like – 0.02%
  covMatrixComponents["CosmicRayBackground"] = (TMatrixD*) covMatrixPtr->Clone();
  corMatrixComponents["CosmicRayBackground"] = (TMatrixD*) covMatrixPtr->Clone();

  for(int iComp = 0 ; iComp < matrixBinCutsList.size() ; iComp++){

    if(     isSubGeV(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isMuLike(getAtmpdEventType(matrixBinCutsList[iComp]))
    ){
      (*covMatrixComponents["CosmicRayBackground"])[iComp][iComp] = 0.02;
    }

    if(     isMultiGeV(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isMuLike(  getAtmpdEventType(matrixBinCutsList[iComp]))
    ){
      (*covMatrixComponents["CosmicRayBackground"])[iComp][iComp] = 0.02;
    }

    if(     isMultiGeV(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isMuLike(  getAtmpdEventType(matrixBinCutsList[iComp]))
    ){
      (*covMatrixComponents["CosmicRayBackground"])[iComp][iComp] = 0.07;
    }

    if(     isPC(getAtmpdEventType(matrixBinCutsList[iComp])) ){
      (*covMatrixComponents["CosmicRayBackground"])[iComp][iComp] = 0.49;
    }

  }

  // Correlations
  for(int iComp = 0 ; iComp < matrixBinCutsList.size() ; iComp++){
    (*corMatrixComponents["CosmicRayBackground"])[iComp][iComp] = 1; // obvious
    for(int jComp = 0 ; jComp < matrixBinCutsList.size() ; jComp++){
      if(jComp == iComp) continue;
      if(
            (*covMatrixComponents["CosmicRayBackground"])[iComp][iComp] != 0
        and (*covMatrixComponents["CosmicRayBackground"])[jComp][jComp] != 0
      ){
        double sign = (*covMatrixComponents["CosmicRayBackground"])[iComp][iComp]*(*covMatrixComponents["CosmicRayBackground"])[jComp][jComp];
        sign /= TMath::Abs(sign);
        (*corMatrixComponents["CosmicRayBackground"])[iComp][jComp] = 1*sign; // 100 % correlated
      }
    }
  }

}
void processFlasherbackgrounds(){

  // Flasher backgrounds
  // Fully correlated normalization change between
  // Sub-GeV e-like – 0.03%
  // Multi-GeV e-like – 0.07%

  // Sub-GeV mu-like – 0.02%
  covMatrixComponents["FlasherBackground"] = (TMatrixD*) covMatrixPtr->Clone();
  corMatrixComponents["FlasherBackground"] = (TMatrixD*) covMatrixPtr->Clone();

  for(int iComp = 0 ; iComp < matrixBinCutsList.size() ; iComp++){

    if(     isSubGeV(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isELike(getAtmpdEventType(matrixBinCutsList[iComp]))
    ){
      (*covMatrixComponents["FlasherBackground"])[iComp][iComp] = 0.03;
    }

    if(     isMultiGeV(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isELike(  getAtmpdEventType(matrixBinCutsList[iComp]))
    ){
      (*covMatrixComponents["FlasherBackground"])[iComp][iComp] = 0.07;
    }

  }

  // Correlations
  for(int iComp = 0 ; iComp < matrixBinCutsList.size() ; iComp++){
    (*corMatrixComponents["FlasherBackground"])[iComp][iComp] = 1; // obvious
    for(int jComp = 0 ; jComp < matrixBinCutsList.size() ; jComp++){
      if(jComp == iComp) continue;
      if(
            (*covMatrixComponents["FlasherBackground"])[iComp][iComp] != 0
        and (*covMatrixComponents["FlasherBackground"])[jComp][jComp] != 0
      ){
        double sign = (*covMatrixComponents["FlasherBackground"])[iComp][iComp]*(*covMatrixComponents["FlasherBackground"])[jComp][jComp];
        sign /= TMath::Abs(sign);
        (*corMatrixComponents["FlasherBackground"])[iComp][jComp] = 1*sign; // 100 % correlated
      }
    }
  }

}
void processFiducialVolume(){

  // 2.0% normalization uncertainty, fully correlated for all FC and PC samples

}

void processRingSeparation(){

  // Separate estimates are made for (50 < dwall < 200 and) 200 < dwall for each 1-ring
  // Sub-GeV e-like , prec < 400 MeV : -0.67 (-1.42) %
  // Sub-GeV µ-like , p < 400 MeV : 0.82 (-0.90) %
  // Sub-GeV e-like , p > 400 MeV : 1.11 (1.89) %
  // Sub-GeV µ-like , p > 400 MeV : 2.32 (0.90) %
  // Multi-GeV e-like : 1.21 (8.61) %
  // Multi-GeV µ-like : -2.3% (2.65) %
  // Here the “-” means that at +1s there are fewer 1-ring events
  covMatrixComponents["RingSeparation"] = (TMatrixD*) covMatrixPtr->Clone();
  corMatrixComponents["RingSeparation"] = (TMatrixD*) covMatrixPtr->Clone();

  for(int iComp = 0 ; iComp < matrixBinCutsList.size() ; iComp++){

    // Sub-GeV e-like , prec < 400 MeV : (-1.42) %
    if(     isSubGeV(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isELike(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isLowEnergy(matrixBinCutsList[iComp])
        and isCloseToWall(matrixBinCutsList[iComp])
    ){
      (*covMatrixComponents["RingSeparation"])[iComp][iComp] = -1.42;
    }

    // Sub-GeV e-like , prec < 400 MeV : -0.67 %
    if(     isSubGeV(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isELike(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isLowEnergy(matrixBinCutsList[iComp])
        and isFarToWall(matrixBinCutsList[iComp])
    ){
      (*covMatrixComponents["RingSeparation"])[iComp][iComp] = -0.67;
    }

    // Sub-GeV µ-like , p < 400 MeV : -0.90 %
    if(     isSubGeV(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isMuLike(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isLowEnergy(matrixBinCutsList[iComp])
        and isCloseToWall(matrixBinCutsList[iComp])
    ){
      (*covMatrixComponents["RingSeparation"])[iComp][iComp] = -0.90;
    }

    // Sub-GeV µ-like , p < 400 MeV : 0.82 %
    if(     isSubGeV(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isMuLike(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isLowEnergy(matrixBinCutsList[iComp])
        and isFarToWall(matrixBinCutsList[iComp])
    ){
      (*covMatrixComponents["RingSeparation"])[iComp][iComp] = 0.82;
    }

    // Sub-GeV e-like , p > 400 MeV : 1.89 %
    if(     isSubGeV(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isELike(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isHighEnergy(matrixBinCutsList[iComp])
        and isCloseToWall(matrixBinCutsList[iComp])
    ){
      (*covMatrixComponents["RingSeparation"])[iComp][iComp] = 1.89;
    }

    // Sub-GeV e-like , p > 400 MeV : 1.11 %
    if(     isSubGeV(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isELike(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isHighEnergy(matrixBinCutsList[iComp])
        and isFarToWall(matrixBinCutsList[iComp])
    ){
      (*covMatrixComponents["RingSeparation"])[iComp][iComp] = 1.11;
    }

    // Sub-GeV µ-like , p > 400 MeV : 0.90 %
    if(     isSubGeV(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isMuLike(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isHighEnergy(matrixBinCutsList[iComp])
        and isCloseToWall(matrixBinCutsList[iComp])
    ){
      (*covMatrixComponents["RingSeparation"])[iComp][iComp] = 0.90;
    }

    // Sub-GeV µ-like , p > 400 MeV : 2.32 %
    if(     isSubGeV(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isMuLike(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isHighEnergy(matrixBinCutsList[iComp])
        and isFarToWall(matrixBinCutsList[iComp])
    ){
      (*covMatrixComponents["RingSeparation"])[iComp][iComp] = 2.32;
    }

    // Multi-GeV e-like : 8.61
    if(     isMultiGeV(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isELike(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isCloseToWall(matrixBinCutsList[iComp])
    ){
      (*covMatrixComponents["RingSeparation"])[iComp][iComp] = 8.61;
    }

    // Multi-GeV e-like : 1.21
    if(     isMultiGeV(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isELike(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isFarToWall(matrixBinCutsList[iComp])
    ){
      (*covMatrixComponents["RingSeparation"])[iComp][iComp] = 1.21;
    }

    // Multi-GeV µ-like : 2.65 %
    if(     isMultiGeV(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isELike(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isCloseToWall(matrixBinCutsList[iComp])
    ){
      (*covMatrixComponents["RingSeparation"])[iComp][iComp] = 2.65;
    }

    // Multi-GeV µ-like : -2.3% %
    if(     isMultiGeV(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isELike(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isFarToWall(matrixBinCutsList[iComp])
    ){
      (*covMatrixComponents["RingSeparation"])[iComp][iComp] = -2.3;
    }

//    // Sub-GeV e-like : -2.1 (-3.1) %
//    if (isSubGeV(getAtmpdEventType(matrixBinCutsList[iComp]))
//        and isELike(getAtmpdEventType(matrixBinCutsList[iComp]))
//        and isCloseToWall(matrixBinCutsList[iComp])
//      ) {
//      (*covMatrixComponents["RingSeparation"])[iComp][iComp] = -3.1;
//    }
//
//    // Sub-GeV e-like : -2.1 (-3.1) %
//    if (isSubGeV(getAtmpdEventType(matrixBinCutsList[iComp]))
//        and isELike(getAtmpdEventType(matrixBinCutsList[iComp]))
//        and isFarToWall(matrixBinCutsList[iComp])
//      ) {
//      (*covMatrixComponents["RingSeparation"])[iComp][iComp] = -2.1;
//    }

//    // Sub-GeV µ-like : -2.0 (-3.3)%
//    if (isSubGeV(getAtmpdEventType(matrixBinCutsList[iComp]))
//        and isMuLike(getAtmpdEventType(matrixBinCutsList[iComp]))
//        and isCloseToWall(matrixBinCutsList[iComp])
//      ) {
//      (*covMatrixComponents["RingSeparation"])[iComp][iComp] = -3.3;
//    }
//
//    // Sub-GeV µ-like : -2.0 (-3.3)%
//    if (isSubGeV(getAtmpdEventType(matrixBinCutsList[iComp]))
//        and isMuLike(getAtmpdEventType(matrixBinCutsList[iComp]))
//        and isFarToWall(matrixBinCutsList[iComp])
//      ) {
//      (*covMatrixComponents["RingSeparation"])[iComp][iComp] = -2.0;
//    }

    // Multi-GeV e-like : -0.6 (-4.2) %
    if (isMultiRing(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isELike(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isCloseToWall(matrixBinCutsList[iComp])
      ) {
      (*covMatrixComponents["RingSeparation"])[iComp][iComp] = -4.2;
    }

    // Multi-GeV e-like : -0.6 (-4.2) %
    if (isMultiRing(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isELike(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isFarToWall(matrixBinCutsList[iComp])
      ) {
      (*covMatrixComponents["RingSeparation"])[iComp][iComp] = -0.6;
    }

    // Multi-GeV µ-like : 0.7% (3.1) %
    if (isMultiRing(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isMuLike(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isCloseToWall(matrixBinCutsList[iComp])
      ) {
      (*covMatrixComponents["RingSeparation"])[iComp][iComp] = 3.1;
    }

    // Multi-GeV µ-like : 0.7% (3.1) %
    if (isMultiRing(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isMuLike(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isFarToWall(matrixBinCutsList[iComp])
      ) {
      (*covMatrixComponents["RingSeparation"])[iComp][iComp] = 0.7;
    }

  }

  // Correlations
  for(int iComp = 0 ; iComp < matrixBinCutsList.size() ; iComp++){
    (*corMatrixComponents["RingSeparation"])[iComp][iComp] = 1; // obvious
    for(int jComp = 0 ; jComp < matrixBinCutsList.size() ; jComp++){
      if(jComp == iComp) continue;
      if(
        (*covMatrixComponents["RingSeparation"])[iComp][iComp] != 0
        and (*covMatrixComponents["RingSeparation"])[jComp][jComp] != 0
        ){
        double sign = (*covMatrixComponents["RingSeparation"])[iComp][iComp]*(*covMatrixComponents["RingSeparation"])[jComp][jComp];
        sign /= TMath::Abs(sign);
        (*corMatrixComponents["RingSeparation"])[iComp][jComp] = 1*sign; // 100 % correlated
      }

    }
  }

}
void processSingleRingPid(){

  // Separate estimates are made for (50 < dwall < 200 and) 200 < dwall for each 1-ring
  // Sub-GeV e-like : 0.36 (0.99) %
  // Sub-GeV µ-like : -0.37 (-0.90)%
  // Multi-GeV e-like : 0.06 (0.23) %
  // Multi-GeV µ-like : -0.06 (-0.21) %
  // At the end there is one fitting parameter for “single-ring PID” and its effect at 1-sigma is
  // encoded with these numbers
  // It is not correlated with the ring-counting error parameter
  covMatrixComponents["SingleRingPid"] = (TMatrixD *) covMatrixPtr->Clone();
  corMatrixComponents["SingleRingPid"] = (TMatrixD *) covMatrixPtr->Clone();

  for (int iComp = 0; iComp < matrixBinCutsList.size(); iComp++) {

    // Sub-GeV e-like : 0.36 (0.99) %
    if (isSubGeV(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isELike(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isCloseToWall(matrixBinCutsList[iComp])
      ) {
      (*covMatrixComponents["SingleRingPid"])[iComp][iComp] = 0.99;
    }

    // Sub-GeV e-like : 0.36
    if (isSubGeV(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isELike(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isFarToWall(matrixBinCutsList[iComp])
      ) {
      (*covMatrixComponents["SingleRingPid"])[iComp][iComp] = 0.36;
    }

    // Sub-GeV µ-like : -0.37 (-0.90)%
    if (isSubGeV(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isMuLike(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isCloseToWall(matrixBinCutsList[iComp])
      ) {
      (*covMatrixComponents["SingleRingPid"])[iComp][iComp] = -0.90;
    }

    // Sub-GeV µ-like : -0.37 (-0.90)%
    if (isSubGeV(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isMuLike(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isFarToWall(matrixBinCutsList[iComp])
      ) {
      (*covMatrixComponents["SingleRingPid"])[iComp][iComp] = -0.37;
    }

    // Multi-GeV e-like : 0.06 (0.23) %
    if (isMultiGeV(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isELike(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isCloseToWall(matrixBinCutsList[iComp])
      ) {
      (*covMatrixComponents["SingleRingPid"])[iComp][iComp] = 0.23;
    }

    // Multi-GeV e-like : 0.06 (0.23) %
    if (isMultiGeV(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isELike(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isFarToWall(matrixBinCutsList[iComp])
      ) {
      (*covMatrixComponents["SingleRingPid"])[iComp][iComp] = 0.06;
    }

    // Multi-GeV µ-like : -0.06 (-0.21) %
    if (isMultiGeV(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isMuLike(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isCloseToWall(matrixBinCutsList[iComp])
      ) {
      (*covMatrixComponents["SingleRingPid"])[iComp][iComp] = -0.21;
    }

    // Multi-GeV µ-like : -0.06 (-0.21) %
    if (isMultiGeV(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isMuLike(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isFarToWall(matrixBinCutsList[iComp])
      ) {
      (*covMatrixComponents["SingleRingPid"])[iComp][iComp] = -0.06;
    }

  }

  // Correlations
  for(int iComp = 0 ; iComp < matrixBinCutsList.size() ; iComp++){
    (*corMatrixComponents["SingleRingPid"])[iComp][iComp] = 1; // obvious
    for(int jComp = 0 ; jComp < matrixBinCutsList.size() ; jComp++){
      if(jComp == iComp) continue;
      if(
        (*covMatrixComponents["SingleRingPid"])[iComp][iComp] != 0
        and (*covMatrixComponents["SingleRingPid"])[jComp][jComp] != 0
        ){
        double sign = (*covMatrixComponents["SingleRingPid"])[iComp][iComp]*(*covMatrixComponents["SingleRingPid"])[jComp][jComp];
        sign /= TMath::Abs(sign);
        (*corMatrixComponents["SingleRingPid"])[iComp][jComp] = 1*sign; // 100 % correlated
      }

    }
  }

}
void processMultiRingPid(){

  cout << __METHOD_NAME__ << endl;

  // Separate estimates are made for (50 < dwall < 200 and) 200 < dwall for each 1-ring
  // Separate estimates are made for (50 < dwall < 200 and) 200 < dwall for each multi-ring
  // Sub-GeV e-like : -0.72 (-3.2) % -> ???
  // Sub-GeV µ-like : 0.31 (1.3)% -> ???
  // Multi-GeV e-like : 1.1 (1.9) %
  // Multi-GeV µ-like : -0.06 (-1.0) %
  covMatrixComponents["MultiRingPid"] = (TMatrixD *) covMatrixPtr->Clone();
  corMatrixComponents["MultiRingPid"] = (TMatrixD *) covMatrixPtr->Clone();

  for (int iComp = 0; iComp < matrixBinCutsList.size(); iComp++) {

    // Multi-GeV e-like : 1.1 (1.9) %
    if (isMultiRing(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isELike(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isCloseToWall(matrixBinCutsList[iComp])
      ) {
      (*covMatrixComponents["MultiRingPid"])[iComp][iComp] = 1.9;
    }

    // Multi-GeV e-like : 1.1 (1.9) %
    if (isMultiRing(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isELike(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isFarToWall(matrixBinCutsList[iComp])
      ) {
      (*covMatrixComponents["MultiRingPid"])[iComp][iComp] = 1.1;
    }

    // Multi-GeV µ-like : -0.06 (-1.0) %
    if (isMultiRing(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isMuLike(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isCloseToWall(matrixBinCutsList[iComp])
      ) {
      (*covMatrixComponents["MultiRingPid"])[iComp][iComp] = -1.0;
    }

    // Multi-GeV µ-like : -0.06 (-1.0) %
    if (isMultiRing(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isMuLike(getAtmpdEventType(matrixBinCutsList[iComp]))
        and isFarToWall(matrixBinCutsList[iComp])
      ) {
      (*covMatrixComponents["MultiRingPid"])[iComp][iComp] = -0.06;
    }

  }

  // Correlations
  for(int iComp = 0 ; iComp < matrixBinCutsList.size() ; iComp++){
    (*corMatrixComponents["MultiRingPid"])[iComp][iComp] = 1; // obvious
    for(int jComp = 0 ; jComp < matrixBinCutsList.size() ; jComp++){

      if(jComp == iComp) continue;

      if(
        (*covMatrixComponents["MultiRingPid"])[iComp][iComp] != 0
        and (*covMatrixComponents["MultiRingPid"])[jComp][jComp] != 0
        ){
        double sign = (*covMatrixComponents["MultiRingPid"])[iComp][iComp]*(*covMatrixComponents["MultiRingPid"])[jComp][jComp];
        sign /= TMath::Abs(sign);
        (*corMatrixComponents["MultiRingPid"])[iComp][jComp] = 1*sign; // 100 % correlated
      }

    }
  }

}

void squareDiag(){

  cout << __METHOD_NAME__ << endl;

  for(auto covElement : covMatrixComponents){
    for(int iComp = 0 ; iComp < matrixBinCutsList.size() ; iComp++){
      (*covElement.second)[iComp][iComp] = (*covElement.second)[iComp][iComp]*0.01; // IN %
      (*covElement.second)[iComp][iComp] *= (*covElement.second)[iComp][iComp]; // SQUARE
    }
  }

}
void propagateOffDiagCov(){

  cout << __METHOD_NAME__ << endl;

  for(auto covElement : covMatrixComponents){
    for(int iComp = 0 ; iComp < matrixBinCutsList.size() ; iComp++){
      for(int jComp = 0 ; jComp < matrixBinCutsList.size() ; jComp++){
        if(iComp == jComp) continue;
        (*covElement.second)[iComp][jComp] =
          TMath::Sqrt((*covElement.second)[iComp][iComp] * (*covElement.second)[jComp][jComp])
          * (*corMatrixComponents[covElement.first])[iComp][jComp];
      }
    }
  }

}
void mergeComponents(){

  cout << __METHOD_NAME__ << endl;

  for(auto covElement : covMatrixComponents){
    for(int iComp = 0 ; iComp < matrixBinCutsList.size() ; iComp++){
      for(int jComp = 0 ; jComp < matrixBinCutsList.size() ; jComp++){
        (*covMatrixPtr)[iComp][jComp] += (*covElement.second)[iComp][jComp];
      }
    }
  }

}
void writeToFile(){

  cout << __METHOD_NAME__ << endl;

  cout << "OUTPUT FILE: " << outFilePath << endl;
  TFile* outFile = TFile::Open(outFilePath.c_str(), "RECREATE");

  covMatrixPtr->Write("covMatrix");
  auto* covHistTemp = GenericToolbox::convertTMatrixDtoTH2D(covMatrixPtr, "covMatrix");
  hookLabels(covHistTemp);
  covHistTemp->Write("covMatrixTH2D");

  auto* corTemp = GenericToolbox::convertToCorrelationMatrix(covMatrixPtr);
  corTemp->Write("corMatrix");

  auto* corHistTemp = GenericToolbox::convertTMatrixDtoTH2D(corTemp, "corMatrix");
  hookLabels(corHistTemp);
  corHistTemp->GetZaxis()->SetRangeUser(-1,1);
  corHistTemp->Write("corMatrixTH2D");

  // https://root-forum.cern.ch/t/save-a-vector-or-array-or-tobjarray-of-strings-or-tstrings-or-tobjstrings-to-a-tfile/19537
  // outFile->WriteObjectAny(&matrixBinCutsList,"vector<string>","myStrings");
  // TObjArray componentsCutsList;
  // TObjString sPtr[matrixBinCutsList.size()];
  // for(size_t iCut = 0 ; iCut < matrixBinCutsList.size() ; iCut++ ){
  //   sPtr[iCut] = TObjString(matrixBinCutsList[iCut].c_str());
  //   // componentsCutsList.Add(sPtr);
  // }
  // TObjArray componentsCutsList(sPtr);
  // componentsCutsList.Write("componentsCutsList");
  outFile->mkdir("components");
  outFile->cd("components");

  for(auto covElement : covMatrixComponents){

    outFile->mkdir(Form("components/%s",covElement.first.c_str()));
    outFile->cd(Form("components/%s",covElement.first.c_str()));

    covElement.second->Write(covElement.first.c_str());
    GenericToolbox::convertTMatrixDtoTH2D(covElement.second, covElement.first)->Write((covElement.first + "TH2D").c_str());
    GenericToolbox::convertTMatrixDtoTH2D(corMatrixComponents[covElement.first], covElement.first)->Write((covElement.first + "CorrTH2D").c_str());

    TMatrixD* sqrtMatrixTemp = (TMatrixD*) &covElement.second->Sqrt();
    GenericToolbox::convertTMatrixDtoTH2D(sqrtMatrixTemp, covElement.first)->Write((covElement.first + "SqrtTH2D").c_str());

  }

  outFile->Close();

}
void hookLabels(TH2D* hist_){

  for(int iBin = 0 ; iBin < int(matrixBinCutsList.size()) ; iBin++){
    auto cutsList = GenericToolbox::splitString(matrixBinCutsList[iBin], " && ");
    std::string parsedName = getAtmpdEventTypeName(getAtmpdEventType(matrixBinCutsList[iBin]));
    parsedName += " && " + GenericToolbox::joinVectorString(cutsList, " && ", 1);
    hist_->GetYaxis()->SetBinLabel(iBin+1, parsedName.c_str());
  }
  hist_->GetZaxis()->SetLabelSize( hist_->GetZaxis()->GetLabelSize()/2. );
  hist_->GetYaxis()->SetLabelSize( hist_->GetYaxis()->GetLabelSize()/2. );

  hist_->GetYaxis()->SetTitle("");

}

void fillBinCutsList(){

  cout << __METHOD_NAME__ << endl;

  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_elike_0dcy) + " && genmom <  400 && 50 < dwall < 200");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_elike_0dcy) + " && genmom <  400 && 200 < dwall");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_elike_0dcy) + " && genmom >= 400 && 50 < dwall < 200");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_elike_0dcy) + " && genmom >= 400 && 200 < dwall");

  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_elike_1dcy) + " && genmom <  400 && 50 < dwall < 200");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_elike_1dcy) + " && genmom <  400 && 200 < dwall");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_elike_1dcy) + " && genmom >= 400 && 50 < dwall < 200");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_elike_1dcy) + " && genmom >= 400 && 200 < dwall");

  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_SingleRing_pi0like));

  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_mulike_0dcy) + " && genmom <  400 && 50 < dwall < 200");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_mulike_0dcy) + " && genmom <  400 && 200 < dwall");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_mulike_0dcy) + " && genmom >= 400 && 50 < dwall < 200");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_mulike_0dcy) + " && genmom >= 400 && 200 < dwall");

  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_mulike_1dcy) + " && genmom <  400 && 50 < dwall < 200");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_mulike_1dcy) + " && genmom <  400 && 200 < dwall");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_mulike_1dcy) + " && genmom >= 400 && 50 < dwall < 200");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_mulike_1dcy) + " && genmom >= 400 && 200 < dwall");

  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_mulike_2dcy) + " && genmom <  400 && 50 < dwall < 200");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_mulike_2dcy) + " && genmom <  400 && 200 < dwall");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_mulike_2dcy) + " && genmom >= 400 && 50 < dwall < 200");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_mulike_2dcy) + " && genmom >= 400 && 200 < dwall");

//  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(SubGeV_pi0like));

  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(MultiGeV_elike_nue) + " && 50 < dwall < 200");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(MultiGeV_elike_nue) + " && 200 < dwall");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(MultiGeV_elike_nuebar) + " && 50 < dwall < 200");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(MultiGeV_elike_nuebar) + " && 200 < dwall");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(MultiGeV_mulike) + " && 50 < dwall < 200");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(MultiGeV_mulike) + " && 200 < dwall");


  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(MultiRing_elike_nue) + " && 50 < dwall < 200");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(MultiRing_elike_nue) + " && 200 < dwall");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(MultiRing_elike_nuebar) + " && 50 < dwall < 200");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(MultiRing_elike_nuebar) + " && 200 < dwall");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(MultiRing_mulike) + " && 50 < dwall < 200");
  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(MultiRing_mulike) + " && 200 < dwall");
//  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(MultiRingOther_1));

//  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(PCStop));
//  matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(PCThru));

  // matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(UpStop_mu));
  // matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(UpThruNonShower_mu));
  // matrixBinCutsList.emplace_back("ATMPDEventType == " + std::to_string(UpThruShower_mu));

  for(size_t iCut = 0 ; iCut < matrixBinCutsList.size() ; iCut++ ){
    cout << "Component #" << iCut << ": \"" << matrixBinCutsList[iCut] << "\"" << endl;
  }

}
void init(){

  cout << __METHOD_NAME__ << endl;

  eventTypesList["SubGeV_elike_0dcy"] = 1;
  eventTypesList["SubGeV_elike_1dcy"] = 2;
  eventTypesList["SubGeV_SingleRing_pi0like"] = 3;
  eventTypesList["SubGeV_mulike_0dcy"] = 4;
  eventTypesList["SubGeV_mulike_1dcy"] = 5;
  eventTypesList["SubGeV_mulike_2dcy"] = 6;
  eventTypesList["SubGeV_pi0like"] = 7;
  eventTypesList["MultiGeV_elike_nue"] = 8;
  eventTypesList["MultiGeV_elike_nuebar"] = 9;
  eventTypesList["MultiGeV_mulike"] = 10;
  eventTypesList["MultiRing_elike_nue"] = 11;
  eventTypesList["MultiRing_elike_nuebar"] = 12;
  eventTypesList["MultiRing_mulike"] = 13;
  eventTypesList["MultiRingOther_1"] = 14;
  eventTypesList["PCStop"] = 15;
  eventTypesList["PCThru"] = 16;
  eventTypesList["UpStop_mu"] = 17;
  eventTypesList["UpThruNonShower_mu"] = 18;
  eventTypesList["UpThruShower_mu"] = 19;

  fillBinCutsList();

  covMatrixPtr = new TMatrixD(matrixBinCutsList.size(), matrixBinCutsList.size());
//  for(int iRow = 0 ; iRow < covMatrixPtr->GetNrows() ; iRow++){
//
//  }

}

string getAtmpdEventTypeName(int atmpdEventTypeId_){
  for(const auto& eventType : eventTypesList){
    if(eventType.second == atmpdEventTypeId_){
      return eventType.first;
    }
  }
  return "";
}
int getAtmpdEventType(std::string cutsStr_){
  // Get ATMPDEventType
  auto tempVec = GenericToolbox::splitString(cutsStr_, " && ");
  int ATMPDEventTypeInt = 0;
  for(const auto& vecElem : tempVec){
    if(GenericToolbox::doesStringStartsWithSubstring(vecElem, "ATMPDEventType")){
      ATMPDEventTypeInt = stoi(GenericToolbox::splitString(vecElem, " == ")[1]);
    }
  }
  return ATMPDEventTypeInt;
}
bool isCloseToWall(std::string cutsStr_){
  if(not GenericToolbox::doesStringContainsSubstring(cutsStr_, "dwall")) return false;
  return GenericToolbox::doesStringContainsSubstring(cutsStr_, "50 < dwall < 200");
}
bool isFarToWall(std::string cutsStr_){
  if(not GenericToolbox::doesStringContainsSubstring(cutsStr_, "dwall")) return false;
  return GenericToolbox::doesStringContainsSubstring(cutsStr_, "200 < dwall");
}
bool isLowEnergy(std::string cutsStr_){
  if(not GenericToolbox::doesStringContainsSubstring(cutsStr_, "genmom")) return false;
  return GenericToolbox::doesStringContainsSubstring(cutsStr_, "genmom <");
}
bool isHighEnergy(std::string cutsStr_){
  if(not GenericToolbox::doesStringContainsSubstring(cutsStr_, "genmom")) return false;
  return GenericToolbox::doesStringContainsSubstring(cutsStr_, "genmom >");
}

bool isFC(int ATMPDEventTypeInt){
  bool isFC = true;
  for(const auto& eventType : eventTypesList){
    if(eventType.second == ATMPDEventTypeInt){
      if(    GenericToolbox::doesStringStartsWithSubstring(eventType.first, "PC")
          or GenericToolbox::doesStringStartsWithSubstring(eventType.first, "Up")){
        isFC = false;
        break;
      }
    }
  }
  return isFC;
}
bool isPC(int ATMPDEventTypeInt){
  bool isPC = false;
  for(const auto& eventType : eventTypesList){
    if(eventType.second == ATMPDEventTypeInt){
      if( GenericToolbox::doesStringStartsWithSubstring(eventType.first, "PC") ){
        isPC = true;
        break;
      }
    }
  }
  return isPC;
}
bool isSubGeV(int ATMPDEventTypeInt){
  bool isSubGeV = false;
  for(const auto& eventType : eventTypesList){
    if(eventType.second == ATMPDEventTypeInt){
      if( GenericToolbox::doesStringStartsWithSubstring(eventType.first, "SubGeV") ){
        isSubGeV = true;
        break;
      }
    }
  }
  return isSubGeV;
}
bool isMultiGeV(int ATMPDEventTypeInt){
  bool isMultiGeV = false;
  for(const auto& eventType : eventTypesList){
    if(eventType.second == ATMPDEventTypeInt){
      if( GenericToolbox::doesStringStartsWithSubstring(eventType.first, "MultiGeV") ){
        isMultiGeV = true;
        break;
      }
    }
  }
  return isMultiGeV;
}
bool isMultiRing(int ATMPDEventTypeInt){
  bool isMultiRing = false;
  for(const auto& eventType : eventTypesList){
    if(eventType.second == ATMPDEventTypeInt){
      if( GenericToolbox::doesStringStartsWithSubstring(eventType.first, "MultiRing") ){
        isMultiRing = true;
        break;
      }
    }
  }
  return isMultiRing;
}

bool isMuLike(int ATMPDEventTypeInt){
  bool isMuLike = false;
  for(const auto& eventType : eventTypesList){
    if(eventType.second == ATMPDEventTypeInt){
      if( GenericToolbox::doesStringContainsSubstring(eventType.first, "mulike") ){
        isMuLike = true;
        break;
      }
    }
  }
  return isMuLike;
}
bool isELike(int ATMPDEventTypeInt){
  bool isELike = false;
  for(const auto& eventType : eventTypesList){
    if(eventType.second == ATMPDEventTypeInt){
      if( GenericToolbox::doesStringContainsSubstring(eventType.first, "elike") ){
        isELike = true;
        break;
      }
    }
  }
  return isELike;
}
