#ifndef SHAREDPARS_C
#define SHAREDPARS_C


#include "sharedPars.h"

using namespace std;

double sharedPars::getParD(const char* parname){
  return kr->getKeyD(parname);
}

TString sharedPars::getParS(const char* parname){
  return kr->getKeyS(parname);
}

int sharedPars::getParI(const char* parname){
  return kr->getKeyI(parname);
}

void sharedPars::readParsFromFile(){

  //create object to read in keys from file
  kr = new keyread(parFileName.Data());

  //read in contents
  kr->readFile();

  //set parameters to values
  physLoBoundAtt0 = kr->getKeyD("physLoBoundAtt0");
  physLoBoundAtt1 = kr->getKeyD("physLoBoundAtt1");
  physLoBoundAtt2 = kr->getKeyD("physLoBoundAtt2");
  physLoBoundAtt3 = kr->getKeyD("physLoBoundAtt3");
  physLoBoundAtt4 = kr->getKeyD("physLoBoundAtt4");
  flgUsePriorsInFit = kr->getKeyI("flgUsePriorsInFit");
  useSplinesFlg = kr->getKeyI("useSplinesFlg");
  MCMCNSteps= kr->getKeyI("MCMCNSteps");
  MCMCNBurnIn = kr->getKeyI("MCMCNBurnIn");
  MCMCTunePar = kr->getKeyD("MCMCTunePar");
  normFactor = kr->getKeyD("normFactor");
  fixAllSmearFlg = kr->getKeyI("fixAllSmearFlg");
  FVBinName0= kr->getKeyS("FVBinName0");
  FVBinName1=kr->getKeyS("FVBinName1");
  FVBinName2=kr->getKeyS("FVBinName2");
  fQAttName0=kr->getKeyS("fQAttName0");
  fQAttName1=kr->getKeyS("fQAttName1");
  fQAttName2=kr->getKeyS("fQAttName2");
  fQAttName3=kr->getKeyS("fQAttName3");
  fQAttName4=kr->getKeyS("fQAttName4");
  fQAttName5=kr->getKeyS("fQAttName5");
  fQAttName6=kr->getKeyS("fQAttName6");
  fQAttName7=kr->getKeyS("fQAttName7");
  fQAttName8=kr->getKeyS("fQAttName8");
  hFactoryDataFiles=kr->getKeyS("hFactoryDataFiles");
  hFactoryMCFiles=kr->getKeyS("hFactoryMCFiles");
  hFactoryOutput=kr->getKeyS("hFactoryOutput");
  sampleName0=kr->getKeyS("sampleName0");
  sampleName1=kr->getKeyS("sampleName1");
  sampleName2=kr->getKeyS("sampleName2");
  sysParName0=kr->getKeyS("sysParName0");
  sysParName1=kr->getKeyS("sysParName1");
  sysParName2=kr->getKeyS("sysParName2");
  sysParName3=kr->getKeyS("sysParName3");
  sysParName4=kr->getKeyS("sysParName4");
  sysParName5=kr->getKeyS("sysParName5");
  sysParName6=kr->getKeyS("sysParName6");
  sysParName7=kr->getKeyS("sysParName7");
  sysParName8=kr->getKeyS("sysParName8");
  nFVBins     = kr->getKeyI("nFVBins");
  nSamples    = kr->getKeyI("nSamples");
  nComponents = kr->getKeyI("nComponents");
  nAttributes  = kr->getKeyI("nAttributes"); 
  nSysPars    = kr->getKeyI("nSysPars");
  smearPriorWidthAtt0 = kr->getKeyD("smearPriorWidthAtt0");
  smearPriorWidthAtt1 = kr->getKeyD("smearPriorWidthAtt1");
  smearPriorWidthAtt2 = kr->getKeyD("smearPriorWidthAtt2");
  smearPriorWidthAtt3 = kr->getKeyD("smearPriorWidthAtt3");
  smearPriorWidthAtt4 = kr->getKeyD("smearPriorWidthAtt4");
  smearPriorWidthAtt5 = kr->getKeyD("smearPriorWidthAtt5");
  biasPriorWidthAtt0 = kr->getKeyD("biasPriorWidthAtt0");
  biasPriorWidthAtt1 = kr->getKeyD("biasPriorWidthAtt1");
  biasPriorWidthAtt2 = kr->getKeyD("biasPriorWidthAtt2");
  biasPriorWidthAtt3 = kr->getKeyD("biasPriorWidthAtt3");
  biasPriorWidthAtt4 = kr->getKeyD("biasPriorWidthAtt4");
  biasPriorWidthAtt5 = kr->getKeyD("biasPriorWidthAtt5");
  flgUsePhysLoBound = kr->getKeyI("flgUsePhysLoBound");
  flgUseFitParFile = kr->getKeyI("flgUseFitParFile");
  nBinBuffer = kr->getKeyI("nBinBuffer");
#ifdef NMODE
  nModes = NMODE;
#endif
  preProcessFilesMC = kr->getKeyS("preProcessFilesMC"); 
  preProcessOutDir = kr->getKeyS("preProcessOutDir"); 
  preProcessFilesData = kr->getKeyS("preProcessFilesData"); 
  preProcessFilesSpline = kr->getKeyS("preProcessFilesSpline");
  preProcessFilesBANFF = kr->getKeyS("preProcessFilesBANFF");
  preProcessMCComponents = kr->getKeyI("preProcessMCComponents");
  preProcessFVBinning = kr->getKeyI("preProcessFVBinning");
  preProcessMCSamples = kr->getKeyI("preProcessMCSamples");
  preProcAddMoreVars = kr->getKeyI("preProcAddMoreVars");
  preProcMaskFile = kr->getKeyS("preProcMaskFile");
  preProcMaskFlg = kr->getKeyI("preProcMaskFlg");
  ntupleType=kr->getKeyS("ntupleType");
  globalRootName = kr->getKeyS("globalRootName");
  splineFactoryOutput = kr->getKeyS("splineFactoryOutput");
  sysParType = kr->getKeyS("sysParType");
  NMCMCPts = kr->getKeyI("NMCMCPts");
  MCMCBurnIn=kr->getKeyI("MCMCBurnIn");
  NMCEvents=kr->getKeyI("NMCEvents");
  MCMCFile=kr->getKeyS("MCMCFile");
  preProcFCCut=kr->getKeyI("preProcFCCut");;
  preProcEVisCut=kr->getKeyD("preProcEVisCut");
  preProcWallMinCut=kr->getKeyD("preProcWallMinCut");
  preProcToWallMinCut=kr->getKeyD("preProcToWallMinCut");
  preProcNseMax0=kr->getKeyI("preProcNseMax");
  preProcNseMin=kr->getKeyI("preProcNseMin");
  preProcInGateCut=kr->getKeyD("preProcInGateCut");
  NDataEvents = kr->getKeyI("NDataEvents");
  flgUseNormPars = kr->getKeyI("flgUseNormPars");
  flgFixAllSmearPars = kr->getKeyI("flgFixAllSmearPars");
}

sharedPars::sharedPars(const char* parfilename){
  parFileName = parfilename;
//  readParsFromFile();
}


#endif
