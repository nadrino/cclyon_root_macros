#ifndef SHAREDPARS_H
#define SHAREDPARS_H

#include "shared.h"
#include "TString.h"
#include "keyread.h"
#include <iostream>


//class to hold and read in shared parameters for fits
class sharedPars{
  public:

  //constructor
  sharedPars(const char* parfilename);

  //name of parameter file
  TString parFileName;
  
  //read in values from parameter file
  void readParsFromFile();

  //shared variables
  TString globalRootName;

  //can get parameters directly by name
  int getParI(const char*);
  double getParD(const char*);
  TString getParS(const char*); 

  //key reading object
  keyread* kr;

  //fit parameters
  int nFVBins;
  int nSamples;
  int nComponents;
  int nAttributes;
  int nSysPars;
  int nModes;
  double normFactor;
  TString preProcessFilesMC;
  TString preProcessFilesSpline;
  TString preProcessFilesBANFF;
  TString preProcessFilesData;
  TString preProcessOutDir;
  int preProcessMCComponents;
  int preProcessMCSamples;
  int preProcessFVBinning; 
  int preProcFCCut;
  double preProcEVisCut;
  double preProcWallMinCut;
  double preProcToWallMinCut;
  int    preProcNseMax0;
  int    preProcNseMin;
  int    preProcMaskFlg;
  TString preProcMaskFile;
  TString ntupleType;
  double  preProcInGateCut;
  int     preProcAddMoreVars;
  TString FVBinName0;
  TString FVBinName1;
  TString FVBinName2;
  TString fQAttName0;
  TString fQAttName1;
  TString fQAttName2;
  TString fQAttName3;
  TString fQAttName4;
  TString fQAttName5;
  TString fQAttName6;
  TString fQAttName7;
  TString fQAttName8;
  TString MCComponentName0;
  TString MCComponentName1;
  TString MCComponentName2;
  TString MCComponentName3;
  TString MCComponentName4;
  TString MCComponentName5;
  TString MCComponentName6;
  TString sampleName0;
  TString sampleName1;
  TString sampleName2;
  TString sysParName0;
  TString sysParName1;
  TString sysParName2;
  TString sysParName3;
  TString sysParName4;
  TString sysParName5;
  TString sysParName6;
  TString sysParName7;
  TString sysParName8;
  TString hFactoryOutput;
  TString hFactoryMCFiles;
  TString hFactoryDataFiles;
  TString splineFactoryOutput;
  int MCMCNSteps;
  int MCMCNBurnIn;
  double MCMCTunePar;
  int useSplinesFlg;
  int fixAllSmearFlg;
  int NMCMCPts;
  int MCMCBurnIn;
  int NMCEvents;
  TString MCMCFile;
  TString sysParType;
  int NDataEvents;
  int flgUseNormPars;
  int flgFixAllSmearPars;
  int flgUsePriorsInFit;
  double smearPriorWidthAtt0;
  double smearPriorWidthAtt1;
  double smearPriorWidthAtt2;
  double smearPriorWidthAtt3;
  double smearPriorWidthAtt4;
  double smearPriorWidthAtt5;
  double biasPriorWidthAtt0;
  double biasPriorWidthAtt1;
  double biasPriorWidthAtt2;
  double biasPriorWidthAtt3;
  double biasPriorWidthAtt4;
  double biasPriorWidthAtt5;
  int flgUsePhysLoBound;
  double physLoBoundAtt0;
  double physLoBoundAtt1;
  double physLoBoundAtt2;
  double physLoBoundAtt3;
  double physLoBoundAtt4;
  int flgUseFitParFile;
  int nBinBuffer;
//  double physBoundAtt0;
//  double physBoundAtt1;
//  double physBoundAtt2;
//  double physBoundAtt3;
 
};


#ifdef CINTMODE
#include "sharedPars.cxx"
#endif

#endif
