#ifndef PREPROCESS_H
#define PREPROCESS_H

#include <iostream>
#include "visRing.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TString.h"
#include "FVCalculators.h"
#include "sharedPars.h"
#include "TGraph.h"
#include "eventSelectors.h"
#include <map>
#include <string>
#include "TH2FV.h"
#include "TObjArray.h"
#include "masktools.h"
//#include "calcEnu.h"
#include "shared.h"

#define NFILEMAX 5000
//#define VERBOSE

using namespace std;


const float SK4AtmLivetime = 2166.5;


//////////////////////////////////////////////////////////////
// Class to take care of preprocessing of all data and MC files
// Usage:
//	1) create instance using preProcess()
//	2) specify parameter file using setParFileName(<name>)
//  3) run using runPreProcessing()
class preProcess{
  public:

  /////////////////////
  // constructors
  preProcess();


  /////////////////////
  // internal variables
  TChain* chmc; //< points to input MC
  TChain* chdat; //< points to input data
  TTree* tr;
  TTree* trout;
  TTree* rfgWeightTree;
  TChain* rfgWeightChain;
  TFile* fout;
  TFile*  rfgFile;
  TString nameTag;
  TString outDir;
  TString fileNames[NFILEMAX];
  TString parFileName;
  TString ntupleType;
  fqEvent* fq;
  visRing*  vis;
  int MCComponents;
  int FVBinning;
  int nFiles;
  int MCSamples; //< flag for sample definitions
  TH1D* hWeight;
  TGraph* gWeight;
  int useWeights;
  map<string,double> attributeMap;
  TString attributeList[50];
  int nAttributes;
  TH2FV* hFVBins;
  TH1D* hmask;
  int fakeShiftFlg; //< puts a fake 50 pt shift in PID likelihood for sample 0
  int fakeNormFlg;  //< puts a fake 15% normalization shift for sub GeV events
  float atmMCNorm;
  float rfgweight;
  int flgUseRFGWgt;

  ////////////////////////
  //for cuts
  int NHITACMax;
  double EVisMin;
  double WallMin;
  double ToWallMin;
  int NSEMax;
  int NSEMin;
  double InGateMin; 

  //////////////////////
  //new branch variables
  float towall;
  float wall;
  float towallv[50];
  float minconev[50];
  float perimv[50];
  float wallv2;
  float fq1rwall[10][7];
  float fq1rtowall[10][7];
  float fq1rmincone[10][7];
  float fq1rperim[10][7];
  float fq1renu[10][7];
  float evtweight;
  float attribute[1000];
  float fqrcpar;
  float fqmrdot;
  float fqrcpmin;
  float fqrclike;
  float fqpi0par;
  float fqpippar;
  int ncomponent;
  int nsample;
  int nbin;
  int best2RID;
  int passmucut;
  int passecut;
  int passe1rpicut;
  
  /////////////////////
  //methods
  void setTree(TTree* trin);
  void setTree(TChain* trin);
  void setParFileName(const char* filename);
  void runPreProcessing();
  void setFVBinHisto();
  void setupNewTree();
  int preProcessIt(); //< process a file and return the number of entries in new tree
  int passCuts();
  void fillAttributes(fqEvent* fqevent);
  void fillAttributeMap(fqEvent* fqevent);
  int absmode;
  int getComponent();
  int getSample();
  int getBin();
  int getBest2RFitID(fqEvent* fqevent); //< returns best 2R fit ID for current event
  float getWeight();
  void processFile(const char* fname,const char* outname="");
  void processAllFiles(TChain* chain);
  void setRFGDir(const char* rfgdir); //< sets RFG weight directory and turns on weights 
  float getAtmWeight(); //< get the atmospheric MC weights
  //sets a histogram to calculate weights for events (use to correct cosmic
  //muon momenum distribution)
  void setWeightHistogram(const char* file, const char * name);
  void makeTestFiles(const char* outdir, int testtype, int nmc, int ndata, int randseed); //< makes test files
  void calcNeutrinoEnergy();
  TH2D* hRCPar;
  float getRCParameter(fqEvent* fqevent); //< calculate the RC parameter for this event
  float getPi0Parameter(fqEvent* fqevent); //< calculate the pi0 parameter for this event
  float getPiPParameter(fqEvent* fqevent); //< calculate the pi0 parameter for this event
//  float getRFGWeight();
  TString getRFGFileName(const char* rootname);
  void applyT2KSelection();

  private:

  int flgAddMoreVars;
  int flgUseSpikeMask;
};

#endif

#ifdef CINTMODE
#include "globalRandom.cxx"
#endif
