#ifndef MASKTOOLS_H
#define MASKTOOLS_H

#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include <iostream>
#include "FVCalculators.h"
#ifndef T2K
#include "fqEvent.h"
#else
#include "t2kfqEvent.h"
#endif

using namespace std;


int passMask(TH1D* hmask, double value);


// a class to mask out spikes in hybrid pi0 data samples
class masktools{
  public:

  // constructor
  masktools();

  // chain with all files
  TChain* chdata;

  // number of bins in histogram
  int nbins;

  // threshold for masking
  double thresh;

  // histgrams
  TH1D* hwall;
  TH1D* hmask;

  // make a mask and save it
  void makethismask(const char* filename);

};

#ifdef CINTMODE
#include "masktools.cxx"
#endif

#endif
