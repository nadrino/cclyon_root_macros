#ifndef DEFINE_H
#define DEFINE_H

#include "TH1D.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TRandom2.h"

//#define T2K // should be uncommented if using T2K xsec parametrization

// uncomment this line if not using xiaoyue's skimmed atmospheric MC
//#define USE_XL_WEIGHTS

// uncomment to use shimpei's MC weights
#define USE_ST_WEIGHTS

// for global random
extern TRandom2* randy;

// uncomment to compile in interactive mode
#define CINTMODE


#define NSAMPMAX 5 // total number of samples (# of sub-events)
#define NBINMAX 10  // total number of FV bins
#define NCOMPMAX 8 // number of MC components based on true info
#define NATTMAX 6  // number of attributes (fQ reconstructed variables)
#define FLGDEBUG 0 // set to 1 to print out some useful things
#ifndef T2K
#define NSYSPARMAX 40 // number of flux and xsec pars from tn186
#else
#define NSYSPARMAX 500
#endif
#define NHBINSMAX 300 // 
#define NPTSMAX 21
#define NMODE 10
#define SCALING 0.1 // 0.1 for 10 years data / 100 years MC 
#define NMCMCPARS 500 // for markovTools

#endif
