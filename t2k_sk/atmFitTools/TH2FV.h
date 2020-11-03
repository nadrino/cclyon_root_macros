#ifndef  TH2FV_H
#define  TH2FV_H
#include "shared.h"
#include "TH2Poly.h"
#include "TVector2.h"
#include "TMath.h"
#include <iostream>

#define NPOLYBINSMAX 9999


using namespace std;



////////////////////////////////////////////////////
// Class to make FV 2D histograms
//    X axis is towall
//    Y axis is wall
class TH2FV:public TH2Poly{

  public:

  TH2FV(const char* name, int bintype=0, int nbinsx =10, double xmin=0, double xmax=5000,
                                       int nbinsy =10, double ymin=0, double ymax=2000.);
//  TH2FV(const char* name, int bintype, int nbinsx, double xmin, double xmax,
                                       //int nbinsy, double ymin, double ymax);





  //vars
  int fBinType;
  int fNBinsX;
  int fNBinsY;
  double fXMin;
  double fXMax;
  double fYMin;
  double fYMax;
  double fCenter[NPOLYBINSMAX][2];
  int fNCells;
  int fAdded;

  //initialize
  void Init();
  
  //add bin and keep track of center
  //(messes up AddBin() for some reason
  void AddBinWithCenter(int n, double *xx, double* yy);

  //add bin center
  void AddBinCenter(int n, double *xx, double* yy);
  void AddSquareDiagonalBins(int nbins, double maxval);

  // set min to be smallest non-zero value
  void SetMinNonZero();

  // get the center of a bin
  double GetBinCenterX(int nbin);
  double GetBinCenterX(int binx, int biny);
  double GetBinCenterY(int nbin);
  double GetBinCenterY(int binx, int biny);
  double GetMaxWall();
  int GetNcells() { return fNCells; }
  // draw standard view
  void DrawStdView(const char* opts);
 
  // test bin consistency
  void TestHisto();

  double fMaxWall;

  protected:

  //for triangle plots:
//  void AddSquareDiagonalBins(int nbins, double maxval);
  int LineIntersects(double* binx, double* biny);
  //double fMaxWall;
  void SetSplitBins(double* splitx, double* splity);
  void InitSplitDiagonal();
  void InitFVBins();
  void InitFVBins2(); 
  void InitFVBins3();
  void InitFVBins4();
  void InitFVBins5();
  void InitFVBins6();
  void InitFVBins7();
  void InitStdBins(double wall1, double wall2, double towall1,
                   double towall2, double towall3, double towall4);
  //ClassDef(TH2FV,1);

};

#ifdef CINTMODE
#include "TH2FV.cxx"
#endif
#endif
