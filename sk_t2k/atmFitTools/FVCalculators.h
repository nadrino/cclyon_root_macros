#ifndef FVCALC_H
#define FVCALC_H

#include "TPolyLine3D.h"
#include <TVector3.h>
#include <TMath.h>
#include <math.h>
#include <iostream>
#include "TRandom2.h"
#include "TCanvas.h"

using namespace std;

double calcWall(TVector3* vpos);
double calcWall2(TVector3* vpos);
double calcToWallCustom(TVector3 *vpos, TVector3 *vdir, double dt);
double calcMinCone(TVector3 *vpos, TVector3* vdir, int npts=50);
double calcPerimeter(TVector3 *vpos, TVector3* vdir, int npts=50, int visflg=0);
double calcToWall(TVector3* vpostmp, TVector3* vdirtmp);
double calcPhiWall(TVector3* vpos, TVector3* vdir);

#ifdef CINTMODE
#include "FVCalculators.cxx"
#endif
#endif
