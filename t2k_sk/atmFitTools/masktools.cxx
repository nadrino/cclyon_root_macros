#ifndef MASKTOOLS_C
#define MASKTOOLS_C

#include "masktools.h"


using namespace std;

int passMask(TH1D* hmask, double value){
  int ibin = hmask->FindBin(value);
  return (int)hmask->GetBinContent(ibin);
}

void masktools::makethismask(const char* filename){
  
  // make histograms
  hwall = new TH1D("hwall","hwall",nbins,0,1000);
  hmask = (TH1D*)hwall->Clone("hmask");

  // fill histogram
#ifndef T2K
  fqEvent* fqevent = new fqEvent(chdata);
#else
  t2kfqEvent *fqevent = new t2kfqEvent(chdata);
#endif
  int nev = chdata->GetEntries();
  TVector3* vpos = new TVector3();
  for (int iev=0; iev<nev; iev++){
    chdata->GetEntry(iev);
    vpos->SetXYZ(fqevent->fq1rpos[0][2][0],
                 fqevent->fq1rpos[0][2][1],
                 fqevent->fq1rpos[0][2][2]);
    double wall = calcWall(vpos);
    hwall->Fill(wall);
  }

  // make mask 
  for (int ibin=1; ibin<=hwall->GetNbinsX(); ibin++){
     double content = hwall->GetBinContent(ibin);
     if (content<thresh) hmask->SetBinContent(ibin,1.);
     else{
       hmask->SetBinContent(ibin,0.);
     }
  }

  // save mask
  hmask->SaveAs(filename);

  return;
}

masktools::masktools(){
  
}

#endif
