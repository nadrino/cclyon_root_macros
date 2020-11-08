#ifndef SELECTORS_H
#define SELECTORS_H

#include "TMath.h"
#include <iostream>


#define USELOGRC

using namespace std;

struct fqcutparams{

  float fqmommu;
  float fqnring;
  float fqmome;
  float fqpid;
  float fqpippar;
  float fqpi0par;
  float fqenumu;
  float fqenue;
  float fqrcpar;
  int   nhitac;
  int   fqnsubev;

};



float calcEnu(int ipid, float pmom, float dirx, float diry, float dirz){

  // masses
  float Mp = 938.3;
  float Mn = 939.5;
  float Me = 0.511;
  float Mmu = 105.65;
  float Ml = 0.;
  // binding
  float V = 27.;
  float Enu = 0.;

  // lepton mass
  if (ipid==1) Ml = Me;
  else if (ipid==2) Ml = Mmu;
  else{ return 0.;}  

  // lepton momenum
  float Pl = pmom;

  // lepton energy
  float El = TMath::Sqrt( (Pl*Pl) + (Ml*Ml) );

  // cosine to beam
  float costh = (0.669764*dirx) + (-0.742179*diry) + (0.024228*dirz);

  Enu = (Mn -V)*El - ((Ml*Ml)/2.);
  Enu += Mn*V - (V*V/2.);
  Enu += ( Mp*Mp - Mn*Mn)/2.;
  Enu /= (Mn-V-El+(Pl*costh));

  //
  return Enu;

}

////////////////////////////////////
// cut using struct
int selectNuMu(fqcutparams fq){


  // nu energy cut
  if (fq.fqenumu>30000.) return 0.;
  if (fq.fqenumu<0.) return 0.;

  // subevent cut
  if (fq.fqnsubev>2) return 0;

  // FC cut
  if (fq.nhitac>16){
    return 0;
  }

  // Evis cut
  if (fq.fqmome<30){
    return 0;
  }

  // muon momentum cut
  if (fq.fqmommu<200.){
    return 0;
  }

  // e/mu PID cut
  if (fq.fqpid>0.2*fq.fqmome){
    return 0;
  }
  
  // mu/pip PID cut
  if (fq.fqpippar>0.){
    return 0;
  }

  // ring-counting cut
//  if (fq.fqnring!=1){
//    return 0;
//  }

  #ifndef USELOGRC
  if (fq.fqnring!=1){
    return 0;
  }
  #else
  if (fq.fqrcpar>0.){
    return 0;
  }
  #endif 

//#ifdef VERBOSE
//  cout<<"PASSED NU MU CUTS"<<endl;
//#endif

  return 1;

}


/*
///////////////////////////////////
int selectNuMu( int nhitac,
                int nsubev,
                float enu,
                float emomentum,
                float mumomentum,
                float fqpid,
                int   fqrc){

  // nu energy cut
  if (enu>30000.) return 0.;
  if (enu<0.) return 0.;

  // subevent cut
  if (nsubev>2) return 0;

  // FC cut
  if (nhitac>16){
    return 0;
  }

  // Evis cut
  if (emomentum<30){
    return 0;
  }

  // muon momentum cut
  if (mumomentum<100.){
    return 0;
  }

  // e/mu PID cut
  if (fqpid>0.2*emomentum){
    return 0;
  }
   
  // ring-counting cut
  if (fqrc!=1){
    return 0;
  }

  return 1;

}
*/              





//////////////////////////////////////////
int selectNuE(fqcutparams fq){

  // neutrino energy cut
  if (fq.fqenue>1250.) return 0.;
  if (fq.fqenue<0.) return 0.;

  // subevent cut
  if (fq.fqnsubev!=1) return 0;

  // FC cut
  if (fq.nhitac>16) return 0;

  // Evis cut
  if (fq.fqmome<100) return 0;

  // e/mu pid cut
  if (fq.fqpid<0.2*fq.fqmome) return 0;

  // RC cut
  #ifndef USELOGRC
  if (fq.fqnring!=1){
    return 0;
  }
  #else
  if (fq.fqrcpar>0.){
    return 0;
  }
  #endif 

  // pi0 cut
  if (fq.fqpi0par>0) return 0;

  
  return 1;

}


//////////////////////////////////////////
int selectNuE1Rpi(fqcutparams fq){

  // neutrino energy cut
  if (fq.fqenue>1250.) return 0.;
  if (fq.fqenue<0.) return 0.;

  // subevent cut
  if (fq.fqnsubev!=2) return 0;

  // FC cut
  if (fq.nhitac>16) return 0;

  // Evis cut
  if (fq.fqmome<100) return 0;

  // e/mu pid cut
  if (fq.fqpid<0.2*fq.fqmome) return 0;

  // RC cut
  #ifndef USELOGRC
  if (fq.fqnring!=1){
    return 0;
  }
  #else
  if (fq.fqrcpar>0.){
    return 0;
  }
  #endif 

  // pi0 cut
  if (fq.fqpi0par>0) return 0;
  
  return 1;

}



///////////////////////////////////////////////
void printCutPars(fqcutparams fq){
  
 cout<<"-----------------------"<<endl;
 cout<<"-- fq cut parameters --"<<endl;
 cout<<"-----------------------"<<endl;
 cout<<" muon momentum: "<<fq.fqmommu<<endl;
 cout<<" electron momentum: "<<fq.fqmome<<endl;
 cout<<" PID param: "<<fq.fqpid<<endl;
 cout<<" Pi0 param: "<<fq.fqpi0par<<endl;
 cout<<" PiP param: "<<fq.fqpippar<<endl;
 cout<<" RC param: "<<fq.fqrcpar<<endl;
 cout<<" nhitac: "<<fq.nhitac<<endl;
 cout<<" Enu electron: "<<fq.fqenue<<endl;
 cout<<" Enu muon: "<<fq.fqenumu<<endl;
 cout<<" Nsubev: "<<fq.fqnsubev<<endl;
 cout<<" Nring:  "<<fq.fqnring<<endl;
 if (selectNuE(fq)){
   cout<<" PASSED NUE SELECTION"<<endl;
 }
 if (selectNuMu(fq)){
   cout<<" PASSED NUMU SELECTION"<<endl;
 }

 return;
}

/*
int selectNuE(  int nhitac,
                int nsubev,
                float enu,
                float emomentum,
                float fqpid,
                int   fqrc,
                float fqpi0){

  // neutrino energy cut
  if (enu>1200.) return 0.;
  if (enu<0.) return 0.;

  // subevent cut
  if (nsubev>1) return 0;

  // FC cut
  if (nhitac>16) return 0;

  // Evis cut
  if (emomentum<100) return 0;

  // e/mu pid cut
  if (fqpid<0.2*emomentum) return 0;

  // ring-counting cut
  if (fqrc!=1) return 0;

  // pi0 cut
  if (fqpi0>0) return 0;

  return 1;

}
*/

float calcpi0par(float fqpi0like, float fqpi0mass){
  float pi0par = fqpi0like- 70. - ((140.-70.)/(40.- 120.))*(fqpi0mass - 120.);
  return pi0par;
}






#endif




