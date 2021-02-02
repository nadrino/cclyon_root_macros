#include "TChain.h"

#include "t2ksk_eventSelections.hh"

void t2ksk_sample(){

  TChain * h1 = new TChain("h1");
  h1->Add("/Users/cvilela/T2K-SK/Data/RawMC/*.root");
  
  t2ksk::setTreeAddresses(h1);

  unsigned int nEntries = h1->GetEntries();
  
  for (unsigned int i = 0; i < nEntries; i++){
    h1->GetEntry(i);

    if (t2ksk::is1Re())         std::cout << "11 " << t2ksk::ComputeErec(0, t2ksk::ELECTRON) << std::endl;
    else if (t2ksk::is1Rmu())   std::cout << "13 " << t2ksk::ComputeErec(0, t2ksk::MUON) << std::endl; 
    else if (t2ksk::is1Re1de()) std::cout << "14 " << t2ksk::ComputeErecCCDel(0, t2ksk::ELECTRON) << std::endl;; 
  }
}
