


void test(){

  vector<double> FHCNumuCCOth_Mom_Bin = {0. , 500. , 700. , 1000., 1250., 1500., 2000., 3000., 10000.};
  vector<double> FHCNumuCCOth_Cos_Bin = {-1.0 , 0.7 , 0.8 , 0.85, 0.90, 0.93, 0.95, 0.96,
                                         0.97, 0.98, 0.99, 1.0};

  for(int iCosBin = 0 ; iCosBin < FHCNumuCCOth_Cos_Bin.size()-1 ; iCosBin++){
    for( int iMomBin = 0 ; iMomBin < FHCNumuCCOth_Mom_Bin.size()-1 ; iMomBin++){
      cout << FHCNumuCCOth_Cos_Bin[iCosBin] << " ";
      cout << FHCNumuCCOth_Cos_Bin[iCosBin+1]  << " ";
      cout << FHCNumuCCOth_Mom_Bin[iMomBin]  << " ";
      cout << FHCNumuCCOth_Mom_Bin[iMomBin+1]  << " ";
      cout << endl;
    }
  }

}
