

ENUM_EXPANDER(Test, 1, test1, test2, test3)


void test(){

  vector<double> FHCNumuCCOth_Mom_Bin = {0. , 600. , 1000., 1250., 2000., 4000., 30000.};
  vector<double> FHCNumuCCOth_Cos_Bin = {-1.0, 0.7, 0.8, 0.86, 0.9 , 0.93, 0.95, 0.97, 0.99, 1.0};

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
