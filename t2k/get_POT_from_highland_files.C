

std::string __irods_pulled_path__ = "/sps/t2k/ablanche/work/results/irods/run*";

double getCumulatedPOT(std::vector<std::string> filesList_);

void get_POT_from_highland_files(){

  std::cout << "Computing FHC Accumulated POT..." << std::endl;
  std::string ls_command = "ls " + __irods_pulled_path__ + "/NumuCCMultiPiAnalysis* &> " + __irods_pulled_path__ + "/temp.txt";
  gSystem->Exec(ls_command.c_str());
  double fhcPOT = getCumulatedPOT(GenericToolBox::dumpFileAsVectorString(__irods_pulled_path__ + "/temp.txt"));
  std::cout << "FHC Accumulated POT: " << fhcPOT << std::endl;


  std::cout << "Computing FHC Accumulated POT..." << std::endl;
  std::string ls_command = "ls " + __irods_pulled_path__ + "/NumuCCMultiPiAnalysis* &> " + __irods_pulled_path__ + "/temp.txt";
  gSystem->Exec(ls_command.c_str());
  double rhcPOT = getCumulatedPOT(GenericToolBox::dumpFileAsVectorString(__irods_pulled_path__ + "/temp.txt"));
  std::cout << "RHC Accumulated POT: " << rhcPOT << std::endl;

  exit(0);

}

double getCumulatedPOT(std::vector<std::string> filesList_){

  DataSample* ds;

  double cumulated_pot = 0;

  for(int i_file = 0 ; i_file < int(file_path_list.size()); i_file++){
    ds = new DataSample( (file_path_list[i_file]).c_str(), kGoodBeamGoodND280);
    GenericToolbox::displayProgressBar(i_file, int(file_path_list.size()), "Accumulating POT...");
    cumulated_pot += ds->GetPOT();
  }

  return cumulated_pot;

}
