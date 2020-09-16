

std::string __irods_pulled_path__ = "/sps/t2k/ablanche/work/results/irods";

double getCumulatedPOT(std::vector<std::string> filesList_);

void get_POT_from_highland_files(){

  std::stringstream lsCommand;
  lsCommand << "ls -d " << __irods_pulled_path__ << "/run*/  &> " + __irods_pulled_path__ + "/temp.txt";;
  gSystem->Exec(lsCommand.str().c_str());
  vector<string> runFolders = TToolBox::read_file(__irods_pulled_path__ + "/temp.txt");

  std::cout << "Computing FHC Accumulated POT..." << std::endl;
  double fhcPOT = 0;
  string runFolder;
  for( int iRun = 0 ; iRun < runFolders.size() ; iRun++ ){
    runFolder = runFolders[iRun];
    std::stringstream lsCommand;
    lsCommand.str("");
    lsCommand << "ls " << __irods_pulled_path__ << "/" << runFolder;
    lsCommand << "/NumuCCMultiPiAnalysis* &> " + __irods_pulled_path__ + "/temp.txt";
    gSystem->Exec(lsCommand.str().c_str());
    double fhcRunPOT = getCumulatedPOT(TToolBox::read_file(__irods_pulled_path__ + "/temp.txt"));
    std::cout << "  FHC Accumulated POT " << runFolder << ": " << fhcRunPOT << std::endl;
    fhcPOT += fhcRunPOT;
  }
  std::cout << "FHC Accumulated POT: " << fhcPOT << std::endl;


  std::cout << "Computing RHC Accumulated POT..." << std::endl;
  double rhcPOT = 0;
  string runFolder;
  for( int iRun = 0 ; iRun < runFolders.size() ; iRun++ ){
    runFolder = runFolders[iRun];
    std::stringstream lsCommand;
    lsCommand.str("");
    lsCommand << "ls " << __irods_pulled_path__ << "/" << runFolder;
    lsCommand << "/AntiNumuCCMultiPiAnalysis* &> " + __irods_pulled_path__ + "/temp.txt";
    gSystem->Exec(lsCommand.str().c_str());
    double rhcRunPOT = getCumulatedPOT(TToolBox::read_file(__irods_pulled_path__ + "/temp.txt"));
    std::cout << "  RHC Accumulated POT " << runFolder << ": " << rhcRunPOT << std::endl;
    rhcPOT += rhcRunPOT;
  }
  std::cout << "RHC Accumulated POT: " << rhcPOT << std::endl;

  exit(0);

}

double getCumulatedPOT(std::vector<std::string> filesList_){

  DataSample* ds;

  double cumulated_pot = 0;

  for(int i_file = 0 ; i_file < int(filesList_.size()); i_file++){
    ds = new DataSample( (filesList_[i_file]).c_str(), kGoodBeamGoodND280);
    TToolBox::display_loading(i_file, int(filesList_.size()), "Accumulating POT...");
    cumulated_pot += ds->GetPOT();
  }

  return cumulated_pot;

}
