

std::string __irods_pulled_path__ = "/sps/t2k/ablanche/work/results/irods";

double getCumulatedPOT(std::vector<std::string> filesList_);

void get_POT_from_highland_files(){

  if(gROOT->GetVersionInt() > 60000){
    cout << "This script has to be ran on ROOT5 with highland libraries linked." << endl;
    exit(1);
  }

  std::stringstream lsCommand;
  lsCommand << "ls -d " << __irods_pulled_path__ << "/run*/  &> " + __irods_pulled_path__ + "/temp.txt";;
  gSystem->Exec(lsCommand.str().c_str());
  vector<string> runFolders = TToolBox::read_file(__irods_pulled_path__ + "/temp.txt");

  vector<double> mapRunFHCPOT;

  std::cout << "Computing FHC Accumulated POT..." << std::endl;
  double fhcPOT = 0;
  string runFolder;
  for( int iRun = 0 ; iRun < runFolders.size() ; iRun++ ){
    runFolder = runFolders[iRun];
    vector<string> pathSlices = TToolBox::split_string(runFolder, "/");
    std::string runName = pathSlices[pathSlices.size()-2];

    std::stringstream lsCommand;
    lsCommand.str("");
    lsCommand << "ls " << runFolder;
    lsCommand << "/NumuCCMultiPiAnalysis* &> " + __irods_pulled_path__ + "/temp.txt";
    gSystem->Exec(lsCommand.str().c_str());

    double fhcRunPOT = getCumulatedPOT(TToolBox::read_file(__irods_pulled_path__ + "/temp.txt"));
    std::cout << "  FHC Accumulated POT " << runName << ": " << fhcRunPOT << std::endl;
    fhcPOT += fhcRunPOT;
    mapRunFHCPOT.push_back(fhcRunPOT);
  }
  std::cout << "FHC Accumulated POT: " << fhcPOT << std::endl;


  std::cout << "Computing RHC Accumulated POT..." << std::endl;
  double rhcPOT = 0;
  vector<double> mapRunRHCPOT;
  string runFolder;
  for( int iRun = 0 ; iRun < runFolders.size() ; iRun++ ){
    runFolder = runFolders[iRun];
    vector<string> pathSlices = TToolBox::split_string(runFolder, "/");
    std::string runName = pathSlices[pathSlices.size()-2];
    std::stringstream lsCommand;
    lsCommand.str("");
    lsCommand << "ls " << runFolder;
    lsCommand << "/AntiNumuCCMultiPiAnalysis* &> " + __irods_pulled_path__ + "/temp.txt";
    gSystem->Exec(lsCommand.str().c_str());
    double rhcRunPOT = getCumulatedPOT(TToolBox::read_file(__irods_pulled_path__ + "/temp.txt"));
    std::cout << "  RHC Accumulated POT " << runName << ": " << rhcRunPOT << std::endl;
    rhcPOT += rhcRunPOT;
    mapRunRHCPOT.push_back(rhcRunPOT);
  }
  std::cout << "RHC Accumulated POT: " << rhcPOT << std::endl;

  // Latex Table:
  cout << "\\begin{table}[hbt]" << endl;
  cout << "\\begin{center}" << endl;
  cout << "\\footnotesize" << endl;
  cout << "\\begin{tabular}{|c|c|c|}" << endl;
  cout << "\\hline" << endl;
  cout << "\\hline" << endl;
  cout << "Run & FHC POT & RHC POT \\\\" << endl;

  for( int iRun = 0 ; iRun < runFolders.size() ; iRun++ ){

    runFolder = runFolders[iRun];
    vector<string> pathSlices = TToolBox::split_string(runFolder, "/");
    std::string runName = pathSlices[pathSlices.size()-2];

    cout << runName << " & ";
    cout << mapRunFHCPOT[iRun] << " & ";
    cout << mapRunRHCPOT[iRun] << " \\\\";
    cout << endl;

    cout << "\\hline" << endl;

  }

  cout << "\\hline" << endl;
  cout << "\\end{tabular}" << endl;
  cout << "\\end{center}" << endl;
  cout << "\\end{table}" << endl;

  exit(0);

}

double getCumulatedPOT(std::vector<std::string> filesList_){

  DataSample* ds;

  double cumulated_pot = 0;

  for(int i_file = 0 ; i_file < int(filesList_.size()); i_file++){
    ds = new DataSample( (filesList_[i_file]).c_str(), kGoodBeamGoodND280);
    TToolBox::display_loading(i_file, int(filesList_.size()), "Accumulating POT...");
    cumulated_pot += ds->GetPOT();
    delete ds;
  }

  return cumulated_pot;

}
