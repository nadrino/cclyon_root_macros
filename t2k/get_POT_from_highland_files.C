


std::vector<std::string> read_file(std::string file_path_);

std::string __irods_pulled_path__ = "/sps/t2k/ablanche/work/results/irods/run2a";

void get_POT_from_highland_files(){

  std::string ls_command = "ls " + __irods_pulled_path__ + "/NumuCCMultiPiAnalysis* &> " + __irods_pulled_path__ + "/temp.txt";
  gSystem->Exec(ls_command.c_str());
  std::vector<std::string> file_path_list = read_file(__irods_pulled_path__ + "/temp.txt");

  DataSample* ds;

  std::sort(file_path_list.begin(), file_path_list.end());

  for(int i_file = 0 ; i_file < int(file_path_list.size()); i_file++){
    ds = new DataSample((__irods_pulled_path__ + "/" + file_path_list[i_file]).c_str(),kGoodBeamGoodND280)
    cerr << ds.GetPOT() << endl;
  }

  exit(EXIT_SUCCESS);

}



std::vector<std::string> read_file(std::string file_path_){

  std::vector<std::string> res;
  std::string line;
  std::ifstream infile(file_path_.c_str());
  while (std::getline(infile, line))
  {
    res.push_back(line);
  }
  infile.close();

  return res;

}
