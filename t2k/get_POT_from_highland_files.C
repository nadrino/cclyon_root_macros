

void display_loading(int current_index_, int end_index_, string title_ = "", bool force_display_ = false);
std::vector<std::string> read_file(std::string file_path_);

std::string __irods_pulled_path__ = "/sps/t2k/ablanche/work/results/irods/run2a";

void get_POT_from_highland_files(){

  std::string ls_command = "ls " + __irods_pulled_path__ + "/NumuCCMultiPiAnalysis* &> " + __irods_pulled_path__ + "/temp.txt";
  gSystem->Exec(ls_command.c_str());
  std::vector<std::string> file_path_list = read_file(__irods_pulled_path__ + "/temp.txt");

  DataSample* ds;
  // std::sort(file_path_list.begin(), file_path_list.end());

  double cumulated_pot = 0;

  for(int i_file = 0 ; i_file < int(file_path_list.size()); i_file++){
    ds = new DataSample( (file_path_list[i_file]).c_str(), kGoodBeamGoodND280);
    display_loading(i_file, int(file_path_list.size()), "Accumulating POT");
    cumulated_pot += ds.GetPOT();
  }

  cerr << "Accumulated POT : " << cumulated_pot << endl;

  return 0;

}

void display_loading(int current_index_, int end_index_, string title_, bool force_display_) {

  int percent = int(round(double(current_index_) / end_index_ * 100.));
  if(force_display_ or current_index_ >= end_index_-1) {
    if(last_displayed_value != -1) clog << "\r" << title_ << " : " << 100 << "%" << endl;
    reset_last_displayed_value();
    return;
  }
  if(last_displayed_value == -1 or last_displayed_value < percent) {
    last_displayed_value = percent;
    clog << "\r" << title_ << " : " << percent << "%" << flush << "\r";
  }

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
