
std::vector<std::string> get_list_of_files_in_subfolders(std::string folder_path_, std::string files_extension_ = ""){
  std::string *folder_path = &folder_path_;
  std::string *files_extension = nullptr;
  if(not files_extension_.empty()) files_extension = &files_extension_;
  return get_list_of_files_in_subfolders(folder_path, files_extension);
}


string __irods_pulled_path__ = "/sps/t2k/ablanche/work/results/irods/run2a";

void get_POT_from_highland_files(){

  DataSample* ds;

  std::vector<std::string> file_path_list = TToolBox::get_list_of_files_in_subfolders(__irods_pulled_path__, ".root");
  std::sort(file_path_list.begin(), file_path_list.end());

  for(int i_file = 0 ; i_file < int(file_path_list.size()); i_file++){
    if(TToolBox::do_string_starts_with_substring(file_path_list[i_file], "NumuCCMultiPiAnalysis")){

      ds = new DataSample((__irods_pulled_path__ + "/" + file_path_list[i_file]).c_str(),kGoodBeamGoodND280)
      cerr << ds.GetPOT() << endl;

    }
  }

  exit(EXIT_SUCCESS);

}
