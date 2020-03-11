

std::vector<std::string> get_list_of_entries_in_folder(std::string *folder_path_){

  if(not do_path_is_folder(*folder_path_)) return std::vector<std::string>();

  std::vector<std::string> entries_list;
  TSystemDirectory dir((*folder_path_).c_str(), (*folder_path_).c_str());
  TList *files = dir.GetListOfFiles();
  for(int i_entry = 0 ; i_entry < dir.GetListOfFiles()->GetEntries() ; i_entry++){
    string entry = dir.GetListOfFiles()->At(i_entry)->GetName();
    if(entry != "." and entry != ".."){
      entries_list.emplace_back(entry);
    }
  }
  return entries_list;

}
std::vector<std::string> get_list_of_subfolders_in_folder(std::string *folder_path_) {

  auto entries_list = get_list_of_entries_in_folder(folder_path_);
  std::vector<std::string> folders_list;
  for(int i_entry = 0 ; i_entry < int(entries_list.size()) ; i_entry++){
    if(do_path_is_folder(*folder_path_ + "/" + entries_list[i_entry]))
      folders_list.emplace_back(entries_list[i_entry]);
  }
  return folders_list;

}
std::vector<std::string> get_list_of_files_in_folder(std::string *folder_path_, std::string *files_extension_ = 0) {

  auto entries_list = get_list_of_entries_in_folder(folder_path_);
  std::vector<std::string> files_list;
  for(int i_entry = 0 ; i_entry < int(entries_list.size()) ; i_entry++){
    if(files_extension_ == 0 or do_string_ends_with_substring(entries_list[i_entry], *files_extension_)){
      if(do_path_is_file(*folder_path_ + "/" + entries_list[i_entry]))
        files_list.emplace_back(entries_list[i_entry]);
    }
  }
  return files_list;

}
std::vector<std::string> get_list_of_files_in_subfolders(std::string *folder_path_, std::string *files_extension_ = 0){

  std::vector<std::string> output_file_paths;

  auto files_list = get_list_of_files_in_folder(folder_path_, files_extension_);
  for(int i_file = 0 ; i_file < int(files_list.size()) ; i_file++){
    output_file_paths.emplace_back(files_list[i_file]);
  }

  auto subfolders_list = get_list_of_subfolders_in_folder(folder_path_);
  for(int i_subfolder = 0 ; i_subfolder < int(subfolders_list.size()) ; i_subfolder++){
    string subfolde_full_path = *folder_path_ + "/" + subfolders_list[i_subfolder];
    auto subfiles_path = get_list_of_files_in_subfolders(&subfolde_full_path, files_extension_); // RECURSIVE
    for(int i_subfile = 0 ; i_subfile < int(subfiles_path.size()) ; i_subfile++){
       output_file_paths.emplace_back(subfolders_list[i_subfolder] + "/" + subfiles_path[i_subfile]);
    }
  }

  return output_file_paths;

}

std::vector<std::string> get_list_of_files_in_subfolders(std::string folder_path_, std::string files_extension_ = ""){
  std::string *folder_path = &folder_path_;
  std::string *files_extension = 0;
  if(not files_extension_.empty()) files_extension = &files_extension_;
  return get_list_of_files_in_subfolders(folder_path, files_extension);
}









std::string __irods_pulled_path__ = "/sps/t2k/ablanche/work/results/irods/run2a";

void get_POT_from_highland_files(){

  DataSample* ds;

  std::vector<std::string> file_path_list = get_list_of_files_in_subfolders(__irods_pulled_path__, ".root");
  std::sort(file_path_list.begin(), file_path_list.end());

  for(int i_file = 0 ; i_file < int(file_path_list.size()); i_file++){
    if(do_string_starts_with_substring(file_path_list[i_file], "NumuCCMultiPiAnalysis")){

      ds = new DataSample((__irods_pulled_path__ + "/" + file_path_list[i_file]).c_str(),kGoodBeamGoodND280)
      cerr << ds.GetPOT() << endl;

    }
  }

  exit(EXIT_SUCCESS);

}
