

string __irods_pulled_path__ = "/sps/t2k/ablanche/work/results/irods/run2a";

void get_POT_from_highland_files(){

  DataSample* ds;

  auto file_path_list = TToolBox::get_list_of_files_in_subfolders(__irods_pulled_path__, ".root");
  std::sort(file_path_list.begin(), file_path_list.end());

  for(int i_file = 0 ; i_file < int(file_path_list.size()); i_file++){
    if(TToolBox::do_string_starts_with_substring(file_path_list[i_file], "NumuCCMultiPiAnalysis")){

      ds = new DataSample((__irods_pulled_path__ + "/" + file_path_list[i_file]).c_str(),kGoodBeamGoodND280)
      cerr << ds.GetPOT() << endl;

    }
  }

  exit(EXIT_SUCCESS);

}
