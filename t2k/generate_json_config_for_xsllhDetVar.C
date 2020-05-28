

// User Parameters
std::string __input_files_list_path__ = "./files_list.txt";
// std::string __tree_converter_json_file_path__ = "./config_TreeConvert_numu_run8.json";
std::string __binning_file_path__ = "/sps/t2k/ablanche/work/results/xsLLhFitter/BANFF_Fit/binning/BANFF_binning.txt";
int __num_samples__ = 3;


void generate_json_config_for_xsllhDetVar(){

  cout << TToolBox::INFO << "Starting json config file generator for xsllhDetVar." << endl;
  cout << TToolBox::ALERT << " > Please check the following parameters:" << endl;
  cout << TToolBox::ALERT << "__input_files_list_path__ = " << __input_files_list_path__ << endl;
  cout << TToolBox::ALERT << "__num_samples__ = " << __num_samples__ << endl;
  cout << TToolBox::ALERT << "__binning_file_path__ = " << __binning_file_path__ << endl;
  // cout << TToolBox::ALERT << "__tree_converter_json_file_path__ = " << __tree_converter_json_file_path__ << endl;

  if(not TToolBox::do_path_is_file(__input_files_list_path__)){
    cout << TToolBox::ERROR << "__input_files_list_path__=" << __input_files_list_path__ << " was not found." << endl;
    exit(1);
  }
  auto input_highland_files_path_list = TToolBox::read_file(__input_files_list_path__);

  cout << TToolBox::WARNING << "Generating json config file for xsllhDetVar..." << endl;

  stringstream json_files_ss;
  TFile* temp_tfile;
  TTree* temp_ttree;
  for(auto &intput_file_path: input_highland_files_path_list){
    if(not json_files_ss.str().empty()){
      json_files_ss << "," << endl;
    }

    // get infos from the file
    temp_tfile = TFile::Open(intput_file_path.c_str(), "READ");
    temp_ttree = (TTree*) temp_tfile->Get("all_syst");
    temp_ttree->GetEntry(0);
    double num_toys = temp_ttree->GetLeaf("NTOYS")->GetValue(0);
    temp_ttree = (TTree*) temp_tfile->Get("config");
    double num_syst = temp_ttree->GetLeaf("NSYST")->GetValue(0);
    temp_tfile->Close();
    delete temp_tfile;

    json_files_ss << "    {" << endl;
    json_files_ss << "      \"fname_input\": \"" << intput_file_path << "\"," << endl;
    json_files_ss << "      \"tree_name\": \"all_syst\"," << endl;
    json_files_ss << "      \"detector\": \"ND280\"," << endl;
    json_files_ss << "      \"num_toys\": " << num_toys << "," << endl;
    json_files_ss << "      \"num_syst\": " << num_syst << "," << endl;
    json_files_ss << "      \"num_samples\": " << 3 << "," << endl;
    json_files_ss << "      \"cut_level\" : [7, 7, 6]," << endl;
    json_files_ss << "      \"samples\" : {\"0\": [0], \"1\": [1], \"2\": [2]}," << endl;
    json_files_ss << "      \"use\": true" << endl;
    json_files_ss << "    }";
  }
  json_files_ss << endl;

  stringstream json_ss;
  json_ss << "{" << endl;
  json_ss << "  \"fname_output\": \"detector_covariance_matrix.root\"," << endl;
  json_ss << "  \"covariance\": true," << endl;
  json_ss << "  \"covariance_name\": \"covariance_matrix\"," << endl;
  json_ss << "  \"correlation_name\": \"correlation_matrix\"," << endl;
  json_ss << "  \"weight_cut\": 10," << endl; // fro discussion with Stephen
  json_ss << "  \"single_syst\": false," << endl;
  json_ss << "  \"syst_idx\": \"\"," << endl;
  json_ss << "  \"var_names\": [\"selmu_costheta\", \"selmu_mom\"]," << endl;
  json_ss << "  \"projection\": false," << endl;
  json_ss << "  \"plot_variable\": \"\"," << endl;
  json_ss << "  \"pdf_print\": false," << endl;
  json_ss << "  \"cov_sample_binning\": {" << endl;
  json_ss << "    \"0\": \"" << __binning_file_path__ << "\"" << endl;
  json_ss << "    \"1\": \"" << __binning_file_path__ << "\"" << endl;
  json_ss << "    \"2\": \"" << __binning_file_path__ << "\"" << endl;
  json_ss << "  }," << endl;
  json_ss << "  \"files\":[" << endl;
  json_ss << json_files_ss.str();
  json_ss << "  ]" << endl;

  string output_json_file_path = "./config_xsllhDetVar_" + TToolBox::get_filename_without_ext_from_file_path(__binning_file_path__) + ".json";
  TToolBox::write_string_in_file(output_json_file_path, json_ss.str());
  cout << INFO << "Config file has been written as : " << output_json_file_path << endl;
  cout << INFO << "xsllhDetVar -j " << output_json_file_path << endl;

}
