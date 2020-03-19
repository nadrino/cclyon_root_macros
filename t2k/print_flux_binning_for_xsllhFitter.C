

TFile* __covariance_matrix_file__;
std::vector<TObject*> __binning_axis_list__;


void print_flux_binning_for_xsllhFitter(){

  if(gROOT->GetListOfFiles()->GetSize() == 0){
    cerr << TToolBox::ERROR << "No file open." << endl;
    exit(1);
  }
  if(gROOT->GetListOfFiles()->GetSize() > 1){
    cerr << TToolBox::ERROR << "More than one file openned." << endl;
    exit(1);
  }

  __covariance_matrix_file__ = (TFile*) gROOT->GetListOfFiles()->At(0);
  __binning_axis_list__ = TToolBox::get_list_of_object_from_directory(
    __covariance_matrix_file__->GetDirectory("/"),
    "TAxis"
  );

  stringstream ss;
  ss << 1 << endl;
  int n_bins = 0;
  for(int i_axis = 0 ; i_axis < int(__binning_axis_list__.size()) ; i_axis++){
    TAxis* the_axis = (TAxis*) __binning_axis_list__[i_axis];
    string axis_name = the_axis->GetName();
    vector<string> splitted_axis_name = TToolBox::split_string(axis_name, "_");

    if(splitted_axis_name.size() < 3) continue; // ignore other objects
    else if(splitted_axis_name[0] != "nd5") continue; // ignore SK

    int beam_mode;
    if(splitted_axis_name[1] == "numode") beam_mode = 1;
    else if(splitted_axis_name[1] == "anumode") beam_mode = -1;
    else continue;

    int particle_id;
    if(splitted_axis_name[2] == "numu") particle_id = 14;
    else if(splitted_axis_name[2] == "numub") particle_id = -14;
    else if(splitted_axis_name[2] == "nue") particle_id = 12;
    else if(splitted_axis_name[2] == "nueb") particle_id = -12;
    else continue;

    for(int i_bin = 0 ; i_bin < the_axis->GetNbins(); i_bin++){
      n_bins++;
      ss << std::setfill('0') << std::setw(4) << int(1000*the_axis->GetBinLowEdge(i_bin+1)) << " ";
      ss << std::setfill('0') << std::setw(4) << int(1000*the_axis->GetBinUpEdge(i_bin+1)) << " ";
      ss << particle_id << " ";
      ss << beam_mode << endl;
    }
  }

  cout << ss.str() << endl;
  TToolBox::write_string_in_file("flux_binning.txt", ss.str());

  TMatrixD* total_flux_cov = (TMatrixD*) __covariance_matrix_file__->Get("total_flux_cov");

  TMatrixD* chopped_flux_cov = new TMatrixD(n_bins,n_bins);
  for(int i_cov = 0 ; i_cov < n_bins ; i_cov++){
    for(int j_cov = 0 ; j_cov < n_bins ; j_cov++){
      (*chopped_flux_cov)[i_cov][j_cov] = (*total_flux_cov)[i_cov][j_cov];
    }
  }

  TFile* output_file = TFile::Open("chopped_flux_cov.root","RECREATE");
  chopped_flux_cov->Write("chopped_flux_cov_TMatrixD");

  TNamed binning_string("binning_xsLLhFitter_TNamed", ss.str().c_str());
  binning_string.Write();

  output_file->Close();

}
