
bool do_tfile_is_valid(TFile *input_tfile_, bool check_if_writable_=false){

    if(input_tfile_ == nullptr){
        return false;
    }

    if(not input_tfile_->IsOpen()){
        return false;
    }

    if(check_if_writable_ && not input_tfile_->IsWritable()){
        return false;
    }

    return true;

}
bool do_tfile_is_valid(std::string input_file_path_){
  bool file_is_valid = false;
  Int_t old_verbosity = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;
  TFile* input_tfile = TFile::Open(input_file_path_.c_str(), "READ");
  if(do_tfile_is_valid(input_tfile)){
      file_is_valid = true;
      input_tfile->Close();
  }
  delete input_tfile;
  gErrorIgnoreLevel = old_verbosity;
  return file_is_valid;
}

void do_root_file_is_valid(string path=""){

  std::remove("do_root_file_is_valid.tmp");

  if(path.empty()){
    cerr << "Usage : root \"do_root_file_is_valid.C(\\\"rootfile.root\\\")\"" << std::endl;
    exit(1);
  }

  bool is_valid = do_tfile_is_valid(path);
  std::ofstream output_file ("do_root_file_is_valid.tmp", std::ofstream::out);
  output_file << is_valid << std::endl;
  output_file.close();
  exit(1);

}
