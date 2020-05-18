


void do_root_file_is_valid(string path=""){

  std::remove("do_root_file_is_valid.tmp");

  if(path.empty()){
    cerr << "Usage : root \"do_root_file_is_valid.C(\\\"rootfile.root\\\")\"" << std::endl;
    exit(1);
  }

  bool is_valid = TToolBox::do_tfile_is_valid(path);
  std::ofstream output_file ("do_root_file_is_valid.tmp", std::ofstream::out);
  output_file << is_valid << std::endl;
  output_file.close();
  exit(1);

}
