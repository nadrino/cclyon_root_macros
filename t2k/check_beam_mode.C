void check_beam_mode(){

  std::remove("check_beam_mode.out");

  int nb_files_openned = 0;
  TFile* flat_file;

  TIter next(gROOT->GetListOfGlobals(1));
  TGlobal *global;
  while ((global=(TGlobal*)next())) {
    TString type = global->GetTypeName();
    if (type=="TFile") {
      TFile *file = (TFile*)gInterpreter->Calc(global->GetName());
      if (file && file->IsOpen()){
        nb_files_openned++;
        flat_file = file;
      }
    }
  }

  if(nb_files_openned != 1){
    cerr << "nb_files_openned != 1" << endl;
    cerr << "Usage : root path/to/root/file.root /path/to/this/script/check_beam_mode.C" << endl;
    exit(1);
  }

  TTree* flattree = (TTree*) flat_file->Get("flattree");

  if(flattree == 0){
    cerr << "flattree could not be found in " << flat_file->GetName() << endl;
    exit(1);
  }

  flattree->GetEntry(0);
  double run = flattree->GetLeaf("sRun")->GetValue(0);

  int beam_mode = 0;
  if(int(run/1E7) == 9){
    beam_mode = 1;
    cout << "Numu beam mode detected." << endl;
  }
  else if(int(run/1E7) == 8){
    beam_mode = -1;
    cout << "AntiNumu beam mode detected." << endl;
  }
  else {
    beam_mode = 0;
    cout << "Unknown beam mode : run=" << run << ". Conventions may have changed." << endl;
    cout << "Please check : https://www.t2k.org/nd280/datacomp/production006/mcp/mcProdSummary" << endl;
  }

  std::ofstream output_file ("check_beam_mode.tmp", std::ofstream::out);
  output_file << beam_mode << std::endl;
  output_file.close();

  exit(0);

}
