#include<vector>

using namespace std;

vector<TFile*> get_list_of_openned_tfiles();

void check_beam_mode(){

  vector<TFile*> openned_tfiles;
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

  if(openned_tfiles.size() != 1){
    cerr << "openned_tfiles.size() != 1" << endl;
    cerr << "Usage : root path/to/root/file.root /path/to/this/script/check_beam_mode.C" << endl;
  }

  TTree* flattree = (TTree*) flat_file->Get("flattree");

  flattree->GetEntry(0);
  double run = flattree->GetLeaf("sRun")->GetValue(0);

  if(int(run/1E7) == 9) cout << "1" << endl;
  else if(int(run/1E7) == 8) cout << "-1" << endl;
  else cout << "0" << endl;

  return;

}
