

void lf()
{
   TIter next(gROOT->GetListOfGlobals(1));
   TGlobal *global;
   while ((global=(TGlobal*)next())) {
     TString type = global->GetTypeName();
     if (type=="TFile") {
       TFile *file = (TFile*)gInterpreter->Calc(global->GetName());
       if (file && file->IsOpen())
         printf("%s: %s\n", global->GetName(),file->GetName());
     }
   }
 }

void check_beam_mode(){

lf();

  // if(argc != 2){
  //   return;
  // }
  //
  // string flat_file_path(argv[1]);
  // TFile* flat_file = TFile::Open(flat_file_path.c_str());
  // TTree* flattree = (TTree*) flat_file->Get("flattree");
  //
  // flattree->GetEntry(0);
  // double run = flattree->GetLeaf("sRun")->GetValue(0);
  //
  // if(int(run/1E7) == 9) cout << "1" << endl;
  // else if(int(run/1E7) == 8) cout << "-1" << endl;
  // else cout << "0" << endl;
  //
  // return;

}
