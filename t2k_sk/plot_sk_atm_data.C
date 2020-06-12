
enum ATMPDEventType{
  SubGeV_elike_0dcy = 1,
  SubGeV_elike_1dcy,
  SubGeV_SingleRing_pi0like,
  SubGeV_mulike_0dcy,
  SubGeV_mulike_1dcy,
  SubGeV_mulike_2dcy,
  SubGeV_pi0like,
  MultiGeV_elike_nue,
  MultiGeV_elike_nuebar,
  MultiGeV_mulike,
  MultiRing_elike_nue,
  MultiRing_elike_nuebar,
  MultiRing_mulike,
  MultiRingOther_1,
  PCStop,
  PCThru,
  UpStop_mu,
  UpThruNonShower_mu,
  UpThruShower_mu
};

void plot_sk_atm_data(){

  TFile* sk_tfile = TFile::Open("$RESOURCES_DIR/T2KSKjoint/sk4_fcmc_18a_fQv6r0_minituple_100yr_01.root");
  TTree* atm_minituple = (TTree*) sk_tfile->Get("atm_minituple");

  map<string, vector<string>> cuts_map;
  cuts_map["FC_Sub-GeV_nue_nuebar"] = vector<string>();
  cuts_map["FC_Sub-GeV_nue_nuebar"].emplace_back("(ATMPDEventType == " + to_string(SubGeV_elike_0dcy) + " || ATMPDEventType == " + to_string(SubGeV_elike_1dcy) + ")"); // SubGeV_elike_0dcy
  // Single-ring
  // cuts_map["FC_Sub-GeV_nue_nuebar"].emplace_back("nring == 1");
  cuts_map["FC_Sub-GeV_nue_nuebar"].emplace_back("evis < 1330"); // visible energy  (MeV/c) :  This is the sum of amome of each rings
  cuts_map["FC_Sub-GeV_nue_nuebar"].emplace_back("evis > 100"); // visible energy  (MeV/c) :  This is the sum of amome of each rings
  // cuts_map["FC_Sub-GeV_nue_nuebar"].emplace_back("amome[0] > 100");
  // cuts_map["FC_Sub-GeV_nue_nuebar"].emplace_back("ip[0] == 2"); // particle type by PID (1:gamma 2:electron 3:muon)
  // FCFV
  cuts_map["FC_Sub-GeV_nue_nuebar"].emplace_back("fqwall > 50.0"); // note
  // cuts_map["FC_Sub-GeV_nue_nuebar"].emplace_back("fqwall > 200"); // atm
  cuts_map["FC_Sub-GeV_nue_nuebar"].emplace_back("fqtowall > 170");

  // string norm_string = "(oscweight3f*solarweight*3244.4/(365.25*100.0))";
  string norm_string = "(oscweight3f*3244.4/(365.25*100.0))";

  string cuts_str;
  cuts_str += "(";
  cuts_str +=  TToolBox::join_vector_string(cuts_map["FC_Sub-GeV_nue_nuebar"], " && ");
  // cuts_str += "nring==1 && evis<1330 && ip[0] == 2 && wall>200";
  cuts_str += ")";

  TH2D* h2 = TToolBox::get_TH2D_log_binning("h2", "FC_Sub-GeV_nue_nuebar", 20, 0.1, 10., 10, -1, 1, "X");
  h2->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
  h2->GetYaxis()->SetTitle("cos zenith");
  h2->GetZaxis()->SetTitle("Events/5000 Days");

  string draw_str = "gencz:evis/1000.";
  draw_str += ">>h2";

  atm_minituple->Draw(
    draw_str.c_str(),
    (cuts_str + "*" + norm_string).c_str(),
    "goff"
  );

  h2->Draw("COLZ");
  gPad->SetLogx();
  TToolBox::fix_TH2D_display(h2);


}
