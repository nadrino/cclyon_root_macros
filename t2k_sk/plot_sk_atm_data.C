
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

string get_Enu_rec_QE(string lepton_name_, string lepton_mom_);

void plot_sk_atm_data(){

  TFile* sk_tfile = TFile::Open("$RESOURCES_DIR/T2KSKjoint/sk4_fcmc_18a_fQv6r0_minituple_100yr_01.root");
  TTree* atm_minituple = (TTree*) sk_tfile->Get("atm_minituple");

  map<string, vector<string>> cuts_map;
  cuts_map["FC_Sub-GeV_nue_nuebar"] = vector<string>();
  cuts_map["FC_Sub-GeV_nue_nuebar"].emplace_back(
    "(ATMPDEventType == " + to_string(SubGeV_elike_0dcy) + " || ATMPDEventType == " + to_string(SubGeV_elike_1dcy) + ")"
  );
  // Single-ring
  cuts_map["FC_Sub-GeV_nue_nuebar"].emplace_back("nring == 1");
  // cuts_map["FC_Sub-GeV_nue_nuebar"].emplace_back("fqnmrfit == 1");
  // cuts_map["FC_Sub-GeV_nue_nuebar"].emplace_back("evis < 1330"); // visible energy  (MeV/c) :  This is the sum of amome of each rings
  // cuts_map["FC_Sub-GeV_nue_nuebar"].emplace_back("evis > 100"); // visible energy  (MeV/c) :  This is the sum of amome of each rings
  // cuts_map["FC_Sub-GeV_nue_nuebar"].emplace_back("amome[0] > 100");
  cuts_map["FC_Sub-GeV_nue_nuebar"].emplace_back("ip[0] == 2"); // particle type by PID (1:gamma 2:electron 3:muon)
  // FCFV
  cuts_map["FC_Sub-GeV_nue_nuebar"].emplace_back("fqwall > 50.0"); // note
  // cuts_map["FC_Sub-GeV_nue_nuebar"].emplace_back("fqwall > 200"); // atm
  cuts_map["FC_Sub-GeV_nue_nuebar"].emplace_back("fqtowall > 170");

  string norm_string = "(oscweight3f*solarweight*5000./(365.25*100.0))"; // 3244.4

  string cuts_str;
  cuts_str += "(";
  cuts_str +=  TToolBox::join_vector_string(cuts_map["FC_Sub-GeV_nue_nuebar"], " && ");
  // cuts_str += "nring==1 && evis<1330 && ip[0] == 2 && wall>200";
  cuts_str += ")";


  TCanvas* c_h2 = new TCanvas("c_h2", "c_h2", 800, 800);
  TH2D* h2 = TToolBox::get_TH2D_log_binning("h2", "FC_Sub-GeV_nue_nuebar", 20, 0.1, 10., 10, -1, 1, "X");
  h2->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
  h2->GetYaxis()->SetTitle("cos zenith");
  h2->GetZaxis()->SetTitle("Events/5000 Days");
  string h2_draw_str = "gencz:genmom/1000.";
  h2_draw_str += ">>h2";

  atm_minituple->Draw(
    h2_draw_str.c_str(),
    (cuts_str + "*" + norm_string).c_str(),
    "goff"
  );

  c_h2->cd();
  h2->Draw("COLZ");
  gPad->SetLogx(1);
  TToolBox::fix_TH2D_display(h2);

  TCanvas* c_h1 = new TCanvas("c_h1", "c_h1", 800, 800);
  TH1D* h1 = TToolBox::get_TH1D_log_binning("h1", "FC_Sub-GeV_nue_nuebar", 20, 0.1, 10.);
  h1->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
  h1->GetYaxis()->SetTitle("Events/5000 Days");
  string h1_draw_str = get_Enu_rec_QE("electron", "genmom") + "/1000.";
  h1_draw_str += ">>h1";

  atm_minituple->Draw(
    h1_draw_str.c_str(),
    (cuts_str + "*" + norm_string).c_str(),
    "goff"
  );

  c_h1->cd();
  h1->Draw("HIST");
  gPad->SetLogx(1);

  cout << "atm_minituple->Draw(\"" << h1_draw_str << "\", \"" << (cuts_str + "*" + norm_string) << "\");" << endl;

}

string get_Enu_rec_QE(string lepton_name_, string lepton_mom_){

  string output;
  string proton_mass = "938.272";
  string neutron_mass = "939.565";
  string muon_mass = "105.658";
  string electron_mass = "105.658";
  string Eb = "25";

  string lepton_energy;
  string lepton_mass;
  if(lepton_name_ == "electron"){
    lepton_mass = electron_mass;
  }
  else if(lepton_name_ == "muon"){
    lepton_mass = muon_mass;
  }
  lepton_energy = "TMath::Sqrt( TMath::Power(" + lepton_mom_ +",2) + TMath::Power(" + lepton_mass + ",2) )";

  string numerator = "(";
  numerator += "2*" + lepton_energy + " * (" + neutron_mass + " - " + Eb + ")";
  numerator += " - ";
  numerator += "TMath::Power(" + lepton_mass + ",2)";
  numerator += " + ";
  numerator += "2*" + neutron_mass + "*" + Eb;
  numerator += " - ";
  numerator += "TMath::Power(" + Eb + ",2)";
  numerator += " + ";
  numerator += "TMath::Power(" + proton_mass + ",2)";
  numerator += " - ";
  numerator += "TMath::Power(" + neutron_mass + ",2)";
  numerator += ")";
  string denominator = "(";
  denominator += "2*( ";
  denominator += neutron_mass + " - " + Eb + " - " + lepton_energy + " + " + lepton_mom_; // assuming forward momentum in average
  denominator += " )";
  denominator += ")";

  output = numerator + "/" + denominator;
  return output;

}
