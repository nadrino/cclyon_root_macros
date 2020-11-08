
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
  PCStop, // not yet in data
  PCThru,
  UpStop_mu,
  UpThruNonShower_mu,
  UpThruShower_mu
};

string get_Enu_rec_QE(string lepton_name_, string lepton_mom_);

void plot_sk_atm_data(){

  cout << "Openning $REPO_DIR/sk-atm-data/atm_minituples/sk4_fcmc_18a_fQv6r0_minituple_100yr_01.root" << endl;
  TFile* sk_tfile = TFile::Open("$REPO_DIR/sk-atm-data/atm_minituples/sk4_fcmc_18a_fQv6r0_minituple_100yr_01.root");
  TTree* atm_minituple = (TTree*) sk_tfile->Get("atm_minituple");

  // string norm_string = "(oscweight3f*solarweight*5000./(365.25*100.0))";
  // string norm_string = "(oscweight3f*solarweight*3244.4/(365.25*100.0))";
  string norm_string = "(solarweight*3244.4/(365.25*100.0))";


  string draw_str;
  string cuts_str;


  ////////////////////////////////////////////
  // CUTS
  ////////////////////////////////////////////
  map<string, vector<string>> cuts_map;
  cuts_map["FC_SubGeV_nue_nuebar"] = vector<string>();
  cuts_map["FC_SubGeV_nue_nuebar"].emplace_back(
    "(ATMPDEventType == " + to_string(SubGeV_elike_0dcy) + " || ATMPDEventType == " + to_string(SubGeV_elike_1dcy) + ")"
  );
  cuts_map["FC_MultiGeV_nue_nuebar"] = vector<string>();
  cuts_map["FC_MultiGeV_nue_nuebar"].emplace_back(
    "( ATMPDEventType == " + to_string(MultiGeV_elike_nue)
    + " || ATMPDEventType == " + to_string(MultiGeV_elike_nuebar)
    + " || ATMPDEventType == " + to_string(MultiRing_elike_nue)
    + " || ATMPDEventType == " + to_string(MultiRing_elike_nuebar)
    + " )"
  );

  cuts_map["FC_SubGeV_numu_numubar"] = vector<string>();
  cuts_map["FC_SubGeV_numu_numubar"].emplace_back(
    "( ATMPDEventType == " + to_string(SubGeV_mulike_0dcy)
    + " || ATMPDEventType == " + to_string(SubGeV_mulike_1dcy)
    + " || ATMPDEventType == " + to_string(SubGeV_mulike_2dcy)
    + " )"
  );
  cuts_map["FC_MultiGeV_numu_numubar"] = vector<string>();
  cuts_map["FC_MultiGeV_numu_numubar"].emplace_back(
    "( ATMPDEventType == " + to_string(MultiGeV_mulike)
    + " || ATMPDEventType == " + to_string(MultiRing_mulike)
    + " )"
  );
  cuts_map["PC_Stop"] = vector<string>();
  cuts_map["PC_Stop"].emplace_back(
    "( ATMPDEventType == " + to_string(PCStop)
    + " )"
  );
  cuts_map["PC_Thru"] = vector<string>();
  cuts_map["PC_Thru"].emplace_back(
    "( ATMPDEventType == " + to_string(PCThru)
    + " )"
  );

  cuts_map["FC_SubGeV_nue_nuebar_true"] = vector<string>();
  cuts_map["FC_SubGeV_nue_nuebar_true"].emplace_back("(ipnu[0] == 12 || ipnu[0] == -12)");


  ////////////////////////////////////////////
  // HISTOGRAMS
  ////////////////////////////////////////////
  TCanvas* c_h1 = new TCanvas("c_h1", "c_h1", 800, 700);
  // draw_str = get_Enu_rec_QE("electron", "genmom") + "/1000.";
  draw_str = "TMath::Sqrt( TMath::Power(pnu[0],2)+TMath::Power(pnu[1],2)+TMath::Power(pnu[2],2) )";

  TH1D* hSubGeV = TToolBox::get_TH1D_log_binning("hSubGeV", "FC_SubGeV_nue_nuebar", 50, 0.1, 10000.);
  hSubGeV->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
  hSubGeV->GetYaxis()->SetTitle("Events/5000 Days");
  hSubGeV->SetLineColor(kBlue);
  hSubGeV->SetLineStyle(kDashed);
  cuts_str = "(";
  cuts_str +=  TToolBox::join_vector_string(cuts_map["FC_SubGeV_nue_nuebar"], " && ");
  cuts_str += ")";
  atm_minituple->Draw( (draw_str + ">>hSubGeV" ).c_str(), (cuts_str + "*" + norm_string).c_str(), "goff");

  TH1D* hMultiGeV = TToolBox::get_TH1D_log_binning("hMultiGeV", "FC_MultiGeV_nue_nuebar", 50, 0.1, 10000.);
  hMultiGeV->SetLineColor(kCyan);
  cuts_str = "(";
  cuts_str +=  TToolBox::join_vector_string(cuts_map["FC_MultiGeV_nue_nuebar"], " && ");
  cuts_str += ")";
  atm_minituple->Draw( (draw_str + ">>hMultiGeV" ).c_str(), (cuts_str + "*" + norm_string).c_str(), "goff");


  TCanvas* c_h2 = new TCanvas("c_h2", "c_h2", 800, 700);
  draw_str = "TMath::Sqrt( TMath::Power(pnu[0],2)+TMath::Power(pnu[1],2)+TMath::Power(pnu[2],2) )";

  TH1D* hmuSubGeV = TToolBox::get_TH1D_log_binning("hmuSubGeV", "FC_SubGeV_numu_numubar", 50, 0.1, 10000.);
  hmuSubGeV->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
  hmuSubGeV->GetYaxis()->SetTitle("Events/5000 Days");
  hmuSubGeV->SetLineColor(kBlue);
  hmuSubGeV->SetLineStyle(kDashed);
  cuts_str = "(";
  cuts_str +=  TToolBox::join_vector_string(cuts_map["FC_SubGeV_numu_numubar"], " && ");
  cuts_str += ")";
  atm_minituple->Draw( (draw_str + ">>hmuSubGeV" ).c_str(), (cuts_str + "*" + norm_string).c_str(), "goff");

  TH1D* hmuMultiGeV = TToolBox::get_TH1D_log_binning("hmuMultiGeV", "FC_MultiGeV_numu_numubar", 50, 0.1, 10000.);
  hmuMultiGeV->SetLineColor(kCyan);
  cuts_str = "(";
  cuts_str +=  TToolBox::join_vector_string(cuts_map["FC_MultiGeV_numu_numubar"], " && ");
  cuts_str += ")";
  atm_minituple->Draw( (draw_str + ">>hmuMultiGeV" ).c_str(), (cuts_str + "*" + norm_string).c_str(), "goff");

  TH1D* hmuPCStop = TToolBox::get_TH1D_log_binning("hmuPCStop", "PC_Stop", 50, 0.1, 10000.);
  hmuPCStop->SetLineColor(kViolet);
  hmuPCStop->SetFillColor(kViolet);
  cuts_str = "(";
  cuts_str +=  TToolBox::join_vector_string(cuts_map["PC_Stop"], " && ");
  cuts_str += ")";
  atm_minituple->Draw( (draw_str + ">>hmuPCStop" ).c_str(), (cuts_str + "*" + norm_string).c_str(), "goff");

  TH1D* hmuPCThru = TToolBox::get_TH1D_log_binning("hmuPCThru", "PC_Thru", 50, 0.1, 10000.);
  hmuPCThru->SetLineColor(kRed);
  hmuPCThru->SetFillColor(kRed);
  hmuPCThru->SetFillStyle(3005);
  cuts_str = "(";
  cuts_str +=  TToolBox::join_vector_string(cuts_map["PC_Thru"], " && ");
  cuts_str += ")";
  atm_minituple->Draw( (draw_str + ">>hmuPCThru" ).c_str(), (cuts_str + "*" + norm_string).c_str(), "goff");


  ////////////////////////////////////////////
  // PLOT
  ////////////////////////////////////////////
  c_h1->cd();
  hSubGeV->Draw("HIST");
  hMultiGeV->Draw("HISTSAME");
  gPad->SetLogx(1);
  gPad->SetGridx();
  gPad->BuildLegend();
  hSubGeV->SetTitle("");
  TToolBox::save_canvas(c_h1, "sk_nue_atm_FC");

  c_h2->cd();
  hmuSubGeV->Draw("HIST");
  hmuMultiGeV->Draw("HISTSAME");
  hmuPCStop->Draw("HISTSAME");
  hmuPCThru->Draw("HISTSAME");
  gPad->SetLogx(1);
  gPad->SetGridx();
  gPad->BuildLegend();
  hmuSubGeV->SetTitle("");
  TToolBox::save_canvas(c_h2, "sk_numu_atm_FC");


  // TCanvas* c_h2 = new TCanvas("c_h2", "c_h2", 800, 800);
  // TH2D* h2 = TToolBox::get_TH2D_log_binning("h2", "FC_SubGeV_nue_nuebar", 30, 0.1, 100., 10, -1, 1, "X");
  // h2->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
  // h2->GetYaxis()->SetTitle("cos zenith");
  // h2->GetZaxis()->SetTitle("Events/5000 Days");
  // string h2_draw_str = "gencz:genmom/1000.";
  // h2_draw_str += ">>h2";

  // atm_minituple->Draw(
  //   h2_draw_str.c_str(),
  //   (cuts_str + "*" + norm_string).c_str(),
  //   "goff"
  // );
  //
  // c_h2->cd();
  // h2->Draw("COLZ");
  // gPad->SetLogx(1);
  // TToolBox::fix_TH2D_display(h2);

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
