


int __nb_events__ = 1000;

double __earth_radius_in_km__ = 6371;
double __neutrino_production_altitude_in_km__ = 15;

std::map<std::string, double> __event_container__;

TFile* __output_TFile__;
TRandom3* __root_PRNG__;
TF1* __earth_density_TF1__;

void build_PREM();

void neutrino_tracks_generator(){

  __root_PRNG__ = new TRandom3(time(NULL));

  __event_container__["X_vertex"] = 0.;
  __event_container__["Y_vertex"] = 0.;
  __event_container__["Z_vertex"] = 0.;
  __event_container__["Z_SK_vertex"] = 0.;

  __event_container__["R_vertex"] = 0.;
  __event_container__["theta_vertex"] = 0.;
  __event_container__["phi_vertex"] = 0.;

  __event_container__["R_SK_vertex"] = 0.;
  __event_container__["theta_SK_vertex"] = 0.;
  __event_container__["phi_SK_vertex"] = 0.;

  __event_container__["SK_solid_angle_vertex"] = 0.;

  __output_TFile__ = TFile::Open("output.root", "RECREATE");

  build_PREM();

  TTree* events_tree = new TTree("neutrino_tracks", "neutrino_tracks");
  for(auto &observable : __event_container__){
    events_tree->Branch(observable.first.c_str(), &observable.second);
  }

  for(int i_event = 0 ; i_event < __nb_events__ ; i_event++ ){

    TToolBox::display_loading(i_event, __nb_events__);

    __event_container__["R_vertex"] = __earth_radius_in_km__ + __neutrino_production_altitude_in_km__;
    __event_container__["theta_vertex"] = __root_PRNG__->Rndm()*TMath::Pi();
    __event_container__["phi_vertex"] = __root_PRNG__->Rndm()*TMath::Pi()*2;

    __event_container__["X_vertex"] = __event_container__["R_vertex"]*TMath::Sin(__event_container__["theta_vertex"])*TMath::Cos(__event_container__["phi_vertex"]);
    __event_container__["Y_vertex"] = __event_container__["R_vertex"]*TMath::Sin(__event_container__["theta_vertex"])*TMath::Sin(__event_container__["phi_vertex"]);
    __event_container__["Z_vertex"] = __event_container__["R_vertex"]*TMath::Cos(__event_container__["theta_vertex"]);
    __event_container__["Z_SK_vertex"] = __event_container__["Z_vertex"] - __earth_radius_in_km__;

    __event_container__["R_SK_vertex"] += __event_container__["X_vertex"]*__event_container__["X_vertex"];
    __event_container__["R_SK_vertex"] += __event_container__["Y_vertex"]*__event_container__["Y_vertex"];
    __event_container__["R_SK_vertex"] += __event_container__["Z_SK_vertex"]*__event_container__["Z_SK_vertex"];
    __event_container__["R_SK_vertex"] = TMath::Sqrt(__event_container__["R_SK_vertex"]);

    __event_container__["theta_SK_vertex"] = TMath::ACos(__event_container__["Z_SK_vertex"]/__event_container__["R_SK_vertex"]);
    __event_container__["phi_SK_vertex"] = TMath::ATan(__event_container__["Y_vertex"]/__event_container__["X_vertex"]);

    __event_container__["SK_solid_angle_vertex"] = 0.040*0.040/(__event_container__["R_SK_vertex"]*__event_container__["R_SK_vertex"]);

    events_tree->Fill();
  }

  events_tree->Write();
  __output_TFile__->Close();
  exit(0);

}


void build_PREM(){

  // Preliminary reference Earth model
  // https://www.sciencedirect.com/science/article/abs/pii/0031920181900467?via%3Dihub

  __output_TFile__->mkdir("TF1");
  __output_TFile__->cd("TF1");


  std::vector<std::string> layer_label;
  std::vector<double> layer_outer_bound;
  std::vector<std::vector<double>> layer_polynomial_coefficients;

  std::vector< std::map< std::string, std::vector<double> > > layers_data;

  layer_label.emplace_back("InnerCore");
  layer_outer_bound.emplace_back(1221.5);
  layer_polynomial_coefficients.emplace_back(std::vector<double>({13.0885, 0, -8.8381}));

  layer_label.emplace_back("OutterCore");
  layer_outer_bound.emplace_back(3480.0);
  layer_polynomial_coefficients.emplace_back(std::vector<double>({12.5815, -1.2638, -3.6426, -5.5281}));

  layer_label.emplace_back("LowerMantle_1");
  layer_outer_bound.emplace_back(3630.0);
  layer_polynomial_coefficients.emplace_back(std::vector<double>({7.9565, -6.4761, 5.5283, -3.0807}));

  layer_label.emplace_back("LowerMantle_2");
  layer_outer_bound.emplace_back(5600.0);
  layer_polynomial_coefficients.emplace_back(std::vector<double>({7.9565, -6.4761, 5.5283, -3.0807}));

  layer_label.emplace_back("LowerMantle_3");
  layer_outer_bound.emplace_back(5701.0);
  layer_polynomial_coefficients.emplace_back(std::vector<double>({7.9565, -6.4761, 5.5283, -3.0807}));

  layer_label.emplace_back("TransitionZone_1");
  layer_outer_bound.emplace_back(5771.0);
  layer_polynomial_coefficients.emplace_back(std::vector<double>({5.3197, -1.4836}));

  layer_label.emplace_back("TransitionZone_2");
  layer_outer_bound.emplace_back(5971.0);
  layer_polynomial_coefficients.emplace_back(std::vector<double>({11.2494, -8.0298}));

  layer_label.emplace_back("TransitionZone_3");
  layer_outer_bound.emplace_back(6151.0);
  layer_polynomial_coefficients.emplace_back(std::vector<double>({7.1089, -3.8045}));

  layer_label.emplace_back("LVZ");
  layer_outer_bound.emplace_back(6291.0);
  layer_polynomial_coefficients.emplace_back(std::vector<double>({2.6910, 0.6924}));

  layer_label.emplace_back("LID");
  layer_outer_bound.emplace_back(6346.6);
  layer_polynomial_coefficients.emplace_back(std::vector<double>({2.6910, 0.6924}));

  layer_label.emplace_back("Crust_1");
  layer_outer_bound.emplace_back(6356.0);
  layer_polynomial_coefficients.emplace_back(std::vector<double>({2.900}));

  layer_label.emplace_back("Crust_2");
  layer_outer_bound.emplace_back(6368.0);
  layer_polynomial_coefficients.emplace_back(std::vector<double>({2.600}));

  std::vector<TF1*> functions_list;
  double last_bound = 0.;
  for(int i_fct = 0 ; i_fct < layer_label.size() ; i_fct++){

    std::stringstream formulae;
    formulae << "(";
    for(int i_exposant = 0 ; i_exposant < layer_polynomial_coefficients[i_fct].size() ; i_exposant++){
      if(i_exposant != 0) formulae << " + ";
      formulae << layer_polynomial_coefficients[i_fct][i_exposant] << "*TMath::Power(x," << i_exposant << ")";
    }
    formulae << ")";
    formulae << "*" << "(x>=" << last_bound << ")*(x<" << layer_outer_bound[i_fct] << ")";

    functions_list.emplace_back(
      new TF1(layer_label[i_fct].c_str(), formulae.str().c_str(), 0., layer_outer_bound.back()+50.)
    );
    functions_list.back()->Write();

  }

  __earth_density_TF1__ = new TF1("PREM", TToolBox::join_vector_string(layer_label, " + ").c_str(), 0., layer_outer_bound.back()+50.);
  __earth_density_TF1__->Write("PREM_TF1");

  __output_TFile__->cd("");

}
