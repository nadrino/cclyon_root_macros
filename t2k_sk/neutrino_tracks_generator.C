


int __nb_events__ = 1000;

double __earth_radius_in_km__ = 6371;
double __neutrino_production_altitude_in_km__ = 15;

std::map<std::string, double> __event_container__;

TFile* __output_TFile__;
TRandom3* __root_PRNG__;


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
  TTree* events_tree = new TTree("neutrino_tracks", "neutrino_tracks");
  for(auto &observable : __event_container__){
    events_tree->Branch(observable.first.c_str(), &observable.second);
  }

  for(int i_event = 0 ; i_event < __nb_events__ ; i_event++ ){

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

    __event_container__["theta_SK_vertex"] = -TMath::ACos(__event_container__["Z_SK_vertex"]/__event_container__["R_SK_vertex"]) + TMath::Pi()/2.;
    __event_container__["phi_SK_vertex"] = TMath::ATan(__event_container__["Y_vertex"]/__event_container__["X_vertex"]);

    __event_container__["SK_solid_angle_vertex"] = 0.040*0.040/(__event_container__["R_SK_vertex"]*__event_container__["R_SK_vertex"]);

    events_tree->Fill();
  }

  events_tree->Write();
  __output_TFile__->Close();
  exit(0);

}
