

std::vector<double> AtmSubGeVRecoEdges;
std::vector<double> CosZenithRecoEdges;

void init();


void detSystReweight(){

  init();



}


void init(){

  // binning edges from Fig. 7.2 in m.jiang's thesis
  const double tenTo2p0 =    100.0;
  const double tenTo2p4 =    251.2;
  const double tenTo2p6 =    398.1;
  const double tenTo2p8 =    631.0;
  const double tenTo3p0 =   1000.;
  const double tenTo3p2 =   1585.;
  const double tenTo3p4 =   2512.;
  const double tenTo3p7 =   5012.;
  const double tenTo4p0 =  10000.;
  const double tenTo5p0 = 100000.;

  AtmSubGeVRecoEdges = { // in MeV
    // 5 bins
    tenTo2p0, tenTo2p4, tenTo2p6, tenTo2p8, tenTo3p0, tenTo3p2
  };

  // 10 bins
  CosZenithRecoEdges = {
    // 0.20 steps all the way
    -1.00, -0.80, -0.60, -0.40, -0.20,
    0.00,  0.20,  0.40,  0.60,  0.80,  1.00
  };

}
