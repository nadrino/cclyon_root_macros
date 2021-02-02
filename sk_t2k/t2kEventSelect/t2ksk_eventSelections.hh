#ifndef T2KSK_EVENTSELECTIONS_HH
#define T2KSK_EVENTSELECTIONS_HH

namespace t2ksk{

  // J-PARC beam direction in SK coordinates
  static const double beamdir[3] = { 0.669764, -0.742179, 0.024223 };
  
  // Particle indices for SK arrays
  typedef enum{
    ELECTRON = 1,
    MUON,
    PION
  } fq_particle;
  
  // FiTQun variables
  // Sub-event finder
  int fqnse;
  
  // Single-ring fits
  float fq1rmom[100][7];
  float fq1rnll[100][7];
  float fq1rpos[100][7][3];
  float fq1rdir[100][7][3];

  // pi0 fit
  float fqpi0nll[2];
  float fqpi0mass[2];

  // Multi=ring fit
  int fqmrnring[100];
  
  // Other variables
  int evclass;

  void setTreeAddresses(TTree * tree){
    tree->SetBranchAddress("fqnse", &fqnse);
    
    tree->SetBranchAddress("fq1rmom", &fq1rmom);
    tree->SetBranchAddress("fq1rnll", &fq1rnll);
    tree->SetBranchAddress("fq1rpos", &fq1rpos);
    tree->SetBranchAddress("fq1rdir", &fq1rdir);

    tree->SetBranchAddress("fqpi0nll", &fqpi0nll);
    tree->SetBranchAddress("fqpi0mass", &fqpi0mass);

    tree->SetBranchAddress("fqmrnring", &fqmrnring);
    
    tree->SetBranchAddress("evclass", &evclass);
  }

  // Some useful quantities
  float electron_momentum(){
    return fq1rmom[0][ELECTRON];
  }

  float muon_momentum(){
    return fq1rmom[0][MUON];
  }

  // PID cuts
  float electron_muon_PID(){
    // Electrons are positive
    return fq1rnll[0][MUON]-fq1rnll[0][ELECTRON]-0.2*electron_momentum();
  }

  float piplus_muon_PID(){
    // Pions are positive
    return fq1rnll[0][MUON]-fq1rnll[0][PION] - 0.15*muon_momentum();
  }

  float pi0_electron_PID(){
    // Pions are positive
    return fq1rnll[0][ELECTRON]-fqpi0nll[0]-175+0.875*fqpi0mass[0];
  }

  // FV Variables
  // Lifted from Minituple code
  float ComputeTowall(int nsubevent, fq_particle i_particle)
  {
    float x = fq1rpos[nsubevent][i_particle][0];
    float y = fq1rpos[nsubevent][i_particle][1];
    float z = fq1rpos[nsubevent][i_particle][2];

    float dx = fq1rdir[nsubevent][i_particle][0];
    float dy = fq1rdir[nsubevent][i_particle][1];
    float dz = fq1rdir[nsubevent][i_particle][2];
    
    Double_t const R(1690);
    Double_t l_b(100000.0), H;
    Double_t l_t(100000.0);
    Double_t A, B, C, RAD;
    if(dx!=0 || dy!=0){
      A = (dx*dx+dy*dy);
      B = 2*(x*dx+y*dy);
      C = (x*x+y*y-R*R);
      RAD = (B*B) - (4*A*C);
      l_b = ((-1*B) + sqrt(RAD))/(2*A);
    }
    if (dz==0){return l_b;}
    else if(dz>0){H=1810;}
    else if(dz<0){H=-1810;}
    l_t=(H - z)/dz;
    return  (l_t > l_b ? l_b:l_t);
  }

  // Lifted from minituple code
  float ComputeWall(int nsubevent, fq_particle i_particle)
  {
    float x = fq1rpos[nsubevent][i_particle][0];
    float y = fq1rpos[nsubevent][i_particle][1];
    float z = fq1rpos[nsubevent][i_particle][2];


    float Rmax = 1690.;
    float Zmax = 1810.;
    float rr   = sqrt(x*x + y*y);
    float absz = TMath::Abs(z);
    //check if vertex is outside tank
    float signflg = 1.;
    if (absz>Zmax) signflg = -1.;
    if (rr>Rmax)   signflg = -1.;
    //find min distance to wall
    float distz = TMath::Abs(Zmax-absz);
    float distr = TMath::Abs(Rmax-rr);
    float dwall = signflg*fmin(distz,distr);
    return dwall;
  }

  // Lifted from minituple code
  double ComputeCosBeam( int nsubevent, fq_particle i_particle  )
  {

    float dirx = fq1rdir[nsubevent][i_particle][0];
    float diry = fq1rdir[nsubevent][i_particle][1];
    float dirz = fq1rdir[nsubevent][i_particle][2];
    
    float cosb = dirx*beamdir[0] +
      diry*beamdir[1] +
      dirz*beamdir[2] ;
    return cosb;
  }
  
  // CCQE Erec
  // Lifted from minituple code
  double ComputeErec(int nsubevent, fq_particle i_particle, float cosBeam)
  {

    float LeptonP = fq1rmom[nsubevent][i_particle];

    int flag = -1;
    if (i_particle == ELECTRON) flag = 1;
    else if (i_particle == MUON) flag = 0;
    else throw std::runtime_error("Cannot compute Erec if particle isn't an electron or a muon");


    
    static const double Vnuc  = 27.0        ; // MeV 
    static const double mn    = 939.565346  ; // MeV
    static const double mp    = 938.272013  ; // MeV
    
    static const double me    = 0.510998   ; // MeV
    static const double mm    = 105.65836  ; // MeV

    double mass ;
    if( flag == 1 ) mass = me ;
    else            mass = mm ;
    
    double E  = 0.;
    double LeptonE = sqrt( mass*mass + LeptonP*LeptonP );
    
    
    E  = ( mn - Vnuc)*LeptonE - mass*mass/2. ;
    E +=   mn*Vnuc  - Vnuc*Vnuc/2.;
    E += ( mp*mp - mn*mn)/2.;
    
    E /= ( mn - Vnuc - LeptonE + LeptonP*cosBeam );
    
    // returns value in MeV 
    return E;
  }

  double ComputeErec(int nsubevent, fq_particle i_particle)
  {

    float cosBeam = ComputeCosBeam(nsubevent, i_particle);
    
    return ComputeErec(nsubevent, i_particle, cosBeam);
  }

  // Delta resonance Erec
  double ComputeErecCCDel(int nsubevent, fq_particle i_particle, float cosBeam)
  {
    
    float LeptonP = fq1rmom[nsubevent][i_particle];
    
    int flag = -1;
    if (i_particle == ELECTRON) flag = 1;
    else if (i_particle == MUON) flag = 0;
    else throw std::runtime_error("Cannot compute Erec if particle isn't an electron or a muon");

    
    static const double mn    = 939.565346  ; // MeV
    static const double mp    = 938.272013  ; // MeV
    static const double md    = 1232.0  ; // MeV
    
    
    static const double me    = 0.510998   ; // MeV
    static const double mm    = 105.65836  ; // MeV
    
    double mass ;
    if( flag == 1 ) mass = me ;
    else            mass = mm ;
    
    double E  = 0.;
    double LeptonE = sqrt( mass*mass + LeptonP*LeptonP );

    
    E  = mp *LeptonE - mass*mass/2. ;
    E += ( md*md - mp*mp)/2.;
    E /= ( mp - LeptonE + LeptonP*cosBeam );
    
    // returns value in MeV 
    return E;
  }
  double ComputeErecCCDel(int nsubevent, fq_particle i_particle)
  {
    float cosBeam = ComputeCosBeam(nsubevent, i_particle);

    return ComputeErecCCDel(nsubevent, i_particle, cosBeam);
  }

  // Event selections
  bool is1Re(){
    
    if (evclass != 1) return false;
    else if (ComputeWall(0, ELECTRON) <= 80) return false;
    else if (ComputeTowall(0, ELECTRON) <= 170) return false;
    else if (electron_momentum() <= 30) return false;
    else if (fqmrnring[0] != 1) return false;
    else if (electron_muon_PID() < 0) return false;
    else if (electron_momentum() <= 100) return false;
    else if (fqnse != 1) return false;
    else if (ComputeErec(0, ELECTRON) >= 1250) return false;
    else if (pi0_electron_PID() >= 0) return false;
    else return true;
  }

  bool is1Re1de(){
    if (evclass != 1) return false;
    else if (ComputeWall(0, ELECTRON) <= 50) return false;
    else if (ComputeTowall(0, ELECTRON) <= 270) return false;
    else if (electron_momentum() <= 30) return false;
    else if (fqmrnring[0] != 1) return false;
    else if (electron_muon_PID() < 0) return false;
    else if (electron_momentum() <= 100) return false;
    else if (fqnse != 2) return false;
    else if (ComputeErecCCDel(0, ELECTRON) >= 1250) return false;
    else if (pi0_electron_PID() >= 0) return false;
    else return true;
  }

  bool is1Rmu(){
    if (evclass != 1) return false;
    else if (ComputeWall(0, MUON) <= 50) return false;
    else if (ComputeTowall(0, MUON) <= 250) return false;
    else if (electron_momentum() <= 30) return false;
    else if (fqmrnring[0] != 1) return false;
    else if (electron_muon_PID() >= 0) return false;
    else if (muon_momentum() <= 200) return false;
    else if ((fqnse < 1) or (fqnse > 2)) return false;
    else if (piplus_muon_PID() >= 0) return false;
    else return true;
  }

} // end t2ksk namespace

#endif
