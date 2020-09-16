{

  cout << "Using ROOT logon : $REPO_DIR/cclyon_root_macros/logon/rootlogon.C " << endl;

  gROOT->SetStyle("Plain");
  gStyle->SetTitleBorderSize(0);
  gStyle->SetStatFont(42);
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(0);

  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetLabelOffset(0.015,"x");
  gStyle->SetLabelOffset(0.015,"y");

  gStyle->SetTitleFont(42);
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetLegendFont(42);

  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetTitleSize(0.06,"xyz");
  gStyle->SetTitleOffset(1.05,"x");
  gStyle->SetTitleOffset(1.20,"y");

  gStyle->SetStripDecimals(kFALSE);

  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadRightMargin(0.07);
  gStyle->SetPadTopMargin(0.09);

  // gStyle->SetStatW(0.35);
  // gStyle->SetStatH(0.25);

  gStyle->SetPadTickX(kTRUE);
  gStyle->SetPadTickY(kTRUE);

  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);

  //  gStyle->SetPalette(1);

  gStyle->SetLineWidth(2);
  gStyle->SetHistLineWidth(3);
  gStyle->SetFuncWidth(3);
  gStyle->SetFrameLineWidth(2);

  gStyle->SetMarkerSize(1.2);

  std::string REPO_DIR(getenv("REPO_DIR"));
  if(gROOT->GetVersionInt() > 60000){
    std::cout << "Loading GenericToolbox namespace..." << std::endl;
    std::string s = ".include " + REPO_DIR + "/cclyon_root_macros/submodules/cpp-generic-toolbox/include";
    // R__ADD_INCLUDE_PATH((REPO_DIR + "/cclyon_root_macros/submodules/cpp-generic-toolbox/include").c_str());
    // gROOT->LoadMacro("/Users/ablanche/Documents/Work/Repositories/cclyon_root_macros/submodules/cpp-generic-toolbox/include/GenericToolbox.h");
    gROOT->ProcessLine(s.c_str());
    gROOT->ProcessLine("#include \"GenericToolbox.h\"");
    gROOT->ProcessLine("#include \"GenericToolboxRootExt.h\"");
  }

  std::cout << "Loading old TToolBox namespace..." << std::endl;
  std::string s = ".include " + REPO_DIR + "/cclyon_root_macros/logon";
  // R__ADD_INCLUDE_PATH((REPO_DIR + "/cclyon_root_macros/submodules/cpp-generic-toolbox/include").c_str());
  // gROOT->LoadMacro("/Users/ablanche/Documents/Work/Repositories/cclyon_root_macros/submodules/cpp-generic-toolbox/include/GenericToolbox.h");
  gROOT->ProcessLine(s.c_str());
  gROOT->ProcessLine("#include \"rootlogon.h\"");

}
