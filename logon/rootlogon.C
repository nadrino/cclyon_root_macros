{
  cout << "Using ROOT logon : $REPO_DIR/cclyon_root_macros/logon/rootlogon.C " << endl;

  cout << "ROOT version is: " << ROOT_RELEASE << std::endl;

  if(gROOT->GetVersionInt() > 60000){
    std::string thisFolderPath(std::string(__FILE__).substr(0, std::string(__FILE__).rfind("/")));

    std::cout << "Loading GenericToolbox namespace..." << std::endl;
    gROOT->ProcessLine( Form(".include %s/../submodules/cpp-generic-toolbox/include", thisFolderPath.c_str()) );
    // std::cout << "Processing \"GenericToolbox.h\"..." << std::endl;
    // gROOT->ProcessLine("#include \"GenericToolbox.h\"");
    std::cout << "Processing \"GenericToolbox.Root.h\"..." << std::endl;
    gROOT->ProcessLine("#include \"GenericToolbox.Root.h\"");
    // std::cout << "Processing \"GenericToolbox.ThreadPool.h\"..." << std::endl;
    // gROOT->ProcessLine("#include \"GenericToolbox.ThreadPool.h\"");
    std::cout << "Processing \"GenericToolbox.TablePrinter.h\"..." << std::endl;
    gROOT->ProcessLine("#include \"GenericToolbox.TablePrinter.h\"");

    // std::cout << "Loading Simple Logger..." << std::endl;
    // #define LOGGER_PREFIX_LEVEL 3
    gROOT->ProcessLine( Form(".include %s/../submodules/simple-cpp-logger/include", thisFolderPath.c_str()) );
    // std::cout << "Processing \"Logger.h\"..." << std::endl;
    // gROOT->ProcessLine("#include \"Logger.h\"");
  }

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

}

void loadOldToolbox(){

  std::cout << "Loading old TToolBox namespace..." << std::endl;
  gROOT->ProcessLine("#include \"rootlogon.h\"");

}
