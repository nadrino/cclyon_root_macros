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

}

namespace TToolBox {

  static double last_displayed_value = -1;
  static int verbosity_level = 0;
  static Color_t colorCells[10] = {kOrange+1, kGreen-3, kTeal+3, kAzure+7, kCyan-2, kBlue-7, kBlue+2, kOrange+9, kRed+2, kPink+9};

  // Matrices/Vector Tools
  TH2D* get_TH2D_from_TMatrixD(TMatrixD *XY_values_, string graph_title_ = "", string Z_title_ = "",string Y_title_ = "Row #", string X_title_ = "Col #") {

    if(graph_title_.empty()) graph_title_ = XY_values_->GetTitle();

    auto* th2_histogram = new TH2D(graph_title_.c_str(), graph_title_.c_str(),
                                   XY_values_->GetNrows(), -0.5, XY_values_->GetNrows()-0.5,
                                   XY_values_->GetNcols(), -0.5, XY_values_->GetNcols()-0.5
    );

    for(int i_col = 0 ; i_col < XY_values_->GetNcols() ; i_col++){
      for(int j_row = 0 ; j_row < XY_values_->GetNrows() ; j_row++){
        th2_histogram->SetBinContent(i_col + 1, j_row + 1, (*XY_values_)[i_col][j_row]);
      }
    }

    th2_histogram->GetXaxis()->SetTitle(X_title_.c_str());
    th2_histogram->GetYaxis()->SetTitle(Y_title_.c_str());
    th2_histogram->GetZaxis()->SetTitle(Z_title_.c_str());

    return th2_histogram;

  }
  TMatrixD* generate_correlation_matrix(TMatrixD *covariance_matrix_) {
    auto* correlation_matrix = (TMatrixD*) covariance_matrix_->Clone();
    for(int i_row = 0 ; i_row < covariance_matrix_->GetNrows() ; i_row++){
      for(int i_col = 0 ; i_col < covariance_matrix_->GetNcols() ; i_col++){
        (*correlation_matrix)[i_row][i_col] /=
            TMath::Sqrt((*covariance_matrix_)[i_row][i_row] * (*covariance_matrix_)[i_col][i_col]);
      }
    }
    return correlation_matrix;
  }

}
