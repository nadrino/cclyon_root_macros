
string baseDirectory = "/sps/t2k/ablanche/work/results/xsLLhFitter/BANFF_Fit";
string treeFile = "/tree_converter/Converted_MC_run1to8_FHC_only.root";
// string treeFile = "/tree_converter/Converted_MC_run1to9.root";
// string treeFile = "/tree_converter/Converted_NumuCCMultiPiAnalysis_EXAMPLE.root";

TCanvas* c;

map<int, string> sampleNames;
map<int, string> reactionNames;
map<int, int> reactionColors;

std::vector<double> D1binning;
std::vector<double> D2binning;

map<string, TH1D*> histogramMap;
map<string, THStack*> histogramStackMap;


void fillBinnings();
void fillParameters();

void generateSamplesPlots(){

  fillBinnings();
  fillParameters();

  TFile* treeConverted = TFile::Open((baseDirectory+treeFile).c_str());
  TTree* selectedEvents = (TTree*) treeConverted->Get("selectedEvents");

  cout << "selectedEvents -> " << selectedEvents << " entries" << endl;

  for(int i_fgd = 0 ; i_fgd < 2 ; i_fgd++){

    for(const auto& sample : sampleNames){

      c->cd(1 + i_fgd*sampleNames.size() + sample.first);
      string stackName = Form("FGD%i_Sample%i", i_fgd+1, sample.first);
      histogramStackMap[ stackName ] = new THStack( stackName.c_str(), stackName.c_str() );

      vector<TH1D*> histManualStack;
      for(const auto& reaction : reactionNames){

        string histName = Form("FGD%i_Sample%i_Reaction%i", i_fgd+1, sample.first, reaction.first);

        auto* histTemp = new TH1D(
          histName.c_str(),
          reaction.second.c_str(),
          D1binning.size() - 1, &D1binning[0]);

        int nbEvents = selectedEvents->Draw(
          ("D1Reco>>" + histName).c_str(),
          Form(
            "fgd == %i && cut_branch == %i && topology == %i && beammode == 1",
            i_fgd, sample.first, reaction.first
          ),
          "goff"
        );
        if(nbEvents == 0) continue;
        cout << "Processing: " << histName << endl;
        for(int iBin = 0 ; iBin < histTemp->GetNbinsX() ; iBin++){
          if(histTemp->GetBinWidth(iBin) == 0) continue;
          histTemp->SetBinContent(iBin,
            histTemp->GetBinContent(iBin)/histTemp->GetBinWidth(iBin)
          );
        }
        histTemp->SetLineWidth(0);
        histTemp->GetXaxis()->SetTitle("p_{#mu} (MeV/c)");
        histTemp->GetYaxis()->SetTitle("Events/(1 MeV/c)");
        histTemp->GetXaxis()->SetRangeUser(0,2000);
        histTemp->GetXaxis()->SetLimits(0,30000);
        histTemp->SetFillColor(reactionColors[reaction.first]);
        histogramMap[ histName ] = histTemp;

        if(histManualStack.size() == 0){
          histManualStack.emplace_back(histogramMap[ histName ]);
          // histManualStack.back()->Draw();
        }
        else{
          for(int iBin = 0 ; iBin < histogramMap[ histName ]->GetNbinsX() ; iBin++){
            histogramMap[ histName ]->SetBinContent(iBin,
              histogramMap[ histName ]->GetBinContent(iBin) +
              histManualStack.back()->GetBinContent(iBin)
            );
          }
          // histogramMap[ histName ]->Merge(histManualStack.back());
          histManualStack.emplace_back(histogramMap[ histName ]);
          // histManualStack.back()->Draw("SAME");
        }

        // histogramStackMap[ stackName ]->Add(histogramMap[ histName ]);

      }

      cout << " >> Drawing..." << endl;
      string nameBuffer;
      for(int iHist = histManualStack.size()-1 ; iHist >= 0 ; iHist--){
        string opt = "SAME";
        if(iHist == histManualStack.size()-1){
          opt= "";
          nameBuffer = histManualStack[iHist]->GetTitle();
          histManualStack[iHist]->SetTitle(
            Form("FGD%i #nu_{#mu} %s", i_fgd+1, sample.second.c_str())
          );
        }
        histManualStack[iHist]->Draw(opt.c_str());
      }

      gPad->BuildLegend();

      histManualStack.back()->SetTitle(
        Form("FGD%i #nu_{#mu} %s", i_fgd+1, sample.second.c_str())
      );
      c->Update();


      // cout << " > Drawing: " << stackName << " -> " << 1 + i_fgd*sampleNames.size() + sample.first << endl;
      // histogramStackMap[stackName]->Draw("stack");
      // c->Update();

    }

  }



}

void fillParameters(){

 c = new TCanvas("c", "c", 1200,700);
 c->Divide(3,2);

  sampleNames[0] = "CC-0#pi";
  sampleNames[1] = "CC-1#pi";
  sampleNames[2] = "CC-Other";

  reactionNames[0] = "CCQE";
  reactionNames[9] = "2p2h";
  reactionNames[1] = "RES";
  reactionNames[2] = "DIS";
  reactionNames[3] = "COH";
  reactionNames[4] = "NC";
  reactionNames[5] = "CC-#bar{#nu}_{#mu}";
  reactionNames[6] = "CC-#nu_{e}, CC-#bar{#nu}_{e}";
  reactionNames[999] = "other";
  reactionNames[7] = "out FV";
  reactionNames[-1] = "no truth";
  reactionNames[777] = "sand #mu";

  reactionColors[0] = 2;
  reactionColors[9] = 874;
  reactionColors[1] = 3;
  reactionColors[2] = 4;
  reactionColors[3] = 7;
  reactionColors[4] = 6;
  reactionColors[5] = 31;
  reactionColors[6] = 65;
  reactionColors[999] = 48;
  reactionColors[7] = 1;
  reactionColors[-1] = 92;
  reactionColors[777] = 51;

}

void fillBinnings(){

  vector<string> binningLines = GenericToolbox::dumpFileAsVectorString(
    baseDirectory + "/binning/BANFF_binning.txt"
  );
  for(const auto& line : binningLines){
    auto lineElements = GenericToolbox::splitString(
      GenericToolbox::removeExtraDoubledCharacters(line, " "),
      " "
    );
    if(lineElements.size() != 4) continue;
    if(D2binning.empty()) D2binning.emplace_back(stod(lineElements[0]));
    if(D1binning.empty()) D1binning.emplace_back(stod(lineElements[2]));
    D2binning.emplace_back(stod(lineElements[1]));
    D1binning.emplace_back(stod(lineElements[3]));
  }
  sort(D1binning.begin(), D1binning.end());
  sort(D2binning.begin(), D2binning.end());

}
