
string baseDirectory = "/sps/t2k/ablanche/work/results/xsLLhFitter/BANFF_Fit";
string treeFile = "/tree_converter/Converted_MC_run1to8_FHC_only.root";

TCanvas* c;

map<int, string> sampleNames;
map<int, string> reactionNames;
map<int, int> reactionColors;


TH1D* hD1;
TH1D* hD2;
map<string, TH1D*> histogramMap;
map<string, THStack*> histogramStackMap;


void createContainerHistograms();
void fillParameters();


void generateSamplesPlots(){

  createContainerHistograms();
  fillParameters();

  TFile* treeConverted = TFile::Open((baseDirectory+treeFile).c_str());
  TTree* selectedEvents = (TTree*) treeConverted->Get("selectedEvents");

  for(int i_fgd = 0 ; i_fgd < 2 ; i_fgd++){

    for(const auto& sample : sampleNames){

      c->cd(1 + i_fgd*sampleNames.size() + sample.first);
      string stackName = Form("FGD%i_Sample%i", i_fgd+1, sample.first);
      histogramStackMap[ stackName ] = new THStack( stackName.c_str(), stackName.c_str() );

      for(const auto& reaction : reactionNames){

        string histName = Form("FGD%i_Sample%i_Reaction%i", i_fgd+1, sample.first, reaction.first);
        cout << "Processing: " << histName << endl;

        int nbEvents = selectedEvents->Draw(
          "D1Reco>>hD1",
          Form(
            "fgd == %i && cut_branch == %i && reaction == %i",
            i_fgd, sample.first, reaction.first
          ),
          "goff"
        );
        if(nbEvents == 0) continue;
        for(int iBin = 0 ; iBin < hD1->GetNbinsX() ; iBin++){
          if(hD1->GetBinWidth(iBin) == 0) continue;
          hD1->SetBinContent(iBin,
            hD1->GetBinContent(iBin)/hD1->GetBinWidth(iBin)
          );
        }
        hD1->GetXaxis()->SetTitle("p_{#mu} (MeV/c)");
        hD1->GetYaxis()->SetTitle("Events/(1 MeV/c)");
        hD1->GetXaxis()->SetRangeUser(0,2000);
        hD1->SetFillColor(reactionColors[reaction.first]);
        histogramMap[ histName ] = hD1;
        histogramStackMap[ stackName ]->Add(histogramMap[ histName ]);

      }

      cout << " > Drawing: " << stackName << " -> " << 1 + i_fgd*sampleNames.size() + sample.first << endl;
      histogramStackMap[stackName]->Draw("stack");
      c->Update();

    }

  }



}

void fillParameters(){

 c = new TCanvas("c", "c", 800,1200);
 c->Divide(2,3);

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

void createContainerHistograms(){

  std::vector<double> D1binning;
  std::vector<double> D2binning;
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
  hD1 = new TH1D("hD1", "hD1", D1binning.size() - 1, &D1binning[0]);
  hD2 = new TH1D("hD2", "hD2", D2binning.size() - 1, &D2binning[0]);

}
