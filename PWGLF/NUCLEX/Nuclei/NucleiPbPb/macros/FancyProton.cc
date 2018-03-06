#include "src/Common.h"
#include "src/Plotting.h"
using namespace plotting;

#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TStyle.h>

const int kKnownMult = 6;

const char* kMultNameInput[7] = {"0to1","1to5","5to10","10to15","15to20","20to30","30to40"};
const int join[kKnownMult][2] = {{0,0},{1,1},{2,2},{3,4},{5,5},{6,6}};
const double pt_bin_limits[16] = {0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6,1.8,2.0,2.2,2.6,3.0,3.4,3.8};//,4.4};
const int n_pt_bins = 15;


void FancyProton(){

  TFile input_file(Form("%sFinal_combined_spectra_pp13TeV.root",kBaseOutputDir.data()));
  TH1F* hInputStat[7];
  TH1F* hInputSyst[7];
  for(int iInput=0; iInput<7; iInput++){
    hInputStat[iInput]=(TH1F*)input_file.Get(Form("hCombinedTPCTOFTOF_ITSsa_Pr_%s_stat",kMultNameInput[iInput]));
    hInputStat[iInput]->SetDirectory(0);
    hInputStat[iInput]->GetXaxis()->SetRangeUser(0.6,3.8);
    hInputSyst[iInput]=(TH1F*)input_file.Get(Form("hCombinedTPCTOFTOF_ITSsa_Pr_%s_syst",kMultNameInput[iInput]));
    hInputSyst[iInput]->SetDirectory(0);
    hInputSyst[iInput]->GetXaxis()->SetRangeUser(0.6,3.8);
  }

  TFile output(Form("%sFancyProton.root",kBaseOutputDir.data()),"recreate");
  TH1F* hOutputStat[kKnownMult];
  TH1F* hOutputSyst[kKnownMult];
  for(int iC=0; iC<kKnownMult; iC++){
    int n_join = join[iC][1]-join[iC][0]+1;
    hOutputStat[iC] = (TH1F*)hInputStat[join[iC][0]]->Clone(Form("hProton_stat_%d",iC));
    hOutputStat[iC]->Reset();
    for(int iB=1; iB<=hOutputStat[iC]->GetNbinsX(); iB++){
      if(hOutputStat[iC]->GetBinCenter(iB) < 0.3 || hOutputStat[iC]->GetBinCenter(iB) > 4.) continue;
      hOutputStat[iC]->SetBinContent(iB,hInputStat[iC]->GetBinContent(iB));
      hOutputStat[iC]->SetBinError(iB,hInputStat[iC]->GetBinError(iB));
    }
    hOutputSyst[iC] = (TH1F*)hInputSyst[join[iC][0]]->Clone(Form("hProton_syst_%d",iC));
    hOutputSyst[iC]->Reset();
    for(int iB=1; iB<=hOutputSyst[iC]->GetNbinsX(); iB++){
      if(hOutputSyst[iC]->GetBinCenter(iB) < 0.3 || hOutputSyst[iC]->GetBinCenter(iB) > 4.) continue;
      hOutputSyst[iC]->SetBinContent(iB,hInputSyst[iC]->GetBinContent(iB));
      hOutputSyst[iC]->SetBinError(iB,hInputSyst[iC]->GetBinError(iB));
    }
    if(n_join>1){
      for(int i=join[iC][0]+1; i<=join[iC][1];i++){
        hOutputSyst[iC]->Add(hInputStat[i]);
        hOutputStat[iC]->Add(hInputSyst[i]);
      }
      for(int iB=1; iB<=hOutputSyst[iC]->GetNbinsX();iB++){
        hOutputStat[iC]->SetBinContent(iB,hOutputStat[iC]->GetBinContent(iB)/n_join);
        hOutputSyst[iC]->SetBinContent(iB,hOutputSyst[iC]->GetBinContent(iB)/n_join);
      }
    }
    SetHistStyle(hOutputStat[iC],kSpectraColors[iC]);
    hOutputStat[iC]->SetTitle(";#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{ev}} #frac{d^{2}#it{N}}{d#it{p}_{T}d#it{y}}");
    SetHistStyle(hOutputSyst[iC],kSpectraColors[iC]);
    hOutputSyst[iC]->SetTitle(";#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{ev}} #frac{d^{2}#it{N}}{d#it{p}_{T}d#it{y}}");
    hOutputStat[iC]->Write();
    hOutputSyst[iC]->Write();
  }
}
