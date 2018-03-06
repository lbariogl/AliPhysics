#include "src/Common.h"
#include "src/Utils.h"
using namespace utils;
#include "src/Plotting.h"
using namespace plotting;

#include <TFile.h>
#include <TDirectory.h>
#include <TLegend.h>
#include <TList.h>
#include <TCanvas.h>
#include <TH1F.h>

const char* kMultLab[kCentLength-1] = {"0-1","1-5","5-10","10-20","20-30","30-40","40-60","60-80","80-100"};
const double pt_bin_limits[16] = {0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6,1.8,2.0,2.2,2.6,3.0,3.4,3.8};//,4.4};
const int n_pt_bins = 15;
const int kNspecies = 9;
const char* kSpeciesName[kNspecies]= {"Pi","Ka","Pr","Phi","Ks","K0","Lambda","Xi","Om"};
const char* kSpeciesLabels[kNspecies] = {"#pi", "k", "p", "#phi", "K_{s}", "K_{0}", "#Lambda", "#Xi", "#Omega"};

void SignalLoss(){

  TFile input(Form("%sSignalLoss_EvCut.root",kBaseOutputDir.data()));
  TFile output(Form("%sSignalLoss_Correction.root",kBaseOutputDir.data()),"recreate");

  TCanvas* cSignalLoss[kCentLength];

  for(int iC=0; iC<kCentLength-1; iC++){

    TDirectory *c_dir = output.mkdir(Form("%d",iC));

    cSignalLoss[iC] = new TCanvas(Form("cSL_%d",iC),Form("cSL_%d",iC));
    TLegend* leg = new TLegend(0.1,0.1,0.5,0.3);
    leg->SetBorderSize(0);

    for(int iS=0; iS<3/*kNspecies*/; iS++){
      TList* lSp = (TList*)input.Get(kSpeciesName[iS]);
      TH1F* hSp_tmp = (TH1F*)lSp->FindObject(Form("sgnLoss_%s_%s",kSpeciesName[iS],kMultLab[iC]));
      Requires(hSp_tmp,Form("Missing %s",Form("%s/sgnLoss_%s_%s",kSpeciesName[iS],kSpeciesName[iS],kMultLab[iC])));
      TH1F* hSp = new TH1F(Form("h%s_%d",kSpeciesName[iS],iC),";#it{p}_{T} (GeV/#it{c});sigLoss #times #epsilon #times Acc",n_pt_bins,pt_bin_limits);
      int iSpBin = 1;
      float content_tmp[n_pt_bins] = {0.};
      float err2_tmp[n_pt_bins] = {0.};
      int count_tmp[n_pt_bins] = {0};
      for(int iB=1; iB<=hSp_tmp->GetNbinsX(); iB++){
        if(hSp_tmp->GetBinCenter(iB) < pt_bin_limits[0]) continue;
        if(hSp_tmp->GetBinCenter(iB) > kCentPtLimits[iC]) break;
        if(hSp_tmp->GetBinCenter(iB) < pt_bin_limits[iSpBin]){
          //printf("iB: %d hUpEdge: %f iCont: %d contLin: %f\n",iB,hSp_tmp->GetBinLowEdge(iB),iSpBin,pt_bin_limits[iSpBin]);
          content_tmp[iSpBin-1] += hSp_tmp->GetBinContent(iB);
          err2_tmp[iSpBin-1] += Sq(hSp_tmp->GetBinError(iB));
          count_tmp[iSpBin-1]++;
        }
        else{
          iSpBin++;
          content_tmp[iSpBin-1] += hSp_tmp->GetBinContent(iB);
          err2_tmp[iSpBin-1] += Sq(hSp_tmp->GetBinError(iB));
          count_tmp[iSpBin-1]++;
        }
      }
      for(int i=0; i<n_pt_bins; i++){
        if(hSp->GetBinCenter(i+1)>kCentPtLimits[iC]) continue;
        hSp->SetBinContent(i+1,content_tmp[i]/count_tmp[i]);
        hSp->SetBinError(i+1,TMath::Sqrt(err2_tmp[i])/count_tmp[i]);
      }
      hSp->SetDirectory(0);
      hSp->SetTitle(Form("%s%%",kMultLab[iC]));
      hSp->GetYaxis()->SetRangeUser(0.7,1.3);
      SetHistStyle(hSp,kSpectraColors[iC],20+iS);
      leg->AddEntry(hSp,kSpeciesLabels[iS],"PE");
      cSignalLoss[iC]->cd();
      if(iS==0) hSp->Draw();
      else hSp->Draw("same");
    }
    leg->Draw();
    c_dir->cd();
    cSignalLoss[iC]->Write();
  }
}
