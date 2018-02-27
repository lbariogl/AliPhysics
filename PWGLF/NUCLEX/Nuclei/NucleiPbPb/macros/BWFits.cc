#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <Riostream.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TMath.h>
#include "src/Common.h"
#include "YieldMean.C"
#include "AliPWGFunc.h"
#include "AliPID.h"
#include "AdditionalFunctions.h"

static TF1 *fBGBlastWave_Integrand = NULL;
static TF1 *fBGBlastWave_Integrand_num = NULL;
static TF1 *fBGBlastWave_Integrand_den = NULL;

constexpr double kParticleMass = 1.87561;
constexpr int kNfitFunctions = 4;
const string kFitFunctionNames[kNfitFunctions] = {"BlastWave", "Boltzmann", "LevyTsallis", "Mt-exp"};

//Levi-Tsallis parameters
const Float_t normal=2e-4, normMin=1e-5, normMax=1.;
const Float_t n=7, nMin=2, nMax=100;
const Float_t C=0.2, CMin=0.01, CMax=0.6;

double BGBlastWave_Integrand(const double *x, const double *p) {

  /*
     x[0] -> r (radius)
     p[0] -> mT (transverse mass)
     p[1] -> pT (transverse momentum)
     p[2] -> beta_max (surface velocity)
     p[3] -> T (freezout temperature)
     p[4] -> n (velocity profile)
     */

  double r = x[0];
  double mt = p[0];
  double pt = p[1];
  double beta_max = p[2];
  double temp_1 = 1. / p[3];
  double n = p[4];

  double beta = beta_max * TMath::Power(r, n);
  if (beta > 0.9999999999999999) beta = 0.9999999999999999;
  double rho = TMath::ATanH(beta);
  double argI0 = pt * TMath::SinH(rho) * temp_1;
  if (argI0 > 700.) argI0 = 700.;
  double argK1 = mt * TMath::CosH(rho) * temp_1;
  //  if (argI0 > 100 || argI0 < -100)
  //    printf("r=%f, pt=%f, beta_max=%f, temp=%f, n=%f, mt=%f, beta=%f, rho=%f, argI0=%f, argK1=%f\n", r, pt, beta_max, 1. / temp_1, n, mt, beta, rho, argI0, argK1);
  return r * mt * TMath::BesselI0(argI0) * TMath::BesselK1(argK1);
}

double BGBlastWave_Func(const double *x, const double *p) {
  /* dN/dpt */

  double pt = x[0];
  double mass = p[0];
  double mt = TMath::Sqrt(pt * pt + mass * mass);
  double beta_max = p[1];
  double temp = p[2];
  double n = p[3];
  double norm = p[4];

  if (!fBGBlastWave_Integrand)
    fBGBlastWave_Integrand = new TF1("fBGBlastWave_Integrand", BGBlastWave_Integrand, 0., 1., 5);
  fBGBlastWave_Integrand->SetParameters(mt, pt, beta_max, temp, n);
  double integral = fBGBlastWave_Integrand->Integral(0., 1.);
  return norm * pt * integral;
}

TF1 * BGBlastWave(const char *name, double mass, double beta_max = 0.9, double temp = 0.1,
    double n = 1., double norm = 1.e6) {

  TF1 *fBGBlastWave = new TF1(name, BGBlastWave_Func, 0., 10., 5);
  fBGBlastWave->SetParameters(mass, beta_max, temp, n, norm);
  fBGBlastWave->SetParNames("mass", "beta_max", "T", "n", "norm");
  fBGBlastWave->FixParameter(0, mass);
  fBGBlastWave->SetParLimits(1, 0.5, 0.94);
  fBGBlastWave->SetParLimits(2, 0.09, 0.250);
  fBGBlastWave->SetParLimits(3, 0.09, 3.1);
  return fBGBlastWave;
}

void Denormalize(TH1 * h) {
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    h->SetBinContent(i,h->GetBinContent(i) * TMath::TwoPi() * h->GetBinCenter(i));
    h->SetBinError(i,h->GetBinError(i) * TMath::TwoPi() * h->GetBinCenter(i));
  }
}

void BWFits() {
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  gStyle->SetTitleXOffset(1.3);
  gStyle->SetTitleYOffset(1.6);
  TFile *mineF = TFile::Open(kFinalOutput.data());
  TCanvas *mResults[kCentLength];
  TCanvas *mCanvM[kCentLength],*mCanvA[kCentLength];
  TH1D * mysystM[kCentLength], *mysystA[kCentLength],*speM[kCentLength],*speA[kCentLength];
  TF1 *bw[3];

  AliPWGFunc pwgfunc;
  pwgfunc.SetVarType(AliPWGFunc::kdNdpt);
  TF1* fit_functions[kNfitFunctions] = {
    BGBlastWave(kFitFunctionNames[0].data(),kParticleMass,0.623,0.084,0.76),
    pwgfunc.GetBoltzmann(kParticleMass, 0.1, 1, kFitFunctionNames[1].data()),
    LevyTsallis(kFitFunctionNames[2].data(), kParticleMass),
    pwgfunc.GetMTExp(kParticleMass, 0.1, 1, kFitFunctionNames[3].data())
  };

  TFile bwfile(kBWfitsOutput.data(),"recreate");
  TDirectory* datadir = bwfile.mkdir("data");
  TDirectory* function_dir[4]{nullptr};
  for (int iF = 0; iF < kNfitFunctions; ++iF)
    function_dir[iF] = bwfile.mkdir(kFitFunctionNames[iF].data());
  for (int iC = 0; iC < kCentLength; ++iC) {
    speM[iC] = (TH1D*)mineF->Get(Form("deuterons/%i/stat_all",iC));
    mysystM[iC] = (TH1D*)mineF->Get(Form("antideuterons/%i/syst_all",iC));
    if (!mysystM[iC]) cout << "Missing " << Form("syst%i",iC) << endl;
    if (!speM[iC]) cout << "Missing " << Form("syst%i",iC) << endl;
    if (!speM[iC] || !mysystM[iC]) return;
    TH1D *m = new TH1D(Form("m%i",iC),";#it{p}_{T} (GeV / #it{c}); 1/#it{N}_{ev} d^{2}#it{N}/d#it{p}_{T}d#it{y} (GeV/#it{c})^{-1}",speM[iC]->GetNbinsX(), speM[iC]->GetXaxis()->GetXbins()->GetArray());
    TH1D *sm = new TH1D(Form("sm%i",iC),";#it{p}_{T} (GeV/#it{c});1/(2#pi#it{p}_{T}) 1/#it{N}_{ev} d^{2}#it{N}/d#it{p}_{T}d#it{y} (GeV/#it{c})^{-2}",speM[iC]->GetNbinsX(), speM[iC]->GetXaxis()->GetXbins()->GetArray());

    for (int iB = 1; iB <= speM[iC]->GetNbinsX(); ++iB) {
      double x = speM[iC]->GetBinCenter(iB);
      if (x < 0.6 || x > kCentPtLimits[iC]) continue;
      float stat = speM[iC]->GetBinError(speM[iC]->FindBin(m->GetBinCenter(iB)));
      float syst = mysystM[iC]->GetBinError(speM[iC]->FindBin(m->GetBinCenter(iB)));
      m->SetBinContent(iB,speM[iC]->GetBinContent(speM[iC]->FindBin(m->GetBinCenter(iB))));
      m->SetBinError(iB,stat);
      sm->SetBinContent(iB,speM[iC]->GetBinContent(speM[iC]->FindBin(m->GetBinCenter(iB))));
      sm->SetBinError(iB,syst);
    }

    for (int iF = 0; iF < kNfitFunctions; ++iF) {
      if (!iF && iC == 1) {
        fit_functions[iF]->SetParLimits(1, 0.04, 0.94);
        fit_functions[iF]->SetParLimits(2, 0.01, 0.250);
        fit_functions[iF]->SetParLimits(3, 0.0, 5.);
      } else if (!iF && iC == 2) {
        fit_functions[iF]->SetParLimits(1, 0.4, 0.94);
        fit_functions[iF]->SetParLimits(2, 0.10, 0.250);
        fit_functions[iF]->SetParLimits(3, 0.09, 2.);
      }
      function_dir[iF]->cd();
		  TH1* h = YieldMean(m,sm,fit_functions[iF],0,10.1);
      cout << "\n*****************************" << endl;
      for (int iP = 0; iP < fit_functions[iF]->GetNpar(); ++iP) {
        cout << fit_functions[iF]->GetParName(iP) << ": " << fit_functions[iF]->GetParameter(iP);
        cout << " +/- " << fit_functions[iF]->GetParError(iP) << endl;
      }
      cout << "*****************************" << endl << endl;
      fit_functions[iF]->Write(Form("%s%i",kFitFunctionNames[iF].data(),iC));
      h->Write(Form("result%i",iC));
    }
    datadir->cd();
    m->Write(Form("stat%i",iC));
    sm->Write(Form("syst%i",iC));
  }

  bwfile.Close();

}
