
#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TStyle.h>


void sixVSone(){
TFile *f_thmu_bs = TFile::Open("mc_NNLOvsLO/100k_gen_theta_mu_gen_full_bigsample.root");
TFile *f_the_bs = TFile::Open("mc_NNLOvsLO/100k_gen_theta_e_gen_full_bigsample.root");

TFile *f_thmu = TFile::Open("mc_NNLOvsLO/100k_gen_theta_mu_gen_full.root");
TFile *f_the = TFile::Open("mc_NNLOvsLO/100k_gen_theta_e_gen_full.root");


TH1::SetDefaultSumw2(kTRUE);


TH1D *h_theta_mu_bs=(TH1D*)f_thmu_bs->Get("theta_mu_gen");
TH1D *h_theta_mu=(TH1D*)f_thmu->Get("theta_mu_gen");

TH1D *h_theta_e_bs=(TH1D*)f_the_bs->Get("theta_e_gen");
TH1D *h_theta_e=(TH1D*)f_the->Get("theta_e_gen");


 auto legend_e = new TLegend();
   legend_e->AddEntry(h_theta_e_bs,"big sample","LEP");
   legend_e->AddEntry(h_theta_e,"six samples","LEP");

 auto legend_mu = new TLegend();
   legend_mu->AddEntry(h_theta_mu_bs,"big sample","LEP");
   legend_mu->AddEntry(h_theta_mu,"six samples","LEP");

TCanvas a1("a1","a1",1000,700);
a1.Divide(1,2);

a1.cd(1);
//h_theta_mu_bs->Scale(1./4.73665e+07);
h_theta_mu_bs->Scale(1./5.29395e+06);
h_theta_mu_bs->SetLineColor(kOrange+10);

//h_theta_mu->Scale(1./4.34487e+07);
h_theta_mu->Scale(1./4.46111e+06);
h_theta_mu->SetLineColor(kAzure+7);
h_theta_mu_bs->Draw("E");
h_theta_mu->Draw("E same");
legend_mu->Draw();
gStyle->SetOptStat(0);
gPad->SetLogy();

a1.cd(2);

//h_theta_e_bs->Scale(1./4.73665e+07);
h_theta_e_bs->Scale(1./5.29395e+06);
h_theta_e_bs->SetLineColor(kOrange+10);

//h_theta_e->Scale(1./4.34487e+07);
h_theta_e->Scale(1./4.46111e+06);
h_theta_e->SetLineColor(kAzure+7);
h_theta_e_bs->Draw("E");
h_theta_e->Draw("E same");
legend_e->Draw();
gStyle->SetOptStat(0);

a1.SaveAs("mc_NNLOvsLO/compare_full_bs_100k.pdf");

cout << "Electron six samples: " << endl;
for(int b=1; b<7; b++){ cout << b << ") " << h_theta_e->GetBinContent(b) << " +- " << h_theta_e->GetBinError(b) << endl;}
cout << "Electron single sample: " << endl;
for(int b=1; b<7; b++){ cout << b << ") " << h_theta_e_bs->GetBinContent(b) << " +- " << h_theta_e_bs->GetBinError(b) << endl;}
}


