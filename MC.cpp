
#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TStyle.h>


void MC(){
TFile *f_thmu_NLO_1hit = TFile::Open("GM_NLOvsLO/NOoutchi2_mu_eff_NLO_02cut_1hit_clone.root");
TFile *f_thmu_NLO_2hit = TFile::Open("GM_NLOvsLO/NOoutchi2_mu_eff_NLO_02cut_2hit_clone.root");

TFile *f_the_NLO_1hit = TFile::Open("GM_NLOvsLO/NOoutchi2_el_eff_NLO_02cut_1hit_clone.root");
TFile *f_the_NLO_2hit = TFile::Open("GM_NLOvsLO/NOoutchi2_el_eff_NLO_02cut_2hit_clone.root");

TFile *f_op_NLO_1hit = TFile::Open("GM_NLOvsLO/NOoutchi2_op_eff_NLO_02cut_1hit_clone.root");
TFile *f_op_NLO_2hit = TFile::Open("GM_NLOvsLO/NOoutchi2_op_eff_NLO_02cut_2hit_clone.root");


TH1::SetDefaultSumw2(kTRUE);


TH1D *h_theta_mu_NLO_1hit=(TH1D*)f_thmu_NLO_1hit->Get("theta_mu_clone");
TH1D *h_theta_mu_NLO_2hit=(TH1D*)f_thmu_NLO_2hit->Get("theta_mu_clone");

TH1D *h_theta_e_NLO_1hit=(TH1D*)f_the_NLO_1hit->Get("theta_e_clone");
TH1D *h_theta_e_NLO_2hit=(TH1D*)f_the_NLO_2hit->Get("theta_e_clone");

TH1D *h_op_NLO_1hit=(TH1D*)f_op_NLO_1hit->Get("h_opening_clone");
TH1D *h_op_NLO_2hit=(TH1D*)f_op_NLO_2hit->Get("h_opening_clone");


 auto legend_e = new TLegend(0.75,0.15,0.9,0.3);
   legend_e->AddEntry(h_theta_e_NLO_1hit,"NLO 1 hit shared","LEP");
   legend_e->AddEntry(h_theta_e_NLO_2hit,"NLO 2 hits shared","LEP");

 auto legend_mu = new TLegend(0.1,0.15,0.3,0.30);
   legend_mu->AddEntry(h_theta_mu_NLO_1hit,"NLO 1 hit shared","LEP");
   legend_mu->AddEntry(h_theta_mu_NLO_2hit,"NLO 2 hits shared","LEP");

 auto legend_op = new TLegend(0.1,0.15,0.3,0.30);
   legend_op->AddEntry(h_op_NLO_1hit,"NLO 1 hit shared","LEP");
   legend_op->AddEntry(h_op_NLO_2hit,"NLO 2 hits shared","LEP");


TCanvas a1("a1","a1",1000,700);
a1.Divide(1,3);

a1.cd(1);

h_theta_mu_NLO_1hit->SetLineColor(kOrange+10);
h_theta_mu_NLO_2hit->SetLineColor(kAzure+7);
h_theta_mu_NLO_2hit->SetMinimum(0.90);
h_theta_mu_NLO_2hit->Draw("E");
h_theta_mu_NLO_1hit->Draw("E same");
legend_mu->Draw();
gStyle->SetOptStat(0);

a1.cd(2);
h_theta_e_NLO_1hit->SetLineColor(kOrange+10);
h_theta_e_NLO_2hit->SetLineColor(kAzure+7);
h_theta_e_NLO_2hit->SetMinimum(0.82);
h_theta_e_NLO_2hit->Draw("E");
h_theta_e_NLO_1hit->Draw("E same");
legend_e->Draw();
gStyle->SetOptStat(0);

a1.cd(3);

h_op_NLO_1hit->SetLineColor(kOrange+10);
h_op_NLO_2hit->SetLineColor(kAzure+7);
h_op_NLO_2hit->Draw("E");
h_op_NLO_1hit->Draw("E same");
legend_op->Draw();
gStyle->SetOptStat(0);
a1.SaveAs("GM_NLOvsLO/NOoutchi2_compare_NLO_2hit_NLO_1hit_02cut_clone.pdf");


}


