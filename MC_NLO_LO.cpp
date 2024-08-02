
#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TStyle.h>


void MC_NLO_LO(){
TFile *f_thmu_NLO = TFile::Open("proposal/NOoutchi2_mu_eff_NLO_1hit.root");
TFile *f_thmu_LO = TFile::Open("proposal/NOoutchi2_mu_eff_LO_1hit.root");

TFile *f_the_NLO = TFile::Open("proposal/NOoutchi2_el_eff_NLO_1hit.root");
TFile *f_the_LO = TFile::Open("proposal/NOoutchi2_el_eff_LO_1hit.root");

TFile *f_op_NLO = TFile::Open("proposal/NOoutchi2_op_eff_NLO_1hit.root");
TFile *f_op_LO = TFile::Open("proposal/NOoutchi2_op_eff_LO_1hit.root");

TFile *f_thmu_single_NLO = TFile::Open("proposal/NOoutchi2_mu_single_eff_NLO_1hit.root");
TFile *f_thmu_single_LO = TFile::Open("proposal/NOoutchi2_mu_single_eff_LO_1hit.root");

TFile *f_the_single_NLO = TFile::Open("proposal/NOoutchi2_el_single_eff_NLO_1hit.root");
TFile *f_the_single_LO = TFile::Open("proposal/NOoutchi2_el_single_eff_LO_1hit.root");



TFile *f_thmu_NLO_clone = TFile::Open("proposal/NOoutchi2_mu_eff_NLO_1hit_clone.root");
TFile *f_thmu_LO_clone = TFile::Open("proposal/NOoutchi2_mu_eff_LO_1hit_clone.root");

TFile *f_the_NLO_clone = TFile::Open("proposal/NOoutchi2_el_eff_NLO_1hit_clone.root");
TFile *f_the_LO_clone = TFile::Open("proposal/NOoutchi2_el_eff_LO_1hit_clone.root");

TFile *f_op_NLO_clone = TFile::Open("proposal/NOoutchi2_op_eff_NLO_1hit_clone.root");
TFile *f_op_LO_clone = TFile::Open("proposal/NOoutchi2_op_eff_LO_1hit_clone.root");

TFile *f_thmu_single_NLO_clone = TFile::Open("proposal/NOoutchi2_mu_single_eff_NLO_1hit_clone.root");
TFile *f_thmu_single_LO_clone = TFile::Open("proposal/NOoutchi2_mu_single_eff_LO_1hit_clone.root");

TFile *f_the_single_NLO_clone = TFile::Open("proposal/NOoutchi2_el_single_eff_NLO_1hit_clone.root");
TFile *f_the_single_LO_clone = TFile::Open("proposal/NOoutchi2_el_single_eff_LO_1hit_clone.root");




TH1::SetDefaultSumw2(kTRUE);


TH1D *h_theta_mu_NLO=(TH1D*)f_thmu_NLO->Get("theta_mu");
TH1D *h_theta_mu_LO=(TH1D*)f_thmu_LO->Get("theta_mu");

TH1D *h_theta_e_NLO=(TH1D*)f_the_NLO->Get("theta_e");
TH1D *h_theta_e_LO=(TH1D*)f_the_LO->Get("theta_e");

TH1D *h_op_NLO=(TH1D*)f_op_NLO->Get("h_opening");
TH1D *h_op_LO=(TH1D*)f_op_LO->Get("h_opening");

TH1D *h_theta_mu_single_NLO=(TH1D*)f_thmu_single_NLO->Get("theta_mu_single");
TH1D *h_theta_mu_single_LO=(TH1D*)f_thmu_single_LO->Get("theta_mu_single");

TH1D *h_theta_e_single_NLO=(TH1D*)f_the_single_NLO->Get("theta_e_single");
TH1D *h_theta_e_single_LO=(TH1D*)f_the_single_LO->Get("theta_e_single");


TH1D *h_theta_mu_NLO_clone=(TH1D*)f_thmu_NLO_clone->Get("theta_mu_clone");
TH1D *h_theta_mu_LO_clone=(TH1D*)f_thmu_LO_clone->Get("theta_mu_clone");

TH1D *h_theta_e_NLO_clone=(TH1D*)f_the_NLO_clone->Get("theta_e_clone");
TH1D *h_theta_e_LO_clone=(TH1D*)f_the_LO_clone->Get("theta_e_clone");

TH1D *h_op_NLO_clone=(TH1D*)f_op_NLO_clone->Get("h_opening_clone");
TH1D *h_op_LO_clone=(TH1D*)f_op_LO_clone->Get("h_opening_clone");

TH1D *h_theta_mu_single_NLO_clone=(TH1D*)f_thmu_single_NLO_clone->Get("theta_mu_single_clone");
TH1D *h_theta_mu_single_LO_clone=(TH1D*)f_thmu_single_LO_clone->Get("theta_mu_single_clone");

TH1D *h_theta_e_single_NLO_clone=(TH1D*)f_the_single_NLO_clone->Get("theta_e_single_clone");
TH1D *h_theta_e_single_LO_clone=(TH1D*)f_the_single_LO_clone->Get("theta_e_single_clone");


 auto legend_e = new TLegend();
   legend_e->AddEntry(h_theta_e_NLO,"NLO","LEP");
   legend_e->AddEntry(h_theta_e_LO,"LO","LEP");

 auto legend_mu = new TLegend();
   legend_mu->AddEntry(h_theta_mu_NLO,"NLO","LEP");
   legend_mu->AddEntry(h_theta_mu_LO,"LO","LEP");

 auto legend_op = new TLegend();
   legend_op->AddEntry(h_op_NLO,"NLO","LEP");
   legend_op->AddEntry(h_op_LO,"LO","LEP");

 auto legend_e_single = new TLegend();
   legend_e_single->AddEntry(h_theta_e_single_NLO,"NLO","LEP");
   legend_e_single->AddEntry(h_theta_e_single_LO,"LO","LEP");

 auto legend_mu_single = new TLegend();
   legend_mu_single->AddEntry(h_theta_mu_single_NLO,"NLO","LEP");
   legend_mu_single->AddEntry(h_theta_mu_single_LO,"LO","LEP");

TCanvas a("a","a",700,700);
a.Divide(2,3);
a.cd(1);
h_theta_mu_single_LO_clone->SetLineColor(kAzure+7);
h_theta_mu_single_LO_clone->SetTitle("Reconstruction efficiency of single mu as a function of muon scattering angle");
h_theta_mu_single_LO_clone->GetXaxis()->SetTitle("Muon angle[rad]");
h_theta_mu_single_LO_clone->GetYaxis()->SetTitle("Efficiency");
h_theta_mu_single_LO_clone->SetMinimum(0.7);
h_theta_mu_single_LO_clone->Draw("E");
gStyle->SetOptStat(0);
a.cd(2);
h_theta_e_single_LO_clone->SetTitle("Reconstruction efficiency of single electron as a function of e- scattering angle");
h_theta_e_single_LO_clone->GetXaxis()->SetTitle("Electron angle[rad]");
h_theta_e_single_LO_clone->GetYaxis()->SetTitle("Efficiency");
h_theta_e_single_LO_clone->SetLineColor(kAzure+7);
h_theta_e_single_LO_clone->SetMinimum(0.7);
h_theta_e_single_LO_clone->Draw("E");
gStyle->SetOptStat(0);
a.cd(3);
h_theta_mu_LO_clone->SetTitle("Reconstruction efficiency of elastic event as a function of mu scattering angle");
h_theta_mu_LO_clone->GetXaxis()->SetTitle("Muon angle[rad]");
h_theta_mu_LO_clone->GetYaxis()->SetTitle("Efficiency");
h_theta_mu_LO_clone->SetLineColor(kAzure+7);
h_theta_mu_LO_clone->SetMinimum(0.7);
h_theta_mu_LO_clone->Draw("E");
gStyle->SetOptStat(0);
a.cd(4);
h_theta_e_LO_clone->SetTitle("Reconstruction efficiency of elastic event as a function of e- scattering angle");
h_theta_e_LO_clone->GetXaxis()->SetTitle("Electron angle[rad]");
h_theta_e_LO_clone->GetYaxis()->SetTitle("Efficiency");
h_theta_e_LO_clone->SetLineColor(kAzure+7);
h_theta_e_LO_clone->SetMinimum(0.7);
h_theta_e_LO_clone->Draw("E");
gStyle->SetOptStat(0);
a.cd(5);
h_op_LO_clone->SetTitle("Reconstruction efficiency of elastic event as a function of e-mu opening angle");
h_op_LO_clone->GetXaxis()->SetTitle("Opening angle[rad]");
h_op_LO_clone->GetYaxis()->SetTitle("Efficiency");
h_op_LO_clone->SetLineColor(kAzure+7);
h_op_LO_clone->SetMinimum(0.7);
h_op_LO_clone->Draw("E");
gStyle->SetOptStat(0);
a.SaveAs("proposal/eff_LO_1hit_clone.pdf");



/*
TCanvas a1("a1","a1",1000,700);
a1.Divide(2,3);

a1.cd(1);

h_theta_mu_single_NLO->SetLineColor(kOrange+10);
h_theta_mu_single_LO->SetLineColor(kAzure+7);
h_theta_mu_single_LO->SetMinimum(0.7);
h_theta_mu_single_LO->Draw("E");
h_theta_mu_single_NLO->Draw("E same");
legend_mu_single->Draw();
gStyle->SetOptStat(0);

a1.cd(2);

h_theta_e_single_NLO->SetLineColor(kOrange+10);
h_theta_e_single_LO->SetLineColor(kAzure+7);
h_theta_e_single_LO->SetMinimum(0.7);
h_theta_e_single_LO->Draw("E");
h_theta_e_single_NLO->Draw("E same");
legend_e_single->Draw();
gStyle->SetOptStat(0);

a1.cd(3);

h_theta_mu_NLO->SetLineColor(kOrange+10);
h_theta_mu_LO->SetLineColor(kAzure+7);
h_theta_mu_NLO->Draw("E");
h_theta_mu_LO->Draw("E same");
legend_mu->Draw();
gStyle->SetOptStat(0);

a1.cd(4);
h_theta_e_NLO->SetLineColor(kOrange+10);
h_theta_e_LO->SetLineColor(kAzure+7);
h_theta_e_LO->SetMinimum(0.7);
h_theta_e_LO->Draw("E");
h_theta_e_NLO->Draw("E same");
legend_e->Draw();
gStyle->SetOptStat(0);

a1.cd(5);

h_op_NLO->SetLineColor(kOrange+10);
h_op_LO->SetLineColor(kAzure+7);
h_op_LO->SetMinimum(0.7);
h_op_LO->Draw("E");
h_op_NLO->Draw("E same");
legend_op->Draw();
gStyle->SetOptStat(0);
a1.SaveAs("proposal/NOoutchi2_compare_LO_NLO_1hit.pdf");
*/

}


