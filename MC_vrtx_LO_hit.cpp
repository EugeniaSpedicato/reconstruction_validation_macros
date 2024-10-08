
#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TStyle.h>


void MC_vrtx_LO_hit(){

TFile *f_thmu_LO_clone1 = TFile::Open("proposal/quality_mu_gen_1hitFirstModules.root");
TFile *f_thmu_LO_clone2 = TFile::Open("proposal/quality_mu_gen_2hitFirstModules.root");
TFile *f_thmu_LO_clone0 = TFile::Open("proposal/quality_mu_gen_0hit.root");

TFile *f_the_LO_clone1 = TFile::Open("proposal/quality_el_gen_1hitFirstModules.root");
TFile *f_the_LO_clone2 = TFile::Open("proposal/quality_el_gen_2hitFirstModules.root");
TFile *f_the_LO_clone0 = TFile::Open("proposal/quality_el_gen_0hit.root");

TFile *f_op_LO_clone1 = TFile::Open("proposal/NOoutchi2_op_VRTX_eff_LO_1hit_clone.root");
TFile *f_op_LO_clone2 = TFile::Open("proposal/NOoutchi2_op_VRTX_eff_LO_2hit_clone.root");
TFile *f_op_LO_clone0 = TFile::Open("proposal/NOoutchi2_op_VRTX_eff_LO_0hit_clone.root");

TFile *f_op_NLO_quality1 = TFile::Open("proposal/VRTX_wrong_NLO_1hitFirstModules_true_clone_quality.root");
TFile *f_op_NLO_quality2 = TFile::Open("proposal/VRTX_wrong_NLO_2hitFirstModules_true_clone_quality.root");
TFile *f_op_NLO_quality0 = TFile::Open("proposal/VRTX_wrong_NLO_0hit_true_clone_quality.root");


TH1::SetDefaultSumw2(kTRUE);




TH1D *h_theta_mu_LO_clone0=(TH1D*)f_thmu_LO_clone0->Get("theta_mu_clone_quality");
TH1D *h_theta_mu_LO_clone1=(TH1D*)f_thmu_LO_clone1->Get("theta_mu_clone_quality");
TH1D *h_theta_mu_LO_clone2=(TH1D*)f_thmu_LO_clone2->Get("theta_mu_clone_quality");

TH1D *h_theta_e_LO_clone0=(TH1D*)f_the_LO_clone0->Get("theta_e_clone_quality");
TH1D *h_theta_e_LO_clone1=(TH1D*)f_the_LO_clone1->Get("theta_e_clone_quality");
TH1D *h_theta_e_LO_clone2=(TH1D*)f_the_LO_clone2->Get("theta_e_clone_quality");

TH1D *h_op_LO_clone0=(TH1D*)f_op_LO_clone0->Get("h_opening_clone_quality");
TH1D *h_op_LO_clone1=(TH1D*)f_op_LO_clone1->Get("h_opening_clone_quality");
TH1D *h_op_LO_clone2=(TH1D*)f_op_LO_clone2->Get("h_opening_clone_quality");

TH1D *h_op_NLO_quality0=(TH1D*)f_op_NLO_quality0->Get("h_opening_wrong_quality");
TH1D *h_op_NLO_quality1=(TH1D*)f_op_NLO_quality1->Get("h_opening_wrong_quality");
TH1D *h_op_NLO_quality2=(TH1D*)f_op_NLO_quality2->Get("h_opening_wrong_quality");



 auto legend_q = new TLegend(0.6,0.65,0.9,0.9);
   legend_q->AddEntry(h_op_NLO_quality0,"0 hit shared","LEP");
   legend_q->AddEntry(h_op_NLO_quality1,"1 hit shared","LEP");
   legend_q->AddEntry(h_op_NLO_quality2,"2 hit shared","LEP");
/*
TCanvas q("q","q",700,700);
h_op_NLO_quality2->SetStats(0);
h_op_NLO_quality1->SetStats(0);
h_op_NLO_quality0->SetStats(0);
h_op_NLO_quality0->SetLineColor(kOrange+10);
h_op_NLO_quality1->SetLineColor(kGreen+1);
h_op_NLO_quality2->SetLineColor(kAzure+7);
h_op_NLO_quality2->SetTitle("#wrong vertexing events / #good elastic reconstruction events");
h_op_NLO_quality2->SetMinimum(0.);
h_op_NLO_quality2->Draw("E");
h_op_NLO_quality1->Draw("E same");
h_op_NLO_quality0->Draw("E same");
legend_q->Draw();

q.SaveAs("proposal/wrong_vrtxing_NLO.pdf");
*/

 auto legend_e = new TLegend(0.6,0.1,0.9,0.25);
   legend_e->AddEntry(h_theta_e_LO_clone0,"LO 0 hit","LEP");
   legend_e->AddEntry(h_theta_e_LO_clone1,"LO 1 hit","LEP");
   legend_e->AddEntry(h_theta_e_LO_clone2,"LO 2 hit","LEP");

 auto legend_mu = new TLegend(0.1,0.1,0.4,0.25);
   legend_mu->AddEntry(h_theta_mu_LO_clone0,"LO 0 hit","LEP");
   legend_mu->AddEntry(h_theta_mu_LO_clone1,"LO 1 hit","LEP");
   legend_mu->AddEntry(h_theta_mu_LO_clone2,"LO 2 hit","LEP");

 auto legend_op = new TLegend(0.6,0.1,0.9,0.25);
   legend_op->AddEntry(h_op_LO_clone0,"LO 0 hit","LEP");
   legend_op->AddEntry(h_op_LO_clone1,"LO 1 hit","LEP");
   legend_op->AddEntry(h_op_LO_clone2,"LO 2 hit","LEP");


TCanvas a1("a1","a1",2100,700);
a1.Divide(3,1);
a1.cd(1);
h_theta_mu_LO_clone2->SetTitle("Reconstruction efficiency of elastic event as a function of mu scattering angle");
h_theta_mu_LO_clone2->GetXaxis()->SetTitle("Muon angle[rad]");
h_theta_mu_LO_clone2->GetYaxis()->SetTitle("Efficiency");
h_theta_mu_LO_clone2->SetLineColor(kAzure+7);
h_theta_mu_LO_clone2->SetMinimum(0.6);
TGaxis::SetMaxDigits(3);
h_theta_mu_LO_clone2->Draw("E");
h_theta_mu_LO_clone1->SetLineColor(kGreen+1);
h_theta_mu_LO_clone1->Draw("E same");
h_theta_mu_LO_clone0->SetLineColor(kOrange+10);
h_theta_mu_LO_clone0->Draw("E same");
legend_mu->Draw();
gStyle->SetOptStat(0);
a1.cd(2);
h_theta_e_LO_clone2->SetTitle("Reconstruction efficiency of elastic event as a function of e- scattering angle");
h_theta_e_LO_clone2->GetXaxis()->SetTitle("Electron angle[rad]");
h_theta_e_LO_clone2->GetYaxis()->SetTitle("Efficiency");
h_theta_e_LO_clone2->SetLineColor(kAzure+7);
h_theta_e_LO_clone2->SetMinimum(0.6);
h_theta_e_LO_clone2->Draw("E");
h_theta_e_LO_clone1->SetLineColor(kGreen+1);
h_theta_e_LO_clone1->Draw("E same");
h_theta_e_LO_clone0->SetLineColor(kOrange+10);
h_theta_e_LO_clone0->Draw("E same");
legend_e->Draw();

gStyle->SetOptStat(0);
/*a1.cd(3);
h_op_LO_clone2->SetTitle("Reconstruction efficiency of elastic event as a function of e-mu opening angle");
h_op_LO_clone2->GetXaxis()->SetTitle("Opening angle[rad]");
h_op_LO_clone2->GetYaxis()->SetTitle("Efficiency");
h_op_LO_clone2->SetLineColor(kAzure+7);
h_op_LO_clone2->SetMinimum(0.);
h_op_LO_clone2->Draw("E");
h_op_LO_clone1->SetLineColor(kGreen+1);
h_op_LO_clone1->Draw("E same");
h_op_LO_clone0->SetLineColor(kOrange+10);
h_op_LO_clone0->Draw("E same");
legend_op->Draw();
gStyle->SetOptStat(0);*/
a1.SaveAs("proposal/vrtx_eff_event_LO_012hit_clone.pdf");


}


