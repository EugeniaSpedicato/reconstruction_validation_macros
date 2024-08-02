
#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TStyle.h>


void quality_NLO_cut(){


TFile *f_thmu_NLO_clone2 = TFile::Open("proposal/quality_mu_gen_2hitFirstModules.root");
TFile *f_thmu_NLO_clone1 = TFile::Open("proposal/quality_mu_gen_1hitFirstModules.root");

TFile *f_the_NLO_clone2 = TFile::Open("proposal/quality_el_gen_2hitFirstModules.root");
TFile *f_the_NLO_clone1 = TFile::Open("proposal/quality_el_gen_1hitFirstModules.root");




TH1::SetDefaultSumw2(kTRUE);

TH1D *h_theta_mu_NLO_clone1=(TH1D*)f_thmu_NLO_clone1->Get("theta_mu_clone_quality");
TH1D *h_theta_mu_NLO_clone2=(TH1D*)f_thmu_NLO_clone2->Get("theta_mu_clone_quality");

TH1D *h_theta_e_NLO_clone1=(TH1D*)f_the_NLO_clone1->Get("theta_e_clone_quality");
TH1D *h_theta_e_NLO_clone2=(TH1D*)f_the_NLO_clone2->Get("theta_e_clone_quality");


 auto legend_e = new TLegend(0.6,0.1,0.9,0.25);
   legend_e->AddEntry(h_theta_e_NLO_clone1,"NLO 1 hit","LEP");
   legend_e->AddEntry(h_theta_e_NLO_clone2,"NLO 2 hit","LEP");

 auto legend_mu = new TLegend(0.1,0.1,0.4,0.25);
   legend_mu->AddEntry(h_theta_mu_NLO_clone1,"NLO 1 hit","LEP");
   legend_mu->AddEntry(h_theta_mu_NLO_clone2,"NLO 2 hit","LEP");


TCanvas a1("a1","a1",2100,700);
a1.Divide(2,1);
a1.cd(1);
h_theta_mu_NLO_clone2->SetTitle("Reconstruction efficiency of elastic event as a function of mu scattering angle");
h_theta_mu_NLO_clone2->GetXaxis()->SetTitle("Muon angle[rad]");
h_theta_mu_NLO_clone2->GetYaxis()->SetTitle("Efficiency");
h_theta_mu_NLO_clone2->SetLineColor(kAzure+7);
h_theta_mu_NLO_clone2->SetMinimum(0.6);
TGaxis::SetMaxDigits(3);
h_theta_mu_NLO_clone2->Draw("E");
h_theta_mu_NLO_clone1->SetLineColor(kGreen+1);
h_theta_mu_NLO_clone1->Draw("E same");
legend_mu->Draw();
gStyle->SetOptStat(0);
a1.cd(2);
h_theta_e_NLO_clone2->SetTitle("Reconstruction efficiency of elastic event as a function of e- scattering angle");
h_theta_e_NLO_clone2->GetXaxis()->SetTitle("Electron angle[rad]");
h_theta_e_NLO_clone2->GetYaxis()->SetTitle("Efficiency");
h_theta_e_NLO_clone2->SetLineColor(kAzure+7);
h_theta_e_NLO_clone2->SetMinimum(0.6);
h_theta_e_NLO_clone2->Draw("E");
h_theta_e_NLO_clone1->SetLineColor(kGreen+1);
h_theta_e_NLO_clone1->Draw("E same");
legend_e->Draw();

gStyle->SetOptStat(0);
a1.SaveAs("proposal/eff_event_NLO_12hitFirstModules_02cut_clone_quality.pdf");






}


