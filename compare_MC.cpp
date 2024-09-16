#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TStyle.h>
#include "TGraphErrors.h"



void compare_MC(int nhits_MC,int nhits_MC2, string MC, string MC2){

TH1::SetDefaultSumw2(kTRUE);
int const NBINS = 6;

TFile *f_thmu_MC=TFile::Open(Form("comparison_RDMC/%s_theta_mu_RD_parallel_pre_basic_el_Bend_%dhit_mc.root",MC.c_str(),nhits_MC));
TFile *f_the_MC=TFile::Open(Form("comparison_RDMC/%s_theta_e_RD_parallel_pre_basic_el_Bend_%dhit_mc.root",MC.c_str(),nhits_MC));


TH1D* h_theta_e_MC=(TH1D*)f_the_MC->Get("theta_e");
TH1D* h_theta_mu_MC=(TH1D*)f_thmu_MC->Get("theta_mu");


TFile *f_thmu_MC2=TFile::Open(Form("comparison_RDMC/%s_theta_mu_RD_parallel_pre_basic_el_Bend_%dhit_mc.root",MC2.c_str(),nhits_MC2));
TFile *f_the_MC2=TFile::Open(Form("comparison_RDMC/%s_theta_e_RD_parallel_pre_basic_el_Bend_%dhit_mc.root",MC2.c_str(),nhits_MC2));


TH1D* h_theta_e_MC2=(TH1D*)f_the_MC2->Get("theta_e");
TH1D* h_theta_mu_MC2=(TH1D*)f_thmu_MC2->Get("theta_mu");


 auto l = new TLegend(0.75,0.15,0.9,0.3);
l->AddEntry(h_theta_e_MC,Form("%d hit sared",nhits_MC),"LEP");
l->AddEntry(h_theta_e_MC2,Form("%d hit sared",nhits_MC2),"LEP");

 auto l2 = new TLegend(0.75,0.15,0.9,0.3);
l2->AddEntry(h_theta_mu_MC,Form("%d hit sared",nhits_MC),"LEP");
l2->AddEntry(h_theta_mu_MC2,Form("%d hit sared",nhits_MC2),"LEP");

TCanvas a1("a1","a1",1000,700);
a1.Divide(2,2);
a1.cd(1);
cout << "theta el. reco:" << endl;
h_theta_e_MC2->GetXaxis()->SetRangeUser(0.003,0.03);
h_theta_e_MC->GetXaxis()->SetRangeUser(0.003,0.03);
h_theta_e_MC2->SetMinimum(1.);
h_theta_e_MC2->Draw("E");
h_theta_e_MC->SetLineColor(kPink+10);
h_theta_e_MC->Draw("E same");
gPad->SetLogy();
gStyle->SetOptStat(0);
l->Draw();
a1.cd(2);

TH1D * h3 = (TH1D*) h_theta_e_MC->Clone();
h3->Divide(h_theta_e_MC2);
h3->SetMinimum(0.);
h3->SetTitle(Form("%d hits/%d hits VS #theta_el",nhits_MC,nhits_MC2));
h3->Draw("E");
a1.cd(3);
h_theta_mu_MC2->SetMinimum(1.);
h_theta_mu_MC2->GetXaxis()->SetRangeUser(0.0003,0.002);
h_theta_mu_MC->GetXaxis()->SetRangeUser(0.0003,0.002);
h_theta_mu_MC2->Draw("E");
h_theta_mu_MC->SetLineColor(kPink+10);
h_theta_mu_MC->Draw("E same");
gStyle->SetOptStat(0);
gPad->SetLogy();
l2->Draw();

a1.cd(4);
TH1D * h4 = (TH1D*) h_theta_mu_MC->Clone();
h4->Divide(h_theta_mu_MC2);
h4->SetMinimum(0.);
h4->SetTitle(Form("%d hits/%d hits VS #theta_mu",nhits_MC,nhits_MC2));
h4->Draw("E");
a1.SaveAs(Form("comparison_RDMC/%s_%s_theta_MC_basic_%dhit_%dhit.pdf",MC.c_str(),MC2.c_str(),nhits_MC,nhits_MC2));

}

