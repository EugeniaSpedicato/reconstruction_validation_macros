#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TStyle.h>
#include "TGraphErrors.h"



void compare_data(int nhits_rd,int nhits_rd2, string rd, string rd2){

TH1::SetDefaultSumw2(kTRUE);
int const NBINS = 6;

TFile *f_thmu_RD=TFile::Open(Form("comparison_RDMC/%s_theta_mu_RD_parallel_pre_basic_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_the_RD=TFile::Open(Form("comparison_RDMC/%s_theta_e_RD_parallel_pre_basic_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));


TH1D* h_theta_e_RD=(TH1D*)f_the_RD->Get("theta_e");
TH1D* h_theta_mu_RD=(TH1D*)f_thmu_RD->Get("theta_mu");


TFile *f_thmu_RD2=TFile::Open(Form("comparison_RDMC/%s_theta_mu_RD_parallel_pre_basic_el_Bend_%dhit_data.root",rd2.c_str(),nhits_rd2));
TFile *f_the_RD2=TFile::Open(Form("comparison_RDMC/%s_theta_e_RD_parallel_pre_basic_el_Bend_%dhit_data.root",rd2.c_str(),nhits_rd2));


TH1D* h_theta_e_RD2=(TH1D*)f_the_RD2->Get("theta_e");
TH1D* h_theta_mu_RD2=(TH1D*)f_thmu_RD2->Get("theta_mu");


 auto l = new TLegend(0.75,0.15,0.9,0.3);
l->AddEntry(h_theta_e_RD,Form("%d hit sared",nhits_rd),"LEP");
l->AddEntry(h_theta_e_RD2,Form("%d hit sared",nhits_rd2),"LEP");

 auto l2 = new TLegend(0.75,0.15,0.9,0.3);
l2->AddEntry(h_theta_mu_RD,Form("%d hit sared",nhits_rd),"LEP");
l2->AddEntry(h_theta_mu_RD2,Form("%d hit sared",nhits_rd2),"LEP");

TCanvas a1("a1","a1",1000,700);
a1.Divide(2,2);
a1.cd(1);
cout << "theta el. reco:" << endl;
h_theta_e_RD2->GetXaxis()->SetRangeUser(0.003,0.03);
h_theta_e_RD->GetXaxis()->SetRangeUser(0.003,0.03);
h_theta_e_RD2->SetMinimum(1.);
h_theta_e_RD2->Draw("E");
h_theta_e_RD->SetLineColor(kPink+10);
h_theta_e_RD->Draw("E same");
gPad->SetLogy();
gStyle->SetOptStat(0);
l->Draw();
a1.cd(2);

TH1D * h3 = (TH1D*) h_theta_e_RD->Clone();
h3->Divide(h_theta_e_RD2);
h3->SetMinimum(0.);
h3->SetTitle(Form("%d hits/%d hits VS #theta_el",nhits_rd,nhits_rd2));
h3->Draw("E");
a1.cd(3);
h_theta_mu_RD2->SetMinimum(1.);
h_theta_mu_RD2->GetXaxis()->SetRangeUser(0.0003,0.002);
h_theta_mu_RD->GetXaxis()->SetRangeUser(0.0003,0.002);
h_theta_mu_RD2->Draw("E");
h_theta_mu_RD->SetLineColor(kPink+10);
h_theta_mu_RD->Draw("E same");
gStyle->SetOptStat(0);
gPad->SetLogy();
l2->Draw();

a1.cd(4);
TH1D * h4 = (TH1D*) h_theta_mu_RD->Clone();
h4->Divide(h_theta_mu_RD2);
h4->SetMinimum(0.);
h4->SetTitle(Form("%d hits/%d hits VS #theta_mu",nhits_rd,nhits_rd2));
h4->Draw("E");
a1.SaveAs(Form("comparison_RDMC/%s_%s_theta_data_basic_%dhit_%dhit.pdf",rd.c_str(),rd2.c_str(),nhits_rd,nhits_rd2));

}

