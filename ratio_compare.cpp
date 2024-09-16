#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TStyle.h>
#include "TGraphErrors.h"



void ratio_compare(int nhits_data,int nhits_MC, string data, string MC, int nhits_data2,int nhits_MC2,string data2,string MC2){

TH1::SetDefaultSumw2(kTRUE);
int const NBINS = 6;


TFile *f_the_MC=TFile::Open(Form("/home/espedica/modifica_fair_install/modificaFairRoot/share/MUonE/macros/comparison_RDMC/%s_%s_ratio_theta_el_%dhit_%dhit.root",data.c_str(),MC.c_str(),nhits_data,nhits_MC));
TFile *f_thmu_MC=TFile::Open(Form("/home/espedica/modifica_fair_install/modificaFairRoot/share/MUonE/macros/comparison_RDMC/%s_%s_ratio_theta_mu_%dhit_%dhit.root",data.c_str(),MC.c_str(),nhits_data,nhits_MC));
TFile *f_opening_MC=TFile::Open(Form("/home/espedica/modifica_fair_install/modificaFairRoot/share/MUonE/macros/comparison_RDMC/%s_%s_ratio_opening_%dhit_%dhit.root",data.c_str(),MC.c_str(),nhits_data,nhits_MC));


TH1D* h_theta_e_MC=(TH1D*)f_the_MC->Get("theta_e");
TH1D* h_theta_mu_MC=(TH1D*)f_thmu_MC->Get("theta_mu");
TH1D* h_opening_MC=(TH1D*)f_opening_MC->Get("h_opening");

TFile *f_the_MC2=TFile::Open(Form("comparison_RDMC/%s_%s_ratio_theta_el_%dhit_%dhit.root",data2.c_str(),MC2.c_str(),nhits_data2,nhits_MC2));
TFile *f_thmu_MC2=TFile::Open(Form("comparison_RDMC/%s_%s_ratio_theta_mu_%dhit_%dhit.root",data2.c_str(),MC2.c_str(),nhits_data2,nhits_MC2));
TFile *f_opening_MC2=TFile::Open(Form("comparison_RDMC/%s_%s_ratio_opening_%dhit_%dhit.root",data2.c_str(),MC2.c_str(),nhits_data2,nhits_MC2));


TH1D* h_theta_e_MC2=(TH1D*)f_the_MC2->Get("theta_e");
TH1D* h_theta_mu_MC2=(TH1D*)f_thmu_MC2->Get("theta_mu");
TH1D* h_opening_MC2=(TH1D*)f_opening_MC2->Get("h_opening");




 auto l = new TLegend(0.55,0.15,0.9,0.3);
l->AddEntry(h_theta_e_MC,Form("MC %s, data %s",MC.c_str(),data.c_str()),"LEP");
l->AddEntry(h_theta_e_MC2,Form("MC %s, data %s",MC2.c_str(),data2.c_str()),"LEP");

 auto l2 = new TLegend(0.55,0.15,0.9,0.3);
l2->AddEntry(h_theta_mu_MC,Form("MC %s, data %s",MC.c_str(),data.c_str()),"LEP");
l2->AddEntry(h_theta_mu_MC2,Form("MC %s, data %s",MC2.c_str(),data2.c_str()),"LEP");

 auto l3 = new TLegend(0.55,0.15,0.9,0.3);
l3->AddEntry(h_opening_MC,Form("MC %s, data %s",MC.c_str(),data.c_str()),"LEP");
l3->AddEntry(h_opening_MC2,Form("MC %s, data %s",MC2.c_str(),data2.c_str()),"LEP");

TCanvas a1("a1","a1",3000,1000);
a1.Divide(3,1);
a1.cd(1);
h_theta_e_MC->SetLineColor(kBlue);
h_theta_e_MC2->SetLineColor(kPink+10);
h_theta_e_MC->Draw();
h_theta_e_MC2->Draw("same");
gStyle->SetOptStat(0);
l->Draw();
a1.cd(2);
h_theta_mu_MC->SetLineColor(kBlue);
h_theta_mu_MC2->SetLineColor(kPink+10);
h_theta_mu_MC->Draw();
h_theta_mu_MC2->Draw("same");
gStyle->SetOptStat(0);
l2->Draw();
a1.cd(3);
h_opening_MC->SetLineColor(kBlue);
h_opening_MC2->SetLineColor(kPink+10);
h_opening_MC->Draw();
h_opening_MC2->Draw("same");
gStyle->SetOptStat(0);
l3->Draw();
a1.SaveAs(Form("comparison_RDMC/%s_%s_ratio_comparison_%dhit_%dhit.pdf",MC.c_str(),MC2.c_str(),nhits_MC,nhits_MC2));

}

