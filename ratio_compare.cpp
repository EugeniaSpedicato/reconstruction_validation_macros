#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TStyle.h>
#include "TGraphErrors.h"



void ratio_compare(string path1, int nhits_data,int nhits_MC, string data, string MC, string path2, int nhits_data2,int nhits_MC2,string data2,string MC2){

TH1::SetDefaultSumw2(kTRUE);
int const NBINS = 6;


TFile *f_the_MC=TFile::Open(Form("%s/comparison_RDMC/%s_%s_ratio_basicAllFiducial_chi150_flatArea_theNA_theNA_thmu04_nstubs15_nochivrtx_zPosNoPeak_theta_el_%dhit_%dhit_smearing.root",path1.c_str(),data.c_str(),MC.c_str(),nhits_data,nhits_MC));
TFile *f_thmu_MC=TFile::Open(Form("%s/comparison_RDMC/%s_%s_ratio_basicAllFiducial_chi150_flatArea_theNA_theNA_thmu04_nstubs15_nochivrtx_zPosNoPeak_theta_mu_%dhit_%dhit_smearing.root",path1.c_str(),data.c_str(),MC.c_str(),nhits_data,nhits_MC));
TFile *f_opening_MC=TFile::Open(Form("%s/comparison_RDMC/%s_%s_ratio_basicAllFiducial_chi150_flatArea_theNA_theNA_thmu04_nstubs15_nochivrtx_zPosNoPeak_opening_%dhit_%dhit_smearing.root",path1.c_str(),data.c_str(),MC.c_str(),nhits_data,nhits_MC));


TH1D* h_theta_e_MC=(TH1D*)f_the_MC->Get("theta_e");
TH1D* h_theta_mu_MC=(TH1D*)f_thmu_MC->Get("theta_mu");
TH1D* h_opening_MC=(TH1D*)f_opening_MC->Get("h_opening");

TFile *f_the_MC2=TFile::Open(Form("%s/comparison_RDMC/%s_%s_ratio_basicAllFiducial_chi150_flatArea_theNA_theNA_thmu04_nstubs15_nochivrtx_zPosNoPeak_theta_el_%dhit_%dhit.root",path2.c_str(),data2.c_str(),MC2.c_str(),nhits_data2,nhits_MC2));
TFile *f_thmu_MC2=TFile::Open(Form("%s/comparison_RDMC/%s_%s_ratio_basicAllFiducial_chi150_flatArea_theNA_theNA_thmu04_nstubs15_nochivrtx_zPosNoPeak_theta_mu_%dhit_%dhit.root",path2.c_str(),data2.c_str(),MC2.c_str(),nhits_data2,nhits_MC2));
TFile *f_opening_MC2=TFile::Open(Form("%s/comparison_RDMC/%s_%s_ratio_basicAllFiducial_chi150_flatArea_theNA_theNA_thmu04_nstubs15_nochivrtx_zPosNoPeak_opening_%dhit_%dhit.root",path2.c_str(),data2.c_str(),MC2.c_str(),nhits_data2,nhits_MC2));


TH1D* h_theta_e_MC2=(TH1D*)f_the_MC2->Get("theta_e");
TH1D* h_theta_mu_MC2=(TH1D*)f_thmu_MC2->Get("theta_mu");
TH1D* h_opening_MC2=(TH1D*)f_opening_MC2->Get("h_opening");


// auto l = new TLegend(0.65,0.75,0.9,0.91);
 auto l = new TLegend(0.25,0.75,0.5,0.91);
l->AddEntry(h_theta_e_MC,"smeared","LEP");
l->AddEntry(h_theta_e_MC2,"not smeared","LEP");

// auto l2 = new TLegend(0.55,0.15,0.9,0.3);
 auto l2 = new TLegend(0.55,0.76,0.9,0.91);
l2->AddEntry(h_theta_mu_MC,"smeared","LEP");
l2->AddEntry(h_theta_mu_MC2,"not smeared","LEP");

 auto l3 = new TLegend(0.55,0.75,0.9,0.91);
l3->AddEntry(h_opening_MC,"smeared","LEP");
l3->AddEntry(h_opening_MC2,"not smeared","LEP");
/*
// auto l = new TLegend(0.55,0.15,0.9,0.3);
 auto l = new TLegend(0.55,0.75,0.9,0.91);
l->AddEntry(h_theta_e_MC,Form("%d MC %s, data %s",nhits_MC,MC.c_str(),data.c_str()),"LEP");
l->AddEntry(h_theta_e_MC2,Form("%d MC %s, data %s",nhits_MC2,MC2.c_str(),data2.c_str()),"LEP");

// auto l2 = new TLegend(0.55,0.15,0.9,0.3);
 auto l2 = new TLegend(0.55,0.76,0.9,0.91);
l2->AddEntry(h_theta_mu_MC,Form("%d MC %s, data %s",nhits_MC,MC.c_str(),data.c_str()),"LEP");
l2->AddEntry(h_theta_mu_MC2,Form("%d MC %s, data %s",nhits_MC2,MC2.c_str(),data2.c_str()),"LEP");

 auto l3 = new TLegend(0.55,0.75,0.9,0.91);
l3->AddEntry(h_opening_MC,Form("%d MC %s, data %s",nhits_MC,MC.c_str(),data.c_str()),"LEP");
l3->AddEntry(h_opening_MC2,Form("%d MC %s, data %s",nhits_MC2,MC2.c_str(),data2.c_str()),"LEP");
*/

/*
 auto l = new TLegend(0.65,0.22,0.9,0.37);
l->AddEntry(h_theta_e_MC,"realistic","LEP");
l->AddEntry(h_theta_e_MC2,"ideal","LEP");

 auto l2 = new TLegend(0.65,0.22,0.9,0.37);
l2->AddEntry(h_theta_mu_MC,"realistic","LEP");
l2->AddEntry(h_theta_mu_MC2,"ideal","LEP");

 auto l3 = new TLegend(0.65,0.22,0.9,0.37);
l3->AddEntry(h_opening_MC,"realistic","LEP");
l3->AddEntry(h_opening_MC2,"ideal","LEP");
*/


/*
 auto l = new TLegend(0.65,0.12,0.9,0.27);
l->AddEntry(h_theta_e_MC,Form("%d hit",nhits_MC),"LEP");
l->AddEntry(h_theta_e_MC2,Form("%d hit",nhits_MC2),"LEP");

// auto l2 = new TLegend(0.35,0.72,0.6,0.87);
 auto l2 = new TLegend(0.65,0.12,0.9,0.27);
l2->AddEntry(h_theta_mu_MC,Form("%d hit",nhits_MC),"LEP");
l2->AddEntry(h_theta_mu_MC2,Form("%d hit",nhits_MC2),"LEP");

 auto l3 = new TLegend(0.65,0.12,0.9,0.27);
l3->AddEntry(h_opening_MC,Form("%d hit",nhits_MC),"LEP");
l3->AddEntry(h_opening_MC2,Form("%d hit",nhits_MC2),"LEP");
*/

/*
auto l = new TLegend(0.65,0.12,0.9,0.27);
l->AddEntry(h_theta_e_MC,Form("%d hit",nhits_MC),"LEP");
l->AddEntry(h_theta_e_MC2,Form("0-lose hit",nhits_MC2),"LEP");

// auto l2 = new TLegend(0.35,0.72,0.6,0.87);
auto l2 = new TLegend(0.65,0.12,0.9,0.27);
l2->AddEntry(h_theta_mu_MC,Form("%d hit",nhits_MC),"LEP");
l2->AddEntry(h_theta_mu_MC2,Form("0-lose hit",nhits_MC2),"LEP");

 auto l3 = new TLegend(0.65,0.12,0.9,0.27);
l3->AddEntry(h_opening_MC,Form("%d hit",nhits_MC),"LEP");
l3->AddEntry(h_opening_MC2,Form("0-lose hit",nhits_MC2),"LEP");
*/

TCanvas a1("a1","a1",3000,1000);
a1.Divide(3,1);
a1.cd(1);
h_theta_e_MC->GetYaxis()->SetTitle("Data/MC ratio");
h_theta_e_MC->GetXaxis()->SetTitle("Electron angle [rad]");
h_theta_e_MC->SetTitle("Electron scattering angle");
h_theta_e_MC->SetLineColor(kBlue);
h_theta_e_MC2->SetLineColor(kPink+10);
h_theta_e_MC->GetXaxis()->SetRangeUser(0.0,0.032);//User(0.005,0.032);
h_theta_e_MC->Draw();
h_theta_e_MC2->Draw("same");
h_theta_e_MC->SetMinimum(0.75);
//h_theta_e_MC->SetMaximum(1.);
h_theta_e_MC->SetMaximum(4.);
h_theta_e_MC->GetXaxis()->SetMaxDigits(3);
gStyle->SetOptStat(0);
l->Draw();
a1.cd(2);
h_theta_mu_MC->SetTitle("Muon scattering angle");
h_theta_mu_MC->GetXaxis()->SetTitle("Muon angle [rad]");
h_theta_mu_MC->GetYaxis()->SetTitle("Data/MC ratio");
h_theta_mu_MC->SetLineColor(kBlue);
h_theta_mu_MC2->SetLineColor(kPink+10);
h_theta_mu_MC->GetXaxis()->SetRangeUser(0.0,0.003);//geUser(0.0002,0.0014);
h_theta_mu_MC->Draw();
h_theta_mu_MC2->Draw("same");
h_theta_mu_MC->SetMinimum(0.7);
h_theta_mu_MC->SetMaximum(1.);
h_theta_mu_MC->GetXaxis()->SetMaxDigits(3);
gStyle->SetOptStat(0);
l2->Draw();
a1.cd(3);
h_opening_MC->SetTitle("#mu-e opening angle");
h_opening_MC->GetXaxis()->SetTitle("Opening angle [rad]");
h_opening_MC->GetYaxis()->SetTitle("Data/MC ratio");
h_opening_MC->SetLineColor(kBlue);
h_opening_MC2->SetLineColor(kPink+10);
h_opening_MC->GetXaxis()->SetRangeUser(0.0,0.032);//User(0.005,0.032);
h_opening_MC->Draw();
h_opening_MC2->Draw("same");
h_opening_MC->SetMinimum(0.75);
h_opening_MC->GetXaxis()->SetMaxDigits(3);
gStyle->SetOptStat(0);
l3->Draw();
a1.SaveAs(Form("comparison_RDMC/%s_%s_ratio_basicAllFiducial_chi150_flatArea_theNA_theNA_thmu04_nstubs15_nochivrtx_zPosNoPeak_comparison_%dhit_%dhit.pdf",MC.c_str(),MC2.c_str(),nhits_MC,nhits_MC2));

}

