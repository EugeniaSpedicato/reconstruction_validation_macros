#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TStyle.h>


void skim_comparison(){

TFile *f_MC = TFile::Open("outputs/histos_fifth_ALLPUmu_eff97.root");//histos_fourth_ALLPUmu_eff97_1hitallowed.root");///histos_third_ALLPUmu_eff97.root");
//TFile *f_RD = TFile::Open("/mnt/raid10/DATA/espedica/merged/filter_GA/job_all/outputs/run_6/histos_run_6_sum.root");//_gio_1hit_singleclean/outputs/run_6/histos_run_6_sumold.root");

TFile *f_RD = TFile::Open("/home/espedica/tr2023_promptanalysis/GiovanniA/job_gio/outputs/run_6/histos_run_6_sumold.root");

//TFile *f_MC = TFile::Open("outputs_mesmer/histos_third_1PUmu.root");


//TFile *f_MC = TFile::Open("outputs/histos_third_ALLPUmu_eff97.root");
//TFile *f_RD = TFile::Open("/mnt/raid10/DATA/espedica/merged/filter_GA/job_all/outputs/run_6/histos_run_6_sum.root");

// /mnt/raid10/DATA/espedica/merged/filter_GA/job_gio_1hit_singleclean/outputs/run_6/histos_run_6_sumold.root");
///home/espedica/tr2023_promptanalysis/GiovanniA/job_/outputs/run_6/histos_run_6_sum.root");
//abbiendi/SKIM/test_LAST_fired12/outputs/run_6/histos_run_6_inputs.root");

TH1::SetDefaultSumw2(kTRUE);


TH1D *h_MC_nTotStubs_trackable=(TH1D*)f_MC->Get("h_nTotStubs_trackable");
TH1D *h_RD_nTotStubs_trackable=(TH1D*)f_RD->Get("h_nTotStubs_trackable");

TH1D *h_MC_h_nStubs_0_trackable=(TH1D*)f_MC->Get("h_nStubs_0_trackable");
TH1D *h_RD_h_nStubs_0_trackable=(TH1D*)f_RD->Get("h_nStubs_0_trackable");

TH1D *h_MC_h_nStubs_1_trackable=(TH1D*)f_MC->Get("h_nStubs_1_trackable");
TH1D *h_RD_h_nStubs_1_trackable=(TH1D*)f_RD->Get("h_nStubs_1_trackable");


TH1D *h_MC_h_nTotStubs_single_clean=(TH1D*)f_MC->Get("h_nTotStubs_single_clean");
TH1D *h_RD_h_nTotStubs_single_clean=(TH1D*)f_RD->Get("h_nTotStubs_single_clean");

TH1D *h_MC_h_nStubs_0_single_clean=(TH1D*)f_MC->Get("h_nStubs_0_single_clean");
TH1D *h_RD_h_nStubs_0_single_clean=(TH1D*)f_RD->Get("h_nStubs_0_single_clean");

TH1D *h_MC_h_nStubs_1_single_clean=(TH1D*)f_MC->Get("h_nStubs_1_single_clean");
TH1D *h_RD_h_nStubs_1_single_clean=(TH1D*)f_RD->Get("h_nStubs_1_single_clean");


TH1D *h_MC_h_nTotStubs_pileup234_skim=(TH1D*)f_MC->Get("h_nTotStubs_pileup234_skim");
TH1D *h_RD_h_nTotStubs_pileup234_skim=(TH1D*)f_RD->Get("h_nTotStubs_pileup234_skim");

TH1D *h_MC_h_nStubs_0_pileup234_skim=(TH1D*)f_MC->Get("h_nStubs_0_pileup234_skim");
TH1D *h_RD_h_nStubs_0_pileup234_skim=(TH1D*)f_RD->Get("h_nStubs_0_pileup234_skim");

TH1D *h_MC_h_nStubs_1_pileup234_skim=(TH1D*)f_MC->Get("h_nStubs_1_pileup234_skim");
TH1D *h_RD_h_nStubs_1_pileup234_skim=(TH1D*)f_RD->Get("h_nStubs_1_pileup234_skim");


TH1D *h_MC_h_nTotStubs_pileup_skim=(TH1D*)f_MC->Get("h_nTotStubs_pileup_skim");
TH1D *h_RD_h_nTotStubs_pileup_skim=(TH1D*)f_RD->Get("h_nTotStubs_pileup_skim");

TH1D *h_MC_h_nStubs_0_pileup_skim=(TH1D*)f_MC->Get("h_nStubs_0_pileup_skim");
TH1D *h_RD_h_nStubs_0_pileup_skim=(TH1D*)f_RD->Get("h_nStubs_0_pileup_skim");

TH1D *h_MC_h_nStubs_1_pileup_skim=(TH1D*)f_MC->Get("h_nStubs_1_pileup_skim");
TH1D *h_RD_h_nStubs_1_pileup_skim=(TH1D*)f_RD->Get("h_nStubs_1_pileup_skim");

TH1D *h_MC_h_nstubsPerModule_trackable=(TH1D*)f_MC->Get("h_nstubsPerModule_trackable");
TH1D *h_RD_h_nstubsPerModule_trackable=(TH1D*)f_RD->Get("h_nstubsPerModule_trackable");

TH1D *h_MC_h_nstubsPerModule_single_clean=(TH1D*)f_MC->Get("h_nstubsPerModule_single_clean");
TH1D *h_RD_h_nstubsPerModule_single_clean=(TH1D*)f_RD->Get("h_nstubsPerModule_single_clean");

TH1D *h_MC_h_nstubsPerModule_pileup234_skim=(TH1D*)f_MC->Get("h_nstubsPerModule_pileup234_skim");
TH1D *h_RD_h_nstubsPerModule_pileup234_skim=(TH1D*)f_RD->Get("h_nstubsPerModule_pileup234_skim");


array<TH1D*,12> h_MC_stub_per_mod_trackable,h_MC_stub_per_mod_single_clean,h_MC_stub_per_mod_pileup234,h_RD_stub_per_mod_trackable,h_RD_stub_per_mod_single_clean,h_RD_stub_per_mod_pileup234;

for(int m=0; m<12; m++){
h_MC_stub_per_mod_trackable.at(m)=(TH1D*)f_MC->Get(("stub_per_mod_trackable"+to_string(m)).c_str());
h_MC_stub_per_mod_single_clean.at(m)=(TH1D*)f_MC->Get(("stub_per_mod_single_clean"+to_string(m)).c_str());
h_MC_stub_per_mod_pileup234.at(m)=(TH1D*)f_MC->Get(("stub_per_mod_pileup234"+to_string(m)).c_str());
h_RD_stub_per_mod_trackable.at(m)=(TH1D*)f_RD->Get(("stub_per_mod_trackable"+to_string(m)).c_str());
h_RD_stub_per_mod_single_clean.at(m)=(TH1D*)f_RD->Get(("stub_per_mod_single_clean"+to_string(m)).c_str());
h_RD_stub_per_mod_pileup234.at(m)=(TH1D*)f_RD->Get(("stub_per_mod_pileup234"+to_string(m)).c_str());
}


TCanvas a0("a0","a0",700,700);
a0.Divide(1,3);
a0.cd(1);
h_MC_h_nstubsPerModule_trackable->SetLineColor(kRed);
h_MC_h_nstubsPerModule_trackable->Scale(1./h_MC_nTotStubs_trackable->GetEntries());//h_MC_h_nstubsPerModule_single_clean->GetEntries());
h_RD_h_nstubsPerModule_trackable->Scale(1./h_RD_nTotStubs_trackable->GetEntries());//h_RD_h_nstubsPerModule_single_clean->GetEntries());
h_RD_h_nstubsPerModule_trackable->SetMinimum(0.);
h_RD_h_nstubsPerModule_trackable->Draw("hist");
h_MC_h_nstubsPerModule_trackable->Draw("hist same");
gStyle->SetOptStat(0);
a0.cd(2);
h_MC_h_nstubsPerModule_single_clean->SetLineColor(kRed);
h_MC_h_nstubsPerModule_single_clean->Scale(1./h_MC_h_nTotStubs_single_clean->GetEntries());//h_MC_h_nstubsPerModule_single_clean->GetEntries());
h_RD_h_nstubsPerModule_single_clean->Scale(1./h_RD_h_nTotStubs_single_clean->GetEntries());//h_RD_h_nstubsPerModule_single_clean->GetEntries());
h_MC_h_nstubsPerModule_single_clean->Draw("hist");
h_RD_h_nstubsPerModule_single_clean->Draw("hist same");
gStyle->SetOptStat(0);
a0.cd(3);
h_MC_h_nstubsPerModule_pileup234_skim->SetLineColor(kRed);
h_MC_h_nstubsPerModule_pileup234_skim->Scale(1./h_MC_h_nTotStubs_pileup234_skim->GetEntries());//h_MC_h_nstubsPerModule_pileup234_skim->GetEntries());
h_RD_h_nstubsPerModule_pileup234_skim->Scale(1./h_RD_h_nTotStubs_pileup234_skim->GetEntries());//h_RD_h_nstubsPerModule_pileup234_skim->GetEntries());
h_MC_h_nstubsPerModule_pileup234_skim->Draw("hist");
h_RD_h_nstubsPerModule_pileup234_skim->Draw("hist same");
gStyle->SetOptStat(0);
a0.SaveAs("compare_skim/compare_local_3cmh_nstubsPerModule.pdf");


TCanvas a1("a1","a1",700,700);
a1.Divide(1,3);

a1.cd(1);
h_MC_nTotStubs_trackable->SetLineWidth(2);
h_MC_nTotStubs_trackable->Scale(1./h_MC_nTotStubs_trackable->GetEntries());
h_MC_nTotStubs_trackable->SetLineColor(kRed);
h_RD_nTotStubs_trackable->SetLineWidth(2);
h_RD_nTotStubs_trackable->Scale(1./h_RD_nTotStubs_trackable->GetEntries());
h_MC_nTotStubs_trackable->Draw("hist");
h_RD_nTotStubs_trackable->Draw("hist same");

a1.cd(2);
//h_MC_h_nStubs_0_trackable->Rebin(2);
//h_RD_h_nStubs_0_trackable->Rebin(2);
h_MC_h_nStubs_0_trackable->SetLineWidth(2);
h_MC_h_nStubs_0_trackable->Scale(1./h_MC_h_nStubs_0_trackable->GetEntries());
h_MC_h_nStubs_0_trackable->SetLineColor(kRed);
h_RD_h_nStubs_0_trackable->SetLineWidth(2);
h_RD_h_nStubs_0_trackable->Scale(1./h_RD_h_nStubs_0_trackable->GetEntries());
h_MC_h_nStubs_0_trackable->Draw("hist");
h_RD_h_nStubs_0_trackable->Draw("hist same");

a1.cd(3);

h_MC_h_nStubs_1_trackable->SetLineWidth(2);
h_MC_h_nStubs_1_trackable->SetLineColor(kRed);
h_MC_h_nStubs_1_trackable->Scale(1./h_MC_h_nStubs_1_trackable->GetEntries());
h_RD_h_nStubs_1_trackable->SetLineWidth(2);
h_RD_h_nStubs_1_trackable->Scale(1./h_RD_h_nStubs_1_trackable->GetEntries());
h_MC_h_nStubs_1_trackable->Draw("hist");
h_RD_h_nStubs_1_trackable->Draw("hist same");

a1.SaveAs("compare_skim/compare_local_3cmh_nStubs_1_trackable_bin2_bin2_eff97.pdf");



auto legend_1 = new TLegend(0.7,0.75,0.9,0.9);
   legend_1->AddEntry(h_MC_h_nStubs_0_single_clean,"MC","LEP");
   legend_1->AddEntry(h_RD_h_nStubs_0_single_clean,"RD","LEP");

auto legend_2 = new TLegend(0.7,0.75,0.9,0.9);
   legend_2->AddEntry(h_MC_h_nStubs_1_single_clean,"MC","LEP");
   legend_2->AddEntry(h_RD_h_nStubs_1_single_clean,"RD","LEP");

auto legend_3 = new TLegend(0.7,0.75,0.9,0.9);
legend_3->AddEntry(h_MC_h_nStubs_0_pileup234_skim,"MC","LEP");
legend_3->AddEntry(h_RD_h_nStubs_0_pileup234_skim,"RD","LEP");

auto legend_4 = new TLegend(0.7,0.75,0.9,0.9);
legend_4->AddEntry(h_MC_h_nStubs_1_pileup234_skim,"MC","LEP");
legend_4->AddEntry(h_RD_h_nStubs_1_pileup234_skim,"RD","LEP");


TCanvas a2("a2","a2",700,700);
a2.Divide(2,2);

a2.cd(1);
//h_MC_h_nStubs_0_single_clean->Rebin(2);
//h_RD_h_nStubs_0_single_clean->Rebin(2);
h_RD_h_nStubs_0_single_clean->GetXaxis()->SetRangeUser(0.,10.);
h_MC_h_nStubs_0_single_clean->GetXaxis()->SetRangeUser(0.,10.);
h_MC_h_nStubs_0_single_clean->SetLineWidth(2);
h_MC_h_nStubs_0_single_clean->Scale(1./h_MC_h_nStubs_0_single_clean->GetEntries());
h_MC_h_nStubs_0_single_clean->SetLineColor(kRed);
h_RD_h_nStubs_0_single_clean->SetLineWidth(2);
h_RD_h_nStubs_0_single_clean->Scale(1./h_RD_h_nStubs_0_single_clean->GetEntries());
h_MC_h_nStubs_0_single_clean->Draw("hist");
h_RD_h_nStubs_0_single_clean->Draw("hist same");
legend_1->Draw();
gStyle->SetOptStat(0);
a2.cd(2);
//h_MC_h_nStubs_1_single_clean->Rebin(2);
//h_RD_h_nStubs_1_single_clean->Rebin(2);
h_MC_h_nStubs_1_single_clean->GetXaxis()->SetRangeUser(0.,30.);
h_RD_h_nStubs_1_single_clean->GetXaxis()->SetRangeUser(0.,30.);
h_MC_h_nStubs_1_single_clean->SetLineWidth(2);
h_MC_h_nStubs_1_single_clean->Scale(1./h_MC_h_nStubs_1_single_clean->GetEntries());
h_MC_h_nStubs_1_single_clean->SetLineColor(kRed);
h_RD_h_nStubs_1_single_clean->SetLineWidth(2);
h_RD_h_nStubs_1_single_clean->Scale(1./h_RD_h_nStubs_1_single_clean->GetEntries());
h_MC_h_nStubs_1_single_clean->Draw("hist");
h_RD_h_nStubs_1_single_clean->Draw("hist same");
legend_2->Draw();
gStyle->SetOptStat(0);

a2.cd(3);
TH1D * h0 = (TH1D*) h_RD_h_nStubs_0_single_clean->Clone();
h0->Divide(h_MC_h_nStubs_0_single_clean);
h0->Draw("E");

a2.cd(4);
TH1D * h1 = (TH1D*) h_RD_h_nStubs_1_single_clean->Clone();
h1->Divide(h_MC_h_nStubs_1_single_clean);
h1->Draw("E");

a2.SaveAs("compare_skim/compare_local_3cmmes_min_h_nStubs_single_clean_bin2_eff97.pdf");


TCanvas a3("a3","a3",700,700);
a3.Divide(1,2);

a3.cd(1);
//h_MC_h_nStubs_0_pileup234_skim->SetMaximum(0.9);
//h_MC_h_nStubs_0_pileup234_skim->Rebin(2);
//h_RD_h_nStubs_0_pileup234_skim->Rebin(2);
h_MC_h_nStubs_0_pileup234_skim->GetXaxis()->SetRangeUser(0.,30.);
h_RD_h_nStubs_0_pileup234_skim->GetXaxis()->SetRangeUser(0.,30.);
h_MC_h_nStubs_0_pileup234_skim->SetLineWidth(2);
h_MC_h_nStubs_0_pileup234_skim->Scale(1./h_MC_h_nStubs_0_pileup234_skim->GetEntries());
h_MC_h_nStubs_0_pileup234_skim->SetLineColor(kRed);
h_RD_h_nStubs_0_pileup234_skim->SetLineWidth(2);
h_RD_h_nStubs_0_pileup234_skim->Scale(1./h_RD_h_nStubs_0_pileup234_skim->GetEntries());
h_MC_h_nStubs_0_pileup234_skim->Draw("hist");
h_RD_h_nStubs_0_pileup234_skim->Draw("hist same");
legend_3->Draw();
gStyle->SetOptStat(0);

a3.cd(2);
//h_MC_h_nStubs_1_pileup234_skim->Rebin(2);
//h_RD_h_nStubs_1_pileup234_skim->Rebin(2);
h_MC_h_nStubs_1_pileup234_skim->GetXaxis()->SetRangeUser(0.,35.);
h_RD_h_nStubs_1_pileup234_skim->GetXaxis()->SetRangeUser(0.,35.);
h_MC_h_nStubs_1_pileup234_skim->SetLineWidth(2);
h_MC_h_nStubs_1_pileup234_skim->Scale(1./h_MC_h_nStubs_1_pileup234_skim->GetEntries());
h_MC_h_nStubs_1_pileup234_skim->SetLineColor(kRed);
h_RD_h_nStubs_1_pileup234_skim->SetLineWidth(2);
h_RD_h_nStubs_1_pileup234_skim->Scale(1./h_RD_h_nStubs_1_pileup234_skim->GetEntries());
h_MC_h_nStubs_1_pileup234_skim->Draw("hist");
h_RD_h_nStubs_1_pileup234_skim->Draw("hist same");
legend_4->Draw();
gStyle->SetOptStat(0);

a3.cd(3);
TH1D * h2 = (TH1D*) h_RD_h_nStubs_0_pileup234_skim->Clone();
h2->Divide(h_MC_h_nStubs_0_pileup234_skim);
h2->Draw("E");

a3.cd(4);
TH1D * h3 = (TH1D*) h_RD_h_nStubs_1_pileup234_skim->Clone();
h3->Divide(h_MC_h_nStubs_1_pileup234_skim);
h3->Draw("E");

a3.SaveAs("compare_skim/compare_local_3cmmes_min_h_nStubs_pileup234_skim_bin2_eff97.pdf");


TCanvas a4("a4","a4",700,700);
a4.Divide(1,3);
a4.cd(1);
h_MC_h_nTotStubs_pileup_skim->SetLineWidth(2);
h_RD_h_nTotStubs_pileup_skim->Scale(1./h_RD_h_nTotStubs_pileup_skim->GetEntries());
h_RD_h_nTotStubs_pileup_skim->SetLineColor(kRed);
h_RD_h_nTotStubs_pileup_skim->SetLineWidth(2);
h_RD_h_nTotStubs_pileup_skim->Scale(1./h_RD_h_nTotStubs_pileup_skim->GetEntries());
h_RD_h_nTotStubs_pileup_skim->Draw("hist");
h_RD_h_nTotStubs_pileup_skim->Draw("hist same");

a4.cd(2);
//h_MC_h_nStubs_0_pileup_skim->Rebin(2);
//h_RD_h_nStubs_0_pileup_skim->Rebin(2);
h_MC_h_nStubs_0_pileup_skim->SetLineWidth(2);
h_MC_h_nStubs_0_pileup_skim->Scale(1./h_MC_h_nStubs_0_pileup_skim->GetEntries());
h_MC_h_nStubs_0_pileup_skim->SetLineColor(kRed);
h_RD_h_nStubs_0_pileup_skim->SetLineWidth(2);
h_RD_h_nStubs_0_pileup_skim->Scale(1./h_RD_h_nStubs_0_pileup_skim->GetEntries());
h_MC_h_nStubs_0_pileup_skim->Draw("hist");
h_RD_h_nStubs_0_pileup_skim->Draw("hist same");

a4.cd(3);
//h_MC_h_nStubs_1_pileup_skim->Rebin(2);
//h_RD_h_nStubs_1_pileup_skim->Rebin(2);
h_MC_h_nStubs_1_pileup_skim->SetLineWidth(2);
h_MC_h_nStubs_1_pileup_skim->Scale(1./h_MC_h_nStubs_1_pileup_skim->GetEntries());
h_MC_h_nStubs_1_pileup_skim->SetLineColor(kRed);
h_RD_h_nStubs_1_pileup_skim->SetLineWidth(2);
h_RD_h_nStubs_1_pileup_skim->Scale(1./h_RD_h_nStubs_1_pileup_skim->GetEntries());
h_RD_h_nStubs_1_pileup_skim->Draw("hist");
h_MC_h_nStubs_1_pileup_skim->Draw("hist same");

a4.SaveAs("compare_skim/compare_local_3cmh_nStubs_pileup_skim_bin2_eff97.pdf");


TCanvas a5("a5","a5",5000,5000);
a5.Divide(3,4);
for(int m=0; m<12; m++){
a5.cd(m+1);
h_MC_stub_per_mod_trackable.at(m)->Scale(1./h_MC_stub_per_mod_trackable.at(m)->GetEntries());
//h_stub_per_mod.at(m)->SetLineWidth(3);
h_MC_stub_per_mod_trackable.at(m)->SetLineColor(kRed);
h_RD_stub_per_mod_trackable.at(m)->Scale(1./h_RD_stub_per_mod_trackable.at(m)->GetEntries());
//h_RD_stub_per_mod_trackable.at(m)->GetXaxis()->SetRangeUser(0.,55.);
h_RD_stub_per_mod_trackable.at(m)->Draw("hist");
h_MC_stub_per_mod_trackable.at(m)->Draw("hist same");
gPad->SetLogy();
}
a5.SaveAs("compare_skim/compare_local_3cm_stub_per_mod_eff97.pdf");


TCanvas a6("a6","a6",1000,1000);
a6.Divide(3,4);
for(int m=0; m<12; m++){
a6.cd(m+1);
h_MC_stub_per_mod_single_clean.at(m)->Scale(1./h_MC_stub_per_mod_single_clean.at(m)->GetEntries());
//h_stub_per_mod.at(m)->SetLineWidth(3);
h_MC_stub_per_mod_single_clean.at(m)->SetLineColor(kRed);
h_MC_stub_per_mod_single_clean.at(m)->Draw("hist");
h_RD_stub_per_mod_single_clean.at(m)->Scale(1./h_RD_stub_per_mod_single_clean.at(m)->GetEntries());
h_RD_stub_per_mod_single_clean.at(m)->Draw("hist same");
gPad->SetLogy();
}
a6.SaveAs("compare_skim/compare_local_3cm_stub_per_mod_single_clean_eff97.pdf");


TCanvas a7("a7","a7",1000,1000);
a7.Divide(3,4);
for(int m=0; m<12; m++){
a7.cd(m+1);
h_MC_stub_per_mod_pileup234.at(m)->Scale(1./h_MC_stub_per_mod_pileup234.at(m)->GetEntries());
//h_stub_per_mod.at(m)->SetLineWidth(3);
h_MC_stub_per_mod_pileup234.at(m)->SetLineColor(kRed);
h_MC_stub_per_mod_pileup234.at(m)->Draw("hist");
h_RD_stub_per_mod_pileup234.at(m)->Scale(1./h_RD_stub_per_mod_pileup234.at(m)->GetEntries());
h_RD_stub_per_mod_pileup234.at(m)->Draw("hist same");
gPad->SetLogy();
}
a7.SaveAs("compare_skim/compare_local_3cm_stub_per_mod_pileup234_eff97.pdf");


}


