
#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TStyle.h>


void fake(){

TFile *f_good_el0=TFile::Open("proposal/good_el_0hit_LO.root");
TFile *f_good_mu0=TFile::Open("proposal/good_mu_0hit_LO.root");
TFile *f_fake0=TFile::Open("proposal/fake_0hit_LO.root");

TFile *f_good_el1=TFile::Open("proposal/good_el_1hitFirstModules_LO.root");
TFile *f_good_mu1=TFile::Open("proposal/good_mu_1hitFirstModules_LO.root");
TFile *f_fake1=TFile::Open("proposal/fake_1hitFirstModules_LO.root");


TFile *f_good_el2=TFile::Open("proposal/good_el_2hitFirstModules_LO.root");
TFile *f_good_mu2=TFile::Open("proposal/good_mu_2hitFirstModules_LO.root");
TFile *f_fake2=TFile::Open("proposal/fake_2hitFirstModules_LO.root");


TH1::SetDefaultSumw2(kTRUE);

TH1D *h_good_el0=(TH1D*)f_good_el0->Get("theta_e_good");
TH1D *h_good_mu0=(TH1D*)f_good_mu0->Get("theta_mu_good");
TH1D *h_fake0=(TH1D*)f_fake0->Get("theta_fake");

TH1D *h_good_el1=(TH1D*)f_good_el1->Get("theta_e_good");
TH1D *h_good_mu1=(TH1D*)f_good_mu1->Get("theta_mu_good");
TH1D *h_fake1=(TH1D*)f_fake1->Get("theta_fake");

TH1D *h_good_el2=(TH1D*)f_good_el2->Get("theta_e_good");
TH1D *h_good_mu2=(TH1D*)f_good_mu2->Get("theta_mu_good");
TH1D *h_fake2=(TH1D*)f_fake2->Get("theta_fake");

 auto legend_e = new TLegend(0.6,0.1,0.9,0.25);
   legend_e->AddEntry(h_good_el0,"0 hit","LEP");
   legend_e->AddEntry(h_good_el1,"1 hit","LEP");
   legend_e->AddEntry(h_good_el2,"2 hit","LEP");

 auto legend_mu = new TLegend(0.1,0.1,0.4,0.25);
   legend_mu->AddEntry(h_good_mu0,"0 hit","LEP");
   legend_mu->AddEntry(h_good_mu1,"1 hit","LEP");
   legend_mu->AddEntry(h_good_mu2,"2 hit","LEP");

 auto legend_f = new TLegend(0.7,0.6,0.85,0.75);
   legend_f->AddEntry(h_fake0,"0 hit","LEP");
   legend_f->AddEntry(h_fake1,"1 hit","LEP");
   legend_f->AddEntry(h_fake2,"2 hit","LEP");

TCanvas a1("a1","a1",700,700);
/*a1.Divide(3,1);
a1.cd(1);
h_good_mu1->SetTitle("Good tracks rate as a function of mu scattering angle");
h_good_mu1->GetXaxis()->SetTitle("Muon angle[rad]");
h_good_mu1->GetYaxis()->SetTitle("Rate");
h_good_mu2->SetLineColor(kAzure+7);
h_good_mu1->SetMinimum(0.6);
TGaxis::SetMaxDigits(3);
h_good_mu1->SetLineColor(kGreen+1);
h_good_mu1->Draw("E");
h_good_mu0->SetLineColor(kOrange+10);
h_good_mu0->Draw("E same");
h_good_mu2->Draw("E same");

legend_mu->Draw();
gStyle->SetOptStat(0);

a1.cd(2);
h_good_el1->SetTitle("Good tracks rate as a function of e- scattering angle");
h_good_el1->GetXaxis()->SetTitle("Electron angle[rad]");
h_good_el1->GetYaxis()->SetTitle("Rate");
h_good_el2->SetLineColor(kAzure+7);
h_good_el1->SetMinimum(0.6);
TGaxis::SetMaxDigits(3);
h_good_el1->SetLineColor(kGreen+1);
h_good_el1->Draw("E");
h_good_el0->SetLineColor(kOrange+10);
h_good_el0->Draw("E same");
h_good_el2->Draw("E same");

legend_e->Draw();
gStyle->SetOptStat(0);

a1.cd(3);
*/
h_fake2->SetTitle("Fake tracks rate as a function of particle scattering angle");
h_fake2->GetXaxis()->SetTitle("Particle angle[rad]");
h_fake2->GetYaxis()->SetTitle("Fake rate");
h_fake2->SetLineColor(kAzure+7);
h_fake2->SetMinimum(0.6);
TGaxis::SetMaxDigits(3);
h_fake2->Draw("E");
h_fake1->SetLineColor(kGreen+1);
h_fake1->Draw("E same");
h_fake0->SetLineColor(kOrange+10);
h_fake0->Draw("E same");

legend_f->Draw();
gStyle->SetOptStat(0);
a1.SaveAs("proposal/good_fakeFirstModules_rate_LO_012hit_02cut.pdf");






}


