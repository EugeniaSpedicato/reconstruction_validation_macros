
#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TStyle.h>


void divide(){
TFile *f_RD = TFile::Open("d_eff_RD.root");
TFile *f_MC = TFile::Open("d_eff_MC.root");

TFile *f_thmu_RD = TFile::Open("theta_mu_RD.root");
TFile *f_thmu_MC = TFile::Open("theta_mu_MC.root");

TFile *f_the_RD = TFile::Open("theta_e_RD.root");
TFile *f_the_MC = TFile::Open("theta_e_MC.root");

TFile *f_thmu_gen_MC = TFile::Open("theta_mu_gen_MC.root");
TFile *f_the_gen_MC = TFile::Open("theta_e_gen_MC.root");

TH1::SetDefaultSumw2(kTRUE);

TH1D *h_RD=(TH1D*)f_RD->Get("d_eff_real");
TH1D *h_MC=(TH1D*)f_MC->Get("d_eff_MC");
//TH1D * h1 = (TH1D*) h_RD->Clone();
//h1->Divide(h_MC);


TH1D *h_theta_mu_RD=(TH1D*)f_thmu_RD->Get("theta_mu");
TH1D *h_theta_mu_MC=(TH1D*)f_thmu_MC->Get("theta_mu");
TH1D *h_theta_mu_gen_MC=(TH1D*)f_thmu_gen_MC->Get("theta_mu_gen");

TH1D *h_theta_e_RD=(TH1D*)f_the_RD->Get("theta_e");
TH1D *h_theta_e_MC=(TH1D*)f_the_MC->Get("theta_e");
TH1D *h_theta_e_gen_MC=(TH1D*)f_the_gen_MC->Get("theta_e_gen");


TH1D * h6 = (TH1D*) h_theta_mu_MC->Clone();
h6->Divide(h_theta_mu_gen_MC);
TH1D * h7 = (TH1D*) h_theta_e_MC->Clone();
h7->Divide(h_theta_e_gen_MC);


TCanvas a("a","a",700,700);
a.Divide(2,3);
a.cd(1);

h_theta_e_gen_MC->SetLineColor(kAzure+7);
h_theta_e_gen_MC->Draw("E");
h_theta_e_MC->Draw("E same");
gPad->SetLogy();
gStyle->SetOptStat(0);

a.cd(2);

h_theta_mu_gen_MC->SetLineColor(kAzure+7);
h_theta_mu_gen_MC->Draw("E");
h_theta_mu_MC->Draw("E same");
gPad->SetLogy();
gStyle->SetOptStat(0);

a.cd(3);
cout << "theta el. reco:" << endl;
for(int i=0; i<10; i++){
double norm=8.70408e+06 *5.5*1E+23*3.*1E-30*1325.*h_theta_e_MC->GetBinContent(i)/121846.;
h_theta_e_MC->SetBinContent(i,norm);
cout << "bin " << i << ") number of norm. event MC: " << norm << " VS number of real data events " << h_theta_e_RD->GetBinContent(i) << endl;
}
//h_theta_e_MC->Scale(8.70408e+06*5.5*1E+23*3.*1E-30*1325.);
h_theta_e_MC->SetMinimum(1.);
h_theta_e_MC->Draw("E");
h_theta_e_RD->SetLineColor(kPink+10);
h_theta_e_RD->Draw("E same");
gPad->SetLogy();
gStyle->SetOptStat(0);

a.cd(4);
/*for(int i=0; i<10; i++){
double norm=8.70408e+06 *5.5*1E+23*3.*1E-30*1325.*h_MC->GetBinContent(i)/121846.;
h_MC->SetBinContent(i,norm);
cout << "bin " << i << ") number of norm. event MC: " << norm << " VS number of real data events " << h_RD->GetBinContent(i) << endl;
}
h_MC->SetMinimum(1.);
h_MC->Draw("E");
h_RD->SetLineColor(kPink+10);
h_RD->Draw("E same");
gPad->SetLogy();
gStyle->SetOptStat(0000);*/

cout << "theta el. gen:" << endl;
for(int i=0; i<10; i++){
double norm=8.70408e+06 *5.5*1E+23*3.*1E-30*1325.*h_theta_e_gen_MC->GetBinContent(i)/121846.;
h_theta_e_gen_MC->SetBinContent(i,norm);
cout << "bin " << i << ") number of norm. event MC: " << norm << " VS number of real data events " << h_theta_e_RD->GetBinContent(i) << endl;
}
//h_theta_e_MC->Scale(8.70408e+06*5.5*1E+23*3.*1E-30*1325.);
h_theta_e_gen_MC->SetMinimum(1.);
h_theta_e_gen_MC->Draw("E");
h_theta_e_RD->SetLineColor(kPink+10);
h_theta_e_RD->Draw("E same");
gPad->SetLogy();
gStyle->SetOptStat(0);

a.cd(5);
cout << "theta mu. reco:" << endl;
for(int i=0; i<20; i++){
double norm=8.70408e+06 *5.5*1E+23*3.*1E-30*1325.*h_theta_mu_MC->GetBinContent(i)/121846.;
h_theta_mu_MC->SetBinContent(i,norm);
cout << "bin " << i << ") number of norm. event MC: " << norm << " VS number of real data events " << h_theta_mu_RD->GetBinContent(i) << endl;
}

//h_theta_mu_MC->Scale(8.70408e+06*5.5*1E+23*3.*1E-30*1325.);
h_theta_mu_MC->SetMinimum(1.);
h_theta_mu_MC->Draw("E");
h_theta_mu_RD->SetLineColor(kPink+10);
//h_theta_mu_RD->GetXaxis()->SetRange(0,10);
h_theta_mu_RD->Draw("E same");
gPad->SetLogy();
gStyle->SetOptStat(0);

a.cd(6);
cout << "theta mu. gen:" << endl;
for(int i=0; i<20; i++){
double norm=8.70408e+06 *5.5*1E+23*3.*1E-30*1325.*h_theta_mu_gen_MC->GetBinContent(i)/121846.;
h_theta_mu_gen_MC->SetBinContent(i,norm);
cout << "bin " << i << ") number of norm. event MC: " << norm << " VS number of real data events " << h_theta_mu_RD->GetBinContent(i) << endl;
}

//h_theta_mu_MC->Scale(8.70408e+06*5.5*1E+23*3.*1E-30*1325.);
h_theta_mu_gen_MC->SetMinimum(1.);
h_theta_mu_gen_MC->Draw("E");
h_theta_mu_RD->SetLineColor(kPink+10);
//h_theta_mu_RD->GetXaxis()->SetRange(0,10);
h_theta_mu_RD->Draw("E same");
gPad->SetLogy();
gStyle->SetOptStat(0);

a.SaveAs("theta_MC_RD.pdf");



TH1D * h2 = (TH1D*) h_theta_mu_RD->Clone();
h2->Divide(h_theta_mu_MC);
TH1D * h3 = (TH1D*) h_theta_e_RD->Clone();
h3->Divide(h_theta_e_MC);

TH1D * h4 = (TH1D*) h_theta_mu_RD->Clone();
h4->Divide(h_theta_mu_gen_MC);
TH1D * h5 = (TH1D*) h_theta_e_RD->Clone();
h5->Divide(h_theta_e_gen_MC);


TCanvas d("d","d",650,700);
d.Divide(2,3);
d.cd(1);
h2->Draw("E");
d.cd(2);
h4->Draw("E");
d.cd(3);
h3->Draw("E");
d.cd(4);
h5->Draw("E");
d.cd(5);
h6->Draw("E");
d.cd(6);
h7->Draw("E");
d.SaveAs("h_d_MC_RD.pdf");

}
