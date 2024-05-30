#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include <TGraphErrors.h>
#include "TCanvas.h"
#include <TStyle.h>


void plots_pre_post_vrtx(){


array<TFile*,14> f_pre_2hit_el, f_post_2hit_el;
array<TFile*,14> f_pre_2hit_mu, f_post_2hit_mu;

for(int m=0; m<14; m++){
f_pre_2hit_el.at(m)=TFile::Open(Form("/home/espedica/macros_fairmu/clean_codes/validation/res_prevrtx_clones_el_%d_2hit_LO_newTuple_MS.root",static_cast<char>(m)));
f_post_2hit_el.at(m)=TFile::Open(Form("/home/espedica/macros_fairmu/clean_codes/validation/res_vrtx_clones_el_%d_2hit_LO_newTuple_MS.root",static_cast<char>(m)));
}


for(int m=0; m<14; m++){
f_pre_2hit_mu.at(m)=TFile::Open(Form("/home/espedica/macros_fairmu/clean_codes/validation/res_prevrtx_clones_AngleMu_%d_2hit_LO_newTuple_MS.root",static_cast<char>(m)));
f_post_2hit_mu.at(m)=TFile::Open(Form("/home/espedica/macros_fairmu/clean_codes/validation/res_vrtx_clones_AngleMu_%d_2hit_LO_newTuple_MS.root",static_cast<char>(m)));
}

TH1::SetDefaultSumw2(kTRUE);


array <TH1D*,14> h_pre_2hit_el, h_post_2hit_el;
array <TH1D*,14> h_pre_2hit_mu, h_post_2hit_mu;

for(int m=0; m<14; m++){

        h_pre_2hit_el.at(m)=(TH1D*)f_pre_2hit_el.at(m)->Get(Form("res%d_el_pre_clones",static_cast<char>(m)));
        h_post_2hit_el.at(m)=(TH1D*)f_post_2hit_el.at(m)->Get(Form("res%d_vrtx_el_clones",static_cast<char>(m)));
}

for(int m=0; m<14; m++){

        h_pre_2hit_mu.at(m)=(TH1D*)f_pre_2hit_mu.at(m)->Get(Form("res%d_mu_pre_clones",static_cast<char>(m)));
        h_post_2hit_mu.at(m)=(TH1D*)f_post_2hit_mu.at(m)->Get(Form("res%d_vrtx_mu_clones",static_cast<char>(m)));

}




   const Int_t NBINS_mu = 14;
   Double_t edges_mu[NBINS_mu + 1] = {0.00005,0.00015,0.00025,0.00035,0.00045,0.00055,0.00065,0.00075,0.00085,0.00095,0.0015,0.0025,0.0035,0.0045};

   const Int_t NBINS = 14;
   Double_t edges_el[NBINS + 1] = {0.0005,0.0015,0.0025,0.0035,0.0045,0.0055,0.0065,0.0075,0.0085,0.0095, 0.0125, 0.0175, 0.0225, 0.0285};



auto h_sigma_mu_pre_2hit_clones = new TGraphErrors();
auto h_sigma_el_pre_2hit_clones = new TGraphErrors();
auto h_sigma_mu_post_2hit_clones = new TGraphErrors();
auto h_sigma_el_post_2hit_clones = new TGraphErrors();

auto sigma_el_diff_2hit = new TGraphErrors();
auto sigma_mu_diff_2hit = new TGraphErrors();

double sigma_post_el2[NBINS]={0.};
double sigma_pre_el2[NBINS]={0.};
double sigma_post_mu2[NBINS_mu]={0.};
double sigma_pre_mu2[NBINS_mu]={0.};



TCanvas n2("n2","n2",2800,2800);
n2.Divide(4,4);
for(int m=0; m<14; m++){
n2.cd(m+1);
TF1 *g1 = new TF1("g1", "gaus");
if(m<4) h_pre_2hit_el.at(m)->Fit("g1","R","",-0.0005,0.0005);
else if(m==4) h_pre_2hit_el.at(m)->Fit("g1","R","",-0.002,0.002);
else h_pre_2hit_el.at(m)->Fit("g1");
h_sigma_el_pre_2hit_clones->SetPoint(m,edges_el[m]*1000,g1->GetParameter(2)*1000);
h_sigma_el_pre_2hit_clones->SetPointError(m,0.,g1->GetParError(2)*1000);
sigma_pre_el2[m]=g1->GetParameter(2)*1000;
gStyle->SetStatH(0.1);
gStyle->SetOptFit(0001);
gStyle->SetOptStat(1101);
}
n2.SaveAs("/home/espedica/macros_fairmu/clean_codes/validation/results/res_el_pre_2hit_newTuple_MS.pdf");


TCanvas n3("n3","n3",1400,1400);
n3.Divide(2,2);
n3.cd(1);
TGaxis::SetMaxDigits(3);
TF1 *g = new TF1("g", "gaus");
h_pre_2hit_el.at(0)->Fit("g","R","",-0.0003,0.0003);
h_pre_2hit_el.at(0)->GetXaxis()->SetTitle("rad");
h_pre_2hit_el.at(0)->SetLineColor(kAzure+7);
h_pre_2hit_el.at(0)->Draw("hist same");
n2.cd(2);
h_pre_2hit_el.at(2)->Fit("g","R","",-0.0004,0.0004);;
h_pre_2hit_el.at(2)->GetXaxis()->SetTitle("rad");
h_pre_2hit_el.at(2)->SetLineColor(kAzure+7);
h_pre_2hit_el.at(2)->Draw("hist same");
n2.cd(3);
h_pre_2hit_el.at(9)->Fit("g");
h_pre_2hit_el.at(9)->GetXaxis()->SetTitle("rad");
h_pre_2hit_el.at(9)->SetLineColor(kAzure+7);
h_pre_2hit_el.at(9)->Draw("hist same");
n2.cd(4);
h_pre_2hit_el.at(13)->Fit("g");
h_pre_2hit_el.at(13)->GetXaxis()->SetTitle("rad");
h_pre_2hit_el.at(13)->SetLineColor(kAzure+7);
h_pre_2hit_el.at(13)->Draw("hist same");
gStyle->SetStatH(0.1);
gStyle->SetOptFit(0001);
gStyle->SetOptStat(1101);

//n3.SaveAs("/home/espedica/macros_fairmu/clean_codes/validation/results/res_el_pre_2hit_newTuple_MS.pdf");





TCanvas n4("n4","n4",2800,2800);
n4.Divide(4,4);
for(int m=0; m<14; m++){
n4.cd(m+1);
TGaxis::SetMaxDigits(3);
TF1 *g1 = new TF1("g1", "gaus");
if(m<4) h_post_2hit_el.at(m)->Fit("g1","R","",-0.0005,0.0005);
else if(m==4) h_post_2hit_el.at(m)->Fit("g1","R","",-0.002,0.002);
else h_post_2hit_el.at(m)->Fit("g1");
h_sigma_el_post_2hit_clones->SetPoint(m,edges_el[m]*1000,g1->GetParameter(2)*1000);
h_sigma_el_post_2hit_clones->SetPointError(m,0.,g1->GetParError(2)*1000);
sigma_post_el2[m]=g1->GetParameter(2)*1000;
gStyle->SetStatH(0.1);
gStyle->SetOptFit(0001);
gStyle->SetOptStat(1101);
}
n4.SaveAs("/home/espedica/macros_fairmu/clean_codes/validation/results/res_el_post_2hit_newTuple_MS.pdf");





TCanvas n7("n7","n7",2800,2800);
n7.Divide(4,4);
for(int m=0; m<14; m++){
n7.cd(m+1);
TF1 *g1 = new TF1("g1", "gaus");
h_pre_2hit_mu.at(m)->Fit("g1");
h_sigma_mu_pre_2hit_clones->SetPoint(m,edges_mu[m]*1000,g1->GetParameter(2)*1000);
h_sigma_mu_pre_2hit_clones->SetPointError(m,0.,g1->GetParError(2)*1000);
sigma_pre_mu2[m]=g1->GetParameter(2)*1000;
}
n7.SaveAs("/home/espedica/macros_fairmu/clean_codes/validation/results/res_mu_pre_2hit_newTuple_MS.pdf");


TCanvas n2_mu("n2_mu","n2_mu",1400,1400);
n2_mu.Divide(2,2);
n2_mu.cd(1);
TF1 *g2 = new TF1("g2", "gaus");
h_pre_2hit_mu.at(2)->Fit("g2");
h_pre_2hit_mu.at(2)->SetLineColor(kAzure+7);
h_pre_2hit_mu.at(2)->GetXaxis()->SetTitle("rad");
h_pre_2hit_mu.at(2)->Draw("hist same");
n2_mu.cd(2);
h_pre_2hit_mu.at(9)->Fit("g2","R","",-0.00012,0.00012);
h_pre_2hit_mu.at(9)->SetLineColor(kAzure+7);
h_pre_2hit_mu.at(9)->GetXaxis()->SetTitle("rad");
//h_pre_2hit_mu.at(9)->Rebin(2);
h_pre_2hit_mu.at(9)->Draw("hist same");
n2_mu.cd(3);
h_pre_2hit_mu.at(11)->Fit("g2","R","",-0.0002,0.0002);
h_pre_2hit_mu.at(11)->SetLineColor(kAzure+7);
h_pre_2hit_mu.at(11)->GetXaxis()->SetTitle("rad");
//h_pre_2hit_mu.at(11)->Rebin(2);
h_pre_2hit_mu.at(11)->Draw("hist same");
n2_mu.cd(4);
h_pre_2hit_mu.at(13)->Rebin(2);
h_pre_2hit_mu.at(13)->SetMinimum(1.);
h_pre_2hit_mu.at(13)->Fit("g2");
h_pre_2hit_mu.at(13)->SetLineColor(kAzure+7);
h_pre_2hit_mu.at(13)->GetXaxis()->SetTitle("rad");
h_pre_2hit_mu.at(13)->Draw("hist same");
gStyle->SetStatH(0.1);
gStyle->SetOptFit(0001);
gStyle->SetOptStat(1101);

//n2_mu.SaveAs("/home/espedica/macros_fairmu/clean_codes/validation/results/res_pre_AngleMu_2hit_newTuple_MS.pdf");


TCanvas n8("n8","n8",2800,2800);
n8.Divide(4,4);
for(int m=0; m<14; m++){
n8.cd(m+1);
TF1 *g1 = new TF1("g1", "gaus");
if(m>6)h_post_2hit_mu.at(m)->Fit("g1","R","",-0.0002,0.0002);
else h_post_2hit_mu.at(m)->Fit("g1");
h_sigma_mu_post_2hit_clones->SetPoint(m,edges_mu[m]*1000,g1->GetParameter(2)*1000);
h_sigma_mu_post_2hit_clones->SetPointError(m,0.,g1->GetParError(2)*1000);
sigma_post_mu2[m]=g1->GetParameter(2)*1000;
}
n8.SaveAs("/home/espedica/macros_fairmu/clean_codes/validation/results/res_mu_post_2hit_newTuple_MS.pdf");


h_sigma_mu_pre_2hit_clones->SetName("mu_pre");
h_sigma_el_pre_2hit_clones->SetName("el_pre");
h_sigma_mu_post_2hit_clones->SetName("mu_post");
h_sigma_el_post_2hit_clones->SetName("el_post");


h_sigma_el_pre_2hit_clones->SaveAs("/home/espedica/macros_fairmu/clean_codes/validation/results/sigma_el_pre_2hit_newTuple_MS.root");
h_sigma_el_post_2hit_clones->SaveAs("/home/espedica/macros_fairmu/clean_codes/validation/results/sigma_el_post_2hit_newTuple_MS.root");
h_sigma_mu_pre_2hit_clones->SaveAs("/home/espedica/macros_fairmu/clean_codes/validation/results/sigma_mu_pre_2hit_newTuple_MS.root");
h_sigma_mu_post_2hit_clones->SaveAs("/home/espedica/macros_fairmu/clean_codes/validation/results/sigma_mu_post_2hit_newTuple_MS.root");



for(int m=0; m<14; m++){
sigma_el_diff_2hit->SetPoint(m,edges_el[m]*1000,(sigma_post_el2[m]-sigma_pre_el2[m])/sigma_pre_el2[m]);
if(m!=0)sigma_mu_diff_2hit->SetPoint(m,edges_mu[m]*1000,(sigma_post_mu2[m]-sigma_pre_mu2[m])/sigma_pre_mu2[m]);
}



sigma_el_diff_2hit->SetName("el_2");
sigma_el_diff_2hit->GetYaxis()->SetTitle("(#sigma_post - #sigma_pre)/#sigma_pre");
sigma_el_diff_2hit->GetXaxis()->SetTitle("#theta_el (mrad)");
sigma_mu_diff_2hit->SetName("mu_2");
sigma_mu_diff_2hit->GetYaxis()->SetTitle("(#sigma_post - #sigma_pre)/#sigma_pre");
sigma_mu_diff_2hit->GetXaxis()->SetTitle("#theta_#mu (mrad)");


sigma_el_diff_2hit->SaveAs("/home/espedica/macros_fairmu/clean_codes/validation/results/el_2_newTuple_MS.root"); 
sigma_mu_diff_2hit->SaveAs("/home/espedica/macros_fairmu/clean_codes/validation/results/mu_2_newTuple_MS.root");

}


