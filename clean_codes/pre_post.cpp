
#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include <TGraphErrors.h>
#include "TCanvas.h"
#include <TStyle.h>


void pre_post(){

array<TFile*,14> f_pre_1hit_el, f_post_1hit_el, f_pre_2hit_el, f_post_2hit_el;
array<TFile*,14> f_pre_1hit_mu, f_post_1hit_mu, f_pre_2hit_mu, f_post_2hit_mu;

for(int m=0; m<14; m++){

f_pre_1hit_el.at(m)=TFile::Open(Form("/home/espedica/macros_fairmu/clean_codes/validation/res_prevrtx_clones_el_%d_1hit_LO.root",static_cast<char>(m)));
f_post_1hit_el.at(m)=TFile::Open(Form("/home/espedica/macros_fairmu/clean_codes/validation/res_vrtx_clones_el_%d_1hit_LO.root",static_cast<char>(m)));
f_pre_2hit_el.at(m)=TFile::Open(Form("/home/espedica/macros_fairmu/clean_codes/validation/res_prevrtx_clones_el_%d_2hit_LO_oldTuple_MS.root",static_cast<char>(m)));
f_post_2hit_el.at(m)=TFile::Open(Form("/home/espedica/macros_fairmu/clean_codes/validation/res_vrtx_clones_el_%d_2hit_LO_oldTuple_MS.root",static_cast<char>(m)));
}

for(int m=0; m<14; m++){
f_pre_1hit_mu.at(m)=TFile::Open(Form("/home/espedica/macros_fairmu/clean_codes/validation/res_prevrtx_clones_AngleMu_%d_1hit_LO.root",static_cast<char>(m)));
f_post_1hit_mu.at(m)=TFile::Open(Form("/home/espedica/macros_fairmu/clean_codes/validation/res_vrtx_clones_AngleMu_%d_1hit_LO.root",static_cast<char>(m)));
f_pre_2hit_mu.at(m)=TFile::Open(Form("/home/espedica/macros_fairmu/clean_codes/validation/res_prevrtx_clones_AngleMu_%d_2hit_LO_oldTuple_MS.root",static_cast<char>(m)));
f_post_2hit_mu.at(m)=TFile::Open(Form("/home/espedica/macros_fairmu/clean_codes/validation/res_vrtx_clones_AngleMu_%d_2hit_LO_oldTuple_MS.root",static_cast<char>(m)));
}


TH1::SetDefaultSumw2(kTRUE);

array <TH1D*,14> h_pre_1hit_el, h_post_1hit_el, h_pre_2hit_el, h_post_2hit_el;
array <TH1D*,14> h_pre_1hit_mu, h_post_1hit_mu, h_pre_2hit_mu, h_post_2hit_mu;


for(int m=0; m<14; m++){

h_pre_1hit_el.at(m)=(TH1D*)f_pre_1hit_el.at(m)->Get(Form("res%d_el_pre_clones",static_cast<char>(m)));
h_post_1hit_el.at(m)=(TH1D*)f_post_1hit_el.at(m)->Get(Form("res%d_vrtx_el_clones",static_cast<char>(m)));
h_pre_2hit_el.at(m)=(TH1D*)f_pre_2hit_el.at(m)->Get(Form("res%d_el_pre_clones",static_cast<char>(m)));
h_post_2hit_el.at(m)=(TH1D*)f_post_2hit_el.at(m)->Get(Form("res%d_vrtx_el_clones",static_cast<char>(m)));
}

for(int m=0; m<14; m++){
h_pre_1hit_mu.at(m)=(TH1D*)f_pre_1hit_mu.at(m)->Get(Form("res%d_mu_pre_clones",static_cast<char>(m)));
h_post_1hit_mu.at(m)=(TH1D*)f_post_1hit_mu.at(m)->Get(Form("res%d_vrtx_mu_clones",static_cast<char>(m)));
h_pre_2hit_mu.at(m)=(TH1D*)f_pre_2hit_mu.at(m)->Get(Form("res%d_mu_pre_clones",static_cast<char>(m)));
h_post_2hit_mu.at(m)=(TH1D*)f_post_2hit_mu.at(m)->Get(Form("res%d_vrtx_mu_clones",static_cast<char>(m)));

}


   const Int_t NBINS_mu = 14;
   Double_t edges_mu[NBINS_mu + 1] = {0.00005,0.00015,0.00025,0.00035,0.00045,0.00055,0.00065,0.00075,0.00085,0.00095,0.0015,0.0025,0.0035,0.0045};

   const Int_t NBINS = 14;
   Double_t edges_el[NBINS + 1] = {0.0005,0.0015,0.0025,0.0035,0.0045,0.0055,0.0065,0.0075,0.0085,0.0095, 0.0125, 0.0175, 0.0225, 0.0285};



auto h_sigma_mu_pre_1hit_clones = new TGraphErrors();
auto h_sigma_el_pre_1hit_clones = new TGraphErrors();
auto h_sigma_mu_post_1hit_clones = new TGraphErrors();
auto h_sigma_el_post_1hit_clones = new TGraphErrors();
auto h_sigma_mu_pre_2hit_clones = new TGraphErrors();
auto h_sigma_el_pre_2hit_clones = new TGraphErrors();
auto h_sigma_mu_post_2hit_clones = new TGraphErrors();
auto h_sigma_el_post_2hit_clones = new TGraphErrors();

auto sigma_el_diff_1hit = new TGraphErrors();
auto sigma_mu_diff_1hit = new TGraphErrors();
auto sigma_el_diff_2hit = new TGraphErrors();
auto sigma_mu_diff_2hit = new TGraphErrors();

double sigma_post_el1[NBINS]={0.};
double sigma_pre_el1[NBINS]={0.};
double sigma_post_mu1[NBINS_mu]={0.};
double sigma_pre_mu1[NBINS_mu]={0.};
double sigma_post_el2[NBINS]={0.};
double sigma_pre_el2[NBINS]={0.};
double sigma_post_mu2[NBINS_mu]={0.};
double sigma_pre_mu2[NBINS_mu]={0.};



TCanvas n("n","n",2800,2800);
n.Divide(4,4);
for(int m=0; m<14; m++){
n.cd(m+1);
TF1 *g1 = new TF1("g1", "gaus");
if(m<4) h_pre_1hit_el.at(m)->Fit("g1","R","",-0.0005,0.0005);
else if(m==4) h_pre_1hit_el.at(m)->Fit("g1","R","",-0.002,0.002);
else h_pre_1hit_el.at(m)->Fit("g1");
h_sigma_el_pre_1hit_clones->SetPoint(m,edges_el[m]*1000,g1->GetParameter(2)*1000);
h_sigma_el_pre_1hit_clones->SetPointError(m,0.,g1->GetParError(2)*1000);
sigma_pre_el1[m]=g1->GetParameter(2)*1000;

gStyle->SetStatH(0.1);
gStyle->SetOptFit(0001);
gStyle->SetOptStat(1101);
}
n.SaveAs("/home/espedica/macros_fairmu/clean_codes/validation/results/res_el_pre_1hit_oldTuple.pdf");


TCanvas n1("n1","n1",2800,2800);
n1.Divide(4,4);
for(int m=0; m<14; m++){
n1.cd(m+1);
TGaxis::SetMaxDigits(3);
TF1 *g1 = new TF1("g1", "gaus");
if(m<4) h_post_1hit_el.at(m)->Fit("g1","R","",-0.0005,0.0005);
else if(m==4) h_post_1hit_el.at(m)->Fit("g1","R","",-0.002,0.002);
else h_post_1hit_el.at(m)->Fit("g1");
h_sigma_el_post_1hit_clones->SetPoint(m,edges_el[m]*1000,g1->GetParameter(2)*1000);
h_sigma_el_post_1hit_clones->SetPointError(m,0.,g1->GetParError(2)*1000);
sigma_post_el1[m]=g1->GetParameter(2)*1000;
gStyle->SetStatH(0.1);
gStyle->SetOptFit(0001);
gStyle->SetOptStat(1101);
}
n1.SaveAs("/home/espedica/macros_fairmu/clean_codes/validation/results/res_el_post_1hit_oldTuple.pdf");


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
n2.SaveAs("/home/espedica/macros_fairmu/clean_codes/validation/results/res_el_pre_2hit_oldTuple.pdf");


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

//n3.SaveAs("/home/espedica/macros_fairmu/clean_codes/validation/results/res_el_pre_2hit.pdf");





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
n4.SaveAs("/home/espedica/macros_fairmu/clean_codes/validation/results/res_el_post_2hit_oldTuple.pdf");









TCanvas n5("n5","n5",2800,2800);
n5.Divide(4,4);
for(int m=0; m<14; m++){
n5.cd(m+1);
if(m==0) h_pre_1hit_mu.at(m)->GetXaxis()->SetRangeUser(-0.001,0.001);
TF1 *g1 = new TF1("g1", "gaus");
h_pre_1hit_mu.at(m)->Fit("g1");
h_sigma_mu_pre_1hit_clones->SetPoint(m,edges_mu[m]*1000,g1->GetParameter(2)*1000);
h_sigma_mu_pre_1hit_clones->SetPointError(m,0.,g1->GetParError(2)*1000);
sigma_pre_mu1[m]=g1->GetParameter(2)*1000;
}
n5.SaveAs("/home/espedica/macros_fairmu/clean_codes/validation/results/res_mu_pre_1hit_oldTuple.pdf");


TCanvas n6("n6","n6",2800,2800);
n6.Divide(4,4);
for(int m=0; m<14; m++){
n6.cd(m+1);
TF1 *g1 = new TF1("g1", "gaus");
if(m<12) h_post_1hit_mu.at(m)->Fit("g1");
else  h_post_1hit_mu.at(m)->Fit("g1","R","",-0.0002,0.0002);
h_sigma_mu_post_1hit_clones->SetPoint(m,edges_mu[m]*1000,g1->GetParameter(2)*1000);
h_sigma_mu_post_1hit_clones->SetPointError(m,0.,g1->GetParError(2)*1000);
sigma_post_mu1[m]=g1->GetParameter(2)*1000;
}
n6.SaveAs("/home/espedica/macros_fairmu/clean_codes/validation/results/res_mu_post_1hit_oldTuple.pdf");



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
n7.SaveAs("/home/espedica/macros_fairmu/clean_codes/validation/results/res_mu_pre_2hit_oldTuple.pdf");


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

//n2_mu.SaveAs("/home/espedica/macros_fairmu/clean_codes/validation/results/res_pre_AngleMu_2hit.pdf");


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
n8.SaveAs("/home/espedica/macros_fairmu/clean_codes/validation/results/res_mu_post_2hit_oldTuple.pdf");





for(int m=0; m<14; m++){
sigma_el_diff_1hit->SetPoint(m,edges_el[m]*1000,(sigma_post_el1[m]-sigma_pre_el1[m])/sigma_pre_el1[m]);
if(m!=0)sigma_mu_diff_1hit->SetPoint(m,edges_mu[m]*1000,(sigma_post_mu1[m]-sigma_pre_mu1[m])/sigma_pre_mu1[m]);
sigma_el_diff_2hit->SetPoint(m,edges_el[m]*1000,(sigma_post_el2[m]-sigma_pre_el2[m])/sigma_pre_el2[m]);
if(m!=0)sigma_mu_diff_2hit->SetPoint(m,edges_mu[m]*1000,(sigma_post_mu2[m]-sigma_pre_mu2[m])/sigma_pre_mu2[m]);
}




h_sigma_mu_pre_2hit_clones->SetName("mu_pre");
h_sigma_el_pre_2hit_clones->SetName("el_pre");
h_sigma_mu_post_2hit_clones->SetName("mu_post");
h_sigma_el_post_2hit_clones->SetName("el_post");


h_sigma_el_pre_2hit_clones->SaveAs("/home/espedica/macros_fairmu/clean_codes/validation/results/sigma_el_pre_2hit_oldTuple.root");
h_sigma_el_post_2hit_clones->SaveAs("/home/espedica/macros_fairmu/clean_codes/validation/results/sigma_el_post_2hit_oldTuple.root");
h_sigma_mu_pre_2hit_clones->SaveAs("/home/espedica/macros_fairmu/clean_codes/validation/results/sigma_mu_pre_2hit_oldTuple.root");
h_sigma_mu_post_2hit_clones->SaveAs("/home/espedica/macros_fairmu/clean_codes/validation/results/sigma_mu_post_2hit_oldTuple.root");


 auto legend_mu2_pt2 = new TLegend(0.75,0.15,0.9,0.3);
legend_mu2_pt2->AddEntry(h_sigma_mu_pre_1hit_clones,"1 hit shared","LEP");
legend_mu2_pt2->AddEntry(h_sigma_mu_pre_2hit_clones,"2 hit shared","LEP");



 auto legend_el2_pt2 = new TLegend(0.75,0.15,0.9,0.3);
legend_el2_pt2->AddEntry(h_sigma_el_pre_1hit_clones,"1 hit shared","LEP");
legend_el2_pt2->AddEntry(h_sigma_el_pre_2hit_clones,"2 hit shared","LEP");




 auto legend_mu1_lin = new TLegend(0.75,0.15,0.9,0.3);
legend_mu1_lin->AddEntry(h_sigma_mu_pre_1hit_clones,"pre-vrtx","LEP");
legend_mu1_lin->AddEntry(h_sigma_mu_post_1hit_clones,"post-vrtx","LEP");

 auto legend_mu2_lin = new TLegend(0.75,0.15,0.9,0.3);
legend_mu2_lin->AddEntry(h_sigma_mu_pre_2hit_clones,"pre-vrtx","LEP");
legend_mu2_lin->AddEntry(h_sigma_mu_post_2hit_clones,"post-vrtx","LEP");





 auto legend_mu1 = new TLegend(0.75,0.15,0.9,0.3);
legend_mu1->AddEntry(h_sigma_mu_pre_1hit_clones,"pre-vrtx","LEP");
legend_mu1->AddEntry(h_sigma_mu_post_1hit_clones,"post-vrtx","LEP");

 auto legend_mu2 = new TLegend(0.75,0.15,0.9,0.3);
legend_mu2->AddEntry(h_sigma_mu_pre_2hit_clones,"pre-vrtx","LEP");
legend_mu2->AddEntry(h_sigma_mu_post_2hit_clones,"post-vrtx","LEP");



 auto legend_e1 = new TLegend(0.75,0.15,0.9,0.3);
legend_e1->AddEntry(h_sigma_mu_pre_1hit_clones,"pre-vrtx","LEP");
legend_e1->AddEntry(h_sigma_mu_post_1hit_clones,"post-vrtx","LEP");

 auto legend_e2 = new TLegend(0.75,0.15,0.9,0.3);
legend_e2->AddEntry(h_sigma_mu_pre_2hit_clones,"pre-vrtx","LEP");
legend_e2->AddEntry(h_sigma_mu_post_2hit_clones,"post-vrtx","LEP");




TCanvas s2_pt2("s2_pt2","s2_pt2",1400,1400);
s2_pt2.Divide(1,2);
s2_pt2.cd(1);

TMultiGraph *mu_1hit_pt2 = new TMultiGraph();

h_sigma_mu_pre_1hit_clones->SetMarkerColor(kGreen+1);
h_sigma_mu_pre_2hit_clones->SetMarkerColor(kAzure+7);

mu_1hit_pt2->Add(h_sigma_mu_pre_1hit_clones,"AP*");
mu_1hit_pt2->Add(h_sigma_mu_pre_2hit_clones,"AP*");

mu_1hit_pt2->Draw("AP*");
legend_mu2_pt2->Draw();
mu_1hit_pt2->SetMinimum(0.03);
mu_1hit_pt2->SetMaximum(0.2);
mu_1hit_pt2->SetTitle("Muon angular resolution");
mu_1hit_pt2->GetXaxis()->SetTitle("Scattering angle (mrad)");
mu_1hit_pt2->GetYaxis()->SetTitle("#sigma(#theta) (mrad)");
mu_1hit_pt2->GetHistogram()->GetXaxis()->SetLimits(0.,5.);
gPad->SetLogy();

s2_pt2.cd(2);

TMultiGraph *el_1hit_pt2 = new TMultiGraph();

h_sigma_el_pre_1hit_clones->SetMarkerColor(kGreen+1);
h_sigma_el_pre_2hit_clones->SetMarkerColor(kAzure+7);


el_1hit_pt2->Add(h_sigma_el_pre_1hit_clones,"AP*");
el_1hit_pt2->Add(h_sigma_el_pre_2hit_clones,"AP*");

el_1hit_pt2->Draw("AP*");
legend_el2_pt2->Draw();
el_1hit_pt2->SetMinimum(0.008);
el_1hit_pt2->SetMaximum(5);
el_1hit_pt2->SetTitle("Electron angular resolution");
el_1hit_pt2->GetXaxis()->SetTitle("Scattering angle (mrad)");
el_1hit_pt2->GetYaxis()->SetTitle("#sigma(#theta) (mrad)");
el_1hit_pt2->GetXaxis()->SetLimits(0.,32.);
gPad->SetLogy();

s2_pt2.SaveAs("/home/espedica/macros_fairmu/clean_codes/validation/results/sigma_oldTuple.pdf");





TCanvas s2("s2","s2",2100,1400);
s2.Divide(2,2);

s2.cd(1);

TMultiGraph *mu_1hit = new TMultiGraph();


h_sigma_mu_pre_1hit_clones->SetMarkerColor(kGreen+1);
h_sigma_mu_post_1hit_clones->SetMarkerColor(kBlack);

mu_1hit->Add(h_sigma_mu_pre_1hit_clones,"AP*");
mu_1hit->Add(h_sigma_mu_post_1hit_clones,"AP*");

mu_1hit->Draw("AP*");
legend_mu1_lin->Draw();
mu_1hit->SetMinimum(0.008);
mu_1hit->SetMaximum(0.2);
mu_1hit->SetTitle("Muon angular resolution 1 hit shared");
mu_1hit->GetXaxis()->SetTitle("Scattering angle (mrad)");
mu_1hit->GetYaxis()->SetTitle("#sigma(#theta) (mrad)");
mu_1hit->GetHistogram()->GetXaxis()->SetLimits(0.,5.);
gPad->SetLogy();

s2.cd(2);

TMultiGraph *mu_2hit = new TMultiGraph();


h_sigma_mu_pre_2hit_clones->SetMarkerColor(kAzure+7);
h_sigma_mu_post_2hit_clones->SetMarkerColor(kBlack);

mu_2hit->Add(h_sigma_mu_pre_2hit_clones,"AP*");
mu_2hit->Add(h_sigma_mu_post_2hit_clones,"AP*");

mu_2hit->Draw("AP*");
legend_mu2_lin->Draw();
mu_2hit->SetMinimum(0.008);
mu_2hit->SetMaximum(0.2);
mu_2hit->SetTitle("Muon angular resolution 2 hit shared");
mu_2hit->GetXaxis()->SetTitle("Scattering angle (mrad)");
mu_2hit->GetYaxis()->SetTitle("#sigma(#theta) (mrad)");
mu_2hit->GetHistogram()->GetXaxis()->SetLimits(0.,5.);
gPad->SetLogy();

s2.cd(3);

TMultiGraph *el_1hit = new TMultiGraph();

h_sigma_el_pre_1hit_clones->SetMarkerColor(kGreen+1);
h_sigma_el_post_1hit_clones->SetMarkerColor(kBlack);


el_1hit->Add(h_sigma_el_pre_1hit_clones,"AP*");
el_1hit->Add(h_sigma_el_post_1hit_clones,"AP*");

el_1hit->Draw("AP*");
legend_e1->Draw();
el_1hit->SetMinimum(0.008);
el_1hit->SetMaximum(5.);
el_1hit->SetTitle("Electron angular resolution 1 hit shared");
el_1hit->GetXaxis()->SetTitle("Scattering angle (mrad)");
el_1hit->GetYaxis()->SetTitle("#sigma(#theta) (mrad)");
el_1hit->GetXaxis()->SetLimits(0.,32.);
gPad->SetLogy();

s2.cd(4);

TMultiGraph *el_2hit = new TMultiGraph();

h_sigma_el_pre_2hit_clones->SetMarkerColor(kAzure+7);
h_sigma_el_post_2hit_clones->SetMarkerColor(kBlack);


el_2hit->Add(h_sigma_el_pre_2hit_clones,"AP*");
el_2hit->Add(h_sigma_el_post_2hit_clones,"AP*");

el_2hit->Draw("AP*");
legend_e2->Draw();
el_2hit->SetMinimum(0.008);
el_2hit->SetMaximum(5.);
el_2hit->SetTitle("Electron angular resolution 2 hit shared");
el_2hit->GetXaxis()->SetTitle("Scattering angle (mrad)");
el_2hit->GetYaxis()->SetTitle("#sigma(#theta) (mrad)");
el_2hit->GetXaxis()->SetLimits(0.,32.);
gPad->SetLogy();
s2.SaveAs("/home/espedica/macros_fairmu/clean_codes/validation/results/sigma_prepost_oldTuple.pdf");


TCanvas d("d","d",2100,1400);
d.Divide(2,2);
d.cd(1);
sigma_el_diff_1hit->SetName("el_1");
sigma_el_diff_1hit->GetYaxis()->SetTitle("(#sigma_post - #sigma_pre)/#sigma_pre");
sigma_el_diff_1hit->GetXaxis()->SetTitle("#theta_el (mrad)");
sigma_el_diff_1hit->Draw("AP*");
d.cd(2);
sigma_mu_diff_1hit->SetName("mu_1");
sigma_mu_diff_1hit->GetYaxis()->SetTitle("(#sigma_post - #sigma_pre)/#sigma_pre");
sigma_mu_diff_1hit->GetXaxis()->SetTitle("#theta_#mu (mrad)");
sigma_mu_diff_1hit->Draw("AP*");
d.cd(3);
sigma_el_diff_2hit->SetName("el_2");
sigma_el_diff_2hit->GetYaxis()->SetTitle("(#sigma_post - #sigma_pre)/#sigma_pre");
sigma_el_diff_2hit->GetXaxis()->SetTitle("#theta_el (mrad)");
sigma_el_diff_2hit->Draw("AP*");
d.cd(4);
sigma_mu_diff_2hit->SetName("mu_2");
sigma_mu_diff_2hit->GetYaxis()->SetTitle("(#sigma_post - #sigma_pre)/#sigma_pre");
sigma_mu_diff_2hit->GetXaxis()->SetTitle("#theta_#mu (mrad)");
sigma_mu_diff_2hit->Draw("AP*");

d.SaveAs("/home/espedica/macros_fairmu/clean_codes/validation/results/difference_oldTuple.pdf");

sigma_el_diff_1hit->SaveAs("/home/espedica/macros_fairmu/clean_codes/validation/results/el_1_oldTuple.root");
sigma_mu_diff_1hit->SaveAs("/home/espedica/macros_fairmu/clean_codes/validation/results/mu_1_oldTuple.root");
sigma_el_diff_2hit->SaveAs("/home/espedica/macros_fairmu/clean_codes/validation/results/el_2_oldTuple.root");
sigma_mu_diff_2hit->SaveAs("/home/espedica/macros_fairmu/clean_codes/validation/results/mu_2_oldTuple.root");

}


