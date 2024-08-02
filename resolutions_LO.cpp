
#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TStyle.h>


void resolutions_LO(){



array<TFile*,14> f_pre_1hitFirstModules_el, f_post_1hitFirstModules_el, f_pre_2hitFirstModules_el, f_post_2hitFirstModules_el;
array<TFile*,14> f_pre_1hitFirstModules_mu, f_post_1hitFirstModules_mu, f_pre_2hitFirstModules_mu, f_post_2hitFirstModules_mu;

for(int m=0; m<14; m++){

f_pre_1hitFirstModules_el.at(m)=TFile::Open(Form("quality_tracks/root/res_prevrtx_clones_el_%d_1hitFirstModules_LO_reassign.root",static_cast<char>(m)));
f_post_1hitFirstModules_el.at(m)=TFile::Open(Form("quality_tracks/root/res_vrtx_clones_el_%d_1hitFirstModules_LO_reassign.root",static_cast<char>(m)));
f_pre_2hitFirstModules_el.at(m)=TFile::Open(Form("quality_tracks/root/res_prevrtx_clones_el_%d_0hit_LO_reassign.root",static_cast<char>(m)));
f_post_2hitFirstModules_el.at(m)=TFile::Open(Form("quality_tracks/root/res_vrtx_clones_el_%d_0hit_LO_reassign.root",static_cast<char>(m)));
}

for(int m=0; m<14; m++){
f_pre_1hitFirstModules_mu.at(m)=TFile::Open(Form("quality_tracks/root/res_prevrtx_clones_AngleMu_%d_1hitFirstModules_LO_reassign.root",static_cast<char>(m)));
f_post_1hitFirstModules_mu.at(m)=TFile::Open(Form("quality_tracks/root/res_vrtx_clones_AngleMu_%d_1hitFirstModules_LO_reassign.root",static_cast<char>(m)));
f_pre_2hitFirstModules_mu.at(m)=TFile::Open(Form("quality_tracks/root/res_prevrtx_clones_AngleMu_%d_0hit_LO_reassign.root",static_cast<char>(m)));
f_post_2hitFirstModules_mu.at(m)=TFile::Open(Form("quality_tracks/root/res_vrtx_clones_AngleMu_%d_0hit_LO_reassign.root",static_cast<char>(m)));
}


TH1::SetDefaultSumw2(kTRUE);

array <TH1D*,14> h_pre_1hitFirstModules_el, h_post_1hitFirstModules_el, h_pre_2hitFirstModules_el, h_post_2hitFirstModules_el;
array <TH1D*,14> h_pre_1hitFirstModules_mu, h_post_1hitFirstModules_mu, h_pre_2hitFirstModules_mu, h_post_2hitFirstModules_mu;


for(int m=0; m<14; m++){

h_pre_1hitFirstModules_el.at(m)=(TH1D*)f_pre_1hitFirstModules_el.at(m)->Get(Form("res%d_el_pre_clones",static_cast<char>(m)));
h_post_1hitFirstModules_el.at(m)=(TH1D*)f_post_1hitFirstModules_el.at(m)->Get(Form("res%d_vrtx_el_clones",static_cast<char>(m)));
h_pre_2hitFirstModules_el.at(m)=(TH1D*)f_pre_2hitFirstModules_el.at(m)->Get(Form("res%d_el_pre_clones",static_cast<char>(m)));
h_post_2hitFirstModules_el.at(m)=(TH1D*)f_post_2hitFirstModules_el.at(m)->Get(Form("res%d_vrtx_el_clones",static_cast<char>(m)));
}

for(int m=0; m<14; m++){
h_pre_1hitFirstModules_mu.at(m)=(TH1D*)f_pre_1hitFirstModules_mu.at(m)->Get(Form("res%d_mu_pre_clones",static_cast<char>(m)));
h_post_1hitFirstModules_mu.at(m)=(TH1D*)f_post_1hitFirstModules_mu.at(m)->Get(Form("res%d_vrtx_mu_clones",static_cast<char>(m)));
h_pre_2hitFirstModules_mu.at(m)=(TH1D*)f_pre_2hitFirstModules_mu.at(m)->Get(Form("res%d_mu_pre_clones",static_cast<char>(m)));
h_post_2hitFirstModules_mu.at(m)=(TH1D*)f_post_2hitFirstModules_mu.at(m)->Get(Form("res%d_vrtx_mu_clones",static_cast<char>(m)));

}


/*
 auto legend_op = new TLegend(0.6,0.1,0.9,0.25);
   legend_op->AddEntry(h_op_clone0,"LO 0 hit","LEP");
   legend_op->AddEntry(h_op_clone1,"LO 1 hit","LEP");
   legend_op->AddEntry(h_op_clone2,"LO 2 hit","LEP");*/


//array <TF1*,14> fit_pre_1hitFirstModules_el, fit_post_1hitFirstModules_el, fit_pre_2hitFirstModules_el, fit_post_2hitFirstModules_el, fit_pre_1hitFirstModules_mu, fit_post_1hitFirstModules_mu, fit_pre_2hitFirstModules_mu, fit_post_2hitFirstModules_mu;
//TF1 *a1 = new TF1("a1", "[2]*TMath::Gaus(x,[0],[1])");
//a1->SetParameters(5.3e-05,110e-05,1);


   const Int_t NBINS_mu = 14;
   Double_t edges_mu[NBINS_mu + 1] = {0.00005,0.00015,0.00025,0.00035,0.00045,0.00055,0.00065,0.00075,0.00085,0.00095,0.0015,0.0025,0.0035,0.0045};

   const Int_t NBINS = 14;
   Double_t edges_el[NBINS + 1] = {0.0005,0.0015,0.0025,0.0035,0.0045,0.0055,0.0065,0.0075,0.0085,0.0095, 0.0125, 0.0175, 0.0225, 0.0285};

/*
TH2D* h_sigma_mu_pre_1hitFirstModules_clones=new TH2D("h_sigma_mu_pre_1hit","sigma VS muon scattering angle 1hitFirstModules",NBINS_mu,edges_mu,80,-0.0001,0.0001);
TH2D* h_sigma_el_pre_1hitFirstModules_clones=new TH2D("h_sigma_el_pre_1hit","sigma VS electron scattering angle 1hitFirstModules",NBINS,edges_el,800,-0.0008,0.0008);
TH2D* h_sigma_mu_post_1hitFirstModules_clones=new TH2D("h_sigma_mu_post_1hit","sigma VS muon scattering angle 1hitFirstModules",NBINS_mu,edges_mu,80,-0.0001,0.0001);
TH2D* h_sigma_el_post_1hitFirstModules_clones=new TH2D("h_sigma_el_post_1hit","sigma VS electrons scattering angle 1hitFirstModules",NBINS,edges_el,800,-0.0008,0.0008);

TH2D* h_sigma_mu_pre_2hitFirstModules_clones=new TH2D("h_sigma_mu_pre_2hit","sigma VS muon scattering angle 2hitFirstModules",NBINS_mu,edges_mu,80,-0.0001,0.0001);
TH2D* h_sigma_el_pre_2hitFirstModules_clones=new TH2D("h_sigma_el_pre_2hit","sigma VS electron scattering angle 2hitFirstModules",NBINS,edges_el,800,-0.0008,0.0008);
TH2D* h_sigma_mu_post_2hitFirstModules_clones=new TH2D("h_sigma_mu_post_2hit","sigma VS muon scattering angle 2hitFirstModules",NBINS_mu,edges_mu,80,-0.0001,0.0001);
TH2D* h_sigma_el_post_2hitFirstModules_clones=new TH2D("h_sigma_el_post_2hit","sigma VS electrons scattering angle 2hitFirstModules",NBINS,edges_el,800,-0.0008,0.0008);*/


auto h_sigma_mu_pre_1hitFirstModules_clones = new TGraph();
auto h_sigma_el_pre_1hitFirstModules_clones = new TGraph();
auto h_sigma_mu_post_1hitFirstModules_clones = new TGraph();
auto h_sigma_el_post_1hitFirstModules_clones = new TGraph();
auto h_sigma_mu_pre_2hitFirstModules_clones = new TGraph();
auto h_sigma_el_pre_2hitFirstModules_clones = new TGraph();
auto h_sigma_mu_post_2hitFirstModules_clones = new TGraph();
auto h_sigma_el_post_2hitFirstModules_clones = new TGraph();



auto h_RMS_mu_pre_1hitFirstModules_clones = new TGraph();
auto h_RMS_el_pre_1hitFirstModules_clones = new TGraph();
auto h_RMS_mu_post_1hitFirstModules_clones = new TGraph();
auto h_RMS_el_post_1hitFirstModules_clones = new TGraph();
auto h_RMS_mu_pre_2hitFirstModules_clones = new TGraph();
auto h_RMS_el_pre_2hitFirstModules_clones = new TGraph();
auto h_RMS_mu_post_2hitFirstModules_clones = new TGraph();
auto h_RMS_el_post_2hitFirstModules_clones = new TGraph();


TCanvas n0("n0","n0",2100,2100);
n0.Divide(3,5);
for(int m=0; m<14; m++){
n0.cd(m+1);
/*if(m<=2) h_pre_1hitFirstModules_el.at(m)->GetXaxis()->SetRangeUser(-0.001,0.001);
else if(m>2 and m<=4) h_pre_1hitFirstModules_el.at(m)->GetXaxis()->SetRangeUser(-0.002,0.002);
else  h_pre_1hitFirstModules_el.at(m)->GetXaxis()->SetRangeUser(-0.02,0.02);
*/TGaxis::SetMaxDigits(3);
h_pre_1hitFirstModules_el.at(m)->SetLineColor(kGreen+1);
h_pre_1hitFirstModules_el.at(m)->GetXaxis()->SetTitle("rad");
h_pre_1hitFirstModules_el.at(m)->SetMinimum(1.);
h_pre_1hitFirstModules_el.at(m)->Rebin(2);
TF1 *g1 = new TF1("g1", "gaus");
if(m<4) h_pre_1hitFirstModules_el.at(m)->Fit("g1","R","",-0.0005,0.0005);
else if(m==4) h_pre_1hitFirstModules_el.at(m)->Fit("g1","R","",-0.002,0.002);
else h_pre_1hitFirstModules_el.at(m)->Fit("g1");
h_pre_1hitFirstModules_el.at(m)->Draw("hist same");
h_sigma_el_pre_1hitFirstModules_clones->SetPoint(h_sigma_el_pre_1hitFirstModules_clones->GetN(),edges_el[m]*1000,g1->GetParameter(2)*1000);
h_RMS_el_pre_1hitFirstModules_clones->SetPoint(h_RMS_el_pre_1hitFirstModules_clones->GetN(),edges_el[m]*1000,h_pre_1hitFirstModules_el.at(m)->GetRMS()*1000);
gStyle->SetStatH(0.1);
gStyle->SetOptFit(0001);
gStyle->SetOptStat(1101);
}
n0.SaveAs("proposal/res_pre_el_1hitFirstModules_clones_LO_reassign.pdf");

TCanvas n1("n1","n1",2100,2100);
n1.Divide(3,5);
for(int m=0; m<14; m++){
n1.cd(m+1);
/*if(m<=2) h_post_1hitFirstModules_el.at(m)->GetXaxis()->SetRangeUser(-0.001,0.001);
else if(m>2 and m<=4) h_post_1hitFirstModules_el.at(m)->GetXaxis()->SetRangeUser(-0.002,0.002);
else  h_post_1hitFirstModules_el.at(m)->GetXaxis()->SetRangeUser(-0.02,0.02);
*/TGaxis::SetMaxDigits(3);
h_post_1hitFirstModules_el.at(m)->GetXaxis()->SetTitle("rad");
h_post_1hitFirstModules_el.at(m)->SetMinimum(1.);
h_post_1hitFirstModules_el.at(m)->SetLineColor(kGreen+1);
h_post_1hitFirstModules_el.at(m)->Rebin(2);
TF1 *g1 = new TF1("g1", "gaus");
if(m<4) h_post_1hitFirstModules_el.at(m)->Fit("g1","R","",-0.0005,0.0005);
else if(m==4) h_post_1hitFirstModules_el.at(m)->Fit("g1","R","",-0.002,0.002);
else h_post_1hitFirstModules_el.at(m)->Fit("g1");
h_sigma_el_post_1hitFirstModules_clones->SetPoint(h_sigma_el_post_1hitFirstModules_clones->GetN(),edges_el[m]*1000,g1->GetParameter(2)*1000);
h_RMS_el_post_1hitFirstModules_clones->SetPoint(h_RMS_el_post_1hitFirstModules_clones->GetN(),edges_el[m]*1000,h_post_1hitFirstModules_el.at(m)->GetRMS()*1000);

h_post_1hitFirstModules_el.at(m)->Draw("hist same");
gStyle->SetStatH(0.1);
gStyle->SetOptFit(0001);
gStyle->SetOptStat(1101);
}
n1.SaveAs("proposal/res_post_el_1hitFirstModules_clones_LO_reassign.pdf");

TCanvas n2("n2","n2",2100,2100);
n2.Divide(3,5);
for(int m=0; m<14; m++){
n2.cd(m+1);
/*if(m<=2) h_pre_2hitFirstModules_el.at(m)->GetXaxis()->SetRangeUser(-0.001,0.001);
else if(m>2 and m<=4) h_pre_2hitFirstModules_el.at(m)->GetXaxis()->SetRangeUser(-0.002,0.002);
else  h_pre_2hitFirstModules_el.at(m)->GetXaxis()->SetRangeUser(-0.02,0.02);
*/TGaxis::SetMaxDigits(3);
h_pre_2hitFirstModules_el.at(m)->SetLineColor(kAzure+7);
h_pre_2hitFirstModules_el.at(m)->GetXaxis()->SetTitle("rad");
h_pre_2hitFirstModules_el.at(m)->SetMinimum(1.);
h_post_2hitFirstModules_el.at(m)->Rebin(2);
TF1 *g1 = new TF1("g1", "gaus");
if(m<4) h_pre_2hitFirstModules_el.at(m)->Fit("g1","R","",-0.0005,0.0005);
else if(m==4) h_pre_2hitFirstModules_el.at(m)->Fit("g1","R","",-0.002,0.002);
else h_pre_2hitFirstModules_el.at(m)->Fit("g1");
h_pre_2hitFirstModules_el.at(m)->Draw("hist same");
h_sigma_el_pre_2hitFirstModules_clones->SetPoint(h_sigma_el_pre_2hitFirstModules_clones->GetN(),edges_el[m]*1000,g1->GetParameter(2)*1000);
h_RMS_el_pre_2hitFirstModules_clones->SetPoint(h_RMS_el_pre_2hitFirstModules_clones->GetN(),edges_el[m]*1000,h_pre_2hitFirstModules_el.at(m)->GetRMS()*1000);

gStyle->SetStatH(0.1);
gStyle->SetOptFit(0001);
gStyle->SetOptStat(1101);
}
n2.SaveAs("proposal/res_pre_el_2hitFirstModules_clones_LO_reassign.pdf");

TCanvas n3("n3","n3",2100,2100);
n3.Divide(3,5);
for(int m=0; m<14; m++){
n3.cd(m+1);
/*if(m<=2) h_post_2hitFirstModules_el.at(m)->GetXaxis()->SetRangeUser(-0.001,0.001);
else if(m>2 and m<=4) h_post_2hitFirstModules_el.at(m)->GetXaxis()->SetRangeUser(-0.002,0.002);
else  h_post_2hitFirstModules_el.at(m)->GetXaxis()->SetRangeUser(-0.02,0.02);
*/TGaxis::SetMaxDigits(3);
h_post_2hitFirstModules_el.at(m)->GetXaxis()->SetTitle("rad");
h_post_2hitFirstModules_el.at(m)->SetMinimum(1.);
h_post_2hitFirstModules_el.at(m)->SetLineColor(kAzure+7);
h_post_2hitFirstModules_el.at(m)->Rebin(2);
TF1 *g1 = new TF1("g1", "gaus");
if(m<4) h_post_2hitFirstModules_el.at(m)->Fit("g1","R","",-0.0005,0.0005);
else if(m==4) h_post_2hitFirstModules_el.at(m)->Fit("g1","R","",-0.002,0.002);
else h_post_2hitFirstModules_el.at(m)->Fit("g1");
h_sigma_el_post_2hitFirstModules_clones->SetPoint(h_sigma_el_post_2hitFirstModules_clones->GetN(),edges_el[m]*1000,g1->GetParameter(2)*1000);
h_RMS_el_post_2hitFirstModules_clones->SetPoint(h_RMS_el_post_2hitFirstModules_clones->GetN(),edges_el[m]*1000,h_post_2hitFirstModules_el.at(m)->GetRMS()*1000);

h_post_2hitFirstModules_el.at(m)->Draw("hist same");
gStyle->SetStatH(0.1);
gStyle->SetOptFit(0001);
gStyle->SetOptStat(1101);
}
n3.SaveAs("proposal/res_post_el_2hitFirstModules_clones_LO_reassign.pdf");



TCanvas n0_mu("n0_mu","n0_mu",2100,2100);
n0_mu.Divide(3,5);
for(int m=0; m<14; m++){
n0_mu.cd(m+1);
if(m==0) h_pre_1hitFirstModules_mu.at(m)->GetXaxis()->SetRangeUser(-0.001,0.001);
h_pre_1hitFirstModules_mu.at(m)->SetLineColor(kGreen+1);
TGaxis::SetMaxDigits(3);
h_pre_1hitFirstModules_mu.at(m)->GetXaxis()->SetTitle("rad");
h_pre_1hitFirstModules_mu.at(m)->SetMinimum(1.);
h_pre_1hitFirstModules_mu.at(m)->Rebin(2);
TF1 *g1 = new TF1("g1", "gaus");
h_pre_1hitFirstModules_mu.at(m)->Fit("g1");
h_pre_1hitFirstModules_mu.at(m)->Draw("hist same");
//h_sigma_mu_pre_1hitFirstModules_clones->Fill(edges_mu[m]*1000,g1->GetParameter(2)*1000);
h_sigma_mu_pre_1hitFirstModules_clones->SetPoint(h_sigma_mu_pre_1hitFirstModules_clones->GetN(),edges_mu[m]*1000,g1->GetParameter(2)*1000);
h_RMS_mu_pre_1hitFirstModules_clones->SetPoint(h_RMS_mu_pre_1hitFirstModules_clones->GetN(),edges_mu[m]*1000,h_pre_1hitFirstModules_mu.at(m)->GetRMS()*1000);

gStyle->SetStatH(0.1);
gStyle->SetOptFit(0001);
gStyle->SetOptStat(1101);
}
n0_mu.SaveAs("proposal/res_pre_AngleMu_1hitFirstModules_clones_LO_reassign.pdf");

TCanvas n1_mu("n1_mu","n1_mu",2100,2100);
n1_mu.Divide(3,5);
for(int m=0; m<14; m++){
n1_mu.cd(m+1);
TGaxis::SetMaxDigits(3);
h_post_1hitFirstModules_mu.at(m)->GetXaxis()->SetTitle("rad");
h_post_1hitFirstModules_mu.at(m)->SetMinimum(1.);
h_post_1hitFirstModules_mu.at(m)->SetLineColor(kGreen+1);
h_post_1hitFirstModules_mu.at(m)->Rebin(2);
TF1 *g1 = new TF1("g1", "gaus");
h_post_1hitFirstModules_mu.at(m)->Fit("g1");
//h_sigma_mu_post_1hitFirstModules_clones->Fill(edges_mu[m]*1000,g1->GetParameter(2)*1000);
h_sigma_mu_post_1hitFirstModules_clones->SetPoint(h_sigma_mu_post_1hitFirstModules_clones->GetN(),edges_mu[m]*1000,g1->GetParameter(2)*1000);
h_RMS_mu_post_1hitFirstModules_clones->SetPoint(h_RMS_mu_post_1hitFirstModules_clones->GetN(),edges_mu[m]*1000,h_post_1hitFirstModules_mu.at(m)->GetRMS()*1000);

h_post_1hitFirstModules_mu.at(m)->Draw("hist same");
gStyle->SetStatH(0.1);
gStyle->SetOptFit(0001);
gStyle->SetOptStat(1101);
}
n1_mu.SaveAs("proposal/res_post_AngleMu_1hitFirstModules_clones_LO_reassign.pdf");

TCanvas n2_mu("n2_mu","n2_mu",2100,2100);
n2_mu.Divide(3,5);
for(int m=0; m<14; m++){
n2_mu.cd(m+1);
//if(m==0) h_pre_2hitFirstModules_mu.at(m)->GetXaxis()->SetRangeUser(-0.001,0.001);
h_pre_2hitFirstModules_mu.at(m)->SetLineColor(kAzure+7);
TGaxis::SetMaxDigits(3);
h_pre_2hitFirstModules_mu.at(m)->GetXaxis()->SetTitle("rad");
h_pre_2hitFirstModules_mu.at(m)->SetMinimum(1.);
h_pre_2hitFirstModules_mu.at(m)->Rebin(2);
TF1 *g1 = new TF1("g1", "gaus");
h_pre_2hitFirstModules_mu.at(m)->Fit("g1");
//h_sigma_mu_pre_2hitFirstModules_clones->Fill(edges_mu[m]*1000,g1->GetParameter(2)*1000);
h_sigma_mu_pre_2hitFirstModules_clones->SetPoint(h_sigma_mu_pre_2hitFirstModules_clones->GetN(),edges_mu[m]*1000,g1->GetParameter(2)*1000);
h_RMS_mu_pre_2hitFirstModules_clones->SetPoint(h_RMS_mu_pre_2hitFirstModules_clones->GetN(),edges_mu[m]*1000,h_pre_2hitFirstModules_mu.at(m)->GetRMS()*1000);

h_pre_2hitFirstModules_mu.at(m)->Draw("hist same");
gStyle->SetStatH(0.1);
gStyle->SetOptFit(0001);
gStyle->SetOptStat(1101);
}
n2_mu.SaveAs("proposal/res_pre_AngleMu_2hitFirstModules_clones_LO_reassign.pdf");

TCanvas n3_mu("n3_mu","n3_mu",2100,2100);
n3_mu.Divide(3,5);
for(int m=0; m<14; m++){
n3_mu.cd(m+1);
TGaxis::SetMaxDigits(3);
h_post_2hitFirstModules_mu.at(m)->GetXaxis()->SetTitle("rad");
h_post_2hitFirstModules_mu.at(m)->SetMinimum(1.);
h_post_2hitFirstModules_mu.at(m)->SetLineColor(kAzure+7);
h_post_2hitFirstModules_mu.at(m)->Rebin(2);
TF1 *g1 = new TF1("g1", "gaus");
if(m>6)h_post_2hitFirstModules_mu.at(m)->Fit("g1","R","",-0.0002,0.0002);
else h_post_2hitFirstModules_mu.at(m)->Fit("g1");
//h_sigma_mu_post_2hitFirstModules_clones->Fill(edges_mu[m]*1000,g1->GetParameter(2)*1000);
h_sigma_mu_post_2hitFirstModules_clones->SetPoint(h_sigma_mu_post_2hitFirstModules_clones->GetN(),edges_mu[m]*1000,g1->GetParameter(2)*1000);
h_RMS_mu_post_2hitFirstModules_clones->SetPoint(h_RMS_mu_post_2hitFirstModules_clones->GetN(),edges_mu[m]*1000,h_post_2hitFirstModules_mu.at(m)->GetRMS()*1000);

h_post_2hitFirstModules_mu.at(m)->Draw("hist same");
gStyle->SetStatH(0.1);
gStyle->SetOptFit(0001);
gStyle->SetOptStat(1101);
}
n3_mu.SaveAs("proposal/res_post_AngleMu_2hitFirstModules_clones_LO_reassign.pdf");

 auto legend_mu1 = new TLegend(0.75,0.15,0.9,0.3);
legend_mu1->AddEntry(h_sigma_mu_pre_1hitFirstModules_clones,"1 hit shared","LEP");
legend_mu1->AddEntry(h_sigma_mu_post_1hitFirstModules_clones,"Post_vrtx","LEP");

 auto legend_mu2 = new TLegend(0.75,0.15,0.9,0.3);
legend_mu2->AddEntry(h_sigma_mu_pre_2hitFirstModules_clones,"Pre-vrtx","LEP");
legend_mu2->AddEntry(h_sigma_mu_post_2hitFirstModules_clones,"Post_vrtx","LEP");


 auto legend_mu_r1 = new TLegend(0.75,0.15,0.9,0.3);
legend_mu_r1->AddEntry(h_RMS_mu_pre_1hitFirstModules_clones,"Pre-vrtx","LEP");
legend_mu_r1->AddEntry(h_RMS_mu_post_1hitFirstModules_clones,"Post_vrtx","LEP");


 auto legend_mu_r2 = new TLegend(0.75,0.15,0.9,0.3);
legend_mu_r2->AddEntry(h_RMS_mu_pre_2hitFirstModules_clones,"Pre-vrtx","LEP");
legend_mu_r2->AddEntry(h_RMS_mu_post_2hitFirstModules_clones,"Post_vrtx","LEP");

TCanvas s1("s1","s1",1600,1600);
s1.Divide(2,2);
s1.cd(1);

TMultiGraph *mu_1hit = new TMultiGraph();
mu_1hit->SetMinimum(0.);
mu_1hit->SetMaximum(3.);


mu_1hit->Add(h_sigma_mu_pre_1hitFirstModules_clones,"A*");
mu_1hit->Add(h_sigma_mu_post_1hitFirstModules_clones,"A*");
mu_1hit->Draw("A*");
legend_mu1->Draw();
gPad->SetLogy();
mu_1hit->SetTitle("Sigma VS muon scattering angle LO 1HitShared");

/*h_sigma_mu_pre_1hitFirstModules_clones->SetTitle("Sigma VS muon scattering angle LO pre-vrtx 1HitShared");
h_sigma_mu_pre_1hitFirstModules_clones->GetXaxis()->SetTitle("Muon scattering angle (mrad)");
h_sigma_mu_pre_1hitFirstModules_clones->GetYaxis()->SetTitle("Sigma (rad)");*/

s1.cd(2);

TMultiGraph *mu_2hit = new TMultiGraph();
mu_2hit->SetMinimum(0.);
mu_2hit->SetMaximum(3.);


mu_2hit->Add(h_sigma_mu_pre_2hitFirstModules_clones,"A*");
mu_2hit->Add(h_sigma_mu_post_2hitFirstModules_clones,"A*");
mu_2hit->Draw("A*");
legend_mu2->Draw();
gPad->SetLogy();
mu_2hit->SetTitle("Sigma VS muon scattering angle LO 2HitShared");

s1.cd(3);
TMultiGraph *mu_1hit_r = new TMultiGraph();
mu_1hit_r->SetMinimum(0.);
mu_1hit_r->SetMaximum(3.);


mu_1hit_r->Add(h_RMS_mu_pre_1hitFirstModules_clones,"A*");
mu_1hit_r->Add(h_RMS_mu_post_1hitFirstModules_clones,"A*");
mu_1hit_r->Draw("A*");
legend_mu_r1->Draw();
gPad->SetLogy(); 
mu_1hit_r->SetTitle("RMS VS muon scattering angle LO 1HitShared");

s1.cd(4);

TMultiGraph *mu_2hit_r = new TMultiGraph();
mu_2hit_r->SetMinimum(0.);
mu_2hit_r->SetMaximum(3.);


mu_2hit_r->Add(h_RMS_mu_pre_2hitFirstModules_clones,"A*");
mu_2hit_r->Add(h_RMS_mu_post_2hitFirstModules_clones,"A*");
mu_2hit_r->Draw("A*");
legend_mu_r2->Draw();
gPad->SetLogy();
mu_2hit_r->SetTitle("RMS VS muon scattering angle LO 1HitShared");


s1.SaveAs("proposal/sigma_RMS_mu_LO_reassign.pdf");

 auto legend_e1 = new TLegend(0.75,0.15,0.9,0.3);
legend_e1->AddEntry(h_sigma_el_pre_1hitFirstModules_clones,"Pre-vrtx","LEP");
legend_e1->AddEntry(h_sigma_el_post_1hitFirstModules_clones,"Post_vrtx","LEP");

 auto legend_e2 = new TLegend(0.75,0.15,0.9,0.3);
legend_e2->AddEntry(h_sigma_el_pre_2hitFirstModules_clones,"Pre-vrtx","LEP");
legend_e2->AddEntry(h_sigma_el_post_2hitFirstModules_clones,"Post_vrtx","LEP");

 auto legend_e_r1 = new TLegend(0.75,0.15,0.9,0.3);
legend_e_r1->AddEntry(h_RMS_el_pre_1hitFirstModules_clones,"Pre-vrtx","LEP");
legend_e_r1->AddEntry(h_RMS_el_post_1hitFirstModules_clones,"Post_vrtx","LEP");

 auto legend_e_r2 = new TLegend(0.75,0.15,0.9,0.3);
legend_e_r2->AddEntry(h_RMS_el_pre_2hitFirstModules_clones,"Pre-vrtx","LEP");
legend_e_r2->AddEntry(h_RMS_el_post_2hitFirstModules_clones,"Post_vrtx","LEP");

TCanvas s("s","s",1600,1600);
s.Divide(2,2);
s.cd(1);

TMultiGraph *el_1hit = new TMultiGraph();
el_1hit->SetMinimum(0.);
el_1hit->SetMaximum(3.);


el_1hit->Add(h_sigma_el_pre_1hitFirstModules_clones,"A*");
el_1hit->Add(h_sigma_el_post_1hitFirstModules_clones,"A*");
/*el_1hit->GetXaxis()->SetTitle("Electron scattering angle (mrad)");
el_1hit->GetYaxis()->SetTitle("Sigma (rad)");
el_1hit->SetTitle("Sigma VS electron scattering angle NLO1HitShared");*/
el_1hit->Draw("A*");
legend_e1->Draw();
gPad->SetLogy();
el_1hit->SetTitle("Sigma VS electron scattering angle LO 1HitShared");


s.cd(2);

TMultiGraph *el_2hit = new TMultiGraph();
el_2hit->SetMinimum(0.);
el_2hit->SetMaximum(3.);


el_2hit->Add(h_sigma_el_pre_1hitFirstModules_clones,"A*");
el_2hit->Add(h_sigma_el_post_1hitFirstModules_clones,"A*");
el_2hit->Draw("A*");
legend_e2->Draw();

gPad->SetLogy();
el_2hit->SetTitle("Sigma VS electron scattering angle LO 2HitShared");

s.cd(3);
TMultiGraph *el_1hit_r = new TMultiGraph();
el_1hit_r->SetMinimum(0.);
el_1hit_r->SetMaximum(3.);


el_1hit_r->Add(h_RMS_el_pre_1hitFirstModules_clones,"A*");
el_1hit_r->Add(h_RMS_el_post_1hitFirstModules_clones,"A*");
el_1hit_r->Draw("A*");
legend_e_r1->Draw();
gPad->SetLogy();
el_1hit_r->SetTitle("RMS VS electron scattering angle LO 1HitShared");


/*h_RMS_el_pre_1hitFirstModules_clones->SetTitle("RMS VS Electron scattering angle LO pre-vrtx 1HitShared");
h_RMS_el_pre_1hitFirstModules_clones->GetXaxis()->SetTitle("Electron scattering angle (mrad)");
h_RMS_el_pre_1hitFirstModules_clones->GetYaxis()->SetTitle("RMS (rad)");*/
s.cd(4);

TMultiGraph *el_2hit_r = new TMultiGraph();
el_2hit_r->SetMinimum(0.);
el_2hit_r->SetMaximum(3.);


el_2hit_r->Add(h_RMS_el_pre_2hitFirstModules_clones,"A*");
el_2hit_r->Add(h_RMS_el_post_2hitFirstModules_clones,"A*");
el_2hit_r->Draw("A*");
legend_e_r2->Draw();
gPad->SetLogy();
el_2hit_r->SetTitle("RMS VS electron scattering angle LO 2HitShared");


/*h_RMS_el_pre_2hitFirstModules_clones->SetTitle("RMS VS Electron scattering angle LO pre-vrtx 2HitShared");
h_RMS_el_pre_2hitFirstModules_clones->GetYaxis()->SetTitle("RMS (rad)");
h_RMS_el_pre_2hitFirstModules_clones->GetXaxis()->SetTitle("Electron scattering angle (mrad)");*/

s.SaveAs("proposal/sigma_RMS_el_LO_reassign.pdf");




 auto legend_mu1_pt2 = new TLegend(0.75,0.15,0.9,0.3);
legend_mu1_pt2->AddEntry(h_sigma_mu_post_1hitFirstModules_clones,"1 hit shared","LEP");
legend_mu1_pt2->AddEntry(h_sigma_mu_post_2hitFirstModules_clones,"2 hit shared","LEP");

 auto legend_mu2_pt2 = new TLegend(0.75,0.15,0.9,0.3);
legend_mu2_pt2->AddEntry(h_sigma_mu_pre_1hitFirstModules_clones,"1 hit shared","LEP");
legend_mu2_pt2->AddEntry(h_sigma_mu_pre_2hitFirstModules_clones,"2 hit shared","LEP");


 auto legend_mu_r1_pt2 = new TLegend(0.75,0.15,0.9,0.3);
legend_mu_r1_pt2->AddEntry(h_RMS_mu_pre_1hitFirstModules_clones,"1 hit shared","LEP");
legend_mu_r1_pt2->AddEntry(h_RMS_mu_pre_2hitFirstModules_clones,"2 hit shared","LEP");


 auto legend_mu_r2_pt2 = new TLegend(0.75,0.15,0.9,0.3);
legend_mu_r2_pt2->AddEntry(h_RMS_mu_post_1hitFirstModules_clones,"1 hit shared","LEP");
legend_mu_r2_pt2->AddEntry(h_RMS_mu_post_2hitFirstModules_clones,"2 hit shared","LEP");


 auto legend_el1_pt2 = new TLegend(0.75,0.15,0.9,0.3);
legend_el1_pt2->AddEntry(h_sigma_el_post_1hitFirstModules_clones,"1 hit shared","LEP");
legend_el1_pt2->AddEntry(h_sigma_el_post_2hitFirstModules_clones,"2 hit shared","LEP");

 auto legend_el2_pt2 = new TLegend(0.75,0.15,0.9,0.3);
legend_el2_pt2->AddEntry(h_sigma_el_pre_1hitFirstModules_clones,"1 hit shared","LEP");
legend_el2_pt2->AddEntry(h_sigma_el_pre_2hitFirstModules_clones,"2 hit shared","LEP");


 auto legend_el_r1_pt2 = new TLegend(0.75,0.15,0.9,0.3);
legend_el_r1_pt2->AddEntry(h_RMS_el_pre_1hitFirstModules_clones,"1 hit shared","LEP");
legend_el_r1_pt2->AddEntry(h_RMS_el_pre_2hitFirstModules_clones,"2 hit shared","LEP");


 auto legend_el_r2_pt2 = new TLegend(0.75,0.15,0.9,0.3);
legend_el_r2_pt2->AddEntry(h_RMS_el_post_1hitFirstModules_clones,"1 hit shared","LEP");
legend_el_r2_pt2->AddEntry(h_RMS_el_post_2hitFirstModules_clones,"2 hit shared","LEP");


h_sigma_mu_pre_1hitFirstModules_clones->SetMarkerColor(kRed);
h_sigma_mu_post_2hitFirstModules_clones->SetMarkerColor(kAzure);

h_RMS_mu_pre_1hitFirstModules_clones->SetMarkerColor(kRed);
h_RMS_mu_post_2hitFirstModules_clones->SetMarkerColor(kAzure);

h_sigma_el_pre_1hitFirstModules_clones->SetMarkerColor(kRed);
h_sigma_el_post_2hitFirstModules_clones->SetMarkerColor(kAzure);

h_RMS_el_pre_1hitFirstModules_clones->SetMarkerColor(kRed);
h_RMS_el_post_2hitFirstModules_clones->SetMarkerColor(kAzure);


TCanvas s1_pt2("s1_pt2","s1_pt2",1600,1600);
s1_pt2.Divide(2,2);
s1_pt2.cd(1);

TMultiGraph *mu_1hit_pt2 = new TMultiGraph();
mu_1hit_pt2->SetMinimum(0.);
mu_1hit_pt2->SetMaximum(3.);


mu_1hit_pt2->Add(h_sigma_mu_pre_1hitFirstModules_clones,"A*");
mu_1hit_pt2->Add(h_sigma_mu_pre_2hitFirstModules_clones,"A*");

mu_1hit_pt2->Draw("A*");
legend_mu1_pt2->Draw();
gPad->SetLogy();
mu_1hit_pt2->SetTitle("Sigma VS muon scattering angle LO pre-vrtx");

/*h_sigma_mu_pre_1hitFirstModules_clones->SetTitle("Sigma VS muon scattering angle LO pre-vrtx 1HitShared");
h_sigma_mu_pre_1hitFirstModules_clones->GetXaxis()->SetTitle("Muon scattering angle (mrad)");
h_sigma_mu_pre_1hitFirstModules_clones->GetYaxis()->SetTitle("Sigma (rad)");*/

s1_pt2.cd(2);

TMultiGraph *mu_2hit_pt2 = new TMultiGraph();
mu_2hit_pt2->SetMinimum(0.);
mu_2hit_pt2->SetMaximum(3.);


mu_2hit_pt2->Add(h_sigma_mu_post_1hitFirstModules_clones,"A*");
mu_2hit_pt2->Add(h_sigma_mu_post_2hitFirstModules_clones,"A*");
mu_2hit_pt2->Draw("A*");
legend_mu2_pt2->Draw();
gPad->SetLogy();
mu_2hit_pt2->SetTitle("Sigma VS muon scattering angle LO post-vrtx");

s1_pt2.cd(3);
TMultiGraph *mu_1hit_r_pt2 = new TMultiGraph();
mu_1hit_r_pt2->SetMinimum(0.);
mu_1hit_r_pt2->SetMaximum(3.);


mu_1hit_r_pt2->Add(h_RMS_mu_pre_1hitFirstModules_clones,"A*");
mu_1hit_r_pt2->Add(h_RMS_mu_pre_2hitFirstModules_clones,"A*");
mu_1hit_r_pt2->Draw("A*");
legend_mu_r1_pt2->Draw();
gPad->SetLogy(); 
mu_1hit_r_pt2->SetTitle("RMS VS muon scattering angle LO pre-vrtx");


s1_pt2.cd(4);

TMultiGraph *mu_2hit_r_pt2 = new TMultiGraph();
mu_2hit_r_pt2->SetMinimum(0.);
mu_2hit_r_pt2->SetMaximum(3.);


mu_2hit_r_pt2->Add(h_RMS_mu_post_1hitFirstModules_clones,"A*");
mu_2hit_r_pt2->Add(h_RMS_mu_post_2hitFirstModules_clones,"A*");
mu_2hit_r_pt2->Draw("A*");
legend_mu_r2_pt2->Draw();
gPad->SetLogy();
mu_2hit_r_pt2->SetTitle("RMS VS muon scattering angle LO post-vrtx");


s1_pt2.SaveAs("proposal/sigma_RMS_mu_LO_12hit_reassign.pdf");


TCanvas s_pt2("s_pt2","s_pt2",1600,1600);
s_pt2.Divide(2,2);
s_pt2.cd(1);

TMultiGraph *el_1hit_pt2 = new TMultiGraph();
el_1hit_pt2->SetMinimum(0.);
el_1hit_pt2->SetMaximum(3.);


el_1hit_pt2->Add(h_sigma_el_pre_1hitFirstModules_clones,"A*");
el_1hit_pt2->Add(h_sigma_el_pre_2hitFirstModules_clones,"A*");

el_1hit_pt2->Draw("A*");
legend_el1_pt2->Draw();
gPad->SetLogy();
el_1hit_pt2->SetTitle("Sigma VS electron scattering angle LO pre-vrtx");


s_pt2.cd(2);

TMultiGraph *el_2hit_pt2 = new TMultiGraph();
el_2hit_pt2->SetMinimum(0.);
el_2hit_pt2->SetMaximum(3.);


el_2hit_pt2->Add(h_sigma_el_post_1hitFirstModules_clones,"A*");
el_2hit_pt2->Add(h_sigma_el_post_2hitFirstModules_clones,"A*");
el_2hit_pt2->Draw("A*");
legend_el2_pt2->Draw();
gPad->SetLogy();
el_2hit_pt2->SetTitle("Sigma VS electron scattering angle LO post-vrtx");

s_pt2.cd(3);
TMultiGraph *el_1hit_r_pt2 = new TMultiGraph();
el_1hit_r_pt2->SetMinimum(0.);
el_1hit_r_pt2->SetMaximum(3.);


el_1hit_r_pt2->Add(h_RMS_el_pre_1hitFirstModules_clones,"A*");
el_1hit_r_pt2->Add(h_RMS_el_pre_2hitFirstModules_clones,"A*");
el_1hit_r_pt2->Draw("A*");
legend_el_r1_pt2->Draw();
gPad->SetLogy(); 
el_1hit_r_pt2->SetTitle("RMS VS electron scattering angle LO pre-vrtx");


s_pt2.cd(4);

TMultiGraph *el_2hit_r_pt2 = new TMultiGraph();
el_2hit_r_pt2->SetMinimum(0.);
el_2hit_r_pt2->SetMaximum(3.);


el_2hit_r_pt2->Add(h_RMS_el_post_1hitFirstModules_clones,"A*");
el_2hit_r_pt2->Add(h_RMS_el_post_2hitFirstModules_clones,"A*");
el_2hit_r_pt2->Draw("A*");
legend_el_r2_pt2->Draw();
gPad->SetLogy();
el_2hit_r_pt2->SetTitle("RMS VS electron scattering angle LO post-vrtx");


s_pt2.SaveAs("proposal/sigma_RMS_el_LO_12hit_reassign.pdf");







}


