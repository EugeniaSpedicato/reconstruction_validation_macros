
#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TStyle.h>


void proposal(){



array<TFile*,14> f_pre_1hitFirstModules_el, f_post_1hitFirstModules_el, f_pre_2hitFirstModules_el, f_post_2hitFirstModules_el, f_pre_0hit_el, f_post_0hit_el;
array<TFile*,14> f_pre_1hitFirstModules_mu, f_post_1hitFirstModules_mu, f_pre_2hitFirstModules_mu, f_post_2hitFirstModules_mu, f_pre_0hit_mu, f_post_0hit_mu;;

for(int m=0; m<14; m++){

f_pre_1hitFirstModules_el.at(m)=TFile::Open(Form("quality_tracks/root/res_prevrtx_clones_el_%d_1hitFirstModules_LO_reassign.root",static_cast<char>(m)));
f_post_1hitFirstModules_el.at(m)=TFile::Open(Form("quality_tracks/root/res_vrtx_clones_el_%d_1hitFirstModules_LO_reassign.root",static_cast<char>(m)));
f_pre_2hitFirstModules_el.at(m)=TFile::Open(Form("quality_tracks/root/res_prevrtx_clones_el_%d_2hitFirstModules_LO_reassign.root",static_cast<char>(m)));
f_post_2hitFirstModules_el.at(m)=TFile::Open(Form("quality_tracks/root/res_vrtx_clones_el_%d_2hitFirstModules_LO_reassign.root",static_cast<char>(m)));
f_pre_0hit_el.at(m)=TFile::Open(Form("quality_tracks/root/res_prevrtx_clones_el_%d_0hit_LO_reassign.root",static_cast<char>(m)));
f_post_0hit_el.at(m)=TFile::Open(Form("quality_tracks/root/res_vrtx_clones_el_%d_0hit_LO_reassign.root",static_cast<char>(m)));
}

for(int m=0; m<14; m++){
f_pre_1hitFirstModules_mu.at(m)=TFile::Open(Form("quality_tracks/root/res_prevrtx_clones_AngleMu_%d_1hitFirstModules_LO_reassign.root",static_cast<char>(m)));
f_post_1hitFirstModules_mu.at(m)=TFile::Open(Form("quality_tracks/root/res_vrtx_clones_AngleMu_%d_1hitFirstModules_LO_reassign.root",static_cast<char>(m)));
f_pre_2hitFirstModules_mu.at(m)=TFile::Open(Form("quality_tracks/root/res_prevrtx_clones_AngleMu_%d_2hitFirstModules_LO_reassign.root",static_cast<char>(m)));
f_post_2hitFirstModules_mu.at(m)=TFile::Open(Form("quality_tracks/root/res_vrtx_clones_AngleMu_%d_2hitFirstModules_LO_reassign.root",static_cast<char>(m)));
f_pre_0hit_mu.at(m)=TFile::Open(Form("quality_tracks/root/res_prevrtx_clones_AngleMu_%d_0hit_LO_reassign.root",static_cast<char>(m)));
f_post_0hit_mu.at(m)=TFile::Open(Form("quality_tracks/root/res_vrtx_clones_AngleMu_%d_0hit_LO_reassign.root",static_cast<char>(m)));
}


TH1::SetDefaultSumw2(kTRUE);

array <TH1D*,14> h_pre_1hitFirstModules_el, h_post_1hitFirstModules_el, h_pre_2hitFirstModules_el, h_post_2hitFirstModules_el, h_pre_0hit_el, h_post_0hit_el;
array <TH1D*,14> h_pre_1hitFirstModules_mu, h_post_1hitFirstModules_mu, h_pre_2hitFirstModules_mu, h_post_2hitFirstModules_mu, h_pre_0hit_mu, h_post_0hit_mu;


for(int m=0; m<14; m++){

h_pre_1hitFirstModules_el.at(m)=(TH1D*)f_pre_1hitFirstModules_el.at(m)->Get(Form("res%d_el_pre_clones",static_cast<char>(m)));
h_post_1hitFirstModules_el.at(m)=(TH1D*)f_post_1hitFirstModules_el.at(m)->Get(Form("res%d_vrtx_el_clones",static_cast<char>(m)));
h_pre_2hitFirstModules_el.at(m)=(TH1D*)f_pre_2hitFirstModules_el.at(m)->Get(Form("res%d_el_pre_clones",static_cast<char>(m)));
h_post_2hitFirstModules_el.at(m)=(TH1D*)f_post_2hitFirstModules_el.at(m)->Get(Form("res%d_vrtx_el_clones",static_cast<char>(m)));
h_pre_0hit_el.at(m)=(TH1D*)f_pre_0hit_el.at(m)->Get(Form("res%d_el_pre_clones",static_cast<char>(m)));
h_post_0hit_el.at(m)=(TH1D*)f_post_0hit_el.at(m)->Get(Form("res%d_vrtx_el_clones",static_cast<char>(m)));
}

for(int m=0; m<14; m++){
h_pre_1hitFirstModules_mu.at(m)=(TH1D*)f_pre_1hitFirstModules_mu.at(m)->Get(Form("res%d_mu_pre_clones",static_cast<char>(m)));
h_post_1hitFirstModules_mu.at(m)=(TH1D*)f_post_1hitFirstModules_mu.at(m)->Get(Form("res%d_vrtx_mu_clones",static_cast<char>(m)));
h_pre_2hitFirstModules_mu.at(m)=(TH1D*)f_pre_2hitFirstModules_mu.at(m)->Get(Form("res%d_mu_pre_clones",static_cast<char>(m)));
h_post_2hitFirstModules_mu.at(m)=(TH1D*)f_post_2hitFirstModules_mu.at(m)->Get(Form("res%d_vrtx_mu_clones",static_cast<char>(m)));
h_pre_0hit_mu.at(m)=(TH1D*)f_pre_0hit_mu.at(m)->Get(Form("res%d_mu_pre_clones",static_cast<char>(m)));
h_post_0hit_mu.at(m)=(TH1D*)f_post_0hit_mu.at(m)->Get(Form("res%d_vrtx_mu_clones",static_cast<char>(m)));

}


   const Int_t NBINS_mu = 14;
   Double_t edges_mu[NBINS_mu + 1] = {0.00005,0.00015,0.00025,0.00035,0.00045,0.00055,0.00065,0.00075,0.00085,0.00095,0.0015,0.0025,0.0035,0.0045};

   const Int_t NBINS = 14;
   Double_t edges_el[NBINS + 1] = {0.0005,0.0015,0.0025,0.0035,0.0045,0.0055,0.0065,0.0075,0.0085,0.0095, 0.0125, 0.0175, 0.0225, 0.0285};



auto h_sigma_mu_pre_1hitFirstModules_clones = new TGraph();
auto h_sigma_el_pre_1hitFirstModules_clones = new TGraph();
auto h_sigma_mu_post_1hitFirstModules_clones = new TGraph();
auto h_sigma_el_post_1hitFirstModules_clones = new TGraph();
auto h_sigma_mu_pre_2hitFirstModules_clones = new TGraph();
auto h_sigma_el_pre_2hitFirstModules_clones = new TGraph();
auto h_sigma_mu_post_2hitFirstModules_clones = new TGraph();
auto h_sigma_el_post_2hitFirstModules_clones = new TGraph();

auto h_sigma_mu_pre_0hit_clones = new TGraph();
auto h_sigma_el_pre_0hit_clones = new TGraph();
auto h_sigma_mu_post_0hit_clones = new TGraph();
auto h_sigma_el_post_0hit_clones = new TGraph();



for(int m=0; m<14; m++){
TF1 *g1 = new TF1("g1", "gaus");
if(m<4) h_pre_1hitFirstModules_el.at(m)->Fit("g1","R","",-0.0005,0.0005);
else if(m==4) h_pre_1hitFirstModules_el.at(m)->Fit("g1","R","",-0.002,0.002);
else h_pre_1hitFirstModules_el.at(m)->Fit("g1");
h_sigma_el_pre_1hitFirstModules_clones->SetPoint(h_sigma_el_pre_1hitFirstModules_clones->GetN(),edges_el[m]*1000,g1->GetParameter(2)*1000);
}

for(int m=0; m<14; m++){
TGaxis::SetMaxDigits(3);
TF1 *g1 = new TF1("g1", "gaus");
if(m<4) h_post_1hitFirstModules_el.at(m)->Fit("g1","R","",-0.0005,0.0005);
else if(m==4) h_post_1hitFirstModules_el.at(m)->Fit("g1","R","",-0.002,0.002);
else h_post_1hitFirstModules_el.at(m)->Fit("g1");
h_sigma_el_post_1hitFirstModules_clones->SetPoint(h_sigma_el_post_1hitFirstModules_clones->GetN(),edges_el[m]*1000,g1->GetParameter(2)*1000);
}


for(int m=0; m<14; m++){
TF1 *g1 = new TF1("g1", "gaus");
if(m<4) h_pre_2hitFirstModules_el.at(m)->Fit("g1","R","",-0.0005,0.0005);
else if(m==4) h_pre_2hitFirstModules_el.at(m)->Fit("g1","R","",-0.002,0.002);
else h_pre_2hitFirstModules_el.at(m)->Fit("g1");
h_sigma_el_pre_2hitFirstModules_clones->SetPoint(h_sigma_el_pre_2hitFirstModules_clones->GetN(),edges_el[m]*1000,g1->GetParameter(2)*1000);
}


TCanvas n2("n2","n2",1400,1400);
n2.Divide(2,2);
n2.cd(1);
TGaxis::SetMaxDigits(3);
TF1 *g = new TF1("g", "gaus");
h_pre_2hitFirstModules_el.at(0)->Fit("g","R","",-0.0003,0.0003);
h_pre_2hitFirstModules_el.at(0)->GetXaxis()->SetTitle("rad");
h_pre_2hitFirstModules_el.at(0)->SetLineColor(kAzure+7);
h_pre_2hitFirstModules_el.at(0)->Draw("hist same");
n2.cd(2);
h_pre_2hitFirstModules_el.at(2)->Fit("g","R","",-0.0004,0.0004);;
h_pre_2hitFirstModules_el.at(2)->GetXaxis()->SetTitle("rad");
h_pre_2hitFirstModules_el.at(2)->SetLineColor(kAzure+7);
h_pre_2hitFirstModules_el.at(2)->Draw("hist same");
n2.cd(3);
h_pre_2hitFirstModules_el.at(9)->Fit("g");
h_pre_2hitFirstModules_el.at(9)->GetXaxis()->SetTitle("rad");
h_pre_2hitFirstModules_el.at(9)->SetLineColor(kAzure+7);
h_pre_2hitFirstModules_el.at(9)->Draw("hist same");
n2.cd(4);
h_pre_2hitFirstModules_el.at(13)->Fit("g");
h_pre_2hitFirstModules_el.at(13)->GetXaxis()->SetTitle("rad");
h_pre_2hitFirstModules_el.at(13)->SetLineColor(kAzure+7);
h_pre_2hitFirstModules_el.at(13)->Draw("hist same");
gStyle->SetStatH(0.1);
gStyle->SetOptFit(0001);
gStyle->SetOptStat(1101);

n2.SaveAs("proposal/prop_res_pre_el_2hitFirstModules.pdf");





for(int m=0; m<14; m++){
TGaxis::SetMaxDigits(3);
TF1 *g1 = new TF1("g1", "gaus");
if(m<4) h_post_2hitFirstModules_el.at(m)->Fit("g1","R","",-0.0005,0.0005);
else if(m==4) h_post_2hitFirstModules_el.at(m)->Fit("g1","R","",-0.002,0.002);
else h_post_2hitFirstModules_el.at(m)->Fit("g1");
h_sigma_el_post_2hitFirstModules_clones->SetPoint(h_sigma_el_post_2hitFirstModules_clones->GetN(),edges_el[m]*1000,g1->GetParameter(2)*1000);
}




for(int m=0; m<14; m++){
TF1 *g1 = new TF1("g1", "gaus");
if(m<4) h_pre_0hit_el.at(m)->Fit("g1","R","",-0.0005,0.0005);
else if(m==4) h_pre_0hit_el.at(m)->Fit("g1","R","",-0.002,0.002);
else h_pre_0hit_el.at(m)->Fit("g1");
h_sigma_el_pre_0hit_clones->SetPoint(h_sigma_el_pre_0hit_clones->GetN(),edges_el[m]*1000,g1->GetParameter(2)*1000);
}


for(int m=0; m<14; m++){
TF1 *g1 = new TF1("g1", "gaus");
if(m<4) h_post_0hit_el.at(m)->Fit("g1","R","",-0.0005,0.0005);
else if(m==4) h_post_0hit_el.at(m)->Fit("g1","R","",-0.002,0.002);
else h_post_0hit_el.at(m)->Fit("g1");
h_sigma_el_post_0hit_clones->SetPoint(h_sigma_el_post_0hit_clones->GetN(),edges_el[m]*1000,g1->GetParameter(2)*1000);

}





for(int m=0; m<14; m++){
if(m==0) h_pre_1hitFirstModules_mu.at(m)->GetXaxis()->SetRangeUser(-0.001,0.001);
TF1 *g1 = new TF1("g1", "gaus");
h_pre_1hitFirstModules_mu.at(m)->Fit("g1");
h_sigma_mu_pre_1hitFirstModules_clones->SetPoint(h_sigma_mu_pre_1hitFirstModules_clones->GetN(),edges_mu[m]*1000,g1->GetParameter(2)*1000);

}


for(int m=0; m<14; m++){
TF1 *g1 = new TF1("g1", "gaus");
h_post_1hitFirstModules_mu.at(m)->Fit("g1");
h_sigma_mu_post_1hitFirstModules_clones->SetPoint(h_sigma_mu_post_1hitFirstModules_clones->GetN(),edges_mu[m]*1000,g1->GetParameter(2)*1000);
}



for(int m=0; m<14; m++){
TF1 *g1 = new TF1("g1", "gaus");
h_pre_2hitFirstModules_mu.at(m)->Fit("g1");
h_sigma_mu_pre_2hitFirstModules_clones->SetPoint(h_sigma_mu_pre_2hitFirstModules_clones->GetN(),edges_mu[m]*1000,g1->GetParameter(2)*1000);
}


TCanvas n2_mu("n2_mu","n2_mu",1400,1400);
n2_mu.Divide(2,2);
n2_mu.cd(1);
TF1 *g2 = new TF1("g2", "gaus");
h_pre_2hitFirstModules_mu.at(2)->Fit("g2");
h_pre_2hitFirstModules_mu.at(2)->SetLineColor(kAzure+7);
h_pre_2hitFirstModules_mu.at(2)->GetXaxis()->SetTitle("rad");
h_pre_2hitFirstModules_mu.at(2)->Draw("hist same");
n2_mu.cd(2);
h_pre_2hitFirstModules_mu.at(9)->Fit("g2","R","",-0.00012,0.00012);
h_pre_2hitFirstModules_mu.at(9)->SetLineColor(kAzure+7);
h_pre_2hitFirstModules_mu.at(9)->GetXaxis()->SetTitle("rad");
//h_pre_2hitFirstModules_mu.at(9)->Rebin(2);
h_pre_2hitFirstModules_mu.at(9)->Draw("hist same");
n2_mu.cd(3);
h_pre_2hitFirstModules_mu.at(11)->Fit("g2","R","",-0.0002,0.0002);
h_pre_2hitFirstModules_mu.at(11)->SetLineColor(kAzure+7);
h_pre_2hitFirstModules_mu.at(11)->GetXaxis()->SetTitle("rad");
//h_pre_2hitFirstModules_mu.at(11)->Rebin(2);
h_pre_2hitFirstModules_mu.at(11)->Draw("hist same");
n2_mu.cd(4);
h_pre_2hitFirstModules_mu.at(13)->Rebin(2);
h_pre_2hitFirstModules_mu.at(13)->SetMinimum(1.);
h_pre_2hitFirstModules_mu.at(13)->Fit("g2");
h_pre_2hitFirstModules_mu.at(13)->SetLineColor(kAzure+7);
h_pre_2hitFirstModules_mu.at(13)->GetXaxis()->SetTitle("rad");
h_pre_2hitFirstModules_mu.at(13)->Draw("hist same");
gStyle->SetStatH(0.1);
gStyle->SetOptFit(0001);
gStyle->SetOptStat(1101);

n2_mu.SaveAs("proposal/prop_res_pre_AngleMu_2hitFirstModules.pdf");


for(int m=0; m<14; m++){
TF1 *g1 = new TF1("g1", "gaus");
if(m>6)h_post_2hitFirstModules_mu.at(m)->Fit("g1","R","",-0.0002,0.0002);
else h_post_2hitFirstModules_mu.at(m)->Fit("g1");
h_sigma_mu_post_2hitFirstModules_clones->SetPoint(h_sigma_mu_post_2hitFirstModules_clones->GetN(),edges_mu[m]*1000,g1->GetParameter(2)*1000);
}


for(int m=0; m<14; m++){
TF1 *g1 = new TF1("g1", "gaus");
h_pre_0hit_mu.at(m)->Fit("g1");
h_sigma_mu_pre_0hit_clones->SetPoint(h_sigma_mu_pre_0hit_clones->GetN(),edges_mu[m]*1000,g1->GetParameter(2)*1000);
}


for(int m=0; m<14; m++){
TF1 *g1 = new TF1("g1", "gaus");
if(m>6)h_post_0hit_mu.at(m)->Fit("g1","R","",-0.0002,0.0002);
else h_post_0hit_mu.at(m)->Fit("g1");
h_sigma_mu_post_0hit_clones->SetPoint(h_sigma_mu_post_0hit_clones->GetN(),edges_mu[m]*1000,g1->GetParameter(2)*1000);
}








 auto legend_mu2_pt2 = new TLegend(0.75,0.15,0.9,0.3);
legend_mu2_pt2->AddEntry(h_sigma_mu_pre_0hit_clones,"0 hit shared","LEP");
legend_mu2_pt2->AddEntry(h_sigma_mu_pre_1hitFirstModules_clones,"1 hit shared","LEP");
legend_mu2_pt2->AddEntry(h_sigma_mu_pre_2hitFirstModules_clones,"2 hit shared","LEP");



 auto legend_el2_pt2 = new TLegend(0.75,0.15,0.9,0.3);
legend_el2_pt2->AddEntry(h_sigma_el_pre_0hit_clones,"0 hit shared","LEP");
legend_el2_pt2->AddEntry(h_sigma_el_pre_1hitFirstModules_clones,"1 hit shared","LEP");
legend_el2_pt2->AddEntry(h_sigma_el_pre_2hitFirstModules_clones,"2 hit shared","LEP");



 auto legend_mu0_lin = new TLegend(0.75,0.75,0.9,0.9);
legend_mu0_lin->AddEntry(h_sigma_mu_pre_0hit_clones,"pre-vrtx","LEP");
legend_mu0_lin->AddEntry(h_sigma_mu_post_0hit_clones,"post-vrtx","LEP");

 auto legend_mu1_lin = new TLegend(0.75,0.75,0.9,0.9);
legend_mu1_lin->AddEntry(h_sigma_mu_pre_1hitFirstModules_clones,"pre-vrtx","LEP");
legend_mu1_lin->AddEntry(h_sigma_mu_post_1hitFirstModules_clones,"post-vrtx","LEP");

 auto legend_mu2_lin = new TLegend(0.75,0.75,0.9,0.9);
legend_mu2_lin->AddEntry(h_sigma_mu_pre_2hitFirstModules_clones,"pre-vrtx","LEP");
legend_mu2_lin->AddEntry(h_sigma_mu_post_2hitFirstModules_clones,"post-vrtx","LEP");




 auto legend_mu0 = new TLegend(0.75,0.15,0.9,0.3);
legend_mu0->AddEntry(h_sigma_mu_pre_0hit_clones,"pre-vrtx","LEP");
legend_mu0->AddEntry(h_sigma_mu_post_0hit_clones,"post-vrtx","LEP");

 auto legend_mu1 = new TLegend(0.75,0.15,0.9,0.3);
legend_mu1->AddEntry(h_sigma_mu_pre_1hitFirstModules_clones,"pre-vrtx","LEP");
legend_mu1->AddEntry(h_sigma_mu_post_1hitFirstModules_clones,"post-vrtx","LEP");

 auto legend_mu2 = new TLegend(0.75,0.15,0.9,0.3);
legend_mu2->AddEntry(h_sigma_mu_pre_2hitFirstModules_clones,"pre-vrtx","LEP");
legend_mu2->AddEntry(h_sigma_mu_post_2hitFirstModules_clones,"post-vrtx","LEP");


 auto legend_e0 = new TLegend(0.75,0.15,0.9,0.3);
legend_e0->AddEntry(h_sigma_mu_pre_0hit_clones,"pre-vrtx","LEP");
legend_e0->AddEntry(h_sigma_mu_post_0hit_clones,"post-vrtx","LEP");

 auto legend_e1 = new TLegend(0.75,0.15,0.9,0.3);
legend_e1->AddEntry(h_sigma_mu_pre_1hitFirstModules_clones,"pre-vrtx","LEP");
legend_e1->AddEntry(h_sigma_mu_post_1hitFirstModules_clones,"post-vrtx","LEP");

 auto legend_e2 = new TLegend(0.75,0.15,0.9,0.3);
legend_e2->AddEntry(h_sigma_mu_pre_2hitFirstModules_clones,"pre-vrtx","LEP");
legend_e2->AddEntry(h_sigma_mu_post_2hitFirstModules_clones,"post-vrtx","LEP");




TCanvas s2_pt2("s2_pt2","s2_pt2",1400,1400);
s2_pt2.Divide(1,2);
s2_pt2.cd(1);

TMultiGraph *mu_1hit_pt2 = new TMultiGraph();


h_sigma_mu_pre_1hitFirstModules_clones->SetMarkerColor(kGreen+1);
h_sigma_mu_pre_2hitFirstModules_clones->SetMarkerColor(kAzure+7);
h_sigma_mu_pre_0hit_clones->SetMarkerColor(kOrange+10);

mu_1hit_pt2->Add(h_sigma_mu_pre_0hit_clones,"A*");
mu_1hit_pt2->Add(h_sigma_mu_pre_1hitFirstModules_clones,"A*");
mu_1hit_pt2->Add(h_sigma_mu_pre_2hitFirstModules_clones,"A*");

mu_1hit_pt2->Draw("A*");
legend_mu2_pt2->Draw();
mu_1hit_pt2->SetMinimum(0.03);
mu_1hit_pt2->SetMaximum(0.14);
mu_1hit_pt2->SetTitle("Muon angular resolution");
mu_1hit_pt2->GetXaxis()->SetTitle("Scattering angle (mrad)");
mu_1hit_pt2->GetYaxis()->SetTitle("#sigma(#theta) (mrad)");
mu_1hit_pt2->GetHistogram()->GetXaxis()->SetLimits(0.,5.);

s2_pt2.cd(2);

TMultiGraph *el_1hit_pt2 = new TMultiGraph();

h_sigma_el_pre_1hitFirstModules_clones->SetMarkerColor(kGreen+1);
h_sigma_el_pre_2hitFirstModules_clones->SetMarkerColor(kAzure+7);
h_sigma_el_pre_0hit_clones->SetMarkerColor(kOrange+10);


el_1hit_pt2->Add(h_sigma_el_pre_0hit_clones,"A*");
el_1hit_pt2->Add(h_sigma_el_pre_1hitFirstModules_clones,"A*");
el_1hit_pt2->Add(h_sigma_el_pre_2hitFirstModules_clones,"A*");

el_1hit_pt2->Draw("A*");
legend_el2_pt2->Draw();
el_1hit_pt2->SetMinimum(0.);
el_1hit_pt2->SetMaximum(3.5);
el_1hit_pt2->SetTitle("Electron angular resolution");
el_1hit_pt2->GetXaxis()->SetTitle("Scattering angle (mrad)");
el_1hit_pt2->GetYaxis()->SetTitle("#sigma(#theta) (mrad)");
el_1hit_pt2->GetXaxis()->SetLimits(0.,32.);

s2_pt2.SaveAs("proposal/sigma_linear.pdf");






TCanvas s1_pt2("s1_pt2","s1_pt2",1400,1400);
s1_pt2.Divide(1,2);
s1_pt2.cd(1);

mu_1hit_pt2->SetMinimum(0.008);
mu_1hit_pt2->SetMaximum(0.2);

mu_1hit_pt2->Draw("A*");
legend_mu2_pt2->Draw();
gPad->SetLogy();
mu_1hit_pt2->SetTitle("Muon angular resolution");


s1_pt2.cd(2);

el_1hit_pt2->SetMinimum(0.008);
el_1hit_pt2->SetMaximum(5.);

el_1hit_pt2->Draw("A*");
legend_el2_pt2->Draw();
gPad->SetLogy();
el_1hit_pt2->SetTitle("Electron angular resolution");


s1_pt2.SaveAs("proposal/sigma_log.pdf");





TCanvas s2("s2","s2",2100,1400);
s2.Divide(3,2);

s2.cd(1);
TMultiGraph *mu_0hit = new TMultiGraph();


h_sigma_mu_pre_0hit_clones->SetMarkerColor(kOrange+10);
h_sigma_mu_post_0hit_clones->SetMarkerColor(kBlack);

mu_0hit->Add(h_sigma_mu_pre_0hit_clones,"A*");
mu_0hit->Add(h_sigma_mu_post_0hit_clones,"A*");

mu_0hit->Draw("A*");
legend_mu0_lin->Draw();
mu_0hit->SetMinimum(0.03);
mu_0hit->SetMaximum(0.2);
mu_0hit->SetTitle("Muon angular resolution 0 hit shared");
mu_0hit->GetXaxis()->SetTitle("Scattering angle (mrad)");
mu_0hit->GetYaxis()->SetTitle("#sigma(#theta) (mrad)");
mu_0hit->GetHistogram()->GetXaxis()->SetLimits(0.,5.);

s2.cd(2);

TMultiGraph *mu_1hit = new TMultiGraph();


h_sigma_mu_pre_1hitFirstModules_clones->SetMarkerColor(kGreen+1);
h_sigma_mu_post_1hitFirstModules_clones->SetMarkerColor(kBlack);

mu_1hit->Add(h_sigma_mu_pre_1hitFirstModules_clones,"A*");
mu_1hit->Add(h_sigma_mu_post_1hitFirstModules_clones,"A*");

mu_1hit->Draw("A*");
legend_mu1_lin->Draw();
mu_1hit->SetMinimum(0.03);
mu_1hit->SetMaximum(0.2);
mu_1hit->SetTitle("Muon angular resolution 1 hit shared");
mu_1hit->GetXaxis()->SetTitle("Scattering angle (mrad)");
mu_1hit->GetYaxis()->SetTitle("#sigma(#theta) (mrad)");
mu_1hit->GetHistogram()->GetXaxis()->SetLimits(0.,5.);

s2.cd(3);

TMultiGraph *mu_2hit = new TMultiGraph();


h_sigma_mu_pre_2hitFirstModules_clones->SetMarkerColor(kAzure+7);
h_sigma_mu_post_2hitFirstModules_clones->SetMarkerColor(kBlack);

mu_2hit->Add(h_sigma_mu_pre_2hitFirstModules_clones,"A*");
mu_2hit->Add(h_sigma_mu_post_2hitFirstModules_clones,"A*");

mu_2hit->Draw("A*");
legend_mu2_lin->Draw();
mu_2hit->SetMinimum(0.03);
mu_2hit->SetMaximum(0.2);
mu_2hit->SetTitle("Muon angular resolution 2 hit shared");
mu_2hit->GetXaxis()->SetTitle("Scattering angle (mrad)");
mu_2hit->GetYaxis()->SetTitle("#sigma(#theta) (mrad)");
mu_2hit->GetHistogram()->GetXaxis()->SetLimits(0.,5.);


s2.cd(4);

TMultiGraph *el_0hit = new TMultiGraph();

h_sigma_el_pre_0hit_clones->SetMarkerColor(kOrange+10);
h_sigma_el_post_0hit_clones->SetMarkerColor(kBlack);


el_0hit->Add(h_sigma_el_pre_0hit_clones,"A*");
el_0hit->Add(h_sigma_el_post_0hit_clones,"A*");

el_0hit->Draw("A*");
legend_e0->Draw();
el_0hit->SetMinimum(0.);
el_0hit->SetMaximum(5.);
el_0hit->SetTitle("Electron angular resolution 0 hit shared");
el_0hit->GetXaxis()->SetTitle("Scattering angle (mrad)");
el_0hit->GetYaxis()->SetTitle("#sigma(#theta) (mrad)");
el_0hit->GetXaxis()->SetLimits(0.,32.);

s2.cd(5);

TMultiGraph *el_1hit = new TMultiGraph();

h_sigma_el_pre_1hitFirstModules_clones->SetMarkerColor(kGreen+1);
h_sigma_el_post_1hitFirstModules_clones->SetMarkerColor(kBlack);


el_1hit->Add(h_sigma_el_pre_1hitFirstModules_clones,"A*");
el_1hit->Add(h_sigma_el_post_1hitFirstModules_clones,"A*");

el_1hit->Draw("A*");
legend_e1->Draw();
el_1hit->SetMinimum(0.);
el_1hit->SetMaximum(5.);
el_1hit->SetTitle("Electron angular resolution 1 hit shared");
el_1hit->GetXaxis()->SetTitle("Scattering angle (mrad)");
el_1hit->GetYaxis()->SetTitle("#sigma(#theta) (mrad)");
el_1hit->GetXaxis()->SetLimits(0.,32.);


s2.cd(6);

TMultiGraph *el_2hit = new TMultiGraph();

h_sigma_el_pre_2hitFirstModules_clones->SetMarkerColor(kAzure+7);
h_sigma_el_post_2hitFirstModules_clones->SetMarkerColor(kBlack);


el_2hit->Add(h_sigma_el_pre_2hitFirstModules_clones,"A*");
el_2hit->Add(h_sigma_el_post_2hitFirstModules_clones,"A*");

el_2hit->Draw("A*");
legend_e2->Draw();
el_2hit->SetMinimum(0.);
el_2hit->SetMaximum(5.);
el_2hit->SetTitle("Electron angular resolution 2 hit shared");
el_2hit->GetXaxis()->SetTitle("Scattering angle (mrad)");
el_2hit->GetYaxis()->SetTitle("#sigma(#theta) (mrad)");
el_2hit->GetXaxis()->SetLimits(0.,32.);

s2.SaveAs("proposal/sigma_linear_prepost.pdf");




TCanvas s1("s1","s1",2100,1400);
s1.Divide(3,2);
s1.cd(1);

mu_0hit->SetMinimum(0.008);
mu_0hit->SetMaximum(0.2);

mu_0hit->Draw("A*");
legend_mu0->Draw();
gPad->SetLogy();

s1.cd(2);

mu_1hit->SetMinimum(0.008);
mu_1hit->SetMaximum(0.2);

mu_1hit->Draw("A*");
legend_mu1->Draw();
gPad->SetLogy();

s1.cd(3);

mu_2hit->SetMinimum(0.008);
mu_2hit->SetMaximum(0.2);

mu_2hit->Draw("A*");
legend_mu2->Draw();
gPad->SetLogy();


s1.cd(4);

el_0hit->SetMinimum(0.008);
el_0hit->SetMaximum(5.);

el_0hit->Draw("A*");
legend_e0->Draw();
gPad->SetLogy();

s1.cd(5);

el_1hit->SetMinimum(0.008);
el_1hit->SetMaximum(5.);

el_1hit->Draw("A*");
legend_e1->Draw();
gPad->SetLogy();

s1.cd(6);

el_2hit->SetMinimum(0.008);
el_2hit->SetMaximum(5.);

el_2hit->Draw("A*");
legend_e2->Draw();
gPad->SetLogy();


s1.SaveAs("proposal/sigma_log_prepost.pdf");



}


