#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TStyle.h>


void compare_version_0hit(string path1, string path2, string info1, string info2){


TFile * f_el_pre = TFile::Open(Form("%ssigma_el_pre_0hit_%s.root",path1.c_str(),info1.c_str()));
TFile * f_el_post = TFile::Open(Form("%ssigma_el_post_0hit_%s.root",path1.c_str(),info1.c_str()));
TFile * f_mu_pre = TFile::Open(Form("%ssigma_mu_pre_0hit_%s.root",path1.c_str(),info1.c_str()));
TFile * f_mu_post = TFile::Open(Form("%ssigma_mu_post_0hit_%s.root",path1.c_str(),info1.c_str()));


TFile *f_el_0 = TFile::Open(Form("%sel_0_%s.root",path1.c_str(),info1.c_str()));
TFile *f_mu_0 = TFile::Open(Form("%smu_0_%s.root",path1.c_str(),info1.c_str()));


TFile * f_el_pre_testing =TFile::Open(Form("%ssigma_el_pre_0hit_%s.root",path2.c_str(),info2.c_str()));
TFile * f_el_post_testing = TFile::Open(Form("%ssigma_el_post_0hit_%s.root",path2.c_str(),info2.c_str()));
TFile * f_mu_pre_testing = TFile::Open(Form("%ssigma_mu_pre_0hit_%s.root",path2.c_str(),info2.c_str()));
TFile * f_mu_post_testing = TFile::Open(Form("%ssigma_mu_post_0hit_%s.root",path2.c_str(),info2.c_str()));


TFile *f_el_0_testing = TFile::Open(Form("%sel_0_%s.root",path2.c_str(),info2.c_str()));
TFile *f_mu_0_testing = TFile::Open(Form("%smu_0_%s.root",path2.c_str(),info2.c_str()));


TFile *f_eff = TFile::Open(Form("%seff_%s.root",path1.c_str(),info1.c_str()));
TFile *f_eff_testing = TFile::Open(Form("%seff_%s.root",path2.c_str(),info2.c_str()));



TH1::SetDefaultSumw2(kTRUE);

TGraph *el_pre_testing=(TGraph*)f_el_pre_testing->Get("el_pre");
TGraph *mu_pre_testing=(TGraph*)f_mu_pre_testing->Get("mu_pre");
TGraph *el_post_testing=(TGraph*)f_el_post_testing->Get("el_post");
TGraph *mu_post_testing=(TGraph*)f_mu_post_testing->Get("mu_post");


TGraph *el_pre=(TGraph*)f_el_pre->Get("el_pre");
TGraph *mu_pre=(TGraph*)f_mu_pre->Get("mu_pre");
TGraph *el_post=(TGraph*)f_el_post->Get("el_post");
TGraph *mu_post=(TGraph*)f_mu_post->Get("mu_post");



TGraph *el_0_testing=(TGraph*)f_el_0_testing->Get("el_0");
TGraph *mu_0_testing=(TGraph*)f_mu_0_testing->Get("mu_0");


TGraph *el_0=(TGraph*)f_el_0->Get("el_0");
TGraph *mu_0=(TGraph*)f_mu_0->Get("mu_0");


TH1D *h_op_LO_clone0=(TH1D*)f_eff->Get("h_op_LO_clone0");
TH1D *h_op_LO_clone1=(TH1D*)f_eff->Get("h_op_LO_clone1");
TH1D *h_op_LO_clone2=(TH1D*)f_eff->Get("h_op_LO_clone2");
TH1D *h_theta_mu_single_LO_clone0=(TH1D*)f_eff->Get("h_theta_mu_single_LO_clone0");
TH1D *h_theta_mu_single_LO_clone1=(TH1D*)f_eff->Get("h_theta_mu_single_LO_clone1");
TH1D *h_theta_mu_single_LO_clone2=(TH1D*)f_eff->Get("h_theta_mu_single_LO_clone2");
TH1D *h_theta_e_single_LO_clone0=(TH1D*)f_eff->Get("h_theta_e_single_LO_clone0");
TH1D *h_theta_e_single_LO_clone1=(TH1D*)f_eff->Get("h_theta_e_single_LO_clone1");
TH1D *h_theta_e_single_LO_clone2=(TH1D*)f_eff->Get("h_theta_e_single_LO_clone2");

TH1D *h_op_LO_clone0_testing=(TH1D*)f_eff_testing->Get("h_op_LO_clone0");
TH1D *h_op_LO_clone1_testing=(TH1D*)f_eff_testing->Get("h_op_LO_clone1");
TH1D *h_op_LO_clone2_testing=(TH1D*)f_eff_testing->Get("h_op_LO_clone2");
TH1D *h_theta_mu_single_LO_clone0_testing=(TH1D*)f_eff_testing->Get("h_theta_mu_single_LO_clone0");
TH1D *h_theta_mu_single_LO_clone1_testing=(TH1D*)f_eff_testing->Get("h_theta_mu_single_LO_clone1");
TH1D *h_theta_mu_single_LO_clone2_testing=(TH1D*)f_eff_testing->Get("h_theta_mu_single_LO_clone2");
TH1D *h_theta_e_single_LO_clone0_testing=(TH1D*)f_eff_testing->Get("h_theta_e_single_LO_clone0");
TH1D *h_theta_e_single_LO_clone1_testing=(TH1D*)f_eff_testing->Get("h_theta_e_single_LO_clone1");
TH1D *h_theta_e_single_LO_clone2_testing=(TH1D*)f_eff_testing->Get("h_theta_e_single_LO_clone2");


 auto legend_mu2 = new TLegend(0.75,0.15,0.9,0.3);
legend_mu2->AddEntry(mu_0,Form("%s",info1.c_str()),"LEP");
legend_mu2->AddEntry(mu_0_testing,Form("%s",info2.c_str()),"LEP");

 auto legend_e2 = new TLegend(0.75,0.75,0.9,0.9);
legend_e2->AddEntry(el_0,Form("%s",info1.c_str()),"LEP");
legend_e2->AddEntry(el_0_testing,Form("%s",info2.c_str()),"LEP");


TCanvas d("d","d",2100,1400);
d.Divide(1,2);
d.cd(1);
TMultiGraph *mg2 = new TMultiGraph();
el_0->SetMarkerColor(kRed);
el_0_testing->SetMarkerColor(kBlue);
mg2->Add(el_0,"A*");
mg2->Add(el_0_testing,"A*");
mg2->Draw("A*");
legend_e2->Draw();
mg2->SetTitle("Sigma difference el 0 hit shared");
mg2->GetYaxis()->SetTitle("(#sigma_post [mrad] - #sigma_pre [mrad])/#sigma_pre [mrad]");
mg2->GetXaxis()->SetTitle("Electron angle [mrad]");

d.cd(2);
TMultiGraph *mg4 = new TMultiGraph();
mu_0->SetMarkerColor(kRed);
mu_0_testing->SetMarkerColor(kBlue);
mg4->Add(mu_0,"A*");
mg4->Add(mu_0_testing,"A*");
mg4->Draw("A*");
legend_mu2->Draw();
mg4->SetTitle("Sigma difference #mu 0 hit shared");
mg4->GetYaxis()->SetTitle("(#sigma_post [mrad] - #sigma_pre [mrad])/#sigma_pre [mrad]");
mg4->GetXaxis()->SetTitle("Muon angle [mrad]");

d.SaveAs(Form("%sdiff_0hit_%s.pdf",path2.c_str(),info2.c_str()));



 auto legend_mu1pre = new TLegend(0.75,0.15,0.9,0.3);
legend_mu1pre->AddEntry(mu_pre,"pre WP_14","LEP");
legend_mu1pre->AddEntry(mu_pre_testing,Form("pre WP_14 %s",info2.c_str()),"LEP");

 auto legend_mu2post = new TLegend(0.75,0.15,0.9,0.3);
legend_mu2post->AddEntry(mu_post,"post WP_14","LEP");
legend_mu2post->AddEntry(mu_post_testing,Form("post WP_14 %s",info2.c_str()),"LEP");


 auto legend_e1pre = new TLegend(0.75,0.15,0.9,0.3);
legend_e1pre->AddEntry(el_pre,"pre WP_14","LEP");
legend_e1pre->AddEntry(el_pre_testing,Form("pre WP_14 %s",info2.c_str()),"LEP");

 auto legend_e2post = new TLegend(0.75,0.15,0.9,0.3);
legend_e2post->AddEntry(el_post,"post WP_14","LEP");
legend_e2post->AddEntry(el_post_testing,Form("post WP_14 %s",info2.c_str()),"LEP");


auto el_pre_diff = new TGraphErrors();
auto el_post_diff = new TGraphErrors();
auto mu_pre_diff = new TGraphErrors();
auto mu_post_diff = new TGraphErrors();

   const Int_t NBINS_mu = 14;
   Double_t edges_mu[NBINS_mu + 1] = {0.00005,0.00015,0.00025,0.00035,0.00045,0.00055,0.00065,0.00075,0.00085,0.00095,0.0015,0.0025,0.0035,0.0045};

   const Int_t NBINS = 14;
   Double_t edges_el[NBINS + 1] = {0.0005,0.0015,0.0025,0.0035,0.0045,0.0055,0.0065,0.0075,0.0085,0.0095, 0.0125, 0.0175, 0.0225, 0.0285};


for(int m=0; m< el_pre->GetN(); m++){
el_pre_diff->SetPoint(m,edges_el[m]*1000,el_pre->GetPointY(m)-el_pre_testing->GetPointY(m));
el_post_diff->SetPoint(m,edges_el[m]*1000,el_post->GetPointY(m)-el_post_testing->GetPointY(m));
}

for(int m=0; m< mu_pre->GetN(); m++){
mu_pre_diff->SetPoint(m,edges_mu[m]*1000,mu_pre->GetPointY(m)-mu_pre_testing->GetPointY(m));
mu_post_diff->SetPoint(m,edges_mu[m]*1000,mu_post->GetPointY(m)-mu_post_testing->GetPointY(m));
}



TCanvas d2("d2","d2",2100,1400);
d2.Divide(2,2);
d2.cd(1);
TMultiGraph *mg1pre = new TMultiGraph();
el_pre->SetMarkerColor(kRed);
el_pre_testing->SetMarkerColor(kBlue);
mg1pre->Add(el_pre,"A*");
mg1pre->Add(el_pre_testing,"A*");
mg1pre->Draw("A*");
legend_e1pre->Draw();
mg1pre->SetTitle("Electron angular resolution - PRE Vertexing");
mg1pre->GetYaxis()->SetTitle("#sigma_pre [mrad]");
mg1pre->GetXaxis()->SetTitle("Electron angle [mrad]");
mg1pre->SetMinimum(0.0001);


d2.cd(2);
TMultiGraph *mg2post = new TMultiGraph();
el_post->SetMarkerColor(kRed);
el_post_testing->SetMarkerColor(kBlue);
mg2post->Add(el_post,"A*");
mg2post->Add(el_post_testing,"A*");
mg2post->Draw("A*");
legend_e2post->Draw();
mg2post->SetTitle("Electron angular resolution - POST Vertexing");
mg2post->GetYaxis()->SetTitle("#sigma_post [mrad]");
mg2post->GetXaxis()->SetTitle("Electron angle [mrad]");
mg2post->SetMinimum(0.0001);


d2.cd(3);
TMultiGraph *mg3pre = new TMultiGraph();
mu_pre->SetMarkerColor(kRed);
mu_pre_testing->SetMarkerColor(kBlue);
mg3pre->Add(mu_pre,"A*");
mg3pre->Add(mu_pre_testing,"A*");
mg3pre->Draw("A*");
legend_mu1pre->Draw();
mg3pre->SetMinimum(0.0001);
mg3pre->SetTitle("Muon angular resolution - PRE Vertexing");
mg3pre->GetYaxis()->SetTitle("#sigma_pre [mrad]");
mg3pre->GetXaxis()->SetTitle("Muon angle [mrad]");

d2.cd(4);
TMultiGraph *mg4post = new TMultiGraph();
mu_post->SetMarkerColor(kRed);
mu_post_testing->SetMarkerColor(kBlue);
mg4post->Add(mu_post,"A*");
mg4post->Add(mu_post_testing,"A*");
mg4post->Draw("A*");
legend_mu2post->Draw();
mg4post->SetTitle("Muon angular resolution - POST Vertexing");
mg4post->GetYaxis()->SetTitle("#sigma_post [mrad]");
mg4post->GetXaxis()->SetTitle("Muon angle [mrad]");
mg4post->SetMinimum(0.0001);

d2.SaveAs(Form("%ssigma_overlap_0hit_%s.pdf",path2.c_str(),info2.c_str()));


h_op_LO_clone0->SetTitle("Eff. elastic event VS opening angle 0 hit shared");
h_op_LO_clone1->SetTitle("Eff. elastic event VS opening angle 1 hit shared");
h_op_LO_clone2->SetTitle("Eff. elastic event VS opening angle 0 hit shared");
h_theta_mu_single_LO_clone0->SetTitle("Eff. single particle VS #theta mu 0 hit shared");
h_theta_mu_single_LO_clone1->SetTitle("Eff. single particle VS #theta mu 1 hit shared");
h_theta_mu_single_LO_clone2->SetTitle("Eff. single particle VS #theta mu 0 hit shared");
h_theta_e_single_LO_clone0->SetTitle("Eff. single particle VS #theta e 0 hit shared");
h_theta_e_single_LO_clone1->SetTitle("Eff. single particle VS #theta e 1 hit shared");
h_theta_e_single_LO_clone2->SetTitle("Eff. single particle VS #theta e 0 hit shared");


}



