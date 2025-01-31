#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TStyle.h>


void compare_version_2hit(string path1, string path2, string info1, string info2){



TFile * f_el_pre = TFile::Open(Form("%ssigma_el_pre_2hit_%s.root",path1.c_str(),info1.c_str()));
TFile * f_el_post = TFile::Open(Form("%ssigma_el_post_2hit_%s.root",path1.c_str(),info1.c_str()));
TFile * f_mu_pre = TFile::Open(Form("%ssigma_mu_pre_2hit_%s.root",path1.c_str(),info1.c_str()));
TFile * f_mu_post = TFile::Open(Form("%ssigma_mu_post_2hit_%s.root",path1.c_str(),info1.c_str()));


TFile *f_el_2 = TFile::Open(Form("%sel_2_%s.root",path1.c_str(),info1.c_str()));
TFile *f_mu_2 = TFile::Open(Form("%smu_2_%s.root",path1.c_str(),info1.c_str()));


TFile * f_el_pre_testing =TFile::Open(Form("%ssigma_el_pre_2hit_%s.root",path2.c_str(),info2.c_str()));
TFile * f_el_post_testing = TFile::Open(Form("%ssigma_el_post_2hit_%s.root",path2.c_str(),info2.c_str()));
TFile * f_mu_pre_testing = TFile::Open(Form("%ssigma_mu_pre_2hit_%s.root",path2.c_str(),info2.c_str()));
TFile * f_mu_post_testing = TFile::Open(Form("%ssigma_mu_post_2hit_%s.root",path2.c_str(),info2.c_str()));


TFile *f_el_2_testing = TFile::Open(Form("%sel_2_%s.root",path2.c_str(),info2.c_str()));
TFile *f_mu_2_testing = TFile::Open(Form("%smu_2_%s.root",path2.c_str(),info2.c_str()));


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



TGraph *el_2_testing=(TGraph*)f_el_2_testing->Get("el_2");
TGraph *mu_2_testing=(TGraph*)f_mu_2_testing->Get("mu_2");


TGraph *el_2=(TGraph*)f_el_2->Get("el_2");
TGraph *mu_2=(TGraph*)f_mu_2->Get("mu_2");


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
legend_mu2->AddEntry(mu_2,Form("%s",info1.c_str()),"LEP");
legend_mu2->AddEntry(mu_2_testing,Form("%s",info2.c_str()),"LEP");

 auto legend_e2 = new TLegend(0.75,0.75,0.9,0.9);
legend_e2->AddEntry(el_2,Form("%s",info1.c_str()),"LEP");
legend_e2->AddEntry(el_2_testing,Form("%s",info2.c_str()),"LEP");


TCanvas d("d","d",2100,1400);
d.Divide(1,2);
d.cd(1);
TMultiGraph *mg2 = new TMultiGraph();
el_2->SetMarkerColor(kRed);
el_2_testing->SetMarkerColor(kBlue);
mg2->Add(el_2,"A*");
mg2->Add(el_2_testing,"A*");
mg2->Draw("A*");
legend_e2->Draw();
mg2->SetTitle("Sigma difference el 2 hit shared");
mg2->GetYaxis()->SetTitle("(#sigma_post - #sigma_pre)/#sigma_pre");
mg2->GetXaxis()->SetTitle("#theta_#el (mrad)");

d.cd(2);
TMultiGraph *mg4 = new TMultiGraph();
mu_2->SetMarkerColor(kRed);
mu_2_testing->SetMarkerColor(kBlue);
mg4->Add(mu_2,"A*");
mg4->Add(mu_2_testing,"A*");
mg4->Draw("A*");
legend_mu2->Draw();
mg4->SetTitle("Sigma difference #mu 2 hit shared");
mg4->GetYaxis()->SetTitle("(#sigma_post - #sigma_pre)/#sigma_pre");
mg4->GetXaxis()->SetTitle("#theta_#mu (mrad)");

d.SaveAs(Form("%sdiff_2hit_%s.pdf",path2.c_str(),info2.c_str()));



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
mg1pre->SetTitle("Pre sigma el 2 hit shared");
mg1pre->GetYaxis()->SetTitle("#sigma_pre");
mg1pre->GetXaxis()->SetTitle("#theta_#el (mrad)");
mg1pre->SetMinimum(0.0001);


d2.cd(2);
TMultiGraph *mg2post = new TMultiGraph();
el_post->SetMarkerColor(kRed);
el_post_testing->SetMarkerColor(kBlue);
mg2post->Add(el_post,"A*");
mg2post->Add(el_post_testing,"A*");
mg2post->Draw("A*");
legend_e2post->Draw();
mg2post->SetTitle("Post sigma el 2 hit shared");
mg2post->GetYaxis()->SetTitle("#sigma_post");
mg2post->GetXaxis()->SetTitle("#theta_#el (mrad)");
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
mg3pre->SetTitle("Pre sigma #mu 2 hit shared");
mg3pre->GetYaxis()->SetTitle("#sigma_pre");
mg3pre->GetXaxis()->SetTitle("#theta_#mu (mrad)");

d2.cd(4);
TMultiGraph *mg4post = new TMultiGraph();
mu_post->SetMarkerColor(kRed);
mu_post_testing->SetMarkerColor(kBlue);
mg4post->Add(mu_post,"A*");
mg4post->Add(mu_post_testing,"A*");
mg4post->Draw("A*");
legend_mu2post->Draw();
mg4post->SetTitle("Post sigma #mu 2 hit shared");
mg4post->GetYaxis()->SetTitle("#sigma_post");
mg4post->GetXaxis()->SetTitle("#theta_#mu (mrad)");
mg4post->SetMinimum(0.0001);

d2.SaveAs(Form("%ssigma_overlap_2hit_%s.pdf",path2.c_str(),info2.c_str()));


h_op_LO_clone0->SetTitle("Eff. elastic event VS opening angle 0 hit shared");
h_op_LO_clone1->SetTitle("Eff. elastic event VS opening angle 1 hit shared");
h_op_LO_clone2->SetTitle("Eff. elastic event VS opening angle 2 hit shared");
h_theta_mu_single_LO_clone0->SetTitle("Eff. single particle VS #theta mu 0 hit shared");
h_theta_mu_single_LO_clone1->SetTitle("Eff. single particle VS #theta mu 1 hit shared");
h_theta_mu_single_LO_clone2->SetTitle("Eff. single particle VS #theta mu 2 hit shared");
h_theta_e_single_LO_clone0->SetTitle("Eff. single particle VS #theta e 0 hit shared");
h_theta_e_single_LO_clone1->SetTitle("Eff. single particle VS #theta e 1 hit shared");
h_theta_e_single_LO_clone2->SetTitle("Eff. single particle VS #theta e 2 hit shared");


//info1="Default_v.14";
//info2="This_work";

TCanvas e("e","e",2300,2100);
e.Divide(3,3);
e.cd(1);
h_op_LO_clone0->SetTitle("0 hit shared");
h_op_LO_clone0->GetXaxis()->SetTitle("Opening angle [rad]");
h_op_LO_clone0->GetYaxis()->SetTitle("Event reconstruction efficiency");

h_op_LO_clone0->SetLineColor(kRed);
h_op_LO_clone0_testing->SetLineColor(kBlue);
h_op_LO_clone0_testing->SetLineStyle(2);
h_op_LO_clone0->SetMinimum(0.);
h_op_LO_clone0->SetMaximum(1.05);
TGaxis::SetMaxDigits(3);
h_op_LO_clone0->Draw();
h_op_LO_clone0_testing->Draw(" same");
 auto l0 = new TLegend(0.50,0.15,0.9,0.4);
l0->AddEntry(h_op_LO_clone0,Form("%s",info1.c_str()),"LEP");
l0->AddEntry(h_op_LO_clone0_testing,Form("%s",info2.c_str()),"LEP");
l0->Draw("same");
gStyle->SetOptStat(0);

e.cd(2);
h_op_LO_clone1->SetTitle("1 hit shared");
h_op_LO_clone1->GetXaxis()->SetTitle("Opening angle [rad]");
h_op_LO_clone1->GetYaxis()->SetTitle("Event reconstruction efficiency");

h_op_LO_clone1->SetLineColor(kRed);
h_op_LO_clone1_testing->SetLineColor(kBlue);
h_op_LO_clone1_testing->SetLineStyle(2);
h_op_LO_clone1->SetMinimum(0.);
h_op_LO_clone1->SetMaximum(1.05);
TGaxis::SetMaxDigits(3);
h_op_LO_clone1->Draw();
h_op_LO_clone1_testing->Draw(" same");
 auto l1 = new TLegend(0.50,0.15,0.9,0.4);
l1->AddEntry(h_op_LO_clone1,Form("%s",info1.c_str()),"LEP");
l1->AddEntry(h_op_LO_clone1_testing,Form("%s",info2.c_str()),"LEP");
l1->Draw("same");
gStyle->SetOptStat(0);

e.cd(3);
h_op_LO_clone2->SetTitle("2 hit shared");
h_op_LO_clone2->GetXaxis()->SetTitle("Opening angle [rad]");
h_op_LO_clone2->GetYaxis()->SetTitle("Event reconstruction efficiency");

h_op_LO_clone2->SetLineColor(kRed);
h_op_LO_clone2_testing->SetLineColor(kBlue);
h_op_LO_clone2_testing->SetLineStyle(2);
h_op_LO_clone2->SetMinimum(0.);
h_op_LO_clone2->SetMaximum(1.05);
TGaxis::SetMaxDigits(3);
h_op_LO_clone2->Draw();
h_op_LO_clone2_testing->Draw(" same");
 auto l2 = new TLegend(0.50,0.15,0.9,0.4);
l2->AddEntry(h_op_LO_clone2,Form("%s",info1.c_str()),"LEP");
l2->AddEntry(h_op_LO_clone2_testing,Form("%s",info2.c_str()),"LEP");
l2->Draw("same");
gStyle->SetOptStat(0);

e.cd(4);
h_theta_mu_single_LO_clone0->SetTitle("0 hit shared");
h_theta_mu_single_LO_clone0->GetXaxis()->SetTitle("Muon angle [rad]");
h_theta_mu_single_LO_clone0->GetYaxis()->SetTitle("Muon reconstruction efficiency");

h_theta_mu_single_LO_clone0->SetLineColor(kRed);
h_theta_mu_single_LO_clone0_testing->SetLineColor(kBlue);
h_theta_mu_single_LO_clone0_testing->SetLineStyle(2);
h_theta_mu_single_LO_clone0->SetMinimum(0.6);
h_theta_mu_single_LO_clone0->SetMaximum(1.05);
TGaxis::SetMaxDigits(3);
h_theta_mu_single_LO_clone0->Draw();
h_theta_mu_single_LO_clone0_testing->Draw(" same");
 auto l0_a = new TLegend(0.50,0.15,0.9,0.4);
l0_a->AddEntry(h_theta_mu_single_LO_clone0,Form("%s",info1.c_str()),"LEP");
l0_a->AddEntry(h_theta_mu_single_LO_clone0_testing,Form("%s",info2.c_str()),"LEP");
l0_a->Draw("same");
gStyle->SetOptStat(0);

e.cd(5);
h_theta_mu_single_LO_clone1->SetTitle("1 hit shared");
h_theta_mu_single_LO_clone1->GetXaxis()->SetTitle("Muon angle [rad]");
h_theta_mu_single_LO_clone1->GetYaxis()->SetTitle("Muon reconstruction efficiency");

h_theta_mu_single_LO_clone1->SetLineColor(kRed);
h_theta_mu_single_LO_clone1_testing->SetLineColor(kBlue);
h_theta_mu_single_LO_clone1_testing->SetLineStyle(2);
h_theta_mu_single_LO_clone1->SetMinimum(0.6);
h_theta_mu_single_LO_clone1->SetMaximum(1.05);
TGaxis::SetMaxDigits(3);
h_theta_mu_single_LO_clone1->Draw();
h_theta_mu_single_LO_clone1_testing->Draw(" same");
 auto l1_a = new TLegend(0.50,0.15,0.9,0.4);
l1_a->AddEntry(h_theta_mu_single_LO_clone1,Form("%s",info1.c_str()),"LEP");
l1_a->AddEntry(h_theta_mu_single_LO_clone1_testing,Form("%s",info2.c_str()),"LEP");
l1_a->Draw("same");
gStyle->SetOptStat(0);

e.cd(6);
h_theta_mu_single_LO_clone2->SetTitle("2 hit shared");
h_theta_mu_single_LO_clone2->GetXaxis()->SetTitle("Muon angle [rad]");
h_theta_mu_single_LO_clone2->GetYaxis()->SetTitle("Muon reconstruction efficiency");

h_theta_mu_single_LO_clone2->SetLineColor(kRed);
h_theta_mu_single_LO_clone2_testing->SetLineColor(kBlue);
h_theta_mu_single_LO_clone2_testing->SetLineStyle(2);
h_theta_mu_single_LO_clone2->SetMinimum(0.6);
h_theta_mu_single_LO_clone2->SetMaximum(1.05);
TGaxis::SetMaxDigits(3);
h_theta_mu_single_LO_clone2->Draw();
h_theta_mu_single_LO_clone2_testing->Draw(" same");
 auto l2_a = new TLegend(0.50,0.15,0.9,0.4);
l2_a->AddEntry(h_theta_mu_single_LO_clone2,Form("%s",info1.c_str()),"LEP");
l2_a->AddEntry(h_theta_mu_single_LO_clone2_testing,Form("%s",info2.c_str()),"LEP");
l2_a->Draw("same");
gStyle->SetOptStat(0);

e.cd(7);
h_theta_e_single_LO_clone0->SetTitle("0 hit shared");
h_theta_e_single_LO_clone0->GetXaxis()->SetTitle("Electron angle [rad]");
h_theta_e_single_LO_clone0->GetYaxis()->SetTitle("Electron reconstruction efficiency");

h_theta_e_single_LO_clone0->SetLineColor(kRed);
h_theta_e_single_LO_clone0_testing->SetLineColor(kBlue);
h_theta_e_single_LO_clone0_testing->SetLineStyle(2);
h_theta_e_single_LO_clone0->SetMinimum(0.6);
h_theta_e_single_LO_clone0->SetMaximum(1.05);
TGaxis::SetMaxDigits(3);

h_theta_e_single_LO_clone0->Draw();
h_theta_e_single_LO_clone0_testing->Draw(" same");
 auto l0_e = new TLegend(0.50,0.15,0.9,0.4);
l0_e->AddEntry(h_theta_e_single_LO_clone0,Form("%s",info1.c_str()),"LEP");
l0_e->AddEntry(h_theta_e_single_LO_clone0_testing,Form("%s",info2.c_str()),"LEP");
l0_e->Draw("same");

gStyle->SetOptStat(0);

e.cd(8);
gStyle->SetOptStat(0);

h_theta_e_single_LO_clone1->SetTitle("1 hit shared");
h_theta_e_single_LO_clone1->GetXaxis()->SetTitle("Electron angle [rad]");
h_theta_e_single_LO_clone1->GetYaxis()->SetTitle("Electron reconstruction efficiency");

h_theta_e_single_LO_clone1->SetLineColor(kRed);
h_theta_e_single_LO_clone1_testing->SetLineColor(kBlue);
h_theta_e_single_LO_clone1_testing->SetLineStyle(2);
h_theta_e_single_LO_clone1->SetMinimum(0.6);
h_theta_e_single_LO_clone1->SetMaximum(1.05);
TGaxis::SetMaxDigits(3);
h_theta_e_single_LO_clone1->Draw();
h_theta_e_single_LO_clone1_testing->Draw(" same");
 auto l1_e = new TLegend(0.50,0.15,0.9,0.4);
l1_e->AddEntry(h_theta_e_single_LO_clone1,Form("%s",info1.c_str()),"LEP");
l1_e->AddEntry(h_theta_e_single_LO_clone1_testing,Form("%s",info2.c_str()),"LEP");
l1_e->Draw("same");


e.cd(9);
h_theta_e_single_LO_clone2->SetTitle("2 hit shared");
h_theta_e_single_LO_clone2->GetXaxis()->SetTitle("Electron angle [rad]");
h_theta_e_single_LO_clone2->GetYaxis()->SetTitle("Electron reconstruction efficiency");

h_theta_e_single_LO_clone2->SetLineColor(kRed);
h_theta_e_single_LO_clone2_testing->SetLineColor(kBlue);
h_theta_e_single_LO_clone2_testing->SetLineStyle(2);
h_theta_e_single_LO_clone2->SetMinimum(0.6);
h_theta_e_single_LO_clone2->SetMaximum(1.05);
TGaxis::SetMaxDigits(3);

h_theta_e_single_LO_clone2->Draw();
h_theta_e_single_LO_clone2_testing->Draw(" same");
 auto l2_e = new TLegend(0.50,0.15,0.9,0.4);
l2_e->AddEntry(h_theta_e_single_LO_clone2,Form("%s",info1.c_str()),"LEP");
l2_e->AddEntry(h_theta_e_single_LO_clone2_testing,Form("%s",info2.c_str()),"LEP");
l2_e->Draw("same");
gStyle->SetOptStat(0);

e.SaveAs(Form("%scompare_eff_%s.pdf",path2.c_str(),info2.c_str()));


}



