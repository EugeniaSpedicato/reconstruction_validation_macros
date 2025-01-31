#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include <TGraphErrors.h>
#include "TCanvas.h"
#include <TStyle.h>


void plots_pre_post_vrtx_2hit(string path,string info){


array<TFile*,84> f_pre_1hit_el, f_post_1hit_el, f_pre_2hit_el, f_post_2hit_el,f_umb_e,f_umb_e_bv;
array<TFile*,84> f_pre_1hit_mu, f_post_1hit_mu, f_pre_2hit_mu, f_post_2hit_mu,f_umb_mu,f_umb_mu_bv;

for(int m=0; m<14; m++){
 for(int id=0; id<6; id++){
string range;
if(id==0)range = "0_5";
else if(id==1)range = "5_10";
else if(id==2)range = "10_15";
else if(id==3)range = "15_20";
else if(id==4)range = "20_25";
else if(id==5)range = "25_32";
//f_pre_1hit_el.at(6*m+id)=TFile::Open(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/res_prevrtx_clones_el_%d_1hit_LO_index%s_258faf6b.root",static_cast<char>(m),range.c_str(),info.c_str()));
//f_post_1hit_el.at(6*m+id)=TFile::Open(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/res_vrtx_clones_el_%d_1hit_LO_index%s_258faf6b.root",static_cast<char>(m),range.c_str(),info.c_str()));
f_pre_2hit_el.at(6*m+id)=TFile::Open(Form("%sres_prevrtx_clones_el_%d_2hit_LO_index%s_%s.root",path.c_str(),static_cast<char>(m),range.c_str(),info.c_str()));
f_post_2hit_el.at(6*m+id)=TFile::Open(Form("%sres_vrtx_clones_el_%d_2hit_LO_index%s_%s.root",path.c_str(),static_cast<char>(m),range.c_str(),info.c_str()));

if(m==0){f_umb_e.at(id)=TFile::Open(Form("%sumb_e_2hit_index%s_%s.root",path.c_str(),range.c_str(),info.c_str()));
f_umb_e_bv.at(id)=TFile::Open(Form("%sumb_e_bv_2hit_index%s_%s.root",path.c_str(),range.c_str(),info.c_str()));}
 }
}

for(int m=0; m<14; m++){
 for(int id=0; id<6; id++){
string range;
if(id==0)range = "0_5";
else if(id==1)range = "5_10";
else if(id==2)range = "10_15";
else if(id==3)range = "15_20";
else if(id==4)range = "20_25";
else if(id==5)range = "25_32";
//f_pre_1hit_mu.at(6*m+id)=TFile::Open(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/res_prevrtx_clones_mu_%d_1hit_LO_index%s_258faf6b_%s.root",static_cast<char>(m),range.c_str(),info.c_str()));
//f_post_1hit_mu.at(6*m+id)=TFile::Open(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/res_vrtx_clones_mu_%d_1hit_LO_index%s_258faf6b_%s.root",static_cast<char>(m),range.c_str(),info.c_str()));
f_pre_2hit_mu.at(6*m+id)=TFile::Open(Form("%sres_prevrtx_clones_mu_%d_2hit_LO_index%s_%s.root",path.c_str(),static_cast<char>(m),range.c_str(),info.c_str()));
f_post_2hit_mu.at(6*m+id)=TFile::Open(Form("%sres_vrtx_clones_mu_%d_2hit_LO_index%s_%s.root",path.c_str(),static_cast<char>(m),range.c_str(),info.c_str()));

if(m==0){f_umb_mu.at(id)=TFile::Open(Form("%sumb_mu_2hit_index%s_%s.root",path.c_str(),range.c_str(),info.c_str()));
f_umb_mu_bv.at(id)=TFile::Open(Form("%sumb_mu_bv_2hit_index%s_%s.root",path.c_str(),range.c_str(),info.c_str()));}
 }
}

TH1::SetDefaultSumw2(kTRUE);


array <TH1D*,14> h_pre_1hit_el, h_post_1hit_el, h_pre_2hit_el, h_post_2hit_el;
array <TH1D*,14> h_pre_1hit_mu, h_post_1hit_mu, h_pre_2hit_mu, h_post_2hit_mu;

array <TH2D*,6> t0e{0};
array <TH2D*,6> t1e{0};
array <TH2D*,6> t0m{0};
array <TH2D*,6> t1m{0};

for(int id=0;id<6;id++){

        t0e.at(id)=(TH2D*)f_umb_e.at(id)->Get("hres2fm");
        t1e.at(id)=(TH2D*)f_umb_e_bv.at(id)->Get("hres2fmbv");
        t0m.at(id)=(TH2D*)f_umb_mu.at(id)->Get("hres1fm");
        t1m.at(id)=(TH2D*)f_umb_mu_bv.at(id)->Get("hres1fmbv");

        if(id!=0){
        t0e.at(0)->Add(t0e.at(id));
        t1e.at(0)->Add(t1e.at(id));
        t0m.at(0)->Add(t0m.at(id));
        t1m.at(0)->Add(t1m.at(id));}
}

 auto legend0 = new TLegend(0.75,0.15,0.9,0.3);
 auto legend1 = new TLegend(0.75,0.15,0.9,0.3);


TCanvas aa("aa","aa",1400,1400);
aa.Divide(2,2);
aa.cd(1);
t0e.at(0)->SetMarkerColor(kRed);
t0e.at(0)->SetMarkerStyle(21);
t0e.at(0)->SetMarkerSize(0.5);
t0e.at(0)->Draw();

t1e.at(0)->SetMarkerColor(kBlack);
t1e.at(0)->SetMarkerStyle(21);
t1e.at(0)->SetMarkerSize(0.5);
t1e.at(0)->Draw("same");
aa.cd(2);

t0m.at(0)->SetMarkerStyle(21);
t0m.at(0)->SetMarkerSize(0.5);
t0m.at(0)->SetMarkerColor(kRed);
t0m.at(0)->Draw();

t1m.at(0)->SetMarkerColor(kBlack);
t1m.at(0)->SetMarkerStyle(21);
t1m.at(0)->SetMarkerSize(0.5);
t1m.at(0)->Draw("same");
aa.cd(3);
auto m0 = t0e.at(0)->ProfileX();
legend0->AddEntry(m0,"after vertex","LEP");

m0->SetMarkerColor(kRed);
m0->SetMaximum(0.001);
m0->SetMinimum(-0.018);
m0->Draw();

auto m1 = t1e.at(0)->ProfileX();
legend0->AddEntry(m1,"befor vertex","LEP");
m1->SetMarkerColor(kBlack);
m1->Draw("same");
legend0->Draw("same");

aa.cd(4);
auto m3 = t0m.at(0)->ProfileX();
legend1->AddEntry(m3,"after vertex","LEP");
m3->SetMarkerColor(kRed);
m3->Draw();
auto m4 = t1m.at(0)->ProfileX();
legend0->AddEntry(m4,"befor vertex","LEP");
m4->SetMarkerColor(kBlack);
m4->Draw("same");
legend1->Draw("same");
aa.SaveAs(Form("%sresults/umb_%s.pdf",path.c_str(),info.c_str()));

for(int m=0; m<14; m++){

array <TH1D*,6> t2{0};
array <TH1D*,6> t3{0};

 for(int id=0; id<6; id++){


        t2.at(id)=(TH1D*)f_pre_2hit_el.at(6*m+id)->Get(Form("res%d_el_pre_clones",static_cast<char>(m)));
        t3.at(id)=(TH1D*)f_post_2hit_el.at(6*m+id)->Get(Form("res%d_vrtx_el_clones",static_cast<char>(m)));

        if(id!=0){//t0.at(0)->Add(t0.at(id));
        //t1.at(0)->Add(t1.at(id));
        t2.at(0)->Add(t2.at(id));
        t3.at(0)->Add(t3.at(id));}
 }
//h_pre_1hit_el.at(m)=t0.at(0);
//h_post_1hit_el.at(m)=t1.at(0);
h_pre_2hit_el.at(m)=t2.at(0);
h_post_2hit_el.at(m)=t3.at(0);
}

for(int m=0; m<14; m++){

array <TH1D*,6> t2{0};
array <TH1D*,6> t3{0};

 for(int id=0; id<6; id++){
//        t0.at(id)=(TH1D*)f_pre_1hit_mu.at(6*m+id)->Get(Form("res%d_mu_pre_clones",static_cast<char>(m)));
//        t1.at(id)=(TH1D*)f_post_1hit_mu.at(6*m+id)->Get(Form("res%d_vrtx_mu_clones",static_cast<char>(m)));
        t2.at(id)=(TH1D*)f_pre_2hit_mu.at(6*m+id)->Get(Form("res%d_mu_pre_clones",static_cast<char>(m)));
        t3.at(id)=(TH1D*)f_post_2hit_mu.at(6*m+id)->Get(Form("res%d_vrtx_mu_clones",static_cast<char>(m)));

        if(id!=0){//t0.at(0)->Add(t0.at(id));
//        t1.at(0)->Add(t1.at(id));
        t2.at(0)->Add(t2.at(id));
        t3.at(0)->Add(t3.at(id));}
 }

//h_pre_1hit_mu.at(m)=t0.at(0);
//h_post_1hit_mu.at(m)=t1.at(0);
h_pre_2hit_mu.at(m)=t2.at(0);
h_post_2hit_mu.at(m)=t3.at(0);
}




   const Int_t NBINS_mu = 14;
   Double_t edges_mu[NBINS_mu + 1] = {0.00005,0.00015,0.00025,0.00035,0.00045,0.00055,0.00065,0.00075,0.00085,0.00095,0.0015,0.0025,0.0035,0.0045};

   const Int_t NBINS = 14;
   Double_t edges_el[NBINS + 1] = {0.0005,0.0015,0.0025,0.0035,0.0045,0.0055,0.0065,0.0075,0.0085,0.0095, 0.0125, 0.0175, 0.0225, 0.0285};


/*
auto h_sigma_mu_pre_1hit_clones = new TGraphErrors();
auto h_sigma_el_pre_1hit_clones = new TGraphErrors();
auto h_sigma_mu_post_1hit_clones = new TGraphErrors();
auto h_sigma_el_post_1hit_clones = new TGraphErrors();
*/auto h_sigma_mu_pre_2hit_clones = new TGraphErrors();
auto h_sigma_el_pre_2hit_clones = new TGraphErrors();
auto h_sigma_mu_post_2hit_clones = new TGraphErrors();
auto h_sigma_el_post_2hit_clones = new TGraphErrors();

//auto sigma_el_diff_1hit = new TGraphErrors();
//auto sigma_mu_diff_1hit = new TGraphErrors();
auto sigma_el_diff_2hit = new TGraphErrors();
auto sigma_mu_diff_2hit = new TGraphErrors();

/*double sigma_post_el1[NBINS]={0.};
double sigma_pre_el1[NBINS]={0.};
double sigma_post_mu1[NBINS_mu]={0.};
double sigma_pre_mu1[NBINS_mu]={0.};
*/double sigma_post_el2[NBINS]={0.};
double sigma_pre_el2[NBINS]={0.};
double sigma_post_mu2[NBINS_mu]={0.};
double sigma_pre_mu2[NBINS_mu]={0.};



/*TCanvas n("n","n",2800,2800);
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
n.SaveAs("/home/espedica/macros_fairmu/clean_codes/separate/validation/results/res_el_pre_1hit_%s.pdf");


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
n1.SaveAs("/home/espedica/macros_fairmu/clean_codes/separate/validation/results/res_el_post_1hit_%s.pdf");
*/

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
n2.SaveAs(Form("%sresults/res_el_pre_2hit_%s.pdf",path.c_str(),info.c_str()));


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

//n3.SaveAs("%sresults/res_el_pre_2hit_%s.pdf");





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
n4.SaveAs(Form("%sresults/res_el_post_2hit_%s.pdf",path.c_str(),info.c_str()));








/*
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
n5.SaveAs("%sresults/res_mu_pre_1hit_%s.pdf");


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
n6.SaveAs("%sresults/res_mu_post_1hit_%s.pdf");
*/


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
n7.SaveAs(Form("%sresults/res_mu_pre_2hit_%s.pdf",path.c_str(),info.c_str()));


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

//n2_mu.SaveAs("%sresults/res_pre_mu_2hit_%s.pdf");


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
n8.SaveAs(Form("%sresults/res_mu_post_2hit_%s.pdf",path.c_str(),info.c_str()));


h_sigma_mu_pre_2hit_clones->SetName("mu_pre");
h_sigma_el_pre_2hit_clones->SetName("el_pre");
h_sigma_mu_post_2hit_clones->SetName("mu_post");
h_sigma_el_post_2hit_clones->SetName("el_post");


h_sigma_el_pre_2hit_clones->SaveAs(Form("%sresults/sigma_el_pre_2hit_%s.root",path.c_str(),info.c_str()));
h_sigma_el_post_2hit_clones->SaveAs(Form("%sresults/sigma_el_post_2hit_%s.root",path.c_str(),info.c_str()));
h_sigma_mu_pre_2hit_clones->SaveAs(Form("%sresults/sigma_mu_pre_2hit_%s.root",path.c_str(),info.c_str()));
h_sigma_mu_post_2hit_clones->SaveAs(Form("%sresults/sigma_mu_post_2hit_%s.root",path.c_str(),info.c_str()));



for(int m=0; m<14; m++){
//sigma_el_diff_1hit->SetPoint(m,edges_el[m]*1000,(sigma_post_el1[m]-sigma_pre_el1[m])/sigma_pre_el1[m]);
//if(m!=0)sigma_mu_diff_1hit->SetPoint(m,edges_mu[m]*1000,(sigma_post_mu1[m]-sigma_pre_mu1[m])/sigma_pre_mu1[m]);
sigma_el_diff_2hit->SetPoint(m,edges_el[m]*1000,(sigma_post_el2[m]-sigma_pre_el2[m])/sigma_pre_el2[m]);
if(m!=0)sigma_mu_diff_2hit->SetPoint(m,edges_mu[m]*1000,(sigma_post_mu2[m]-sigma_pre_mu2[m])/sigma_pre_mu2[m]);
}



/*sigma_el_diff_1hit->SetName("el_1");
sigma_el_diff_1hit->GetYaxis()->SetTitle("(#sigma_post - #sigma_pre)/#sigma_pre");
sigma_el_diff_1hit->GetXaxis()->SetTitle("#theta_el (mrad)");
sigma_mu_diff_1hit->SetName("mu_1");
sigma_mu_diff_1hit->GetYaxis()->SetTitle("(#sigma_post - #sigma_pre)/#sigma_pre");
sigma_mu_diff_1hit->GetXaxis()->SetTitle("#theta_#mu (mrad)");
*/sigma_el_diff_2hit->SetName("el_2");
sigma_el_diff_2hit->GetYaxis()->SetTitle("(#sigma_post - #sigma_pre)/#sigma_pre");
sigma_el_diff_2hit->GetXaxis()->SetTitle("#theta_el (mrad)");
sigma_mu_diff_2hit->SetName("mu_2");
sigma_mu_diff_2hit->GetYaxis()->SetTitle("(#sigma_post - #sigma_pre)/#sigma_pre");
sigma_mu_diff_2hit->GetXaxis()->SetTitle("#theta_#mu (mrad)");


//sigma_el_diff_1hit->SaveAs(Form("%sresults/el_1.root"); 
//sigma_mu_diff_1hit->SaveAs(Form("%sresults/mu_1.root"); 
sigma_el_diff_2hit->SaveAs(Form("%sresults/el_2_%s.root",path.c_str(),info.c_str())); 
sigma_mu_diff_2hit->SaveAs(Form("%sresults/mu_2_%s.root",path.c_str(),info.c_str()));

}


