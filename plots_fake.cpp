
#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TStyle.h>


void plots_fake(string path,string info){

array<TFile*,6> f_good_el0,f_good_mu0,f_fake0,f_good_el1,f_good_mu1,f_fake1,f_good_el2,f_good_mu2,f_fake2;
array<TFile*,6> f_reco_el0,f_reco_mu0,f_reco0,f_reco_el1,f_reco_mu1,f_reco1,f_reco_el2,f_reco_mu2,f_reco2;

for(int id=0; id<6;id++){

string range;
if(id==0)range = "0_5";
else if(id==1)range = "5_10";
else if(id==2)range = "10_15";
else if(id==3)range = "15_20";
else if(id==4)range = "20_25";
else if(id==5)range = "25_32";

f_reco_el0.at(id)=TFile::Open((path+"reco_el_0hit_LO_"+range+"_"+info+".root").c_str());
f_reco_mu0.at(id)=TFile::Open((path+"reco_mu_0hit_LO_"+range+"_"+info+".root").c_str());
f_reco0.at(id)=TFile::Open((path+"reco_0hit_LO_"+range+"_"+info+".root").c_str());

f_reco_el1.at(id)=TFile::Open((path+"reco_el_1hit_LO_"+range+"_"+info+".root").c_str());
f_reco_mu1.at(id)=TFile::Open((path+"reco_mu_1hit_LO_"+range+"_"+info+".root").c_str());
f_reco1.at(id)=TFile::Open((path+"reco_1hit_LO_"+range+"_"+info+".root").c_str());

f_reco_el2.at(id)=TFile::Open((path+"reco_el_2hit_LO_"+range+"_"+info+".root").c_str());
f_reco_mu2.at(id)=TFile::Open((path+"reco_mu_2hit_LO_"+range+"_"+info+".root").c_str());
f_reco2.at(id)=TFile::Open((path+"reco_2hit_LO_"+range+"_"+info+".root").c_str());


f_good_el0.at(id)=TFile::Open((path+"good_el_0hit_LO_"+range+"_"+info+".root").c_str());
f_good_mu0.at(id)=TFile::Open((path+"good_mu_0hit_LO_"+range+"_"+info+".root").c_str());
f_fake0.at(id)=TFile::Open((path+"fake_0hit_LO_"+range+"_"+info+".root").c_str());

f_good_el1.at(id)=TFile::Open((path+"good_el_1hit_LO_"+range+"_"+info+".root").c_str());
f_good_mu1.at(id)=TFile::Open((path+"good_mu_1hit_LO_"+range+"_"+info+".root").c_str());
f_fake1.at(id)=TFile::Open((path+"fake_1hit_LO_"+range+"_"+info+".root").c_str());

f_good_el2.at(id)=TFile::Open((path+"good_el_2hit_LO_"+range+"_"+info+".root").c_str());
f_good_mu2.at(id)=TFile::Open((path+"good_mu_2hit_LO_"+range+"_"+info+".root").c_str());
f_fake2.at(id)=TFile::Open((path+"fake_2hit_LO_"+range+"_"+info+".root").c_str());
}

TH1::SetDefaultSumw2(kTRUE);

array <TH1D*,6> tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8;
array <TH1D*,6> r_tmp0,r_tmp1,r_tmp2,r_tmp3,r_tmp4,r_tmp5,r_tmp6,r_tmp7,r_tmp8;


 for(int id=0; id<6; id++){

tmp0.at(id)=(TH1D*)f_good_el0.at(id)->Get("theta_e_good");
tmp1.at(id)=(TH1D*)f_good_mu0.at(id)->Get("theta_mu_good");
tmp2.at(id)=(TH1D*)f_fake0.at(id)->Get("theta_fake");
tmp3.at(id)=(TH1D*)f_good_el1.at(id)->Get("theta_e_good");
tmp4.at(id)=(TH1D*)f_good_mu1.at(id)->Get("theta_mu_good");
tmp5.at(id)=(TH1D*)f_fake1.at(id)->Get("theta_fake");
tmp6.at(id)=(TH1D*)f_good_el2.at(id)->Get("theta_e_good");
tmp7.at(id)=(TH1D*)f_good_mu2.at(id)->Get("theta_mu_good");
tmp8.at(id)=(TH1D*)f_fake2.at(id)->Get("theta_fake");

r_tmp0.at(id)=(TH1D*)f_reco_el0.at(id)->Get("theta_e_reco");
r_tmp1.at(id)=(TH1D*)f_reco_mu0.at(id)->Get("theta_mu_reco");
r_tmp2.at(id)=(TH1D*)f_reco0.at(id)->Get("theta_reco");
r_tmp3.at(id)=(TH1D*)f_reco_el1.at(id)->Get("theta_e_reco");
r_tmp4.at(id)=(TH1D*)f_reco_mu1.at(id)->Get("theta_mu_reco");
r_tmp5.at(id)=(TH1D*)f_reco1.at(id)->Get("theta_reco");
r_tmp6.at(id)=(TH1D*)f_reco_el2.at(id)->Get("theta_e_reco");
r_tmp7.at(id)=(TH1D*)f_reco_mu2.at(id)->Get("theta_mu_reco");
r_tmp8.at(id)=(TH1D*)f_reco2.at(id)->Get("theta_reco");

        if(id!=0){
	tmp0.at(0)->Add(tmp0.at(id));
        tmp1.at(0)->Add(tmp1.at(id));
        tmp2.at(0)->Add(tmp2.at(id));
        tmp3.at(0)->Add(tmp3.at(id));
        tmp4.at(0)->Add(tmp4.at(id));
        tmp5.at(0)->Add(tmp5.at(id));
        tmp6.at(0)->Add(tmp6.at(id));
        tmp7.at(0)->Add(tmp7.at(id));
        tmp8.at(0)->Add(tmp8.at(id));

        r_tmp0.at(0)->Add(r_tmp0.at(id));
        r_tmp1.at(0)->Add(r_tmp1.at(id));
        r_tmp2.at(0)->Add(r_tmp2.at(id));
        r_tmp3.at(0)->Add(r_tmp3.at(id));
        r_tmp4.at(0)->Add(r_tmp4.at(id));
        r_tmp5.at(0)->Add(r_tmp5.at(id));
        r_tmp6.at(0)->Add(r_tmp6.at(id));
        r_tmp7.at(0)->Add(r_tmp7.at(id));
        r_tmp8.at(0)->Add(r_tmp8.at(id));


	}
}

TH1D * h0_0 = (TH1D*) r_tmp0.at(0)->Clone();
TH1D * h1_0 = (TH1D*) r_tmp1.at(0)->Clone();
TH1D * h2_0 = (TH1D*) r_tmp2.at(0)->Clone();

TH1D *h_good_el0=(TH1D*)tmp0.at(0)->Clone();
h_good_el0->Divide(h_good_el0,h0_0,1,1,"B");

TH1D *h_good_mu0=(TH1D*)tmp1.at(0)->Clone();
h_good_mu0->Divide(h_good_mu0,h1_0,1,1,"B");

TH1D *h_fake0=(TH1D*)tmp2.at(0)->Clone();
h_fake0->Divide(h_fake0,h2_0,1,1,"B");


TH1D * h3_1 = (TH1D*) r_tmp3.at(0)->Clone();
TH1D * h4_1 = (TH1D*) r_tmp4.at(0)->Clone();
TH1D * h5_1 = (TH1D*) r_tmp5.at(0)->Clone();

TH1D *h_good_el1=(TH1D*)tmp3.at(0)->Clone();
h_good_el1->Divide(h_good_el1,h3_1,1,1,"B");

TH1D *h_good_mu1=(TH1D*)tmp4.at(0)->Clone();
h_good_mu1->Divide(h_good_mu1,h4_1,1,1,"B");

TH1D *h_fake1=(TH1D*)tmp5.at(0)->Clone();
h_fake1->Divide(h_fake1,h5_1,1,1,"B");


TH1D * h6_2 = (TH1D*) r_tmp6.at(0)->Clone();
TH1D * h7_2 = (TH1D*) r_tmp7.at(0)->Clone();
TH1D * h8_2 = (TH1D*) r_tmp8.at(0)->Clone();

TH1D *h_good_el2=(TH1D*)tmp6.at(0)->Clone();
h_good_el2->Divide(h_good_el2,h6_2,1,1,"B");

TH1D *h_good_mu2=(TH1D*)tmp7.at(0)->Clone();
h_good_mu2->Divide(h_good_mu2,h7_2,1,1,"B");

TH1D *h_fake2=(TH1D*)tmp8.at(0)->Clone();
h_fake2->Divide(h_fake2,h8_2,1,1,"B");



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
h_fake2->SetTitle("");//Fake tracks rate as a function of particle scattering angle");
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
a1.SaveAs(Form("%sresults/good_fake_rate_LO_012hit_%s.pdf",path.c_str(),info.c_str()));






}


