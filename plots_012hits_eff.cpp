
#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TStyle.h>


void plots_012hits_eff(string path,string info){

array<TFile*,6> f_thmu_LO_clone1,f_thmu_LO_clone2,f_thmu_LO_clone0,f_the_LO_clone1,f_the_LO_clone2,f_the_LO_clone0,f_op_LO_clone1,f_op_LO_clone2,f_op_LO_clone0,f_thmu_single_LO_clone1,f_thmu_single_LO_clone2,f_thmu_single_LO_clone0,f_the_single_LO_clone1,f_the_single_LO_clone2,f_the_single_LO_clone0,f_thmu_LO_gen1,f_thmu_LO_gen2,f_thmu_LO_gen0,f_the_LO_gen1,f_the_LO_gen2,f_the_LO_gen0,f_op_LO_gen1,f_op_LO_gen2,f_op_LO_gen0;

/*
for(int id=0; id<6;id++){
string range;
if(id==0)range = "0_5";
else if(id==1)range = "5_10";
else if(id==2)range = "10_15";
else if(id==3)range = "15_20";
else if(id==4)range = "20_25";
else if(id==5)range = "25_32";

f_thmu_LO_clone1.at(id) = TFile::Open((path+"mu_eff_LO_1hit_"+range+"_"+info+".root").c_str());
f_thmu_LO_clone2.at(id) = TFile::Open((path+"mu_eff_LO_2hit_"+range+"_"+info+".root").c_str());
f_thmu_LO_clone0.at(id) = TFile::Open((path+"mu_eff_LO_0hit_"+range+"_"+info+".root").c_str());

f_the_LO_clone1.at(id) = TFile::Open((path+"el_eff_LO_1hit_"+range+"_"+info+".root").c_str());
f_the_LO_clone2.at(id) = TFile::Open((path+"el_eff_LO_2hit_"+range+"_"+info+".root").c_str());
f_the_LO_clone0.at(id) = TFile::Open((path+"el_eff_LO_0hit_"+range+"_"+info+".root").c_str());

f_op_LO_clone1.at(id) = TFile::Open((path+"op_eff_LO_1hit_"+range+"_"+info+".root").c_str());
f_op_LO_clone2.at(id) = TFile::Open((path+"op_eff_LO_2hit_"+range+"_"+info+".root").c_str());
f_op_LO_clone0.at(id) = TFile::Open((path+"op_eff_LO_0hit_"+range+"_"+info+".root").c_str());

f_thmu_single_LO_clone1.at(id) = TFile::Open((path+"mu_single_eff_LO_1hit_"+range+"_"+info+".root").c_str());
f_thmu_single_LO_clone2.at(id) = TFile::Open((path+"mu_single_eff_LO_2hit_"+range+"_"+info+".root").c_str());
f_thmu_single_LO_clone0.at(id) = TFile::Open((path+"mu_single_eff_LO_0hit_"+range+"_"+info+".root").c_str());

f_the_single_LO_clone1.at(id) = TFile::Open((path+"el_single_eff_LO_1hit_"+range+"_"+info+".root").c_str());
f_the_single_LO_clone2.at(id) = TFile::Open((path+"el_single_eff_LO_2hit_"+range+"_"+info+".root").c_str());
f_the_single_LO_clone0.at(id) = TFile::Open((path+"el_single_eff_LO_0hit_"+range+"_"+info+".root").c_str());
}
*/


for(int id=0; id<6;id++){

string range;
if(id==0)range = "0_5";
else if(id==1)range = "5_10";
else if(id==2)range = "10_15";
else if(id==3)range = "15_20";
else if(id==4)range = "20_25";
else if(id==5)range = "25_32";

f_thmu_LO_clone1.at(id) = TFile::Open((path+"theta_mu_clone_1hit_"+range+"_"+info+".root").c_str());
f_thmu_LO_clone2.at(id) = TFile::Open((path+"theta_mu_clone_2hit_"+range+"_"+info+".root").c_str());
f_thmu_LO_clone0.at(id) = TFile::Open((path+"theta_mu_clone_0hit_"+range+"_"+info+".root").c_str());

f_the_LO_clone1.at(id) = TFile::Open((path+"theta_e_clone_1hit_"+range+"_"+info+".root").c_str());
f_the_LO_clone2.at(id) = TFile::Open((path+"theta_e_clone_2hit_"+range+"_"+info+".root").c_str());
f_the_LO_clone0.at(id) = TFile::Open((path+"theta_e_clone_0hit_"+range+"_"+info+".root").c_str());

f_op_LO_clone1.at(id) = TFile::Open((path+"h_opening_clone_1hit_"+range+"_"+info+".root").c_str());
f_op_LO_clone2.at(id) = TFile::Open((path+"h_opening_clone_2hit_"+range+"_"+info+".root").c_str());
f_op_LO_clone0.at(id) = TFile::Open((path+"h_opening_clone_0hit_"+range+"_"+info+".root").c_str());


f_thmu_LO_gen1.at(id) = TFile::Open((path+"theta_mu_gen_1hit_"+range+"_"+info+".root").c_str());
f_thmu_LO_gen2.at(id) = TFile::Open((path+"theta_mu_gen_2hit_"+range+"_"+info+".root").c_str());
f_thmu_LO_gen0.at(id) = TFile::Open((path+"theta_mu_gen_0hit_"+range+"_"+info+".root").c_str());

f_the_LO_gen1.at(id) = TFile::Open((path+"theta_e_gen_1hit_"+range+"_"+info+".root").c_str());
f_the_LO_gen2.at(id) = TFile::Open((path+"theta_e_gen_2hit_"+range+"_"+info+".root").c_str());
f_the_LO_gen0.at(id) = TFile::Open((path+"theta_e_gen_0hit_"+range+"_"+info+".root").c_str());

f_op_LO_gen1.at(id) = TFile::Open((path+"h_opening_gen_1hit_"+range+"_"+info+".root").c_str());
f_op_LO_gen2.at(id) = TFile::Open((path+"h_opening_gen_2hit_"+range+"_"+info+".root").c_str());
f_op_LO_gen0.at(id) = TFile::Open((path+"h_opening_gen_0hit_"+range+"_"+info+".root").c_str());


f_thmu_single_LO_clone1.at(id) = TFile::Open((path+"theta_mu_single_clone_1hit_"+range+"_"+info+".root").c_str());
f_thmu_single_LO_clone2.at(id) = TFile::Open((path+"theta_mu_single_clone_2hit_"+range+"_"+info+".root").c_str());
f_thmu_single_LO_clone0.at(id) = TFile::Open((path+"theta_mu_single_clone_0hit_"+range+"_"+info+".root").c_str());

f_the_single_LO_clone1.at(id) = TFile::Open((path+"theta_e_single_clone_1hit_"+range+"_"+info+".root").c_str());
f_the_single_LO_clone2.at(id) = TFile::Open((path+"theta_e_single_clone_2hit_"+range+"_"+info+".root").c_str());
f_the_single_LO_clone0.at(id) = TFile::Open((path+"theta_e_single_clone_0hit_"+range+"_"+info+".root").c_str());
}


TH1::SetDefaultSumw2(kTRUE);

//TH1D* h_theta_mu_LO_clone0;TH1D* h_theta_mu_LO_clone1;TH1D* h_theta_mu_LO_clone2;TH1D* h_theta_e_LO_clone0;TH1D* h_theta_e_LO_clone1;TH1D* h_theta_e_LO_clone2;TH1D* h_op_LO_clone0;TH1D* h_op_LO_clone1;TH1D* h_op_LO_clone2;TH1D* h_theta_mu_single_LO_clone0;TH1D* h_theta_mu_single_LO_clone1;TH1D* h_theta_mu_single_LO_clone2;TH1D* h_theta_e_single_LO_clone0;TH1D* h_theta_e_single_LO_clone1;TH1D* h_theta_e_single_LO_clone2;
array <TH1D*,6> tmp_h_theta_mu_LO_clone0,tmp_h_theta_mu_LO_clone1,tmp_h_theta_mu_LO_clone2,tmp_h_theta_e_LO_clone0,tmp_h_theta_e_LO_clone1,tmp_h_theta_e_LO_clone2,tmp_h_op_LO_clone0,tmp_h_op_LO_clone1,tmp_h_op_LO_clone2,tmp_h_theta_mu_single_LO_clone0,tmp_h_theta_mu_single_LO_clone1,tmp_h_theta_mu_single_LO_clone2,tmp_h_theta_e_single_LO_clone0,tmp_h_theta_e_single_LO_clone1,tmp_h_theta_e_single_LO_clone2,tmp_h_theta_mu_LO_gen0,tmp_h_theta_mu_LO_gen1,tmp_h_theta_mu_LO_gen2,tmp_h_theta_e_LO_gen0,tmp_h_theta_e_LO_gen1,tmp_h_theta_e_LO_gen2,tmp_h_op_LO_gen0,tmp_h_op_LO_gen1,tmp_h_op_LO_gen2;

 for(int id=0; id<6; id++){

tmp_h_theta_mu_LO_clone0.at(id)=(TH1D*)f_thmu_LO_clone0.at(id)->Get("theta_mu_clone");
tmp_h_theta_mu_LO_clone1.at(id)=(TH1D*)f_thmu_LO_clone1.at(id)->Get("theta_mu_clone");
tmp_h_theta_mu_LO_clone2.at(id)=(TH1D*)f_thmu_LO_clone2.at(id)->Get("theta_mu_clone");

tmp_h_theta_e_LO_clone0.at(id)=(TH1D*)f_the_LO_clone0.at(id)->Get("theta_e_clone");
tmp_h_theta_e_LO_clone1.at(id)=(TH1D*)f_the_LO_clone1.at(id)->Get("theta_e_clone");
tmp_h_theta_e_LO_clone2.at(id)=(TH1D*)f_the_LO_clone2.at(id)->Get("theta_e_clone");

tmp_h_op_LO_clone0.at(id)=(TH1D*)f_op_LO_clone0.at(id)->Get("h_opening_clone");
tmp_h_op_LO_clone1.at(id)=(TH1D*)f_op_LO_clone1.at(id)->Get("h_opening_clone");
tmp_h_op_LO_clone2.at(id)=(TH1D*)f_op_LO_clone2.at(id)->Get("h_opening_clone");

tmp_h_theta_mu_LO_gen0.at(id)=(TH1D*)f_thmu_LO_gen0.at(id)->Get("theta_mu_gen");
tmp_h_theta_mu_LO_gen1.at(id)=(TH1D*)f_thmu_LO_gen1.at(id)->Get("theta_mu_gen");
tmp_h_theta_mu_LO_gen2.at(id)=(TH1D*)f_thmu_LO_gen2.at(id)->Get("theta_mu_gen");

tmp_h_theta_e_LO_gen0.at(id)=(TH1D*)f_the_LO_gen0.at(id)->Get("theta_e_gen");
tmp_h_theta_e_LO_gen1.at(id)=(TH1D*)f_the_LO_gen1.at(id)->Get("theta_e_gen");
tmp_h_theta_e_LO_gen2.at(id)=(TH1D*)f_the_LO_gen2.at(id)->Get("theta_e_gen");

tmp_h_op_LO_gen0.at(id)=(TH1D*)f_op_LO_gen0.at(id)->Get("h_opening_gen");
tmp_h_op_LO_gen1.at(id)=(TH1D*)f_op_LO_gen1.at(id)->Get("h_opening_gen");
tmp_h_op_LO_gen2.at(id)=(TH1D*)f_op_LO_gen2.at(id)->Get("h_opening_gen");


tmp_h_theta_mu_single_LO_clone0.at(id)=(TH1D*)f_thmu_single_LO_clone0.at(id)->Get("theta_mu_single_clone");
tmp_h_theta_mu_single_LO_clone1.at(id)=(TH1D*)f_thmu_single_LO_clone1.at(id)->Get("theta_mu_single_clone");
tmp_h_theta_mu_single_LO_clone2.at(id)=(TH1D*)f_thmu_single_LO_clone2.at(id)->Get("theta_mu_single_clone");

tmp_h_theta_e_single_LO_clone0.at(id)=(TH1D*)f_the_single_LO_clone0.at(id)->Get("theta_e_single_clone");
tmp_h_theta_e_single_LO_clone1.at(id)=(TH1D*)f_the_single_LO_clone1.at(id)->Get("theta_e_single_clone");
tmp_h_theta_e_single_LO_clone2.at(id)=(TH1D*)f_the_single_LO_clone2.at(id)->Get("theta_e_single_clone");

        if(id!=0){
tmp_h_theta_mu_LO_clone0.at(0)->Add(tmp_h_theta_mu_LO_clone0.at(id));
tmp_h_theta_mu_LO_clone1.at(0)->Add(tmp_h_theta_mu_LO_clone1.at(id));
tmp_h_theta_mu_LO_clone2.at(0)->Add(tmp_h_theta_mu_LO_clone2.at(id));
tmp_h_theta_e_LO_clone0.at(0)->Add(tmp_h_theta_e_LO_clone0.at(id));
tmp_h_theta_e_LO_clone1.at(0)->Add(tmp_h_theta_e_LO_clone1.at(id));
tmp_h_theta_e_LO_clone2.at(0)->Add(tmp_h_theta_e_LO_clone2.at(id));
tmp_h_op_LO_clone0.at(0)->Add(tmp_h_op_LO_clone0.at(id));
tmp_h_op_LO_clone1.at(0)->Add(tmp_h_op_LO_clone1.at(id));
tmp_h_op_LO_clone2.at(0)->Add(tmp_h_op_LO_clone2.at(id));

tmp_h_theta_mu_LO_gen0.at(0)->Add(tmp_h_theta_mu_LO_gen0.at(id));
tmp_h_theta_mu_LO_gen1.at(0)->Add(tmp_h_theta_mu_LO_gen1.at(id));
tmp_h_theta_mu_LO_gen2.at(0)->Add(tmp_h_theta_mu_LO_gen2.at(id));
tmp_h_theta_e_LO_gen0.at(0)->Add(tmp_h_theta_e_LO_gen0.at(id));
tmp_h_theta_e_LO_gen1.at(0)->Add(tmp_h_theta_e_LO_gen1.at(id));
tmp_h_theta_e_LO_gen2.at(0)->Add(tmp_h_theta_e_LO_gen2.at(id));
tmp_h_op_LO_gen0.at(0)->Add(tmp_h_op_LO_gen0.at(id));
tmp_h_op_LO_gen1.at(0)->Add(tmp_h_op_LO_gen1.at(id));
tmp_h_op_LO_gen2.at(0)->Add(tmp_h_op_LO_gen2.at(id));

tmp_h_theta_mu_single_LO_clone0.at(0)->Add(tmp_h_theta_mu_single_LO_clone0.at(id));
tmp_h_theta_mu_single_LO_clone1.at(0)->Add(tmp_h_theta_mu_single_LO_clone1.at(id));
tmp_h_theta_mu_single_LO_clone2.at(0)->Add(tmp_h_theta_mu_single_LO_clone2.at(id));
tmp_h_theta_e_single_LO_clone0.at(0)->Add(tmp_h_theta_e_single_LO_clone0.at(id));
tmp_h_theta_e_single_LO_clone1.at(0)->Add(tmp_h_theta_e_single_LO_clone1.at(id));
tmp_h_theta_e_single_LO_clone2.at(0)->Add(tmp_h_theta_e_single_LO_clone2.at(id));
	}

}

TH1D * h0gen_0 = (TH1D*) tmp_h_theta_e_LO_gen0.at(0)->Clone();
TH1D * h1gen_0 = (TH1D*) tmp_h_theta_mu_LO_gen0.at(0)->Clone();
TH1D * h4gen_0 = (TH1D*) tmp_h_op_LO_gen0.at(0)->Clone();

TH1D * h_theta_e_single_LO_clone0 = (TH1D*) tmp_h_theta_e_single_LO_clone0.at(0)->Clone();
h_theta_e_single_LO_clone0->Divide(h_theta_e_single_LO_clone0,h0gen_0,1,1,"B");

TH1D * h_theta_mu_single_LO_clone0 = (TH1D*) tmp_h_theta_mu_single_LO_clone0.at(0)->Clone();
h_theta_mu_single_LO_clone0->Divide(h_theta_mu_single_LO_clone0,h1gen_0,1,1,"B");

TH1D * h_theta_e_LO_clone0 = (TH1D*) tmp_h_theta_e_LO_clone0.at(0)->Clone();
h_theta_e_LO_clone0->Divide(h_theta_e_LO_clone0,h0gen_0,1,1,"B");

TH1D * h_theta_mu_LO_clone0 = (TH1D*) tmp_h_theta_mu_LO_clone0.at(0)->Clone();
h_theta_mu_LO_clone0->Divide(h_theta_mu_LO_clone0,h1gen_0,1,1,"B");

TH1D * h_op_LO_clone0 = (TH1D*) tmp_h_op_LO_clone0.at(0)->Clone();
h_op_LO_clone0->Divide(h_op_LO_clone0,h4gen_0,1,1,"B");




TH1D * h0gen_1 = (TH1D*) tmp_h_theta_e_LO_gen1.at(0)->Clone();
TH1D * h1gen_1 = (TH1D*) tmp_h_theta_mu_LO_gen1.at(0)->Clone();
TH1D * h4gen_1 = (TH1D*) tmp_h_op_LO_gen1.at(0)->Clone();

TH1D * h_theta_e_single_LO_clone1 = (TH1D*) tmp_h_theta_e_single_LO_clone1.at(0)->Clone();
h_theta_e_single_LO_clone1->Divide(h_theta_e_single_LO_clone1,h0gen_1,1,1,"B");

TH1D * h_theta_mu_single_LO_clone1 = (TH1D*) tmp_h_theta_mu_single_LO_clone1.at(0)->Clone();
h_theta_mu_single_LO_clone1->Divide(h_theta_mu_single_LO_clone1,h1gen_1,1,1,"B");

TH1D * h_theta_e_LO_clone1 = (TH1D*) tmp_h_theta_e_LO_clone1.at(0)->Clone();
h_theta_e_LO_clone1->Divide(h_theta_e_LO_clone1,h0gen_0,1,1,"B");

TH1D * h_theta_mu_LO_clone1 = (TH1D*) tmp_h_theta_mu_LO_clone1.at(0)->Clone();
h_theta_mu_LO_clone1->Divide(h_theta_mu_LO_clone1,h1gen_1,1,1,"B");

TH1D * h_op_LO_clone1 = (TH1D*) tmp_h_op_LO_clone1.at(0)->Clone();
h_op_LO_clone1->Divide(h_op_LO_clone1,h4gen_1,1,1,"B");



TH1D * h0gen_2 = (TH1D*) tmp_h_theta_e_LO_gen2.at(0)->Clone();
TH1D * h1gen_2 = (TH1D*) tmp_h_theta_mu_LO_gen2.at(0)->Clone();
TH1D * h4gen_2 = (TH1D*) tmp_h_op_LO_gen2.at(0)->Clone();

TH1D * h_theta_e_single_LO_clone2 = (TH1D*) tmp_h_theta_e_single_LO_clone2.at(0)->Clone();
h_theta_e_single_LO_clone2->Divide(h_theta_e_single_LO_clone2,h0gen_2,1,1,"B");

TH1D * h_theta_mu_single_LO_clone2 = (TH1D*) tmp_h_theta_mu_single_LO_clone2.at(0)->Clone();
h_theta_mu_single_LO_clone2->Divide(h_theta_mu_single_LO_clone2,h1gen_2,1,1,"B");

TH1D * h_theta_e_LO_clone2 = (TH1D*) tmp_h_theta_e_LO_clone2.at(0)->Clone();
h_theta_e_LO_clone2->Divide(h_theta_e_LO_clone2,h0gen_2,1,1,"B");

TH1D * h_theta_mu_LO_clone2 = (TH1D*) tmp_h_theta_mu_LO_clone2.at(0)->Clone();
h_theta_mu_LO_clone2->Divide(h_theta_mu_LO_clone2,h1gen_2,1,1,"B");

TH1D * h_op_LO_clone2 = (TH1D*) tmp_h_op_LO_clone2.at(0)->Clone();
h_op_LO_clone2->Divide(h_op_LO_clone2,h4gen_2,1,1,"B");




 auto legend_e = new TLegend(0.6,0.1,0.9,0.25);
   legend_e->AddEntry(h_theta_e_LO_clone0,"0 hit shared","LEP");
   legend_e->AddEntry(h_theta_e_LO_clone1,"1 hit shared","LEP");
   legend_e->AddEntry(h_theta_e_LO_clone2,"2 hit shared","LEP");

 auto legend_mu = new TLegend(0.1,0.1,0.4,0.25);
   legend_mu->AddEntry(h_theta_mu_LO_clone0,"0 hit shared","LEP");
   legend_mu->AddEntry(h_theta_mu_LO_clone1,"1 hit shared","LEP");
   legend_mu->AddEntry(h_theta_mu_LO_clone2,"2 hit shared","LEP");

 auto legend_op = new TLegend(0.6,0.1,0.9,0.25);
   legend_op->AddEntry(h_op_LO_clone0,"0 hit shared","LEP");
   legend_op->AddEntry(h_op_LO_clone1,"1 hit shared","LEP");
   legend_op->AddEntry(h_op_LO_clone2,"2 hit shared","LEP");

 auto legend_e_single = new TLegend();
   legend_e_single->AddEntry(h_theta_e_single_LO_clone0,"0 hit shared","LEP");
   legend_e_single->AddEntry(h_theta_e_single_LO_clone1,"1 hit shared","LEP");
   legend_e_single->AddEntry(h_theta_e_single_LO_clone2,"2 hit shared","LEP");

 auto legend_mu_single = new TLegend();
   legend_mu_single->AddEntry(h_theta_mu_single_LO_clone0,"0 hit shared","LEP");
   legend_mu_single->AddEntry(h_theta_mu_single_LO_clone1,"1 hit shared","LEP");
   legend_mu_single->AddEntry(h_theta_mu_single_LO_clone2,"2 hit shared","LEP");

TCanvas a("a","a",1800,700);
a.Divide(2,1);
a.cd(1);
h_theta_mu_single_LO_clone2->SetLineColor(kAzure+7);
h_theta_mu_single_LO_clone2->SetTitle("");//Reconstruction efficiency of single mu as a function of muon scattering angle");
h_theta_mu_single_LO_clone2->GetXaxis()->SetTitle("Muon angle[rad]");
h_theta_mu_single_LO_clone2->GetYaxis()->SetTitle("Efficiency");
h_theta_mu_single_LO_clone2->SetMinimum(0.7);
TGaxis::SetMaxDigits(3);
h_theta_mu_single_LO_clone2->Draw("E");
h_theta_mu_single_LO_clone1->SetLineColor(kGreen+1);
h_theta_mu_single_LO_clone1->Draw("E same");
h_theta_mu_single_LO_clone0->SetLineColor(kOrange+10);
h_theta_mu_single_LO_clone0->Draw("E same");
legend_mu_single->Draw();
gStyle->SetOptStat(0);
a.cd(2);
h_theta_e_single_LO_clone2->SetTitle("");//Reconstruction efficiency of single electron as a function of e- scattering angle");
h_theta_e_single_LO_clone2->GetXaxis()->SetTitle("Electron angle[rad]");
h_theta_e_single_LO_clone2->GetYaxis()->SetTitle("Efficiency");
h_theta_e_single_LO_clone2->SetLineColor(kAzure+7);
h_theta_e_single_LO_clone2->SetMinimum(0.7);
h_theta_e_single_LO_clone2->Draw("E");
h_theta_e_single_LO_clone1->SetLineColor(kGreen+1);
h_theta_e_single_LO_clone1->Draw("E same");
h_theta_e_single_LO_clone0->SetLineColor(kOrange+10);
h_theta_e_single_LO_clone0->Draw("E same");
legend_e_single->Draw();

gStyle->SetOptStat(0);

a.SaveAs(Form("%s/results/eff_single_LO_012hit_%s.pdf",path.c_str(),info.c_str()));

//TCanvas a1("a1","a1",1400,1400);
TCanvas a1("a1","a1",2400,700);
a1.Divide(3,1);
a1.cd(1);
h_theta_mu_LO_clone2->SetTitle("");//Reconstruction efficiency of elastic event as a function of mu scattering angle");
h_theta_mu_LO_clone2->GetXaxis()->SetTitle("Muon angle[rad]");
h_theta_mu_LO_clone2->GetYaxis()->SetTitle("Efficiency");
h_theta_mu_LO_clone2->SetLineColor(kAzure+7);
h_theta_mu_LO_clone2->SetMinimum(0.6);
TGaxis::SetMaxDigits(3);
h_theta_mu_LO_clone2->Draw("E");
h_theta_mu_LO_clone1->SetLineColor(kGreen+1);
h_theta_mu_LO_clone1->Draw("E same");
h_theta_mu_LO_clone0->SetLineColor(kOrange+10);
h_theta_mu_LO_clone0->Draw("E same");
legend_mu->Draw();
gStyle->SetOptStat(0);
a1.cd(2);
h_theta_e_LO_clone2->SetTitle("");//Reconstruction efficiency of elastic event as a function of e- scattering angle");
h_theta_e_LO_clone2->GetXaxis()->SetTitle("Electron angle[rad]");
h_theta_e_LO_clone2->GetYaxis()->SetTitle("Efficiency");
h_theta_e_LO_clone2->SetLineColor(kAzure+7);
h_theta_e_LO_clone2->SetMinimum(0.6);
h_theta_e_LO_clone2->Draw("E");
h_theta_e_LO_clone1->SetLineColor(kGreen+1);
h_theta_e_LO_clone1->Draw("E same");
h_theta_e_LO_clone0->SetLineColor(kOrange+10);
h_theta_e_LO_clone0->Draw("E same");
legend_e->Draw();

gStyle->SetOptStat(0);
a1.cd(3);
h_op_LO_clone2->SetTitle("");//Reconstruction efficiency of elastic event as a function of e-mu opening angle");
h_op_LO_clone2->GetXaxis()->SetTitle("Opening angle[rad]");
h_op_LO_clone2->GetYaxis()->SetTitle("Efficiency");
h_op_LO_clone2->GetXaxis()->SetRangeUser(0.,0.032);
h_op_LO_clone2->SetLineColor(kAzure+7);
h_op_LO_clone2->SetMinimum(0.);
h_op_LO_clone2->Draw("E");
h_op_LO_clone1->SetLineColor(kGreen+1);
h_op_LO_clone1->Draw("E same");
h_op_LO_clone0->SetLineColor(kOrange+10);
h_op_LO_clone0->Draw("E same");
legend_op->Draw();
gStyle->SetOptStat(0);
a1.SaveAs(Form("%s/results/eff_event_LO_012hit_%s.pdf",path.c_str(),info.c_str()));



TFile *o_file=new TFile(Form("%s/results/eff_%s.root",path.c_str(),info.c_str()),"RECREATE");

h_op_LO_clone0->SetName("h_op_LO_clone0");
h_op_LO_clone1->SetName("h_op_LO_clone1");
h_op_LO_clone2->SetName("h_op_LO_clone2");

h_op_LO_clone0->Write();
h_op_LO_clone1->Write();
h_op_LO_clone2->Write();

h_theta_mu_single_LO_clone0->SetName("h_theta_mu_single_LO_clone0");
h_theta_mu_single_LO_clone1->SetName("h_theta_mu_single_LO_clone1");
h_theta_mu_single_LO_clone2->SetName("h_theta_mu_single_LO_clone2");

h_theta_mu_single_LO_clone0->Write();
h_theta_mu_single_LO_clone1->Write();
h_theta_mu_single_LO_clone2->Write();

h_theta_e_single_LO_clone0->SetName("h_theta_e_single_LO_clone0");
h_theta_e_single_LO_clone1->SetName("h_theta_e_single_LO_clone1");
h_theta_e_single_LO_clone2->SetName("h_theta_e_single_LO_clone2");

h_theta_e_single_LO_clone0->Write();
h_theta_e_single_LO_clone1->Write();
h_theta_e_single_LO_clone2->Write();

}


