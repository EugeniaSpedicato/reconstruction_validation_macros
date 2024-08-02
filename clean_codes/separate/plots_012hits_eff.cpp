
#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TStyle.h>


void plots_012hits_eff(){

array<TFile*,6> f_thmu_LO_clone1,f_thmu_LO_clone2,f_thmu_LO_clone0,f_the_LO_clone1,f_the_LO_clone2,f_the_LO_clone0,f_op_LO_clone1,f_op_LO_clone2,f_op_LO_clone0,f_thmu_single_LO_clone1,f_thmu_single_LO_clone2,f_thmu_single_LO_clone0,f_the_single_LO_clone1,f_the_single_LO_clone2,f_the_single_LO_clone0;

 for(int id=0; id<6; id++){

string range;
if(id==0)range = "0_5";
else if(id==1)range = "5_10";
else if(id==2)range = "10_15";
else if(id==3)range = "15_20";
else if(id==4)range = "20_25";
else if(id==5)range = "25_32";

f_thmu_LO_clone1.at(id) = TFile::Open(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/mu_eff_LO_1hit_%d_first2.root",static_cast<char>(id)));
f_thmu_LO_clone2.at(id) = TFile::Open(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/mu_eff_LO_2hit_%d_first2.root",static_cast<char>(id)));
f_thmu_LO_clone0.at(id) = TFile::Open(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/mu_eff_LO_0hit_%d_first2.root",static_cast<char>(id)));

f_the_LO_clone1.at(id) = TFile::Open(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/el_eff_LO_1hit_%d_first2.root",static_cast<char>(id)));
f_the_LO_clone2.at(id) = TFile::Open(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/el_eff_LO_2hit_%d_first2.root",static_cast<char>(id)));
f_the_LO_clone0.at(id) = TFile::Open(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/el_eff_LO_0hit_%d_first2.root",static_cast<char>(id)));

f_op_LO_clone1.at(id) = TFile::Open(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/op_eff_LO_1hit_%d_first2.root",static_cast<char>(id)));
f_op_LO_clone2.at(id) = TFile::Open(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/op_eff_LO_2hit_%d_first2.root",static_cast<char>(id)));
f_op_LO_clone0.at(id) = TFile::Open(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/op_eff_LO_0hit_%d_first2.root",static_cast<char>(id)));

f_thmu_single_LO_clone1.at(id) = TFile::Open(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/mu_single_eff_LO_1hit_%d_first2.root",static_cast<char>(id)));
f_thmu_single_LO_clone2.at(id) = TFile::Open(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/mu_single_eff_LO_2hit_%d_first2.root",static_cast<char>(id)));
f_thmu_single_LO_clone0.at(id) = TFile::Open(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/mu_single_eff_LO_0hit_%d_first2.root",static_cast<char>(id)));

f_the_single_LO_clone1.at(id) = TFile::Open(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/el_single_eff_LO_1hit_%d_first2.root",static_cast<char>(id)));
f_the_single_LO_clone2.at(id) = TFile::Open(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/el_single_eff_LO_2hit_%d_first2.root",static_cast<char>(id)));
f_the_single_LO_clone0.at(id) = TFile::Open(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/el_single_eff_LO_0hit_%d_first2.root",static_cast<char>(id)));
}

TH1::SetDefaultSumw2(kTRUE);

TH1D* h_theta_mu_LO_clone0;TH1D* h_theta_mu_LO_clone1;TH1D* h_theta_mu_LO_clone2;TH1D* h_theta_e_LO_clone0;TH1D* h_theta_e_LO_clone1;TH1D* h_theta_e_LO_clone2;TH1D* h_op_LO_clone0;TH1D* h_op_LO_clone1;TH1D* h_op_LO_clone2;TH1D* h_theta_mu_single_LO_clone0;TH1D* h_theta_mu_single_LO_clone1;TH1D* h_theta_mu_single_LO_clone2;TH1D* h_theta_e_single_LO_clone0;TH1D* h_theta_e_single_LO_clone1;TH1D* h_theta_e_single_LO_clone2;
array <TH1D*,6> tmp_h_theta_mu_LO_clone0,tmp_h_theta_mu_LO_clone1,tmp_h_theta_mu_LO_clone2,tmp_h_theta_e_LO_clone0,tmp_h_theta_e_LO_clone1,tmp_h_theta_e_LO_clone2,tmp_h_op_LO_clone0,tmp_h_op_LO_clone1,tmp_h_op_LO_clone2,tmp_h_theta_mu_single_LO_clone0,tmp_h_theta_mu_single_LO_clone1,tmp_h_theta_mu_single_LO_clone2,tmp_h_theta_e_single_LO_clone0,tmp_h_theta_e_single_LO_clone1,tmp_h_theta_e_single_LO_clone2;

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
tmp_h_theta_mu_single_LO_clone0.at(0)->Add(tmp_h_theta_mu_single_LO_clone0.at(id));
tmp_h_theta_mu_single_LO_clone1.at(0)->Add(tmp_h_theta_mu_single_LO_clone1.at(id));
tmp_h_theta_mu_single_LO_clone2.at(0)->Add(tmp_h_theta_mu_single_LO_clone2.at(id));
tmp_h_theta_e_single_LO_clone0.at(0)->Add(tmp_h_theta_e_single_LO_clone0.at(id));
tmp_h_theta_e_single_LO_clone1.at(0)->Add(tmp_h_theta_e_single_LO_clone1.at(id));
tmp_h_theta_e_single_LO_clone2.at(0)->Add(tmp_h_theta_e_single_LO_clone2.at(id));
	}

}


h_theta_mu_LO_clone0=tmp_h_theta_mu_LO_clone0.at(0);
h_theta_mu_LO_clone1=tmp_h_theta_mu_LO_clone1.at(0);
h_theta_mu_LO_clone2=tmp_h_theta_mu_LO_clone2.at(0);
h_theta_e_LO_clone0=tmp_h_theta_e_LO_clone0.at(0);
h_theta_e_LO_clone1=tmp_h_theta_e_LO_clone1.at(0);
h_theta_e_LO_clone2=tmp_h_theta_e_LO_clone2.at(0);
h_op_LO_clone0=tmp_h_op_LO_clone0.at(0);
h_op_LO_clone1=tmp_h_op_LO_clone1.at(0);
h_op_LO_clone2=tmp_h_op_LO_clone2.at(0);
h_theta_mu_single_LO_clone0=tmp_h_theta_mu_single_LO_clone0.at(0);
h_theta_mu_single_LO_clone1=tmp_h_theta_mu_single_LO_clone1.at(0);
h_theta_mu_single_LO_clone2=tmp_h_theta_mu_single_LO_clone2.at(0);
h_theta_e_single_LO_clone0=tmp_h_theta_e_single_LO_clone0.at(0);
h_theta_e_single_LO_clone1=tmp_h_theta_e_single_LO_clone1.at(0);
h_theta_e_single_LO_clone2=tmp_h_theta_e_single_LO_clone2.at(0);


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

a.SaveAs("/home/espedica/macros_fairmu/clean_codes/separate/validation/results/eff_single_LO_012hit.pdf");

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
h_op_LO_clone2->SetMinimum(0.8);
h_op_LO_clone2->Draw("E");
h_op_LO_clone1->SetLineColor(kGreen+1);
h_op_LO_clone1->Draw("E same");
h_op_LO_clone0->SetLineColor(kOrange+10);
h_op_LO_clone0->Draw("E same");
legend_op->Draw();
gStyle->SetOptStat(0);
a1.SaveAs("/home/espedica/macros_fairmu/clean_codes/separate/validation/results/eff_event_LO_012hit.pdf");


}


