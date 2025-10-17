#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TStyle.h>
#include "TGraphErrors.h"


//qui normalizzo al numero di eventi nei dati

void same_comparison(int nhits_MC2,int nhits_mc, string rd, string mc, string type_d, string type_m, string bend_MC2, string bend_mc){

TH1::SetDefaultSumw2(kTRUE);
int const NBINS = 6;



TFile *f_MC2=TFile::Open(Form("comparison_RDMC/%s_d_aco_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",rd.c_str(),bend_MC2.c_str(),nhits_MC2,type_d.c_str()));
TFile *f_thmu_MC2=TFile::Open(Form("comparison_RDMC/%s_theta_mu_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",rd.c_str(),bend_MC2.c_str(),nhits_MC2,type_d.c_str()));
TFile *f_the_MC2=TFile::Open(Form("comparison_RDMC/%s_theta_e_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",rd.c_str(),bend_MC2.c_str(),nhits_MC2,type_d.c_str()));
TFile *f_2D_MC2=TFile::Open(Form("comparison_RDMC/%s_2D_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",rd.c_str(),bend_MC2.c_str(),nhits_MC2,type_d.c_str()));
TFile *f_op_MC2=TFile::Open(Form("comparison_RDMC/%s_opening_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",rd.c_str(),bend_MC2.c_str(),nhits_MC2,type_d.c_str()));

TFile* f_the_MC=TFile::Open(Form("comparison_RDMC/%s_theta_e_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",mc.c_str(),bend_mc.c_str(),nhits_mc,type_m.c_str()));
TFile* f_thmu_MC=TFile::Open(Form("comparison_RDMC/%s_theta_mu_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",mc.c_str(),bend_mc.c_str(),nhits_mc,type_m.c_str()));
TFile* f_2D_MC=TFile::Open(Form("comparison_RDMC/%s_2D_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",mc.c_str(),bend_mc.c_str(),nhits_mc,type_m.c_str()));
TFile* f_op_MC=TFile::Open(Form("comparison_RDMC/%s_opening_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",mc.c_str(),bend_mc.c_str(),nhits_mc,type_m.c_str()));
TFile* f_aco_MC=TFile::Open(Form("comparison_RDMC/%s_d_aco_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",mc.c_str(),bend_mc.c_str(),nhits_mc,type_m.c_str())); 


TFile *f_vrtx_chi2_MC2=TFile::Open(Form("comparison_RDMC/%s_vrtx_chi2_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_%s_%dhit_%s.root",rd.c_str(),bend_MC2.c_str(),nhits_MC2,type_d.c_str()));
TFile *f_track_el_MC2=TFile::Open(Form("comparison_RDMC/%s_track_el_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_%s_%dhit_%s.root",rd.c_str(),bend_MC2.c_str(),nhits_MC2,type_d.c_str()));
TFile *f_track_mu_MC2=TFile::Open(Form("comparison_RDMC/%s_track_mu_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_%s_%dhit_%s.root",rd.c_str(),bend_MC2.c_str(),nhits_MC2,type_d.c_str()));

TFile *f_vrtx_chi2_MC=TFile::Open(Form("comparison_RDMC/%s_vrtx_chi2_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_%s_%dhit_%s.root",mc.c_str(),bend_mc.c_str(),nhits_mc,type_m.c_str()));
TFile *f_track_el_MC=TFile::Open(Form("comparison_RDMC/%s_track_el_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_%s_%dhit_%s.root",mc.c_str(),bend_mc.c_str(),nhits_mc,type_m.c_str()));
TFile *f_track_mu_MC=TFile::Open(Form("comparison_RDMC/%s_track_mu_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_%s_%dhit_%s.root",mc.c_str(),bend_mc.c_str(),nhits_mc,type_m.c_str()));

TFile * f_posZ_MC2=TFile::Open(Form("comparison_RDMC/%s_z_post_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",rd.c_str(),bend_MC2.c_str(),nhits_MC2,type_d.c_str()));
TFile * f_posZ_MC=TFile::Open(Form("comparison_RDMC/%s_z_post_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",mc.c_str(),bend_mc.c_str(),nhits_mc,type_m.c_str()));

TFile * f_ns1_MC=TFile::Open(Form("comparison_RDMC/%s_nstubs_s1_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",mc.c_str(),bend_mc.c_str(),nhits_mc,type_m.c_str()));
TFile * f_ns1_MC2=TFile::Open(Form("comparison_RDMC/%s_nstubs_s1_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",rd.c_str(),bend_MC2.c_str(),nhits_MC2,type_d.c_str()));

TFile* f_aco_pre_MC=TFile::Open(Form("comparison_RDMC/%s_preCuts_d_aco_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",mc.c_str(),bend_mc.c_str(),nhits_mc,type_m.c_str())); 
TFile* f_aco_pre_MC2=TFile::Open(Form("comparison_RDMC/%s_preCuts_d_aco_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",rd.c_str(),bend_MC2.c_str(),nhits_MC2,type_d.c_str())); 

TFile * f_ns1_pre_MC=TFile::Open(Form("comparison_RDMC/%s_nstubs_pre_s1_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",mc.c_str(),bend_mc.c_str(),nhits_mc,type_m.c_str()));
TFile * f_ns1_pre_MC2=TFile::Open(Form("comparison_RDMC/%s_nstubs_pre_s1_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",rd.c_str(),bend_MC2.c_str(),nhits_MC2,type_d.c_str()));

TFile * f_posZ_pre_MC2=TFile::Open(Form("comparison_RDMC/%s_z_pre_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",rd.c_str(),bend_MC2.c_str(),nhits_MC2,type_d.c_str()));
TFile * f_posZ_pre_MC=TFile::Open(Form("comparison_RDMC/%s_z_pre_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",mc.c_str(),bend_mc.c_str(),nhits_mc,type_m.c_str()));





TH1D* h_posZ_MC2=(TH1D*)f_posZ_MC2->Get("h_z_pos");
TH1D* h_posZ_MC=(TH1D*)f_posZ_MC->Get("h_z_pos");

TH1D* h_theta_e_MC=(TH1D*)f_the_MC->Get("theta_e");
TH1D* h_theta_mu_MC=(TH1D*)f_thmu_MC->Get("theta_mu");
TH2D* h_2D_MC=(TH2D*)f_2D_MC->Get("h2D");
TH1D* h_opening_MC=(TH1D*)f_op_MC->Get("h_opening");
TH1D* h_aco_MC=(TH1D*)f_aco_MC->Get("d_aco_real");

TH1D *h_theta_e_MC2=(TH1D*)f_the_MC2->Get("theta_e");
TH1D *h_theta_mu_MC2=(TH1D*)f_thmu_MC2->Get("theta_mu");
TH2D *h_2D_MC2=(TH2D*)f_2D_MC2->Get("h2D");
TH1D *h_opening_MC2=(TH1D*)f_op_MC2->Get("h_opening");
TH1D *h_aco_MC2=(TH1D*)f_MC2->Get("d_aco_real");



TH1D *h_vrtx_chi2_MC=(TH1D*)f_vrtx_chi2_MC->Get("vrtx_chi2");
TH1D *h_track_el_MC=(TH1D*)f_track_el_MC->Get("track_el");
TH1D *h_track_mu_MC=(TH1D*)f_track_mu_MC->Get("track_mu");


TH1D *h_vrtx_chi2_MC2=(TH1D*)f_vrtx_chi2_MC2->Get("vrtx_chi2");
TH1D *h_track_el_MC2=(TH1D*)f_track_el_MC2->Get("track_el");
TH1D *h_track_mu_MC2=(TH1D*)f_track_mu_MC2->Get("track_mu");


TH1D *h_ns1_MC2=(TH1D*)f_ns1_MC2->Get("h_nstubs_s1");
TH1D *h_ns1_MC=(TH1D*)f_ns1_MC->Get("h_nstubs_s1");

TH1D *h_aco_pre_MC=(TH1D*)f_aco_pre_MC->Get("d_aco_real_pre");
TH1D *h_aco_pre_MC2=(TH1D*)f_aco_pre_MC2->Get("d_aco_real_pre");
TH1D *h_ns1_pre_MC=(TH1D*)f_ns1_pre_MC->Get("h_nstubs_pre_s1");
TH1D *h_ns1_pre_MC2=(TH1D*)f_ns1_pre_MC2->Get("h_nstubs_pre_s1");
TH1D *h_posZ_pre_MC2=(TH1D*)f_posZ_pre_MC2->Get("h_z_pos_pre");
TH1D *h_posZ_pre_MC=(TH1D*)f_posZ_pre_MC->Get("h_z_pos_pre");


cout << "before scaling entries MC: " << h_opening_MC->Integral() << endl;
cout << "before scaling entries RD: " << h_opening_MC2->Integral() << endl;

//	 h_theta_e_MC->Scale(1./h_theta_e_MC->Integral());//9.3792854e+08,7.9729075e+08
//         h_theta_mu_MC->Scale(1./h_theta_mu_MC->Integral());
//	 h_opening_MC->Scale(1./h_opening_MC->Integral());


//         h_theta_e_MC2->Scale(1./h_theta_e_MC2->Integral());//9.3792854e+08,7.9729075e+08
//         h_theta_mu_MC2->Scale(1./h_theta_mu_MC2->Integral());
//         h_opening_MC2->Scale(1./h_opening_MC2->Integral());


 auto l = new TLegend(0.75,0.15,0.9,0.3);
l->AddEntry(h_theta_e_MC,Form("%s",bend_mc.c_str()),"LEP");;
l->AddEntry(h_theta_e_MC2,Form("%s",bend_MC2.c_str()),"LEP");

 auto l1 = new TLegend(0.55,0.15,0.7,0.3);
l1->AddEntry(h_theta_mu_MC,"MC","LEP");;
l1->AddEntry(h_theta_mu_MC2,"data","LEP");

 auto l2 = new TLegend(0.75,0.15,0.9,0.3);
l2->AddEntry(h_opening_MC,"MC","LEP");;
l2->AddEntry(h_opening_MC2,"data","LEP");




TCanvas a("a","a",700,700);
a.Divide(1,2);
a.cd(1);
h_opening_MC->SetMinimum(1.);
h_opening_MC->GetXaxis()->SetRangeUser(0.005,0.02);
h_opening_MC2->GetXaxis()->SetRangeUser(0.005,0.02);
h_opening_MC2->SetMinimum(1.);
h_opening_MC->Draw("hist");
h_opening_MC2->SetLineColor(kPink+10);
h_opening_MC2->Draw("hist same");
cout << "after scaling h_opening_MC entries: " << h_opening_MC->Integral() << endl;
cout << "after scaling h_opening_MC2 entries: " << h_opening_MC2->Integral() << endl;
a.cd(2);

TH1D * h = (TH1D*) h_opening_MC2->Clone();
h->Divide(h_opening_MC);
//h->SetMaximum(3.);
h->SetMinimum(0.9);
h->SetMaximum(1.15);
h->GetXaxis()->SetTitle("Opening angle [rad]");
h->Draw("E");
gStyle->SetOptStat(0);
h->SaveAs(Form("comparison_RDMC/%s_%s_ratio_MCMC2normalization_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_opening_%dhit_%dhit.root",rd.c_str(),mc.c_str(),nhits_MC2,nhits_mc));
a.SaveAs(Form("comparison_RDMC/%s_%s_opening_reb_MCMC2normalization_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%dhit.pdf",rd.c_str(),mc.c_str(),nhits_MC2,nhits_mc));


TCanvas a1("a1","a1",1000,700);
a1.Divide(2,2);
a1.cd(1);
/*
h_theta_e_MC->GetXaxis()->SetRangeUser(0.003,0.01);
h_theta_e_MC2->GetXaxis()->SetRangeUser(0.003,0.01);
*/
h_theta_e_MC->GetXaxis()->SetRangeUser(0.005,0.02);
h_theta_e_MC2->GetXaxis()->SetRangeUser(0.005,0.02);
h_theta_e_MC->GetYaxis()->SetTitle("Entries");
h_theta_e_MC->SetTitle("Electron scattering angle");
h_theta_e_MC->GetXaxis()->SetTitle("Electron angle [rad]");
h_theta_e_MC->SetMinimum(100.);
h_theta_e_MC->Draw("E");
h_theta_e_MC2->SetLineColor(kPink+10);
h_theta_e_MC2->Draw("E same");
l->Draw();
gPad->SetLogy();
gStyle->SetOptStat(0);
cout << "after scaling h_theta_e_MC entries: " << h_theta_e_MC->Integral() << endl;
cout << "after scaling h_theta_e_MC2 entries: " << h_theta_e_MC2->Integral() << endl;

a1.cd(2);

TH1D * h3 = (TH1D*) h_theta_e_MC2->Clone();
h3->Divide(h_theta_e_MC);
h3->SetMinimum(0.9);
h3->SetMaximum(1.15);
h3->SetTitle("Electron scattering angle");
h3->GetYaxis()->SetTitle("Data/MC ratio");
h3->GetXaxis()->SetTitle("Electron angle [rad]");
TLine *line_ = new TLine(0.00,1.,0.02,1.);
TLine *line_p = new TLine(0.00,1.03,0.02,1.03);
TLine *line_m = new TLine(0.00,0.97,0.02,0.97);

TLine *line_p05 = new TLine(0.00,1.005,0.02,1.005);
TLine *line_m05 = new TLine(0.00,0.995,0.02,0.995);

h3->Draw("E");
line_->SetLineStyle(7);
line_p->SetLineStyle(9);
line_m->SetLineStyle(9);
line_->SetLineColor(39);
line_p->SetLineColor(40);
line_m->SetLineColor(40);
line_->Draw("same");
line_p->Draw("same");
line_m->Draw("same");
line_p05->SetLineColor(45);
line_m05->SetLineColor(45);
//line_p05->Draw("same");
//line_m05->Draw("same");

h3->SaveAs(Form("comparison_RDMC/%s_%s_ratio_MCMC2normalization_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_theta_el_%dhit_%dhit.root",rd.c_str(),mc.c_str(),nhits_MC2,nhits_mc));
a1.cd(3);
//h_theta_e_MC->Scale(7.9729075e+08*5.5*1E+23*d_tar*1E-30*315.44638/6.79333e+07);//26563993,605675    137720.00    1.3242886e+09   6.79333e+07 9.3792854e+08
h_theta_mu_MC->SetMinimum(100.);
/*h_theta_mu_MC->GetXaxis()->SetRangeUser(0.0005,0.002);
h_theta_mu_MC2->GetXaxis()->SetRangeUser(0.0005,0.002);*/
h_theta_mu_MC->GetXaxis()->SetRangeUser(0.0003,0.0013);
h_theta_mu_MC2->GetXaxis()->SetRangeUser(0.0003,0.0013);
h_theta_mu_MC->GetYaxis()->SetTitle("Entries");
h_theta_mu_MC->SetTitle("Muon scattering angle");
h_theta_mu_MC->GetXaxis()->SetTitle("Muon angle [rad]");
h_theta_mu_MC->Draw("E");
h_theta_mu_MC2->SetLineColor(kPink+10);
h_theta_mu_MC2->Draw("E same");
l1->Draw();
gStyle->SetOptStat(0);
gPad->SetLogy();
cout << "after scaling h_theta_mu_MC entries: " << h_theta_mu_MC->Integral() << endl;
cout << "after scaling h_theta_mu_MC2 entries: " << h_theta_mu_MC2->Integral() << endl;

a1.cd(4);
TH1D * h4 = (TH1D*) h_theta_mu_MC2->Clone();
h4->Divide(h_theta_mu_MC);
h4->SetTitle("Muon scattering angle");
h4->SetMinimum(0.9);
h4->SetMaximum(1.15);
h4->GetYaxis()->SetTitle("Data/MC ratio");
h4->GetXaxis()->SetTitle("Muon angle [rad]");
TLine *mu_line_ = new TLine(0.0003,1.,0.0013,1.);
TLine *mu_line_p = new TLine(0.0003,1.03,0.0013,1.03);
TLine *mu_line_m = new TLine(0.0003,0.97,0.0013,0.97);
h4->Draw("E");
TLine *mu_line_p05 = new TLine(0.00,1.005,0.02,1.005);
TLine *mu_line_m05 = new TLine(0.00,0.995,0.02,0.995);
mu_line_p05->SetLineColor(45);
mu_line_m05->SetLineColor(45);
//mu_line_p05->Draw("same");
//mu_line_m05->Draw("same");
mu_line_->SetLineStyle(7);
mu_line_p->SetLineStyle(9);
mu_line_m->SetLineStyle(9);
mu_line_->SetLineColor(39);
mu_line_p->SetLineColor(40);
mu_line_m->SetLineColor(40);
mu_line_->Draw();
mu_line_p->Draw("same");
mu_line_m->Draw("same");
h4->SaveAs(Form("comparison_RDMC/%s_%s_ratio_MCMC2normalization_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_theta_mu_%dhit_%dhit.root",rd.c_str(),mc.c_str(),nhits_MC2,nhits_mc));
a1.SaveAs(Form("comparison_RDMC/%s_%s_theta_MC_reco_MC2_add_MCMC2normalization_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%dhit.pdf",rd.c_str(),mc.c_str(),nhits_MC2,nhits_mc));

TF1 *Elastic = new TF1("Elastic","0.5109989461*0.001*((1+(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*cos(x)*cos(x))/(1-(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*cos(x)*cos(x)))",0,0.030); 
TF1 *Elastic2 = new TF1("Elastic2","asin( (sin(x)*sqrt(Elastic(x)*Elastic(x)-0.5109989461*0.001*0.5109989461*0.001))/sqrt( (160+0.5109989461*0.001-Elastic(x))*(160+0.5109989461*0.001-Elastic(x))-105.6583745 *0.001*105.6583745 *0.001 ) )",0,0.030);




TCanvas t("t","t",3000,3000);
t.Divide(3,3);
t.cd(1);
h3->SetLineColor(kRed+1);
h3->Rebin(1);
h3->GetXaxis()->SetRangeUser(0.005,0.02);
h3->Draw();
line_->Draw("same");
line_p->Draw("same");
line_m->Draw("same");
t.cd(2);
h4->SetLineColor(kOrange-3);
h4->Rebin(1);
h4->GetXaxis()->SetRangeUser(0.0002,0.0016);
h4->Draw();
mu_line_->Draw("same");
mu_line_p->Draw("same");
mu_line_m->Draw("same");
t.cd(3);
h->SetLineColor(kGreen-2);
h->Rebin(1);
h->GetXaxis()->SetRangeUser(0.005,0.02);
h->Draw();
line_->Draw("same");
line_p->Draw("same");
line_m->Draw("same");

t.cd(4);
h_theta_e_MC->SetMinimum(100.);
h_theta_e_MC->GetXaxis()->SetTitle("Electron angle [rad]");
h_theta_e_MC->Rebin(1);
h_theta_e_MC2->Rebin(1);
h_theta_e_MC->GetXaxis()->SetRangeUser(0.005,0.02);
h_theta_e_MC2->GetXaxis()->SetRangeUser(0.005,0.02);
h_theta_e_MC->Draw("E");
h_theta_e_MC2->SetLineColor(kRed+1);
h_theta_e_MC2->Draw("E same");
gStyle->SetOptStat(0);
//gPad->SetLogy();
l->Draw();
t.cd(5);
h_theta_mu_MC->SetMinimum(100.);
h_theta_mu_MC->GetXaxis()->SetTitle("Muon angle [rad]");
h_theta_mu_MC->Rebin(1);
h_theta_mu_MC2->Rebin(1);
h_theta_mu_MC->GetXaxis()->SetRangeUser(0.0002,0.0016);
h_theta_mu_MC2->GetXaxis()->SetRangeUser(0.0002,0.0016);
h_theta_mu_MC->Draw("E");
h_theta_mu_MC2->SetLineColor(kOrange-3);
h_theta_mu_MC2->Draw("E same");
gStyle->SetOptStat(0);
//gPad->SetLogy();
l1->Draw();
t.cd(6);
h_opening_MC->SetMinimum(100.);
h_opening_MC->GetXaxis()->SetTitle("Opening angle [rad]");
h_opening_MC->Rebin(1);
h_opening_MC2->Rebin(1);
h_opening_MC->GetXaxis()->SetRangeUser(0.005,0.02);
h_opening_MC2->GetXaxis()->SetRangeUser(0.005,0.02);
h_opening_MC->Draw("E");
h_opening_MC2->SetLineColor(kGreen-2);
h_opening_MC2->Draw("E same");
gStyle->SetOptStat(0);
//gPad->SetLogy();
l2->Draw();
t.cd(7);
h_theta_e_MC->Draw("E");
h_theta_e_MC2->Draw("E same");
h_theta_e_MC->SetMinimum(100.);
gStyle->SetOptStat(0);
gPad->SetLogy();
l->Draw();
t.cd(8);
h_theta_mu_MC->Draw("E");
h_theta_mu_MC2->Draw("E same");
h_theta_mu_MC->SetMinimum(100.);
gStyle->SetOptStat(0);
gPad->SetLogy();
l1->Draw();
t.cd(9);
h_opening_MC->Draw("E");
h_opening_MC2->Draw("E same");
h_opening_MC->SetMinimum(100.);
gStyle->SetOptStat(0);
gPad->SetLogy();
l2->Draw();
t.SaveAs(Form("comparison_RDMC/%s_%s_tesi_comparison_MCMC2normalization_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_%dhit_%dhit.pdf",rd.c_str(),mc.c_str(),nhits_MC2,nhits_mc));


h_aco_MC->Scale(1./h_aco_MC->Integral());
h_vrtx_chi2_MC->Scale(1./h_vrtx_chi2_MC->Integral());
h_posZ_MC->Scale(1./h_posZ_MC->Integral());

h_aco_MC2->Scale(1./h_aco_MC2->Integral());
h_vrtx_chi2_MC2->Scale(1./h_vrtx_chi2_MC2->Integral());
h_posZ_MC2->Scale(1./h_posZ_MC2->Integral());

h_ns1_MC->Scale(1./h_ns1_MC->Integral());
h_ns1_MC2->Scale(1./h_ns1_MC2->Integral());


TCanvas v("v","v",2000,2000);
v.Divide(2,2);
v.cd(1);
h_aco_MC->GetXaxis()->SetRangeUser(-0.4,0.4);
h_aco_MC->Draw("E hist");
h_aco_MC2->SetLineColor(kPink);
h_aco_MC2->Draw("E hist same");
v.cd(2);
h_vrtx_chi2_MC->Draw("E hist");
h_vrtx_chi2_MC2->SetLineColor(kPink);
h_vrtx_chi2_MC2->Draw("E hist same");
gPad->SetLogy();
v.cd(3);
h_posZ_MC->GetXaxis()->SetRangeUser(800.,920);
h_posZ_MC2->SetLineColor(kPink);
h_posZ_MC->Draw("E hist");
h_posZ_MC2->Draw("E hist same");
v.cd(4);
h_ns1_MC->Draw("E hist");
h_ns1_MC2->SetLineColor(kPink);
h_ns1_MC2->Draw("E hist same");
v.SaveAs(Form("comparison_RDMC/%s_%s_variables_comparison_MCMC2normalization_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_%dhit_%dhit.pdf",rd.c_str(),mc.c_str(),nhits_MC2,nhits_mc));


h_aco_pre_MC->Scale(h_aco_pre_MC2->Integral()/h_aco_pre_MC->Integral());
//h_ns1_pre_MC->Scale(h_ns1_pre_MC2->Integral()/h_ns1_pre_MC->Integral());
h_posZ_pre_MC->Scale(h_posZ_pre_MC2->Integral()/h_posZ_pre_MC->Integral());

TCanvas v_pre("v_pre","v_pre",2000,2000);
v_pre.Divide(2,3);
v_pre.cd(1);
h_aco_pre_MC->SetTitle("Acoplanarity - Preselection");
h_aco_pre_MC->GetXaxis()->SetRangeUser(-3.2,3.2);
h_aco_pre_MC->GetYaxis()->SetTitle("Entries");
h_aco_pre_MC->GetXaxis()->SetTitle("Acoplanarity [rad]");
h_aco_pre_MC->Draw("E hist");
h_aco_pre_MC2->SetLineColor(kPink);
h_aco_pre_MC2->Draw("E hist same");
v_pre.cd(2);
h_aco_MC->SetTitle("Acoplanarity - Elastic selection");
h_aco_MC->GetYaxis()->SetTitle("Entries");
h_aco_MC->GetXaxis()->SetTitle("Acoplanarity [rad]");
h_aco_MC->GetXaxis()->SetRangeUser(-0.4,0.4);
h_aco_MC->Draw("E hist");
h_aco_MC2->SetLineColor(kPink);
h_aco_MC2->Draw("E hist same");
v_pre.cd(3);
h_posZ_pre_MC->SetTitle("Vertex Z - Preselection");
h_posZ_pre_MC->GetYaxis()->SetTitle("Entries");
h_posZ_pre_MC->GetXaxis()->SetTitle("Vertex Z [cm]");
h_posZ_pre_MC->GetXaxis()->SetRangeUser(800.,920);
h_posZ_pre_MC2->SetLineColor(kPink);
h_posZ_pre_MC->Draw("E hist");
h_posZ_pre_MC2->Draw("E hist same");
v_pre.cd(4);
h_posZ_MC->SetTitle("Vertex Z - Elastic selection");
h_posZ_MC->GetYaxis()->SetTitle("Entries");
h_posZ_MC->GetXaxis()->SetTitle("Vertex Z [cm]");
h_posZ_MC->GetXaxis()->SetRangeUser(800.,920);
h_posZ_MC2->SetLineColor(kPink);
h_posZ_MC->Draw("E hist");
h_posZ_MC2->Draw("E hist same");
v_pre.cd(5);
h_ns1_pre_MC->SetTitle("Number of stubs Station1 - Preselection");
h_ns1_pre_MC->GetYaxis()->SetTitle("Entries");
h_ns1_pre_MC->GetXaxis()->SetTitle("Number of stubs");
h_ns1_pre_MC2->SetLineColor(kPink);
h_ns1_pre_MC->Draw("E hist ");
h_ns1_pre_MC2->Draw("E hist same");
h_ns1_pre_MC->SetMinimum(1.);
gPad->SetLogy();
v_pre.cd(6);
h_ns1_MC->SetTitle("Number of stubs Station1 - Elastic selection");
h_ns1_MC->GetYaxis()->SetTitle("Entries");
h_ns1_MC->GetXaxis()->SetTitle("Number of stubs");
h_ns1_MC->Draw("E hist");
h_ns1_MC2->SetLineColor(kPink);
h_ns1_MC2->Draw("E hist  same");
h_ns1_MC->SetMinimum(1.);


gPad->SetLogy();
v_pre.SaveAs(Form("comparison_RDMC/%s_%s_variables_preCuts_comparison_MCMC2normalization_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_%dhit_%dhit.pdf",rd.c_str(),mc.c_str(),nhits_MC2,nhits_mc));


}

