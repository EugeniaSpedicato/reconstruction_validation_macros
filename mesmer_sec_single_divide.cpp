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

void mesmer_sec_single_divide(int nhits_rd,int nhits_mc, string rd, string mc, string type_d, string type_m){

TH1::SetDefaultSumw2(kTRUE);
int const NBINS = 6;



TFile *f_RD=TFile::Open(Form("comparison_RDMC/%s_d_aco_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%s.root",rd.c_str(),nhits_rd,type_d.c_str()));
TFile *f_thmu_RD=TFile::Open(Form("comparison_RDMC/%s_theta_mu_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%s.root",rd.c_str(),nhits_rd,type_d.c_str()));
TFile *f_the_RD=TFile::Open(Form("comparison_RDMC/%s_theta_e_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%s.root",rd.c_str(),nhits_rd,type_d.c_str()));
TFile *f_2D_RD=TFile::Open(Form("comparison_RDMC/%s_2D_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%s.root",rd.c_str(),nhits_rd,type_d.c_str()));
TFile *f_op_RD=TFile::Open(Form("comparison_RDMC/%s_opening_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%s.root",rd.c_str(),nhits_rd,type_d.c_str()));

TFile* f_the_MC=TFile::Open(Form("comparison_RDMC/%s_theta_e_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%s_smearing.root",mc.c_str(),nhits_mc,type_m.c_str()));
TFile* f_thmu_MC=TFile::Open(Form("comparison_RDMC/%s_theta_mu_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%s_smearing.root",mc.c_str(),nhits_mc,type_m.c_str()));
TFile* f_2D_MC=TFile::Open(Form("comparison_RDMC/%s_2D_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%s_smearing.root",mc.c_str(),nhits_mc,type_m.c_str()));
TFile* f_op_MC=TFile::Open(Form("comparison_RDMC/%s_opening_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%s_smearing.root",mc.c_str(),nhits_mc,type_m.c_str()));
TFile* f_aco_MC=TFile::Open(Form("comparison_RDMC/%s_d_aco_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%s_smearing.root",mc.c_str(),nhits_mc,type_m.c_str())); 


TFile *f_vrtx_chi2_RD=TFile::Open(Form("comparison_RDMC/%s_vrtx_chi2_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_Bend_%dhit_%s.root",rd.c_str(),nhits_rd,type_d.c_str()));
TFile *f_track_el_RD=TFile::Open(Form("comparison_RDMC/%s_track_el_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_Bend_%dhit_%s.root",rd.c_str(),nhits_rd,type_d.c_str()));
TFile *f_track_mu_RD=TFile::Open(Form("comparison_RDMC/%s_track_mu_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_Bend_%dhit_%s.root",rd.c_str(),nhits_rd,type_d.c_str()));

TFile *f_vrtx_chi2_MC=TFile::Open(Form("comparison_RDMC/%s_vrtx_chi2_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_Bend_%dhit_%s_smearing.root",mc.c_str(),nhits_mc,type_m.c_str()));
TFile *f_track_el_MC=TFile::Open(Form("comparison_RDMC/%s_track_el_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_Bend_%dhit_%s_smearing.root",mc.c_str(),nhits_mc,type_m.c_str()));
TFile *f_track_mu_MC=TFile::Open(Form("comparison_RDMC/%s_track_mu_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_Bend_%dhit_%s_smearing.root",mc.c_str(),nhits_mc,type_m.c_str()));

TFile * f_posZ_RD=TFile::Open(Form("comparison_RDMC/%s_z_post_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%s.root",rd.c_str(),nhits_rd,type_d.c_str()));
TFile * f_posZ_MC=TFile::Open(Form("comparison_RDMC/%s_z_post_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%s_smearing.root",mc.c_str(),nhits_mc,type_m.c_str()));

TFile * f_ns1_MC=TFile::Open(Form("comparison_RDMC/%s_nstubs_s1_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%s_smearing.root",mc.c_str(),nhits_mc,type_m.c_str()));
TFile * f_ns1_RD=TFile::Open(Form("comparison_RDMC/%s_nstubs_s1_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%s.root",rd.c_str(),nhits_rd,type_d.c_str()));

TFile* f_aco_pre_MC=TFile::Open(Form("comparison_RDMC/%s_preCuts_d_aco_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%s_smearing.root",mc.c_str(),nhits_mc,type_m.c_str())); 
TFile* f_aco_pre_RD=TFile::Open(Form("comparison_RDMC/%s_preCuts_d_aco_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%s.root",rd.c_str(),nhits_rd,type_d.c_str())); 

TFile * f_ns1_pre_MC=TFile::Open(Form("comparison_RDMC/%s_nstubs_pre_s1_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%s_smearing.root",mc.c_str(),nhits_mc,type_m.c_str()));
TFile * f_ns1_pre_RD=TFile::Open(Form("comparison_RDMC/%s_nstubs_pre_s1_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%s.root",rd.c_str(),nhits_rd,type_d.c_str()));

TFile * f_posZ_pre_RD=TFile::Open(Form("comparison_RDMC/%s_z_pre_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%s.root",rd.c_str(),nhits_rd,type_d.c_str()));
TFile * f_posZ_pre_MC=TFile::Open(Form("comparison_RDMC/%s_z_pre_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%s_smearing.root",mc.c_str(),nhits_mc,type_m.c_str()));

TFile * f_thetaL_RD=TFile::Open(Form("comparison_RDMC/%s_thetaLeft_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%s.root",rd.c_str(),nhits_rd,type_d.c_str()));
TFile * f_thetaL_MC=TFile::Open(Form("comparison_RDMC/%s_thetaLeft_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%s_smearing.root",mc.c_str(),nhits_mc,type_m.c_str()));
TFile * f_thetaR_RD=TFile::Open(Form("comparison_RDMC/%s_thetaRight_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%s.root",rd.c_str(),nhits_rd,type_d.c_str()));
TFile * f_thetaR_MC=TFile::Open(Form("comparison_RDMC/%s_thetaRight_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%s_smearing.root",mc.c_str(),nhits_mc,type_m.c_str()));



TH1D* h_posZ_RD=(TH1D*)f_posZ_RD->Get("h_z_pos");
TH1D* h_posZ_MC=(TH1D*)f_posZ_MC->Get("h_z_pos");

TH1D* h_theta_e_MC=(TH1D*)f_the_MC->Get("theta_e");
TH1D* h_theta_mu_MC=(TH1D*)f_thmu_MC->Get("theta_mu");
TH2D* h_2D_MC=(TH2D*)f_2D_MC->Get("h2D");
TH1D* h_opening_MC=(TH1D*)f_op_MC->Get("h_opening");
TH1D* h_aco_MC=(TH1D*)f_aco_MC->Get("d_aco_real");

TH1D *h_theta_e_RD=(TH1D*)f_the_RD->Get("theta_e");
TH1D *h_theta_mu_RD=(TH1D*)f_thmu_RD->Get("theta_mu");
TH2D *h_2D_RD=(TH2D*)f_2D_RD->Get("h2D");
TH1D *h_opening_RD=(TH1D*)f_op_RD->Get("h_opening");
TH1D *h_aco_RD=(TH1D*)f_RD->Get("d_aco_real");



TH1D *h_vrtx_chi2_MC=(TH1D*)f_vrtx_chi2_MC->Get("vrtx_chi2");
TH1D *h_track_el_MC=(TH1D*)f_track_el_MC->Get("track_el");
TH1D *h_track_mu_MC=(TH1D*)f_track_mu_MC->Get("track_mu");


TH1D *h_vrtx_chi2_RD=(TH1D*)f_vrtx_chi2_RD->Get("vrtx_chi2");
TH1D *h_track_el_RD=(TH1D*)f_track_el_RD->Get("track_el");
TH1D *h_track_mu_RD=(TH1D*)f_track_mu_RD->Get("track_mu");


TH1D *h_ns1_RD=(TH1D*)f_ns1_RD->Get("h_nstubs_s1");
TH1D *h_ns1_MC=(TH1D*)f_ns1_MC->Get("h_nstubs_s1");

TH1D *h_aco_pre_MC=(TH1D*)f_aco_pre_MC->Get("d_aco_real_pre");
TH1D *h_aco_pre_RD=(TH1D*)f_aco_pre_RD->Get("d_aco_real_pre");
TH1D *h_ns1_pre_MC=(TH1D*)f_ns1_pre_MC->Get("h_nstubs_pre_s1");
TH1D *h_ns1_pre_RD=(TH1D*)f_ns1_pre_RD->Get("h_nstubs_pre_s1");
TH1D *h_posZ_pre_RD=(TH1D*)f_posZ_pre_RD->Get("h_z_pos_pre");
TH1D *h_posZ_pre_MC=(TH1D*)f_posZ_pre_MC->Get("h_z_pos_pre");

TH1D *h_thetaL_RD=(TH1D*)f_thetaL_RD->Get("theta_L");
TH1D *h_thetaL_MC=(TH1D*)f_thetaL_MC->Get("theta_L");

TH1D *h_thetaR_RD=(TH1D*)f_thetaR_RD->Get("theta_R");
TH1D *h_thetaR_MC=(TH1D*)f_thetaR_MC->Get("theta_R");

cout << "before scaling entries MC: " << h_opening_MC->Integral() << endl;
cout << "before scaling entries RD: " << h_opening_RD->Integral() << endl;

if(rd=="my_modifica_golden_all" or rd=="run8") {h_theta_e_MC->Rebin(2); h_theta_mu_MC->Rebin(2);}
else{h_theta_e_MC->Rebin(2); h_theta_mu_MC->Rebin(2);
h_theta_e_RD->Rebin(2); h_theta_mu_RD->Rebin(2);}


	 h_theta_e_MC->Scale(h_theta_e_RD->Integral()/h_theta_e_MC->Integral());//9.3792854e+08,7.9729075e+08
         h_theta_mu_MC->Scale(h_theta_mu_RD->Integral()/h_theta_mu_MC->Integral());
	 h_opening_MC->Scale(h_opening_RD->Integral()/h_opening_MC->Integral());

	h_thetaL_MC->Scale(h_thetaL_RD->Integral()/h_thetaL_MC->Integral());
	h_thetaR_MC->Scale(h_thetaR_RD->Integral()/h_thetaR_MC->Integral());

/*	 h_theta_e_MC->Scale(1./h_theta_e_MC->Integral());//9.3792854e+08,7.9729075e+08
         h_theta_mu_MC->Scale(1./h_theta_mu_MC->Integral());
	 h_opening_MC->Scale(1./h_opening_MC->Integral());

	 h_theta_e_RD->Scale(1./h_theta_e_RD->Integral());//9.3792854e+08,7.9729075e+08
         h_theta_mu_RD->Scale(1./h_theta_mu_RD->Integral());
	 h_opening_RD->Scale(1./h_opening_RD->Integral());
*/

 auto l = new TLegend(0.75,0.15,0.9,0.3);
l->AddEntry(h_theta_e_MC,"MC","LEP");;
l->AddEntry(h_theta_e_RD,"data","LEP");

 auto l1 = new TLegend(0.55,0.15,0.7,0.3);
l1->AddEntry(h_theta_mu_MC,"MC","LEP");;
l1->AddEntry(h_theta_mu_RD,"data","LEP");

 auto l2 = new TLegend(0.75,0.15,0.9,0.3);
l2->AddEntry(h_opening_MC,"MC","LEP");;
l2->AddEntry(h_opening_RD,"data","LEP");



TCanvas a("a","a", 1400,700);
a.Divide(1,2);
a.cd(1);
h_opening_MC->SetMinimum(1.);
h_opening_MC->GetXaxis()->SetRangeUser(0.00,0.032);//0.005,0.02);
h_opening_RD->GetXaxis()->SetRangeUser(0.00,0.032);//0.005,0.02);
h_opening_RD->SetMinimum(1.);
h_opening_MC->Draw("hist");
h_opening_RD->SetLineColor(kPink+10);
h_opening_RD->Draw("hist same");
cout << "after scaling h_opening_MC entries: " << h_opening_MC->Integral() << endl;
cout << "after scaling h_opening_RD entries: " << h_opening_RD->Integral() << endl;
a.cd(2);

TH1D * h = (TH1D*) h_opening_RD->Clone();
h->Divide(h_opening_MC);
//h->SetMaximum(3.);
h->SetMinimum(0.9);
h->SetMaximum(1.15);
h->GetXaxis()->SetTitle("Opening angle [rad]");
h->Draw("E");
gStyle->SetOptStat(0);
h->SaveAs(Form("comparison_RDMC/%s_%s_ratio_RDnormalization_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_opening_%dhit_%dhit_smearing.root",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));
a.SaveAs(Form("comparison_RDMC/%s_%s_opening_reb_RDnormalization_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%dhit_smearing.pdf",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));


TCanvas a1("a1","a1",1000,700);
a1.Divide(2,2);
a1.cd(1);
/*
h_theta_e_MC->GetXaxis()->SetRangeUser(0.003,0.01);
h_theta_e_RD->GetXaxis()->SetRangeUser(0.003,0.01);
*/
h_theta_e_MC->GetXaxis()->SetRangeUser(0.00,0.032);//0.005,0.02);
h_theta_e_RD->GetXaxis()->SetRangeUser(0.00,0.032);//0.005,0.02);
h_theta_e_MC->GetYaxis()->SetTitle("Entries");
h_theta_e_MC->SetTitle("Electron scattering angle");
h_theta_e_MC->GetXaxis()->SetTitle("Electron angle [rad]");
h_theta_e_MC->SetMinimum(100.);
h_theta_e_MC->Draw("E");
h_theta_e_RD->SetLineColor(kPink+10);
h_theta_e_RD->Draw("E same");
l->Draw();
gPad->SetLogy();
gStyle->SetOptStat(0);
cout << "after scaling h_theta_e_MC entries: " << h_theta_e_MC->Integral() << endl;
cout << "after scaling h_theta_e_RD entries: " << h_theta_e_RD->Integral() << endl;

a1.cd(2);

TH1D * h3 = (TH1D*) h_theta_e_RD->Clone();
h3->Divide(h_theta_e_MC);
h3->SetMinimum(0.9);
h3->SetMaximum(1.15);
h3->SetTitle("Electron scattering angle");
h3->GetYaxis()->SetTitle("Data/MC ratio");
h3->GetXaxis()->SetTitle("Electron angle [rad]");
TLine *line_ = new TLine(0.005,1.,0.02,1.);
TLine *line_p = new TLine(0.005,1.03,0.02,1.03);
TLine *line_m = new TLine(0.005,0.97,0.02,0.97);

TLine *line_p05 = new TLine(0.005,1.005,0.02,1.005);
TLine *line_m05 = new TLine(0.005,0.995,0.02,0.995);

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

h3->SaveAs(Form("comparison_RDMC/%s_%s_ratio_RDnormalization_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_theta_el_%dhit_%dhit_smearing.root",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));
a1.cd(3);
//h_theta_e_MC->Scale(7.9729075e+08*5.5*1E+23*d_tar*1E-30*315.44638/6.79333e+07);//26563993,605675    137720.00    1.3242886e+09   6.79333e+07 9.3792854e+08
h_theta_mu_MC->SetMinimum(100.);
/*h_theta_mu_MC->GetXaxis()->SetRangeUser(0.0005,0.002);
h_theta_mu_RD->GetXaxis()->SetRangeUser(0.0005,0.002);*/
h_theta_mu_MC->GetXaxis()->SetRangeUser(0.0002,0.0015);//0.0003,0.0013);
h_theta_mu_RD->GetXaxis()->SetRangeUser(0.0002,0.0015);//0.0003,0.0013);
h_theta_mu_MC->GetYaxis()->SetTitle("Entries");
h_theta_mu_MC->SetTitle("Muon scattering angle");
h_theta_mu_MC->GetXaxis()->SetTitle("Muon angle [rad]");
h_theta_mu_MC->Draw("E");
h_theta_mu_RD->SetLineColor(kPink+10);
h_theta_mu_RD->Draw("E same");
l1->Draw();
gStyle->SetOptStat(0);
gPad->SetLogy();
cout << "after scaling h_theta_mu_MC entries: " << h_theta_mu_MC->Integral() << endl;
cout << "after scaling h_theta_mu_RD entries: " << h_theta_mu_RD->Integral() << endl;

a1.cd(4);
TH1D * h4 = (TH1D*) h_theta_mu_RD->Clone();
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
TLine *mu_line_p05 = new TLine(0.005,1.005,0.02,1.005);
TLine *mu_line_m05 = new TLine(0.005,0.995,0.02,0.995);
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
h4->SaveAs(Form("comparison_RDMC/%s_%s_ratio_RDnormalization_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_theta_mu_%dhit_%dhit_smearing.root",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));
a1.SaveAs(Form("comparison_RDMC/%s_%s_theta_MC_reco_RD_add_RDnormalization_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%dhit_smearing.pdf",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));

TF1 *Elastic = new TF1("Elastic","0.5109989461*0.001*((1+(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*cos(x)*cos(x))/(1-(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*cos(x)*cos(x)))",0,0.030); 
TF1 *Elastic2 = new TF1("Elastic2","asin( (sin(x)*sqrt(Elastic(x)*Elastic(x)-0.5109989461*0.001*0.5109989461*0.001))/sqrt( (160+0.5109989461*0.001-Elastic(x))*(160+0.5109989461*0.001-Elastic(x))-105.6583745 *0.001*105.6583745 *0.001 ) )",0,0.030);




TCanvas t("t","t",3000,3000);
t.Divide(3,3);
t.cd(1);
h3->SetLineColor(kRed+1);
h3->Rebin(1);
h3->GetXaxis()->SetRangeUser(0.00,0.032);//0.005,0.02);
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
h->GetXaxis()->SetRangeUser(0.00,0.032);//0.005,0.02);
h->Draw();
line_->Draw("same");
line_p->Draw("same");
line_m->Draw("same");

t.cd(4);
h_theta_e_MC->SetMinimum(100.);
h_theta_e_MC->GetXaxis()->SetTitle("Electron angle [rad]");
h_theta_e_MC->Rebin(1);
h_theta_e_RD->Rebin(1);
h_theta_e_MC->GetXaxis()->SetRangeUser(0.00,0.032);//0.005,0.02);
h_theta_e_RD->GetXaxis()->SetRangeUser(0.00,0.032);//0.005,0.02);
h_theta_e_MC->Draw("E");
h_theta_e_RD->SetLineColor(kRed+1);
h_theta_e_RD->Draw("E same");
gStyle->SetOptStat(0);
//gPad->SetLogy();
l->Draw();
t.cd(5);
h_theta_mu_MC->SetMinimum(100.);
h_theta_mu_MC->GetXaxis()->SetTitle("Muon angle [rad]");
h_theta_mu_MC->Rebin(1);
h_theta_mu_RD->Rebin(1);
h_theta_mu_MC->GetXaxis()->SetRangeUser(0.0002,0.0016);
h_theta_mu_RD->GetXaxis()->SetRangeUser(0.0002,0.0016);
h_theta_mu_MC->Draw("E");
h_theta_mu_RD->SetLineColor(kOrange-3);
h_theta_mu_RD->Draw("E same");
gStyle->SetOptStat(0);
//gPad->SetLogy();
l1->Draw();
t.cd(6);
h_opening_MC->SetMinimum(100.);
h_opening_MC->GetXaxis()->SetTitle("Opening angle [rad]");
h_opening_MC->Rebin(1);
h_opening_RD->Rebin(1);
h_opening_MC->GetXaxis()->SetRangeUser(0.00,0.032);//0.005,0.02);
h_opening_RD->GetXaxis()->SetRangeUser(0.00,0.032);//0.005,0.02);
h_opening_MC->Draw("E");
h_opening_RD->SetLineColor(kGreen-2);
h_opening_RD->Draw("E same");
gStyle->SetOptStat(0);
//gPad->SetLogy();
l2->Draw();
t.cd(7);
h_theta_e_MC->Draw("E");
h_theta_e_RD->Draw("E same");
h_theta_e_MC->SetMinimum(100.);
gStyle->SetOptStat(0);
gPad->SetLogy();
l->Draw();
t.cd(8);
h_theta_mu_MC->Draw("E");
h_theta_mu_RD->Draw("E same");
h_theta_mu_MC->SetMinimum(100.);
gStyle->SetOptStat(0);
gPad->SetLogy();
l1->Draw();
t.cd(9);
h_opening_MC->Draw("E");
h_opening_RD->Draw("E same");
h_opening_MC->SetMinimum(100.);
gStyle->SetOptStat(0);
gPad->SetLogy();
l2->Draw();
t.SaveAs(Form("comparison_RDMC/%s_%s_tesi_comparison_RDnormalization_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_%dhit_%dhit_smearing.pdf",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));


h_aco_MC->Scale(h_aco_RD->Integral()/h_aco_MC->Integral());
h_vrtx_chi2_MC->Scale(h_vrtx_chi2_RD->Integral()/h_vrtx_chi2_MC->Integral());
h_posZ_MC->Scale(h_posZ_RD->Integral()/h_posZ_MC->Integral());

h_ns1_MC->Scale(h_ns1_RD->Integral()/h_ns1_MC->Integral());


TCanvas v("v","v",2000,2000);
v.Divide(2,2);
v.cd(1);
h_aco_MC->GetXaxis()->SetRangeUser(-0.4,0.4);
h_aco_MC->Draw("E hist");
h_aco_RD->SetLineColor(kPink);
h_aco_RD->Draw("E hist same");
v.cd(2);
h_vrtx_chi2_MC->Draw("E hist");
h_vrtx_chi2_RD->SetLineColor(kPink);
h_vrtx_chi2_RD->Draw("E hist same");
gPad->SetLogy();
v.cd(3);
h_posZ_MC->GetXaxis()->SetRangeUser(800.,920);
h_posZ_RD->SetLineColor(kPink);
h_posZ_MC->Draw("E hist");
h_posZ_RD->Draw("E hist same");
v.cd(4);
h_ns1_MC->Draw("E hist");
h_ns1_RD->SetLineColor(kPink);
h_ns1_RD->Draw("E hist same");
v.SaveAs(Form("comparison_RDMC/%s_%s_variables_comparison_RDnormalization_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_%dhit_%dhit_smearing.pdf",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));


h_aco_pre_MC->Scale(h_aco_pre_RD->Integral()/h_aco_pre_MC->Integral());
//h_ns1_pre_MC->Scale(h_ns1_pre_RD->Integral()/h_ns1_pre_MC->Integral());
h_posZ_pre_MC->Scale(h_posZ_pre_RD->Integral()/h_posZ_pre_MC->Integral());

TCanvas v_pre("v_pre","v_pre",2000,2000);
v_pre.Divide(2,3);
v_pre.cd(1);
h_aco_pre_MC->SetTitle("Acoplanarity - Preselection");
h_aco_pre_MC->GetXaxis()->SetRangeUser(-3.2,3.2);
h_aco_pre_MC->GetYaxis()->SetTitle("Entries");
h_aco_pre_MC->GetXaxis()->SetTitle("Acoplanarity [rad]");
h_aco_pre_MC->Draw("E hist");
h_aco_pre_RD->SetLineColor(kPink);
h_aco_pre_RD->Draw("E hist same");
v_pre.cd(2);
h_aco_MC->SetTitle("Acoplanarity - Elastic selection");
h_aco_MC->GetYaxis()->SetTitle("Entries");
h_aco_MC->GetXaxis()->SetTitle("Acoplanarity [rad]");
h_aco_MC->GetXaxis()->SetRangeUser(-0.4,0.4);
h_aco_MC->Draw("E hist");
h_aco_RD->SetLineColor(kPink);
h_aco_RD->Draw("E hist same");
v_pre.cd(3);
h_posZ_pre_MC->SetTitle("Vertex Z - Preselection");
h_posZ_pre_MC->GetYaxis()->SetTitle("Entries");
h_posZ_pre_MC->GetXaxis()->SetTitle("Vertex Z [cm]");
h_posZ_pre_MC->GetXaxis()->SetRangeUser(800.,920);
h_posZ_pre_RD->SetLineColor(kPink);
h_posZ_pre_MC->Draw("E hist");
h_posZ_pre_RD->Draw("E hist same");
v_pre.cd(4);
h_posZ_MC->SetTitle("Vertex Z - Elastic selection");
h_posZ_MC->GetYaxis()->SetTitle("Entries");
h_posZ_MC->GetXaxis()->SetTitle("Vertex Z [cm]");
h_posZ_MC->GetXaxis()->SetRangeUser(800.,920);
h_posZ_RD->SetLineColor(kPink);
h_posZ_MC->Draw("E hist");
h_posZ_RD->Draw("E hist same");
v_pre.cd(5);
h_ns1_pre_MC->SetTitle("Number of stubs Station1 - Preselection");
h_ns1_pre_MC->GetYaxis()->SetTitle("Entries");
h_ns1_pre_MC->GetXaxis()->SetTitle("Number of stubs");
h_ns1_pre_RD->SetLineColor(kPink);
h_ns1_pre_MC->Draw("E hist ");
h_ns1_pre_RD->Draw("E hist same");
h_ns1_pre_MC->SetMinimum(1.);
gPad->SetLogy();
v_pre.cd(6);
h_ns1_MC->SetTitle("Number of stubs Station1 - Elastic selection");
h_ns1_MC->GetYaxis()->SetTitle("Entries");
h_ns1_MC->GetXaxis()->SetTitle("Number of stubs");
h_ns1_MC->Draw("E hist");
h_ns1_RD->SetLineColor(kPink);
h_ns1_RD->Draw("E hist  same");
h_ns1_MC->SetMinimum(1.);


gPad->SetLogy();
v_pre.SaveAs(Form("comparison_RDMC/%s_%s_variables_preCuts_comparison_RDnormalization_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_%dhit_%dhit_smearing.pdf",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));



double w2_elastic=0.;

    if(mc=="my_modifica")w2_elastic = 1.29605e+06;
    else if(mc=="my_modifica_norun")w2_elastic = 1.26266e+06;
    else if(mc=="my_modifica_LO")w2_elastic = 997218;
    else if(mc=="my_modifica_LO_pure")w2_elastic = 971478;
    else if(mc=="my_modifica_NLO")w2_elastic = 1.91465e+06;


cout << "chi2 per bin" << endl;
cout << "ELECTRON ANGLE: " << endl;

h_theta_e_RD->Chi2Test(h_theta_e_MC,"UW P");

for(int b=0; b < h_theta_e_MC->GetNbinsX(); b++){

Double_t errormc;
Double_t integral = h_theta_e_MC->IntegralAndError(b,b,errormc,"");
cout << "errormc " << errormc << endl;
cout << "errorrd " << sqrt(h_theta_e_RD->GetBinContent(b)) << endl;
cout << "bin " << b << ") chi2 = " << pow((h_theta_e_RD->GetBinContent(b)-h_theta_e_MC->GetBinContent(b))/(sqrt( errormc*errormc  + h_theta_e_RD->GetBinContent(b) )),2) << endl;

}

cout << "MUON ANGLE: " << endl;

h_theta_mu_RD->Chi2Test(h_theta_mu_MC,"UW P");

for(int b=0; b < h_theta_mu_MC->GetNbinsX(); b++){

Double_t errormc;
Double_t integral = h_theta_mu_MC->IntegralAndError(b,b,errormc,"");
cout << "errormc " << errormc << endl;
cout << "errorrd " << sqrt(h_theta_mu_RD->GetBinContent(b)) << endl;
cout << "bin " << b << ") chi2 = " << pow((h_theta_mu_RD->GetBinContent(b)-h_theta_mu_MC->GetBinContent(b))/(sqrt( errormc*errormc + h_theta_mu_RD->GetBinContent(b) )),2) << endl;

}


TCanvas rl("rl","rl",1000,500);
rl.Divide(2,1);
rl.cd(1);

h_thetaL_RD->Rebin(4);
h_thetaL_MC->Rebin(4);
h_thetaR_RD->Rebin(4);
h_thetaR_MC->Rebin(4);

TLine *line_2 = new TLine(0.0,1.,0.02,1.);
TLine *line_p2 = new TLine(0.0,1.03,0.02,1.03);
TLine *line_m2 = new TLine(0.0,0.97,0.02,0.97);

TH1D * hrl = (TH1D*) h_thetaL_RD->Clone();
hrl->Divide(h_thetaL_MC);
//h->SetMaximum(3.);
hrl->SetMinimum(0.8);
hrl->SetMaximum(1.2);
hrl->GetXaxis()->SetTitle("#theta L [rad]");
hrl->Draw("E");
line_2->Draw("same");
line_p2->Draw("same");
line_m2->Draw("same");
gStyle->SetOptStat(0);
rl.cd(2);
TH1D * hrl1 = (TH1D*) h_thetaR_RD->Clone();
hrl1->Divide(h_thetaR_MC);
//h->SetMaximum(3.);
hrl1->SetMinimum(0.8);
hrl1->SetMaximum(1.2);
hrl1->GetXaxis()->SetTitle("#theta R [rad]");
hrl1->Draw("E");
line_2->Draw("same");
line_p2->Draw("same");
line_m2->Draw("same");
gStyle->SetOptStat(0);
rl.SaveAs(Form("comparison_RDMC/%s_%s_thetaRL_comparison_RDnormalization_basicAllFiducial_chi150_flatArea_the20_theNA_thmu02_nstubs15_nochivrtx_zPosNoPeak_%dhit_%dhit_smearing.pdf",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));







}
