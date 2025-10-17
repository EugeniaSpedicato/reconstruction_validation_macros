#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TStyle.h>
#include "TGraphErrors.h"



void norun_ratio(int nhits_rd,int nhits_mc, string rd, string mc){

TH1::SetDefaultSumw2(kTRUE);
int const NBINS = 6;


TFile *f_RD=TFile::Open(Form("comparison_RDMC/%s_d_aco_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_thmu_RD=TFile::Open(Form("comparison_RDMC/%s_theta_mu_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_the_RD=TFile::Open(Form("comparison_RDMC/%s_theta_e_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_2D_RD=TFile::Open(Form("comparison_RDMC/%s_2D_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_op_RD=TFile::Open(Form("comparison_RDMC/%s_opening_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));

TFile* f_the_MC=TFile::Open(Form("comparison_RDMC/%s_theta_e_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile* f_thmu_MC=TFile::Open(Form("comparison_RDMC/%s_theta_mu_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile* f_2D_MC=TFile::Open(Form("comparison_RDMC/%s_2D_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile* f_op_MC=TFile::Open(Form("comparison_RDMC/%s_opening_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile* f_aco_MC=TFile::Open(Form("comparison_RDMC/%s_d_aco_RD_parallel_pre_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc)); 

TFile *f_d0_MC=TFile::Open(Form("comparison_RDMC/%s_distance_MOD0_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile *f_d1_MC=TFile::Open(Form("comparison_RDMC/%s_distance_MOD1_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile *f_d2_MC=TFile::Open(Form("comparison_RDMC/%s_distance_MOD2_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile *f_d3_MC=TFile::Open(Form("comparison_RDMC/%s_distance_MOD3_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile *f_d4_MC=TFile::Open(Form("comparison_RDMC/%s_distance_MOD4_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile *f_d5_MC=TFile::Open(Form("comparison_RDMC/%s_distance_MOD5_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));

TFile *f_d0_RD=TFile::Open(Form("comparison_RDMC/%s_distance_MOD0_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_d1_RD=TFile::Open(Form("comparison_RDMC/%s_distance_MOD1_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_d2_RD=TFile::Open(Form("comparison_RDMC/%s_distance_MOD2_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_d3_RD=TFile::Open(Form("comparison_RDMC/%s_distance_MOD3_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_d4_RD=TFile::Open(Form("comparison_RDMC/%s_distance_MOD4_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_d5_RD=TFile::Open(Form("comparison_RDMC/%s_distance_MOD5_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));


TFile *f_vrtx_chi2_RD=TFile::Open(Form("comparison_RDMC/%s_vrtx_chi2_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_track_el_RD=TFile::Open(Form("comparison_RDMC/%s_track_el_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_track_mu_RD=TFile::Open(Form("comparison_RDMC/%s_track_mu_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_Bend_%dhit_data.root",rd.c_str(),nhits_rd));

TFile *f_vrtx_chi2_MC=TFile::Open(Form("comparison_RDMC/%s_vrtx_chi2_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile *f_track_el_MC=TFile::Open(Form("comparison_RDMC/%s_track_el_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile *f_track_mu_MC=TFile::Open(Form("comparison_RDMC/%s_track_mu_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));


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

TH1D * h_d0_MC=(TH1D*)f_d0_MC->Get("h_dist0");
TH1D * h_d1_MC=(TH1D*)f_d1_MC->Get("h_dist1");
TH1D * h_d2_MC=(TH1D*)f_d2_MC->Get("h_dist2");
TH1D * h_d3_MC=(TH1D*)f_d3_MC->Get("h_dist3");
TH1D * h_d4_MC=(TH1D*)f_d4_MC->Get("h_dist4");
TH1D * h_d5_MC=(TH1D*)f_d5_MC->Get("h_dist5");

TH1D * h_d0_RD=(TH1D*)f_d0_RD->Get("h_dist0");
TH1D * h_d1_RD=(TH1D*)f_d1_RD->Get("h_dist1");
TH1D * h_d2_RD=(TH1D*)f_d2_RD->Get("h_dist2");
TH1D * h_d3_RD=(TH1D*)f_d3_RD->Get("h_dist3");
TH1D * h_d4_RD=(TH1D*)f_d4_RD->Get("h_dist4");
TH1D * h_d5_RD=(TH1D*)f_d5_RD->Get("h_dist5");


TH1D *h_vrtx_chi2_MC=(TH1D*)f_vrtx_chi2_MC->Get("vrtx_chi2");
TH1D *h_track_el_MC=(TH1D*)f_track_el_MC->Get("track_el");
TH1D *h_track_mu_MC=(TH1D*)f_track_mu_MC->Get("track_mu");


TH1D *h_vrtx_chi2_RD=(TH1D*)f_vrtx_chi2_RD->Get("vrtx_chi2");
TH1D *h_track_el_RD=(TH1D*)f_track_el_RD->Get("track_el");
TH1D *h_track_mu_RD=(TH1D*)f_track_mu_RD->Get("track_mu");


if(rd=="my_modifica_golden_all" or rd=="run8") {h_theta_e_MC->Rebin(2); h_theta_mu_MC->Rebin(2);}


   double cross_section = 0.;
   double entries = 0.;
   double w2_fiducial = 0.;
   double n_mu_ot = 0.;
   double d_tar = 0.;
   double w2_elastic = 0.;
   double error_CS=0.;

   double wnorm=1335.01;

   double a_sum_Nw[11]={1492877.,1492853.,1492597.,1493119.,1492810.,1492642.,1493036.,1492805.,1492731.,1492824.,1492634.};

   int N = 11;

   double sum_Nw = 0;
   sum_Nw = accumulate(a_sum_Nw, a_sum_Nw+N, sum_Nw);


   double sum_w_wnorm = wnorm*1.63541e+07;

   double sum_w2_wnorm = wnorm*wnorm*5.25094e+07;


  cout<< "sum_Nw " << sum_Nw << ", sum_w_wnorm " << sum_w_wnorm << ", sum_w2_wnorm " << sum_w2_wnorm << endl;

	error_CS=sqrt(( (sum_w2_wnorm/sum_Nw)  -  (sum_w_wnorm/sum_Nw)*(sum_w_wnorm/sum_Nw) )/sum_Nw);

  cout<< "error_CS " << error_CS << endl;

  cout << "cs " << sum_w_wnorm/sum_Nw << endl;

    entries = 1.11042e+07;
    cross_section = wnorm * 1.63541e+07/sum_Nw;

  cout << "cs " << sum_w_wnorm/sum_Nw << endl;

if(rd=="my_modifica_golden_all") n_mu_ot=8.98816e+08;
if(rd=="my_modifica_eos") n_mu_ot=8.8016548e+08*0.8567;//n_mu_ot=7.9729075e+08;
if(rd=="run8") n_mu_ot=4.2365303e+10*0.8567;//n_mu_ot=7.9729075e+08;

    d_tar=3.;
    w2_fiducial=3.5394e+07;
    w2_elastic= 1.26266e+06;

 double lumi_rd=n_mu_ot* 5.5*1E+23*d_tar*1E-30 * 0.973;
 double lumi_mc=entries/cross_section;
 double ratio_lumi=lumi_rd/lumi_mc;

cout << "lumi_mc: " << lumi_mc << ", lumi_rd: " << lumi_rd << ", ratio_lumi: " << ratio_lumi << endl;


double error_lumi_rd= 6/12.011 * 6.02214076*1e+23 * 1e-30 * sqrt(
			 (sqrt(n_mu_ot)*0.973*d_tar*1.83)*(sqrt(n_mu_ot)*0.973*d_tar*1.83)
			 + (n_mu_ot*0.001*d_tar*1.83)*(n_mu_ot*0.001*d_tar*1.83)
			 + (n_mu_ot*0.973*0.001*1.83)*(n_mu_ot*0.973*0.001*1.83)
			 + (n_mu_ot*0.973*d_tar*0.01)*(n_mu_ot*0.973*d_tar*0.01)
			);


double error_lumi_mc= 1/sqrt(w2_fiducial)*lumi_mc;

double error_lumi_ratio=sqrt( (error_lumi_rd/lumi_mc)*(error_lumi_rd/lumi_mc) + (lumi_rd/(lumi_mc*lumi_mc)*error_lumi_mc)*(lumi_rd/(lumi_mc*lumi_mc)*error_lumi_mc) );

cout<<"error_lumi_mc: " << error_lumi_mc << ", error_lumi_rd: " << error_lumi_rd << ", error_lumi_ratio: " << error_lumi_ratio << endl;

double error_integral_ratio= sqrt( (sqrt(h_theta_e_RD->Integral())/(h_theta_e_MC->Integral()*ratio_lumi))*(sqrt(h_theta_e_RD->Integral())/(h_theta_e_MC->Integral()*ratio_lumi))
				  +(h_theta_e_RD->Integral()/(h_theta_e_MC->Integral()*h_theta_e_MC->Integral()*ratio_lumi)*sqrt(w2_elastic))*(h_theta_e_RD->Integral()/(h_theta_e_MC->Integral()*h_theta_e_MC->Integral()*ratio_lumi)*sqrt(w2_elastic))
				  +(h_theta_e_RD->Integral()/(h_theta_e_MC->Integral()*ratio_lumi*ratio_lumi)*error_lumi_ratio)*(h_theta_e_RD->Integral()/(h_theta_e_MC->Integral()*ratio_lumi*ratio_lumi)*error_lumi_ratio)
				 );


cout << "Ndata/NMC before scaling " << h_theta_e_RD->Integral()/h_theta_e_MC->Integral() << endl;
cout << "error_integral_ratio Ndata/NMC " << error_integral_ratio << endl;

cout << endl;
cout << "---------------------------" << endl;


	 h_theta_e_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section*0.973/entries);//9.3792854e+08,7.9729075e+08
         h_theta_mu_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section*0.973/entries);
	 h_2D_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section*0.973/entries);
	 h_opening_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section*0.973/entries);
         h_aco_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section*0.973/entries);
h_d0_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section*0.973/entries);
h_d1_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section*0.973/entries);
h_d2_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section*0.973/entries);
h_d3_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section*0.973/entries);
h_d4_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section*0.973/entries);
h_d5_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section*0.973/entries);

h_vrtx_chi2_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section*0.973/entries);
h_track_el_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section*0.973/entries);
h_track_mu_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section*0.973/entries);


cout << "Ndata/NMC after scaling " << h_theta_e_RD->Integral()/h_theta_e_MC->Integral() << endl;


TCanvas chi("chi","chi",1500,500);
chi.Divide(3,1);
chi.cd(1);

h_vrtx_chi2_MC->GetXaxis()->SetRangeUser(0.,6.);

h_vrtx_chi2_MC->Draw("hist");
h_vrtx_chi2_RD->SetLineColor(kPink);
h_vrtx_chi2_RD->Draw("hist same");
chi.cd(2);

h_track_el_MC->GetXaxis()->SetRangeUser(0.,6.);
h_track_el_MC->Draw("hist");
h_track_el_RD->SetLineColor(kPink);
h_track_el_RD->Draw("hist same");
chi.cd(3);

h_track_mu_MC->GetXaxis()->SetRangeUser(0.,6.);
h_track_mu_MC->Draw("hist");
h_track_mu_RD->SetLineColor(kPink);
h_track_mu_RD->Draw("hist same");

chi.SaveAs(Form("comparison_RDMC/%s_%s_chi2_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_%dhit_%dhit.pdf",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));


TCanvas di("di","di",1400,2100);
di.Divide(2,3);
di.cd(1);
h_d0_MC->SetTitle("Muon-electron hit separation - Module X0");
h_d0_MC->GetXaxis()->SetTitle("#mu-e hit distance [cm]");
h_d0_MC->GetYaxis()->SetTitle("Entries");
h_d0_MC->GetXaxis()->SetRangeUser(-1.,1.);
h_d0_MC->SetMinimum(0.);
h_d0_MC->Draw("hist");
h_d0_RD->SetLineColor(kPink);
h_d0_RD->Draw("hist same");
di.cd(2);
h_d1_MC->SetTitle("Muon-electron hit separation - Module Y0");
h_d1_MC->GetXaxis()->SetTitle("#mu-e hit distance [cm]");
h_d1_MC->GetYaxis()->SetTitle("Entries");
h_d1_MC->SetMinimum(0.);
h_d1_MC->GetXaxis()->SetRangeUser(-1.,1.);
h_d1_MC->Draw("hist");
h_d1_RD->SetLineColor(kPink);
h_d1_RD->Draw("hist same");
di.cd(3);
h_d2_MC->SetTitle("Muon-electron hit separation - Module U");
h_d2_MC->GetXaxis()->SetTitle("#mu-e hit distance [cm]");
h_d2_MC->GetYaxis()->SetTitle("Entries");
h_d2_MC->GetXaxis()->SetRangeUser(-1.,1.);
h_d2_MC->SetMinimum(0.);
h_d2_MC->Draw("hist");
h_d2_RD->SetLineColor(kPink);
h_d2_RD->Draw("hist same");
di.cd(4);
h_d3_MC->SetTitle("Muon-electron hit separation - Module V");
h_d3_MC->GetXaxis()->SetTitle("#mu-e hit distance [cm]");
h_d3_MC->GetYaxis()->SetTitle("Entries");
h_d3_MC->GetXaxis()->SetRangeUser(-1.,1.);
h_d3_MC->SetMinimum(0.);
h_d3_MC->Draw("hist");
h_d3_RD->SetLineColor(kPink);
h_d3_RD->Draw("hist same");
di.cd(5);
h_d4_MC->SetTitle("Muon-electron hit separation - Module X1");
h_d4_MC->GetXaxis()->SetTitle("#mu-e hit distance [cm]");
h_d4_MC->GetYaxis()->SetTitle("Entries");
h_d4_MC->SetMinimum(0.);
h_d4_MC->Draw("hist");
h_d4_RD->SetLineColor(kPink);
h_d4_RD->Draw("hist same");
di.cd(6);
h_d5_MC->SetTitle("Muon-electron hit separation - Module Y1");
h_d5_MC->GetXaxis()->SetTitle("#mu-e hit distance [cm]");
h_d5_MC->GetYaxis()->SetTitle("Entries");
h_d5_MC->SetMinimum(0.);
h_d5_MC->Draw("hist");
h_d5_RD->SetLineColor(kPink);
h_d5_RD->Draw("hist same");
di.SaveAs(Form("comparison_RDMC/%s_%s_h_dist_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%dhit.pdf",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));

TCanvas aco("aco","aco",2100,700);
aco.Divide(2,1);
aco.cd(1);
h_aco_MC->Rebin(8);
h_aco_RD->Rebin(8);
h_aco_RD->SetLineColor(kPink);
h_aco_MC->SetMinimum(1.);
h_aco_MC->Draw("hist");
h_aco_RD->Draw("hist same");
gPad->SetLogy();
aco.cd(2);
TH1D * h_aco = (TH1D*) h_aco_RD->Clone();
h_aco->Divide(h_aco_MC);
h_aco->SetMaximum(1.);
h_aco->GetXaxis()->SetRangeUser(-1,1);
h_aco->Draw();
aco.SaveAs(Form("comparison_RDMC/%s_%s_h_aco_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%dhit.pdf",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));


TCanvas a("a","a",700,700);
a.Divide(1,2);
a.cd(1);

h_opening_MC->SetTitle("#mu-e opening angle");
h_opening_MC->SetMinimum(1.);
h_opening_MC->GetXaxis()->SetRangeUser(0.005,0.032);
h_opening_RD->GetXaxis()->SetRangeUser(0.005,0.032);
h_opening_RD->SetMinimum(1.);
h_opening_MC->Draw("hist");
h_opening_RD->SetLineColor(kPink+10);
h_opening_RD->Draw("hist same");
cout << "h_opening_MC entries: " << h_opening_MC->Integral() << endl;
cout << "h_opening_RD entries: " << h_opening_RD->Integral() << endl;
a.cd(2);

TH1D * h = (TH1D*) h_opening_RD->Clone();
h->Divide(h_opening_MC);
//h->SetMaximum(3.);
h->SetMinimum(0.);

double x[h->GetNbinsX()];
double y[h->GetNbinsX()];
double ex[h->GetNbinsX()];
double ey[h->GetNbinsX()];
double x_sum[h->GetNbinsX()];
double y_sum[h->GetNbinsX()];
double ex_sum[h->GetNbinsX()];
double ey_sum[h->GetNbinsX()];

double rd_ =0.;
double mc_ =0.;

for(int n=0; n < h->GetNbinsX(); n++){ x[n]=h->GetBinCenter(n); y[n]=h->GetBinContent(n); ex[n]=0.; ey[n]=1./sqrt(h_opening_RD->GetBinContent(n))*h->GetBinContent(n);}//h->GetBinError(n);}
TGraphErrors *op_er=new TGraphErrors(h->GetNbinsX(),x,y,ex,ey);

for(int n=0; n < h->GetNbinsX(); n++){
rd_ = h_opening_RD->GetBinContent(n);
mc_ = h_opening_MC->GetBinContent(n);
 x_sum[n]=h->GetBinCenter(n); y_sum[n]=h->GetBinContent(n); ex_sum[n]=0.; ey[n]=sqrt( (rd_)/(mc_*mc_)+(rd_*rd_)/(mc_*mc_*mc_) );}//h->GetBinError(n);}
TGraphErrors *op_er_sum=new TGraphErrors(h->GetNbinsX(),x,y,ex,ey);

h->SetMaximum(1.);
op_er_sum->SetMaximum(1.);
op_er->SetMaximum(1.);
op_er->SetMinimum(0.);
op_er->GetXaxis()->SetRangeUser(0.,0.032);
//op_er->Draw("AP");
op_er_sum->SetMarkerColor(kOrange);
//op_er_sum->Draw("AP");

h->SetTitle("#mu-e opening angle");
h->GetXaxis()->SetTitle("Opening angle [rad]");
h->Draw("E");
gStyle->SetOptStat(0);
h->SaveAs(Form("comparison_RDMC/%s_%s_ratio_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_opening_%dhit_%dhit.root",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));
a.SaveAs(Form("comparison_RDMC/%s_%s_opening_reb_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%dhit.pdf",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));

TCanvas a1("a1","a1",1000,700);
a1.Divide(2,2);
a1.cd(1);

h_theta_e_MC->GetXaxis()->SetRangeUser(0.005,0.032);
h_theta_e_RD->GetXaxis()->SetRangeUser(0.005,0.032);
h_theta_e_MC->SetTitle("Electron scattering angle");
h_theta_e_MC->SetMinimum(1.);
h_theta_e_MC->Draw("E");
h_theta_e_RD->SetLineColor(kPink+10);
h_theta_e_RD->Draw("E same");
gPad->SetLogy();
gStyle->SetOptStat(0);
cout << "h_theta_e_MC entries: " << h_theta_e_MC->Integral() << endl;
cout << "h_theta_e_RD entries: " << h_theta_e_RD->Integral() << endl;


a1.cd(2);

TH1D * h3 = (TH1D*) h_theta_e_RD->Clone();
h3->Divide(h_theta_e_MC);
h3->SetTitle("Electron scattering angle");
h3->SetMinimum(0.);
h3->GetXaxis()->SetTitle("Electron angle [rad]");
h3->Draw("E");
h3->SaveAs(Form("comparison_RDMC/%s_%s_ratio_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_theta_el_%dhit_%dhit.root",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));
a1.cd(3);
h_theta_mu_MC->SetMinimum(1.);
h_theta_mu_MC->GetXaxis()->SetRangeUser(0.0002,0.0014);
h_theta_mu_RD->GetXaxis()->SetRangeUser(0.0002,0.0014);
h_theta_mu_MC->SetTitle("Muon scattering angle");
h_theta_mu_MC->Draw("E");
h_theta_mu_RD->SetLineColor(kPink+10);
h_theta_mu_RD->Draw("E same");
gStyle->SetOptStat(0);
gPad->SetLogy();
cout << "h_theta_mu_MC entries: " << h_theta_mu_MC->Integral() << endl;
cout << "h_theta_mu_RD entries: " << h_theta_mu_RD->Integral() << endl;

a1.cd(4);
TH1D * h4 = (TH1D*) h_theta_mu_RD->Clone();
h4->Divide(h_theta_mu_MC);
h4->SetTitle("Muon scattering angle");
h4->SetMinimum(0.);
h4->GetYaxis()->SetTitle("Data/MC ratio");
h4->GetXaxis()->SetTitle("Muon angle [rad]");
h4->Draw("E");
h4->SaveAs(Form("comparison_RDMC/%s_%s_ratio_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_theta_mu_%dhit_%dhit.root",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));
a1.SaveAs(Form("comparison_RDMC/%s_%s_theta_MC_reco_RD_add_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%dhit.pdf",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));

TF1 *Elastic = new TF1("Elastic","0.5109989461*0.001*((1+(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*cos(x)*cos(x))/(1-(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*cos(x)*cos(x)))",0,0.032); 
TF1 *Elastic2 = new TF1("Elastic2","asin( (sin(x)*sqrt(Elastic(x)*Elastic(x)-0.5109989461*0.001*0.5109989461*0.001))/sqrt( (160+0.5109989461*0.001-Elastic(x))*(160+0.5109989461*0.001-Elastic(x))-105.6583745 *0.001*105.6583745 *0.001 ) )",0,0.032);


TCanvas a2("a2","a2",6000,3000);
a2.Divide(2,1);

a2.cd(1);
h_2D_MC->SetTitle("MESMER MC - Selected sample");

h_2D_MC->GetXaxis()->SetTitle("Electron angle [rad]");
h_2D_MC->GetYaxis()->SetTitle("Muon angle [rad]");
//h_2D_MC->RebinY(2);
//h_2D_MC->RebinX(2);
h_2D_MC->Draw("COLZ");
Elastic2->Draw("same");
a2.cd(2);
h_2D_RD->SetTitle("Real Data - Selected sample");
h_2D_RD->GetXaxis()->SetTitle("Electron angle [rad]");
h_2D_RD->GetYaxis()->SetTitle("Muon angle [rad]");
//h_2D_RD->RebinX(2);
//h_2D_RD->RebinY(2);
h_2D_RD->Draw("COLZ");
Elastic2->Draw("same");
TGaxis::SetMaxDigits(3);

a2.SaveAs(Form("comparison_RDMC/%s_%s_comparison_2D_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%dhit.pdf",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));



 auto l = new TLegend(0.55,0.15,0.9,0.3);
l->AddEntry(h_theta_e_MC,"MC","LEP");;
l->AddEntry(h_theta_e_RD,"data","LEP");

 auto l1 = new TLegend(0.55,0.15,0.9,0.3);
l1->AddEntry(h_theta_mu_MC,"MC","LEP");;
l1->AddEntry(h_theta_mu_RD,"data","LEP");

 auto l2 = new TLegend(0.55,0.15,0.9,0.3);
l2->AddEntry(h_opening_MC,"MC","LEP");;
l2->AddEntry(h_opening_RD,"data","LEP");


TCanvas t("t","t",3100,2000);
t.Divide(3,2);
t.cd(4);
h3->GetYaxis()->SetTitle("Data/MC ratio");
h3->SetLineColor(kRed+1);
h3->Draw();
t.cd(5);
h4->GetYaxis()->SetTitle("Data/MC ratio");
h4->SetLineColor(kOrange-3);
h4->Draw();
t.cd(6);
h->GetYaxis()->SetTitle("Data/MC ratio");
h->SetLineColor(kGreen-2);
h->Draw();

t.cd(1);
h_theta_e_MC->SetMinimum(1.);
h_theta_e_MC->GetYaxis()->SetTitle("Entries");
h_theta_e_MC->GetXaxis()->SetTitle("Electron angle [rad]");
h_theta_e_MC->Draw("E");
h_theta_e_RD->SetLineColor(kRed+1);
h_theta_e_RD->Draw("E same");
gStyle->SetOptStat(0);
gPad->SetLogy();
l->Draw();
t.cd(2);
h_theta_mu_MC->SetMinimum(1.);
h_theta_mu_MC->GetYaxis()->SetTitle("Entries");
h_theta_mu_MC->GetXaxis()->SetTitle("Muon angle [rad]");
h_theta_mu_MC->Draw("E");
h_theta_mu_RD->SetLineColor(kOrange-3);
h_theta_mu_RD->Draw("E same");
gStyle->SetOptStat(0);
gPad->SetLogy();
l1->Draw();
t.cd(3);
h_opening_MC->GetYaxis()->SetTitle("Entries");
h_opening_MC->SetMinimum(1.);
h_opening_MC->GetXaxis()->SetTitle("Opening angle [rad]");
h_opening_MC->Draw("E");
h_opening_RD->SetLineColor(kGreen-2);
h_opening_RD->Draw("E same");
gStyle->SetOptStat(0);
gPad->SetLogy();
l2->Draw();
t.SaveAs(Form("comparison_RDMC/%s_%s_tesi_comparison_basicAllFiducial_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_%dhit_%dhit.pdf",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));

}

