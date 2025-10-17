#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TStyle.h>
#include "TGraphErrors.h"



void different_hitsos_sig(int nhits_rd,int nhits_mc, string rd, string rd1, string rd2, string rd3, string mc){

TH1::SetDefaultSumw2(kTRUE);
int const NBINS = 6;



TFile *f_thmu_RD0=TFile::Open(Form("comparison_RDMC/%s_theta_mu_RD_parallel_pre_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_the_RD0=TFile::Open(Form("comparison_RDMC/%s_theta_e_RD_parallel_pre_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_op_RD0=TFile::Open(Form("comparison_RDMC/%s_opening_RD_parallel_pre_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));


TFile *f_thmu_RD1=TFile::Open(Form("comparison_RDMC/%s_theta_mu_RD_parallel_pre_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd1.c_str(),nhits_rd));
TFile *f_the_RD1=TFile::Open(Form("comparison_RDMC/%s_theta_e_RD_parallel_pre_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd1.c_str(),nhits_rd));
TFile *f_op_RD1=TFile::Open(Form("comparison_RDMC/%s_opening_RD_parallel_pre_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd1.c_str(),nhits_rd));


TFile *f_thmu_RD2=TFile::Open(Form("comparison_RDMC/%s_theta_mu_RD_parallel_pre_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd2.c_str(),nhits_rd));
TFile *f_the_RD2=TFile::Open(Form("comparison_RDMC/%s_theta_e_RD_parallel_pre_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd2.c_str(),nhits_rd));
TFile *f_op_RD2=TFile::Open(Form("comparison_RDMC/%s_opening_RD_parallel_pre_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd2.c_str(),nhits_rd));


TFile *f_thmu_RD3=TFile::Open(Form("comparison_RDMC/%s_theta_mu_RD_parallel_pre_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd3.c_str(),nhits_rd));
TFile *f_the_RD3=TFile::Open(Form("comparison_RDMC/%s_theta_e_RD_parallel_pre_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd3.c_str(),nhits_rd));
TFile *f_op_RD3=TFile::Open(Form("comparison_RDMC/%s_opening_RD_parallel_pre_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd3.c_str(),nhits_rd));


TFile* f_the_MC=TFile::Open(Form("comparison_RDMC/%s_theta_e_RD_parallel_pre_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile* f_thmu_MC=TFile::Open(Form("comparison_RDMC/%s_theta_mu_RD_parallel_pre_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile* f_op_MC=TFile::Open(Form("comparison_RDMC/%s_opening_RD_parallel_pre_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));



TH1D* h_theta_e_MC=(TH1D*)f_the_MC->Get("theta_e");
TH1D* h_theta_mu_MC=(TH1D*)f_thmu_MC->Get("theta_mu");
TH1D* h_opening_MC=(TH1D*)f_op_MC->Get("h_opening");

TH1D *h_theta_e_RD0=(TH1D*)f_the_RD0->Get("theta_e");
TH1D *h_theta_mu_RD0=(TH1D*)f_thmu_RD0->Get("theta_mu");
TH1D *h_opening_RD0=(TH1D*)f_op_RD0->Get("h_opening");

TH1D *h_theta_e_RD1=(TH1D*)f_the_RD1->Get("theta_e");
TH1D *h_theta_mu_RD1=(TH1D*)f_thmu_RD1->Get("theta_mu");
TH1D *h_opening_RD1=(TH1D*)f_op_RD1->Get("h_opening");

TH1D *h_theta_e_RD2=(TH1D*)f_the_RD2->Get("theta_e");
TH1D *h_theta_mu_RD2=(TH1D*)f_thmu_RD2->Get("theta_mu");
TH1D *h_opening_RD2=(TH1D*)f_op_RD2->Get("h_opening");

TH1D *h_theta_e_RD3=(TH1D*)f_the_RD3->Get("theta_e");
TH1D *h_theta_mu_RD3=(TH1D*)f_thmu_RD3->Get("theta_mu");
TH1D *h_opening_RD3=(TH1D*)f_op_RD3->Get("h_opening");


TH1D * h_theta_e_RD = (TH1D*) h_theta_e_RD0->Clone();
h_theta_e_RD->Add(h_theta_e_RD1);
h_theta_e_RD->Add(h_theta_e_RD2);
h_theta_e_RD->Add(h_theta_e_RD3);

h_theta_e_RD->SaveAs("comparison_RDMC/my_modifica_golden_theta_e_RD_parallel_pre_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_0hit_data.root");

TH1D * h_theta_mu_RD = (TH1D*) h_theta_mu_RD0->Clone();
h_theta_mu_RD->Add(h_theta_mu_RD1);
h_theta_mu_RD->Add(h_theta_mu_RD2);
h_theta_mu_RD->Add(h_theta_mu_RD3);

h_theta_e_RD->SaveAs("comparison_RDMC/my_modifica_golden_theta_mu_RD_parallel_pre_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_0hit_data.root");

TH1D * h_opening_RD = (TH1D*) h_opening_RD0->Clone();
h_opening_RD->Add(h_opening_RD1);
h_opening_RD->Add(h_opening_RD2);
h_opening_RD->Add(h_opening_RD3);

h_opening_RD->SaveAs("comparison_RDMC/my_modifica_golden_opening_RD_parallel_pre_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_0hit_data.root");


   Double_t cross_section = 0.;
   Double_t entries = 0.;
   Double_t w2_fiducial = 0.;
   Double_t n_mu_ot = 0.;
   Double_t d_tar = 0.;
   Double_t w2_elastic = 0.;

if(mc=="my_modifica" or mc=="my_modifica_eos"){

double error_CS=0.;
double wnorm=1335.01;
double a_sum_Nw[7]={1492877.,1492853.,1492597.,1493119.,1492810.,1492642.,1493036.};
double a_sum_w_wnorm[7]={2.0034687e+09,2.0068266e+09,1.998391e+09,2.0046823e+09,2.0071567e+09,1.997697e+09,2.0063223e+09};
double a_sum_w2_wnorm[7]={8.5077092e+12,9.1146482e+12,6.6389419e+12,8.6412965e+12,9.0973111e+12,7.7243284e+12,1.0926206e+13};
   int n = 7;
   double sum_Nw = 0;
   sum_Nw = accumulate(a_sum_Nw, a_sum_Nw+n, sum_Nw);

   double sum_w_wnorm = 0;
   sum_w_wnorm = accumulate(a_sum_w_wnorm, a_sum_w_wnorm+n, sum_w_wnorm);

   double sum_w2_wnorm = 0;
   sum_w2_wnorm = accumulate(a_sum_w2_wnorm, a_sum_w2_wnorm+n, sum_w2_wnorm);

cout << endl;
cout << "---------------------------" << endl;
  cout<< "sum_Nw " << sum_Nw << ", sum_w_wnorm " << sum_w_wnorm << ", sum_w2_wnorm " << sum_w2_wnorm << endl;

error_CS=sqrt(( (sum_w2_wnorm/sum_Nw)  -  (sum_w_wnorm/sum_Nw)*(sum_w_wnorm/sum_Nw) )/sum_Nw);

  cout<< "error_CS " << error_CS << endl;

    cross_section = 1342.0701;
    entries = 6.91701e+06;
    n_mu_ot=7.9729075e+08;
    //n_mu_ot=8.6157600e+08;
    d_tar=3.;
    w2_fiducial=2.19867e+07;
    w2_elastic=1.68164e+06;

 double lumi_rd=n_mu_ot* 5.5*1E+23*d_tar*1E-30 * 0.973;
 double lumi_mc=6.91701e+06/cross_section;
 double ratio_lumi=lumi_rd/lumi_mc;

cout << "lumi_mc: " << lumi_mc << ", lumi_rd: " << lumi_rd << ", ratio_lumi: " << ratio_lumi << endl;

//double error_lumi_mc0=sqrt( (sqrt(6.91701e+06)/1342.0701)*(sqrt(6.91701e+06)/1342.0701) + (6.91701e+06*error_CS/(1342.0701*1342.0701))*(6.91701e+06*error_CS/(1342.0701*1342.0701)));
//double error_lumi_mc=sqrt( (sqrt(w2_fiducial)/cross_section)*(sqrt(w2_fiducial)/cross_section) + (entries*error_CS/(cross_section*cross_section))*(entries*error_CS/(cross_section*cross_section)));

double error_lumi_rd= 6/12.011 * 6.02214076*1e+23 * 1e-30 * sqrt(
			 (sqrt(n_mu_ot)*0.973*d_tar*1.83)*(sqrt(n_mu_ot)*0.973*d_tar*1.83)
			 + (n_mu_ot*0.001*d_tar*1.83)*(n_mu_ot*0.001*d_tar*1.83)
			 + (n_mu_ot*0.973*0.001*1.83)*(n_mu_ot*0.973*0.001*1.83)
			 + (n_mu_ot*0.973*d_tar*0.01)*(n_mu_ot*0.973*d_tar*0.01)
			);


double error_lumi_mc= 1/sqrt(w2_fiducial)*lumi_mc;

double error_lumi_ratio=sqrt( (error_lumi_rd/lumi_mc)*(error_lumi_rd/lumi_mc) + (lumi_rd/(lumi_mc*lumi_mc)*error_lumi_mc)*(lumi_rd/(lumi_mc*lumi_mc)*error_lumi_mc) );

cout<<"error_lumi_mc: " << error_lumi_mc << ", error_lumi_rd: " << error_lumi_rd << ", error_lumi_ratio: " << error_lumi_ratio << endl;

/*double error_integral_ratio= sqrt( (sqrt(h_theta_e_RD->Integral())/h_theta_e_MC->Integral())*(sqrt(h_theta_e_RD->Integral())/h_theta_e_MC->Integral())
				+ (h_theta_e_RD->Integral()/(h_theta_e_MC->Integral()*h_theta_e_MC->Integral()) *sqrt(h_theta_e_MC->Integral()) )*
				(h_theta_e_RD->Integral()/(h_theta_e_MC->Integral()*h_theta_e_MC->Integral()) *sqrt(h_theta_e_MC->Integral())) );*/

/*
double error_integral_ratio= sqrt( (sqrt(h_theta_e_RD->Integral())/(w2_elastic*ratio_lumi))*(sqrt(h_theta_e_RD->Integral())/(w2_elastic*ratio_lumi))
				  +(h_theta_e_RD->Integral()/(w2_elastic*w2_elastic*ratio_lumi)*sqrt(w2_elastic))*(h_theta_e_RD->Integral()/(w2_elastic*w2_elastic*ratio_lumi)*sqrt(w2_elastic))
				  +(h_theta_e_RD->Integral()/(w2_elastic*ratio_lumi*ratio_lumi)*error_lumi_ratio)*(h_theta_e_RD->Integral()/(w2_elastic*ratio_lumi*ratio_lumi)*error_lumi_ratio)
				 );
*/
double error_integral_ratio= sqrt( (sqrt(h_theta_e_RD->Integral())/(h_theta_e_MC->Integral()*ratio_lumi))*(sqrt(h_theta_e_RD->Integral())/(h_theta_e_MC->Integral()*ratio_lumi))
				  +(h_theta_e_RD->Integral()/(h_theta_e_MC->Integral()*h_theta_e_MC->Integral()*ratio_lumi)*sqrt(w2_elastic))*(h_theta_e_RD->Integral()/(h_theta_e_MC->Integral()*h_theta_e_MC->Integral()*ratio_lumi)*sqrt(w2_elastic))
				  +(h_theta_e_RD->Integral()/(h_theta_e_MC->Integral()*ratio_lumi*ratio_lumi)*error_lumi_ratio)*(h_theta_e_RD->Integral()/(h_theta_e_MC->Integral()*ratio_lumi*ratio_lumi)*error_lumi_ratio)
				 );


cout << "Ndata/NMC before scaling " << h_theta_e_RD->Integral()/h_theta_e_MC->Integral() << endl;
cout << "error_integral_ratio Ndata/NMC " << error_integral_ratio << endl;

cout << endl;
cout << "---------------------------" << endl;

}


	 h_theta_e_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section*0.973/entries);//9.3792854e+08,7.9729075e+08
         h_theta_mu_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section*0.973/entries);
	 h_opening_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section*0.973/entries);

cout << "Ndata/NMC after scaling " << h_theta_e_RD->Integral()/h_theta_e_MC->Integral() << endl;



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
h->SaveAs(Form("comparison_RDMC/%s_%s_ratio_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_opening_%dhit_%dhit.root",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));
a.SaveAs(Form("comparison_RDMC/%s_%s_opening_reb_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%dhit.pdf",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));

/*
for(int n=0; n < h->GetNbinsX(); n++){
cout << n << ") RD: 1./sqrt(N) " << 1./sqrt(h_opening_RD->GetBinContent(n)) << " binError " << 1./h_opening_RD->GetBinError(n) << endl;
cout << n << ") MC: 1./sqrt(N) " << 1./sqrt(h_opening_MC->GetBinContent(n)) << " binError " << 1./h_opening_MC->GetBinError(n) << endl;
cout << n << ") GetBinError " << h->GetBinError(n) << endl;
cout << n << ") manual " <<sqrt( (rd_)/(mc_*mc_)+(rd_*rd_)/(mc_*mc_*mc_) ) << endl;
cout << "perc " << h->GetBinError(n)/h->GetBinContent(n)*100 << "%" <<endl;
}
*/

TCanvas a1("a1","a1",1000,700);
a1.Divide(2,2);
a1.cd(1);
/*
h_theta_e_MC->GetXaxis()->SetRangeUser(0.003,0.01);
h_theta_e_RD->GetXaxis()->SetRangeUser(0.003,0.01);
*/

/*h_theta_e_MC->GetXaxis()->SetRangeUser(0.003,0.03);
h_theta_e_RD->GetXaxis()->SetRangeUser(0.003,0.03);*/

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
/*double x3[h3->GetNbinsX()];
double y3[h3->GetNbinsX()];
double ex3[h3->GetNbinsX()];
double ey3[h3->GetNbinsX()];
double x_sum3[h3->GetNbinsX()];
double y_sum3[h3->GetNbinsX()];
double ex_sum3[h3->GetNbinsX()];
double ey_sum3[h3->GetNbinsX()];

double rd3 =0.;
double mc3 =0.;

for(int n=0; n < h3->GetNbinsX(); n++){ x3[n]=h3->GetBinCenter(n); y3[n]=h3->GetBinContent(n); ex3[n]=0.; ey3[n]=1./sqrt(h_theta_e_RD->GetBinContent(n))*h3->GetBinContent(n);}//h->GetBinError(n);}
TGraphErrors *op_er3=new TGraphErrors(h3->GetNbinsX(),x3,y3,ex3,ey3);

for(int n=0; n < h3->GetNbinsX(); n++){
rd3 = h_theta_e_RD->GetBinContent(n);
mc3 = h_theta_e_MC->GetBinContent(n);
 x_sum3[n]=h3->GetBinCenter(n); y_sum3[n]=h3->GetBinContent(n); ex_sum3[n]=0.; ey3[n]=sqrt( (rd3)/(mc3*mc3)+(rd3*rd3)/(mc3*mc3*mc3) );}//h->GetBinError(n);}
TGraphErrors *op_er_sum3=new TGraphErrors(h3->GetNbinsX(),x3,y3,ex3,ey3);

op_er3->SetMaximum(1.);
op_er3->SetMinimum(0.);
op_er3->GetXaxis()->SetRangeUser(0.,0.032);
op_er3->Draw("AP");
op_er_sum3->SetMarkerColor(kOrange);
op_er_sum3->Draw("AP");
*/
h3->SetMinimum(0.);
h3->GetXaxis()->SetTitle("Electron angle [rad]");
h3->Draw("E");
h3->SaveAs(Form("comparison_RDMC/%s_%s_ratio_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_theta_el_%dhit_%dhit.root",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));
a1.cd(3);
//h_theta_e_MC->Scale(7.9729075e+08*5.5*1E+23*d_tar*1E-30*315.44638/6.79333e+07);//26563993,605675    137720.00    1.3242886e+09   6.79333e+07 9.3792854e+08
h_theta_mu_MC->SetMinimum(1.);
/*h_theta_mu_MC->GetXaxis()->SetRangeUser(0.0005,0.002);
h_theta_mu_RD->GetXaxis()->SetRangeUser(0.0005,0.002);*/

/*h_theta_mu_MC->GetXaxis()->SetRangeUser(0.0003,0.002);
h_theta_mu_RD->GetXaxis()->SetRangeUser(0.0003,0.002);*/

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
h4->SaveAs(Form("comparison_RDMC/%s_%s_ratio_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_theta_mu_%dhit_%dhit.root",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));
a1.SaveAs(Form("comparison_RDMC/%s_%s_theta_MC_reco_RD_add_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%dhit.pdf",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));

TF1 *Elastic = new TF1("Elastic","0.5109989461*0.001*((1+(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*cos(x)*cos(x))/(1-(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*cos(x)*cos(x)))",0,0.032); 
TF1 *Elastic2 = new TF1("Elastic2","asin( (sin(x)*sqrt(Elastic(x)*Elastic(x)-0.5109989461*0.001*0.5109989461*0.001))/sqrt( (160+0.5109989461*0.001-Elastic(x))*(160+0.5109989461*0.001-Elastic(x))-105.6583745 *0.001*105.6583745 *0.001 ) )",0,0.032);


}
