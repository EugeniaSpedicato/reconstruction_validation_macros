#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TStyle.h>
#include "TGraphErrors.h"



void sec_single_divide_sig(int nhits_rd,int nhits_mc, string rd, string mc){

TH1::SetDefaultSumw2(kTRUE);
int const NBINS = 6;


/*
TFile *f_RD=TFile::Open(Form("comparison_RDMC/%s_d_aco_RD_parallel_pre_basic_restrict_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_thmu_RD=TFile::Open(Form("comparison_RDMC/%s_theta_mu_RD_parallel_pre_basic_restrict_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_the_RD=TFile::Open(Form("comparison_RDMC/%s_theta_e_RD_parallel_pre_basic_restrict_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_2D_RD=TFile::Open(Form("comparison_RDMC/%s_2D_RD_parallel_pre_basic_restrict_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_op_RD=TFile::Open(Form("comparison_RDMC/%s_opening_RD_parallel_pre_basic_restrict_el_Bend_%dhit.root",rd.c_str(),nhits_rd));

TFile* f_the_MC=TFile::Open(Form("comparison_RDMC/%s_theta_e_MC_parallel_sig_basic_restrict_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
TFile* f_thmu_MC=TFile::Open(Form("comparison_RDMC/%s_theta_mu_MC_parallel_sig_basic_restrict_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
TFile* f_2D_MC=TFile::Open(Form("comparison_RDMC/%s_2D_MC_parallel_sig_basic_restrict_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
TFile* f_op_MC=TFile::Open(Form("comparison_RDMC/%s_opening_MC_parallel_sig_basic_restrict_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
TFile* f_aco_MC=TFile::Open(Form("comparison_RDMC/%s_d_aco_MC_parallel_sig_basic_restrict_el_Bend_%dhit.root",mc.c_str(),nhits_mc)); 

TFile *f_d0_MC=TFile::Open(Form("comparison_RDMC/%s_distance_MOD0_sig_basic_restrict_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_d1_MC=TFile::Open(Form("comparison_RDMC/%s_distance_MOD1_sig_basic_restrict_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
TFile *f_d2_MC=TFile::Open(Form("comparison_RDMC/%s_distance_MOD2_sig_basic_restrict_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
TFile *f_d3_MC=TFile::Open(Form("comparison_RDMC/%s_distance_MOD3_sig_basic_restrict_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
TFile *f_d4_MC=TFile::Open(Form("comparison_RDMC/%s_distance_MOD4_sig_basic_restrict_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
TFile *f_d5_MC=TFile::Open(Form("comparison_RDMC/%s_distance_MOD5_sig_basic_restrict_el_Bend_%dhit.root",mc.c_str(),nhits_mc));

TFile *f_d0_RD=TFile::Open(Form("comparison_RDMC/%s_distance_MOD0_basic_restrict_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_d1_RD=TFile::Open(Form("comparison_RDMC/%s_distance_MOD1_basic_restrict_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_d2_RD=TFile::Open(Form("comparison_RDMC/%s_distance_MOD2_basic_restrict_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_d3_RD=TFile::Open(Form("comparison_RDMC/%s_distance_MOD3_basic_restrict_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_d4_RD=TFile::Open(Form("comparison_RDMC/%s_distance_MOD4_basic_restrict_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_d5_RD=TFile::Open(Form("comparison_RDMC/%s_distance_MOD5_basic_restrict_el_Bend_%dhit.root",rd.c_str(),nhits_rd));

TFile *f_mod0_r0_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_0_sig_basic_restrict_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
TFile *f_mod0_r1_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_1_sig_basic_restrict_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
TFile *f_mod0_r2_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_2_sig_basic_restrict_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
TFile *f_mod0_r3_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_3_sig_basic_restrict_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
TFile *f_mod0_r4_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_4_sig_basic_restrict_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
//TFile *f_mod0_r5_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_5_sig_basic_restrict_el_Bend_%dhit.root",mc.c_str(),nhits_mc));

TFile *f_mod0_r0_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_0_basic_restrict_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_mod0_r1_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_1_basic_restrict_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_mod0_r2_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_2_basic_restrict_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_mod0_r3_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_3_basic_restrict_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_mod0_r4_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_4_basic_restrict_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
//TFile *f_mod0_r5_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_5_basic_restrict_el_Bend_%dhit.root",rd.c_str(),nhits_rd));


TFile *f_mod2_r0_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_0_sig_basic_restrict_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
TFile *f_mod2_r1_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_1_sig_basic_restrict_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
TFile *f_mod2_r2_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_2_sig_basic_restrict_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
TFile *f_mod2_r3_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_3_sig_basic_restrict_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
TFile *f_mod2_r4_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_4_sig_basic_restrict_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
//TFile *f_mod2_r5_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_5_sig_basic_restrict_el_Bend_%dhit.root",mc.c_str(),nhits_mc));

TFile *f_mod2_r0_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_0_basic_restrict_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_mod2_r1_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_1_basic_restrict_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_mod2_r2_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_2_basic_restrict_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_mod2_r3_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_3_basic_restrict_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_mod2_r4_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_4_basic_restrict_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
//TFile *f_mod2_r5_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_5_basic_restrict_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
*/



TFile *f_RD=TFile::Open(Form("comparison_RDMC/%s_d_aco_RD_parallel_pre_basic_restrict_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_thmu_RD=TFile::Open(Form("comparison_RDMC/%s_theta_mu_RD_parallel_pre_basic_restrict_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_the_RD=TFile::Open(Form("comparison_RDMC/%s_theta_e_RD_parallel_pre_basic_restrict_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_2D_RD=TFile::Open(Form("comparison_RDMC/%s_2D_RD_parallel_pre_basic_restrict_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_op_RD=TFile::Open(Form("comparison_RDMC/%s_opening_RD_parallel_pre_basic_restrict_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));

TFile* f_the_MC=TFile::Open(Form("comparison_RDMC/%s_theta_e_RD_parallel_pre_basic_restrict_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile* f_thmu_MC=TFile::Open(Form("comparison_RDMC/%s_theta_mu_RD_parallel_pre_basic_restrict_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile* f_2D_MC=TFile::Open(Form("comparison_RDMC/%s_2D_RD_parallel_pre_basic_restrict_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile* f_op_MC=TFile::Open(Form("comparison_RDMC/%s_opening_RD_parallel_pre_basic_restrict_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile* f_aco_MC=TFile::Open(Form("comparison_RDMC/%s_d_aco_RD_parallel_pre_basic_restrict_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc)); 

TFile *f_d0_MC=TFile::Open(Form("comparison_RDMC/%s_distance_MOD0_basic_restrict_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile *f_d1_MC=TFile::Open(Form("comparison_RDMC/%s_distance_MOD1_basic_restrict_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile *f_d2_MC=TFile::Open(Form("comparison_RDMC/%s_distance_MOD2_basic_restrict_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile *f_d3_MC=TFile::Open(Form("comparison_RDMC/%s_distance_MOD3_basic_restrict_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile *f_d4_MC=TFile::Open(Form("comparison_RDMC/%s_distance_MOD4_basic_restrict_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile *f_d5_MC=TFile::Open(Form("comparison_RDMC/%s_distance_MOD5_basic_restrict_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));

TFile *f_d0_RD=TFile::Open(Form("comparison_RDMC/%s_distance_MOD0_basic_restrict_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_d1_RD=TFile::Open(Form("comparison_RDMC/%s_distance_MOD1_basic_restrict_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_d2_RD=TFile::Open(Form("comparison_RDMC/%s_distance_MOD2_basic_restrict_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_d3_RD=TFile::Open(Form("comparison_RDMC/%s_distance_MOD3_basic_restrict_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_d4_RD=TFile::Open(Form("comparison_RDMC/%s_distance_MOD4_basic_restrict_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_d5_RD=TFile::Open(Form("comparison_RDMC/%s_distance_MOD5_basic_restrict_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));

TFile *f_mod0_r0_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_0_basic_restrict_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile *f_mod0_r1_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_1_basic_restrict_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile *f_mod0_r2_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_2_basic_restrict_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile *f_mod0_r3_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_3_basic_restrict_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile *f_mod0_r4_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_4_basic_restrict_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
//TFile *f_mod0_r5_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_5_basic_restrict_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));

TFile *f_mod0_r0_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_0_basic_restrict_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_mod0_r1_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_1_basic_restrict_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_mod0_r2_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_2_basic_restrict_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_mod0_r3_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_3_basic_restrict_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_mod0_r4_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_4_basic_restrict_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
//TFile *f_mod0_r5_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_5_basic_restrict_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));


TFile *f_mod2_r0_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_0_basic_restrict_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile *f_mod2_r1_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_1_basic_restrict_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile *f_mod2_r2_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_2_basic_restrict_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile *f_mod2_r3_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_3_basic_restrict_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile *f_mod2_r4_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_4_basic_restrict_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
//TFile *f_mod2_r5_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_5_basic_restrict_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));

TFile *f_mod2_r0_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_0_basic_restrict_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_mod2_r1_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_1_basic_restrict_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_mod2_r2_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_2_basic_restrict_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_mod2_r3_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_3_basic_restrict_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_mod2_r4_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_4_basic_restrict_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
//TFile *f_mod2_r5_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_5_basic_restrict_el_Bend_%dhit.root",rd.c_str(),nhits_rd));



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

TH1D *h_mod0_r0_MC=(TH1D*)f_mod0_r0_MC->Get("h_dist_mod0_range0");
TH1D *h_mod0_r1_MC=(TH1D*)f_mod0_r1_MC->Get("h_dist_mod0_range1");
TH1D *h_mod0_r2_MC=(TH1D*)f_mod0_r2_MC->Get("h_dist_mod0_range2");
//TH1D *h_mod0_r3_MC=(TH1D*)f_mod0_r3_MC->Get("h_dist_mod0_range3");
//TH1D *h_mod0_r4_MC=(TH1D*)f_mod0_r4_MC->Get("h_dist_mod0_range4");
//TH1D *h_mod0_r5_MC=(TH1D*)f_mod0_r5_MC->Get("h_dist_mod0_range5");

TH1D *h_mod0_r0_RD=(TH1D*)f_mod0_r0_RD->Get("h_dist_mod0_range0");
TH1D *h_mod0_r1_RD=(TH1D*)f_mod0_r1_RD->Get("h_dist_mod0_range1");
TH1D *h_mod0_r2_RD=(TH1D*)f_mod0_r2_RD->Get("h_dist_mod0_range2");
//TH1D *h_mod0_r3_RD=(TH1D*)f_mod0_r3_RD->Get("h_dist_mod0_range3");
//TH1D *h_mod0_r4_RD=(TH1D*)f_mod0_r4_RD->Get("h_dist_mod0_range4");
//TH1D *h_mod0_r5_RD=(TH1D*)f_mod0_r5_RD->Get("h_dist_mod0_range5");

TH1D *h_mod2_r0_MC=(TH1D*)f_mod2_r0_MC->Get("h_dist_mod2_range0");
TH1D *h_mod2_r1_MC=(TH1D*)f_mod2_r1_MC->Get("h_dist_mod2_range1");
TH1D *h_mod2_r2_MC=(TH1D*)f_mod2_r2_MC->Get("h_dist_mod2_range2");
//TH1D *h_mod2_r3_MC=(TH1D*)f_mod2_r3_MC->Get("h_dist_mod2_range3");
//TH1D *h_mod2_r4_MC=(TH1D*)f_mod2_r4_MC->Get("h_dist_mod2_range4");
//TH1D *h_mod2_r5_MC=(TH1D*)f_mod2_r5_MC->Get("h_dist_mod2_range5");

TH1D *h_mod2_r0_RD=(TH1D*)f_mod2_r0_RD->Get("h_dist_mod2_range0");
TH1D *h_mod2_r1_RD=(TH1D*)f_mod2_r1_RD->Get("h_dist_mod2_range1");
TH1D *h_mod2_r2_RD=(TH1D*)f_mod2_r2_RD->Get("h_dist_mod2_range2");
//TH1D *h_mod2_r3_RD=(TH1D*)f_mod2_r3_RD->Get("h_dist_mod2_range3");
//TH1D *h_mod2_r4_RD=(TH1D*)f_mod2_r4_RD->Get("h_dist_mod2_range4");
//TH1D *h_mod2_r5_RD=(TH1D*)f_mod2_r5_RD->Get("h_dist_mod2_range5");


/*h_theta_e_RD->SetBins(8,0.,0.032);
for(int m=0; m<NBINS; m++){h_theta_e_MC->SetBins(8,0.,0.032);}
*/
//h_theta_mu_RD->SetBins(9,0.0002,0.003);
//h_theta_mu_MC->SetBins(9,0.0002,0.003);

   Double_t cross_section = 0.;
   Double_t entries = 0.;
   Double_t n_mu_ot = 0.;
   Double_t d_tar = 0.;
if(mc=="default"){
    cross_section = 1341.02;//wip master
    entries = 1.96203e+06;
    n_mu_ot=7.9729075e+08;
    d_tar=3.;}
else if(mc=="wip11"){
    cross_section =1340.59787016883;//wip11
    entries =322382;
    n_mu_ot=7.9729075e+08;
    d_tar=3.;}
if(mc=="25micron"){
    cross_section = 1341.02;//wip master
    entries = 676509;
    n_mu_ot=7.9729075e+08;
    d_tar=3.;}
if(mc=="run17"){
    cross_section = 1340.26661292368;
    entries = 973336;
    n_mu_ot=6.1317676e+08;
    d_tar=2.;}
if(mc=="default_chi2out50"){
    cross_section = 1341.6111;
    entries = 2.61979e+06;
    n_mu_ot=7.9729075e+08;
    d_tar=3.;}
if(mc=="default_bb35b5de"){
    cross_section = 1341.6111;
    entries = 2.61979e+06;
    n_mu_ot=7.9729075e+08;
    d_tar=3.;}
/*else if(mc=="ricmis"){
    cross_section = 1343.8760; //ricmis
    entries =; //ricmis
}*/


//   Double_t entries =651153;


//   Double_t cross_section = 12.9221;
//   Double_t entries = 4.43846e+07;

for(int n=0; n < h_opening_MC->GetNbinsX(); n++){cout << n << ") MC: pre " << 1./sqrt(h_opening_MC->GetBinContent(n)) << " binError " << 1./h_opening_MC->GetBinError(n) << endl;}

	 h_theta_e_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section/entries);//9.3792854e+08,7.9729075e+08
         h_theta_mu_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section/entries);
	 h_2D_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section/entries);
	 h_opening_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section/entries);
         h_aco_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section/entries);
h_d0_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section/entries);
h_d1_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section/entries);
h_d2_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section/entries);
h_d3_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section/entries);
h_d4_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section/entries);
h_d5_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section/entries);

h_mod0_r0_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section/entries);
h_mod0_r1_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section/entries);
h_mod0_r2_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section/entries);

h_mod2_r0_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section/entries);
h_mod2_r1_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section/entries);
h_mod2_r2_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section/entries);

//h_mod0_r3_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section/entries);
//h_mod0_r4_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section/entries);
//h_mod0_r5_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section/entries);

//h_mod2_r3_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section/entries);
//h_mod2_r4_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section/entries);
//h_mod2_r5_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section/entries);


TCanvas dr("dr","dr",1400,2100);
dr.Divide(2,3);
dr.cd(1);
h_mod0_r0_MC->Draw("hist");
h_mod0_r0_RD->SetLineColor(kPink);
h_mod0_r0_RD->Draw("hist same");
dr.cd(2);
h_mod0_r0_RD->SetLineColor(kPink);
h_mod0_r0_RD->Draw("hist");
h_mod0_r1_MC->Draw("hist");
h_mod0_r1_RD->SetLineColor(kPink);
h_mod0_r1_RD->Draw("hist same");
dr.cd(3);
h_mod0_r2_MC->Draw("hist");
h_mod0_r2_RD->SetLineColor(kPink);
h_mod0_r2_RD->Draw("hist same");

/*dr.cd(4);
h_mod0_r3_MC->Draw("hist");
h_mod0_r3_RD->SetLineColor(kPink);
h_mod0_r3_RD->Draw("hist same");
dr.cd(5);
h_mod0_r4_MC->Draw("hist");
h_mod0_r4_RD->SetLineColor(kPink);
h_mod0_r4_RD->Draw("hist same");
*/dr.SaveAs(Form("comparison_RDMC/%s_%s_drange_mod0_basic_restrict_el_Bend_%dhit_%dhit.pdf",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));

TCanvas dr2("dr2","dr2",1400,2100);
dr2.Divide(2,3);
dr2.cd(1);
h_mod2_r0_MC->Draw("hist");
h_mod2_r0_RD->SetLineColor(kPink);
h_mod2_r0_RD->Draw("hist same");
dr2.cd(2);
h_mod2_r1_MC->Draw("hist");
h_mod2_r1_RD->SetLineColor(kPink);
h_mod2_r1_RD->Draw("hist same");
dr2.cd(3);
h_mod2_r2_MC->Draw("hist");
h_mod2_r2_RD->SetLineColor(kPink);
h_mod2_r2_RD->Draw("hist same");

/*dr2.cd(4);
h_mod2_r3_MC->Draw("hist");
h_mod2_r3_RD->SetLineColor(kPink);
h_mod2_r3_RD->Draw("hist same");
dr2.cd(5);
h_mod2_r4_MC->Draw("hist");
h_mod2_r4_RD->SetLineColor(kPink);
h_mod2_r4_RD->Draw("hist same");
*/
dr2.SaveAs(Form("comparison_RDMC/%s_%s_drange_mod2_basic_restrict_el_Bend_%dhit_%dhit.pdf",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));

TCanvas di("di","di",1400,2100);
di.Divide(2,3);
di.cd(1);
h_d0_MC->GetXaxis()->SetRangeUser(-1.,1.);
h_d0_MC->SetMinimum(0.);
h_d0_MC->Draw("hist");
h_d0_RD->SetLineColor(kPink);
h_d0_RD->Draw("hist same");
di.cd(2);
h_d1_MC->GetXaxis()->SetRangeUser(-1.,1.);
h_d1_MC->SetMinimum(0.);
h_d1_MC->Draw("hist");
h_d1_RD->SetLineColor(kPink);
h_d1_RD->Draw("hist same");

di.cd(3);
h_d2_MC->GetXaxis()->SetRangeUser(-1.,1.);
h_d2_MC->SetMinimum(0.);
h_d2_MC->Draw("hist");
h_d2_RD->SetLineColor(kPink);
h_d2_RD->Draw("hist same");
di.cd(4);
h_d3_MC->GetXaxis()->SetRangeUser(-1.,1.);
h_d3_MC->SetMinimum(0.);
h_d3_MC->Draw("hist");
h_d3_RD->SetLineColor(kPink);
h_d3_RD->Draw("hist same");
di.cd(5);
h_d4_MC->SetMinimum(0.);
h_d4_MC->Draw("hist");
h_d4_RD->SetLineColor(kPink);
h_d4_RD->Draw("hist same");
di.cd(6);
h_d5_MC->SetMinimum(0.);
h_d5_MC->Draw("hist");
h_d5_RD->SetLineColor(kPink);
h_d5_RD->Draw("hist same");
di.SaveAs(Form("comparison_RDMC/%s_%s_h_dist_basic_restrict_el_Bend_%dhit_%dhit.pdf",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));

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
aco.SaveAs(Form("comparison_RDMC/%s_%s_h_aco_basic_restrict_el_Bend_%dhit_%dhit.pdf",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));


TCanvas a("a","a",700,700);
a.Divide(1,2);
a.cd(1);
h_opening_MC->SetMinimum(1.);
h_opening_MC->GetXaxis()->SetRangeUser(0.004,0.03);
h_opening_RD->GetXaxis()->SetRangeUser(0.004,0.03);
h_opening_RD->SetMinimum(1.);
h_opening_MC->Draw("hist");
h_opening_RD->SetLineColor(kPink+10);
h_opening_RD->Draw("hist same");
a.cd(2);

TH1D * h = (TH1D*) h_opening_RD->Clone();
h->Divide(h_opening_MC);
//h->SetMaximum(3.);
h->SetMinimum(0.);
/*
double x[h->GetNbinsX()];
double y[h->GetNbinsX()];
double ex[h->GetNbinsX()];
double ey[h->GetNbinsX()];
double x_sum[h->GetNbinsX()];
double y_sum[h->GetNbinsX()];
double ex_sum[h->GetNbinsX()];
double ey_sum[h->GetNbinsX()];

double rd =0.;
double mc =0.;

for(int n=0; n < h->GetNbinsX(); n++){ x[n]=h->GetBinCenter(n); y[n]=h->GetBinContent(n); ex[n]=0.; ey[n]=1./sqrt(h_opening_RD->GetBinContent(n))*h->GetBinContent(n);}//h->GetBinError(n);}
TGraphErrors *op_er=new TGraphErrors(h->GetNbinsX(),x,y,ex,ey);

for(int n=0; n < h->GetNbinsX(); n++){
rd = h_opening_RD->GetBinContent(n);
mc = h_opening_MC->GetBinContent(n);
 x_sum[n]=h->GetBinCenter(n); y_sum[n]=h->GetBinContent(n); ex_sum[n]=0.; ey[n]=sqrt( (rd)/(mc*mc)+(rd*rd)/(mc*mc*mc) );}//h->GetBinError(n);}
TGraphErrors *op_er_sum=new TGraphErrors(h->GetNbinsX(),x,y,ex,ey);

h->SetMaximum(1.);
op_er_sum->SetMaximum(1.);
op_er->SetMaximum(1.);
op_er->SetMinimum(0.);
op_er->GetXaxis()->SetRangeUser(0.,0.032);
//op_er->Draw("AP");
op_er_sum->SetMarkerColor(kOrange);
//op_er_sum->Draw("AP");
*/
h->GetXaxis()->SetTitle("opening #theta(rad)");
h->Draw("E");
gStyle->SetOptStat(0);
h->SaveAs(Form("comparison_RDMC/%s_%s_ratio_opening_%dhit_%dhit.root",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));
a.SaveAs(Form("comparison_RDMC/%s_%s_opening_reb_basic_restrict_el_Bend_%dhit_%dhit.pdf",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));

/*for(int n=0; n < h->GetNbinsX(); n++){
cout << n << ") RD: 1./sqrt(N) " << 1./sqrt(h_opening_RD->GetBinContent(n)) << " binError " << 1./h_opening_RD->GetBinError(n) << endl;
cout << n << ") MC: 1./sqrt(N) " << 1./sqrt(h_opening_MC->GetBinContent(n)) << " binError " << 1./h_opening_MC->GetBinError(n) << endl;
cout << n << ") GetBinError " << h->GetBinError(n) << endl;
cout << n << ") manual " <<sqrt( (rd)/(mc*mc)+(rd*rd)/(mc*mc*mc) ) << endl;
cout << "perc " << h->GetBinError(n)/h->GetBinContent(n)*100 << "%" <<endl;
}
*/

TCanvas a1("a1","a1",1000,700);
a1.Divide(2,2);
a1.cd(1);
cout << "theta el. reco:" << endl;
/*
h_theta_e_MC->GetXaxis()->SetRangeUser(0.003,0.01);
h_theta_e_RD->GetXaxis()->SetRangeUser(0.003,0.01);
*/
h_theta_e_MC->GetXaxis()->SetRangeUser(0.003,0.03);
h_theta_e_RD->GetXaxis()->SetRangeUser(0.003,0.03);
h_theta_e_MC->SetMinimum(1.);
h_theta_e_MC->Draw("E");
h_theta_e_RD->SetLineColor(kPink+10);
h_theta_e_RD->Draw("E same");
gPad->SetLogy();
gStyle->SetOptStat(0);
a1.cd(2);

TH1D * h3 = (TH1D*) h_theta_e_RD->Clone();
h3->Divide(h_theta_e_MC);
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
h3->GetXaxis()->SetTitle("#theta_e(rad)");
h3->Draw("E");
h3->SaveAs(Form("comparison_RDMC/%s_%s_ratio_theta_el_%dhit_%dhit.root",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));
a1.cd(3);
cout << "theta mu. reco:" << endl;
//h_theta_e_MC->Scale(7.9729075e+08*5.5*1E+23*d_tar*1E-30*315.44638/6.79333e+07);//26563993,605675    137720.00    1.3242886e+09   6.79333e+07 9.3792854e+08
h_theta_mu_MC->SetMinimum(1.);
/*h_theta_mu_MC->GetXaxis()->SetRangeUser(0.0005,0.002);
h_theta_mu_RD->GetXaxis()->SetRangeUser(0.0005,0.002);*/
h_theta_mu_MC->GetXaxis()->SetRangeUser(0.0003,0.002);
h_theta_mu_RD->GetXaxis()->SetRangeUser(0.0003,0.002);
h_theta_mu_MC->Draw("E");
h_theta_mu_RD->SetLineColor(kPink+10);
h_theta_mu_RD->Draw("E same");
gStyle->SetOptStat(0);
gPad->SetLogy();

a1.cd(4);
TH1D * h4 = (TH1D*) h_theta_mu_RD->Clone();
h4->Divide(h_theta_mu_MC);

/*
double x4[h4->GetNbinsX()];
double y4[h4->GetNbinsX()];
double ex4[h4->GetNbinsX()];
double ey4[h4->GetNbinsX()];
double x_sum4[h4->GetNbinsX()];
double y_sum4[h4->GetNbinsX()];
double ex_sum4[h4->GetNbinsX()];
double ey_sum4[h4->GetNbinsX()];

double rd4 =0.;
double mc4 =0.;

for(int n=0; n < h4->GetNbinsX(); n++){ x4[n]=h4->GetBinCenter(n); y4[n]=h4->GetBinContent(n); ex4[n]=0.; ey4[n]=1./sqrt(h_theta_mu_RD->GetBinContent(n))*h4->GetBinContent(n);}//h->GetBinError(n);}
TGraphErrors *op_er4=new TGraphErrors(h4->GetNbinsX(),x4,y4,ex4,ey4);

for(int n=0; n < h4->GetNbinsX(); n++){
rd4 = h_theta_mu_RD->GetBinContent(n);
mc4 = h_theta_mu_MC->GetBinContent(n);
 x_sum4[n]=h4->GetBinCenter(n); y_sum4[n]=h4->GetBinContent(n); ex_sum4[n]=0.; ey4[n]=sqrt( (rd4)/(mc4*mc4)+(rd4*rd4)/(mc4*mc4*mc4) );}//h->GetBinError(n);}
TGraphErrors *op_er_sum4=new TGraphErrors(h4->GetNbinsX(),x4,y4,ex4,ey4);

op_er4->SetMaximum(1.);
op_er4->SetMinimum(0.);
//op_er4->GetXaxis()->SetRangeUser(0.,0.032);
op_er4->Draw("AP");
op_er_sum4->SetMarkerColor(kOrange);
op_er_sum4->Draw("AP");


*/
h4->SetMinimum(0.);
h4->GetXaxis()->SetTitle("#theta_mu(rad)");
h4->Draw("E");
h4->SaveAs(Form("comparison_RDMC/%s_%s_ratio_theta_mu_%dhit_%dhit.root",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));
a1.SaveAs(Form("comparison_RDMC/%s_%s_theta_MC_reco_RD_add_basic_restrict_el_Bend_%dhit_%dhit.pdf",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));

TF1 *Elastic = new TF1("Elastic","0.5109989461*0.001*((1+(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*cos(x)*cos(x))/(1-(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*cos(x)*cos(x)))",0,0.030); 
TF1 *Elastic2 = new TF1("Elastic2","asin( (sin(x)*sqrt(Elastic(x)*Elastic(x)-0.5109989461*0.001*0.5109989461*0.001))/sqrt( (160+0.5109989461*0.001-Elastic(x))*(160+0.5109989461*0.001-Elastic(x))-105.6583745 *0.001*105.6583745 *0.001 ) )",0,0.030);


TCanvas a2("a2","a2",6000,3000);
a2.Divide(2,1);

a2.cd(1);
h_2D_MC->SetTitle("#theta_el VS #theta_mu MC");

h_2D_MC->GetXaxis()->SetTitle("#theta_el (rad)");
h_2D_MC->GetYaxis()->SetTitle("#theta_mu (rad)");
//h_2D_MC->RebinY(2);
//h_2D_MC->RebinX(2);
h_2D_MC->Draw("COLZ");
Elastic2->Draw("same");
a2.cd(2);
h_2D_RD->SetTitle("#theta_el VS #theta_mu RD");
h_2D_RD->GetXaxis()->SetTitle("#theta_el (rad)");
h_2D_RD->GetYaxis()->SetTitle("#theta_mu (rad)");
//h_2D_RD->RebinX(2);
//h_2D_RD->RebinY(2);
h_2D_RD->Draw("COLZ");
Elastic2->Draw("same");
TGaxis::SetMaxDigits(3);
/*a2.cd(3);


TH2D * h5 = (TH2D*) h_2D_RD->Clone();
h5->Divide(h_2D_MC);


Int_t nx = h5->GetNbinsX();
Int_t ny = h5->GetNbinsY();
for (Int_t i=1; i<nx+1; i++) {
for (Int_t j=1; j<ny+1; j++) {
if (h5->GetBinError(i,j)/h5->GetBinContent(i,j)*100>=5) h5->SetBinContent(i,j,0);}}
*/

/*
for(int a=1; a<h5->GetNbinsX(); a++){
 for(int b=1; b<h5->GetNbinsY(); b++){
	cout << a<<"," <<b << " ) " << h5->GetBinContent(a,b) << " +- " << h5->GetBinError(a,b)<< " ----> " << h5->GetBinError(a,b)/h5->GetBinContent(a,b)*100 << "%" << endl;
 }
}*/

  //Int_t colors[] = {0, 1, 2, 3, 4, 5, 6}; // #colors >= #levels - 1
  //gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);



/*  Int_t n=13;
  Int_t *colors = new Int_t[n];
//  for (Int_t i = 0; i < n; i++){if(i<5) colors[i] = kRed+i; else if(i>=5 and i<10) colors[i] = kOrange+i-5; else colors[i] = kYellow+i;}
  colors[0] = kRed-10;
  colors[1] = kRed-10;
  colors[2] = kRed-10;
  colors[3] = kOrange+6;
  colors[4] = kOrange+6;
  colors[5] = kOrange+6;
  colors[6] = kOrange+10;
  colors[7] = kOrange+10;
  colors[8] = kOrange+10;
  colors[9] = kViolet-9;
  colors[10] = kViolet-9;
  colors[11] = kAzure+6;
  colors[12] = kAzure;

  gStyle->SetPalette(n, colors);

  // #levels <= #colors + 1 (notes: +-3.4e38 = +-FLT_MAX; +1.17e-38 = +FLT_MIN)
  Double_t levels[] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.00, 10.0, 10., 30., 60};

h5->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
// gStyle->SetPalette(kIsland);
h5->SetTitle("#theta el -vs- #theta mu ratio RD/MC");
gStyle->SetPaintTextFormat(".2f");
h5->Draw("COLZ TEXT");
*/a2.SaveAs(Form("comparison_RDMC/%s_%s_comparison_2D_basic_restrict_el_Bend_%dhit_%dhit.pdf",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));

/*TCanvas a3("a3","a3",1400,700);
a3.Divide(2,1);
a3.cd(1);
auto p1 =h5->ProfileX();
p1->SetTitle("Efficiency as a function of #theta el");
p1->Draw();
a3.cd(2);
auto p2 =h5->ProfileY();
p2->SetTitle("Efficiency as a function of #theta mu");
p2->Draw();
a3.SaveAs(Form("comparison_RDMC/%s_%s_profile1hit_basic_restrict_el_Bend_%dhit_%dhit.pdf",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));
*/
TCanvas t("t","t",3000,1000);
t.Divide(3,1);
t.cd(1);
h3->SetLineColor(kRed+1);
h3->Draw();
t.cd(2);
h4->SetLineColor(kOrange-3);
h4->Draw();
t.cd(3);
h->SetLineColor(kGreen-2);
h->Draw();

t.SaveAs(Form("comparison_RDMC/%s_%s_tesi_comparison_basic_restrict_%dhit_%dhit.pdf",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));
}

