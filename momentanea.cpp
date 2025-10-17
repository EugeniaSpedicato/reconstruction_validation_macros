#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TStyle.h>
#include "TGraphErrors.h"



void momentanea(int nhits_rd,int nhits_mc, string rd, string mc){

TH1::SetDefaultSumw2(kTRUE);
int const NBINS = 6;

//chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak

/*
TFile *f_RD=TFile::Open(Form("comparison_RDMC/%s_d_aco_RD_parallel_pre_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_thmu_RD=TFile::Open(Form("comparison_RDMC/%s_theta_mu_RD_parallel_pre_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_the_RD=TFile::Open(Form("comparison_RDMC/%s_theta_e_RD_parallel_pre_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_2D_RD=TFile::Open(Form("comparison_RDMC/%s_2D_RD_parallel_pre_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_op_RD=TFile::Open(Form("comparison_RDMC/%s_opening_RD_parallel_pre_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",rd.c_str(),nhits_rd));

TFile* f_the_MC=TFile::Open(Form("comparison_RDMC/%s_theta_e_MC_parallel_sig_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
TFile* f_thmu_MC=TFile::Open(Form("comparison_RDMC/%s_theta_mu_MC_parallel_sig_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
TFile* f_2D_MC=TFile::Open(Form("comparison_RDMC/%s_2D_MC_parallel_sig_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
TFile* f_op_MC=TFile::Open(Form("comparison_RDMC/%s_opening_MC_parallel_sig_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
TFile* f_aco_MC=TFile::Open(Form("comparison_RDMC/%s_d_aco_MC_parallel_sig_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",mc.c_str(),nhits_mc)); 

TFile *f_d0_MC=TFile::Open(Form("comparison_RDMC/%s_distance_MOD0_sig_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_d1_MC=TFile::Open(Form("comparison_RDMC/%s_distance_MOD1_sig_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
TFile *f_d2_MC=TFile::Open(Form("comparison_RDMC/%s_distance_MOD2_sig_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
TFile *f_d3_MC=TFile::Open(Form("comparison_RDMC/%s_distance_MOD3_sig_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
TFile *f_d4_MC=TFile::Open(Form("comparison_RDMC/%s_distance_MOD4_sig_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
TFile *f_d5_MC=TFile::Open(Form("comparison_RDMC/%s_distance_MOD5_sig_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",mc.c_str(),nhits_mc));

TFile *f_d0_RD=TFile::Open(Form("comparison_RDMC/%s_distance_MOD0_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_d1_RD=TFile::Open(Form("comparison_RDMC/%s_distance_MOD1_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_d2_RD=TFile::Open(Form("comparison_RDMC/%s_distance_MOD2_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_d3_RD=TFile::Open(Form("comparison_RDMC/%s_distance_MOD3_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_d4_RD=TFile::Open(Form("comparison_RDMC/%s_distance_MOD4_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_d5_RD=TFile::Open(Form("comparison_RDMC/%s_distance_MOD5_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",rd.c_str(),nhits_rd));

TFile *f_mod0_r0_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_0_sig_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
TFile *f_mod0_r1_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_1_sig_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
TFile *f_mod0_r2_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_2_sig_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
TFile *f_mod0_r3_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_3_sig_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
TFile *f_mod0_r4_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_4_sig_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
//TFile *f_mod0_r5_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_5_sig_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",mc.c_str(),nhits_mc));

TFile *f_mod0_r0_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_0_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_mod0_r1_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_1_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_mod0_r2_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_2_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_mod0_r3_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_3_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_mod0_r4_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_4_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
//TFile *f_mod0_r5_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_5_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",rd.c_str(),nhits_rd));


TFile *f_mod2_r0_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_0_sig_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
TFile *f_mod2_r1_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_1_sig_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
TFile *f_mod2_r2_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_2_sig_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
TFile *f_mod2_r3_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_3_sig_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
TFile *f_mod2_r4_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_4_sig_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",mc.c_str(),nhits_mc));
//TFile *f_mod2_r5_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_5_sig_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",mc.c_str(),nhits_mc));

TFile *f_mod2_r0_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_0_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_mod2_r1_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_1_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_mod2_r2_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_2_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_mod2_r3_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_3_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
TFile *f_mod2_r4_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_4_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
//TFile *f_mod2_r5_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_5_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
*/




TFile *f_d0_MC=TFile::Open(Form("comparison_RDMC/%s_distance_MOD0_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile *f_d1_MC=TFile::Open(Form("comparison_RDMC/%s_distance_MOD1_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile *f_d2_MC=TFile::Open(Form("comparison_RDMC/%s_distance_MOD2_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile *f_d3_MC=TFile::Open(Form("comparison_RDMC/%s_distance_MOD3_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile *f_d4_MC=TFile::Open(Form("comparison_RDMC/%s_distance_MOD4_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile *f_d5_MC=TFile::Open(Form("comparison_RDMC/%s_distance_MOD5_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));

TFile *f_d0_RD=TFile::Open(Form("comparison_RDMC/%s_distance_MOD0_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_d1_RD=TFile::Open(Form("comparison_RDMC/%s_distance_MOD1_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_d2_RD=TFile::Open(Form("comparison_RDMC/%s_distance_MOD2_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_d3_RD=TFile::Open(Form("comparison_RDMC/%s_distance_MOD3_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_d4_RD=TFile::Open(Form("comparison_RDMC/%s_distance_MOD4_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_d5_RD=TFile::Open(Form("comparison_RDMC/%s_distance_MOD5_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));

/*
TFile *f_mod0_r0_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_0_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile *f_mod0_r1_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_1_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile *f_mod0_r2_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_2_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile *f_mod0_r3_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_3_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile *f_mod0_r4_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_4_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
//TFile *f_mod0_r5_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_5_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));

TFile *f_mod0_r0_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_0_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_mod0_r1_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_1_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_mod0_r2_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_2_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_mod0_r3_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_3_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_mod0_r4_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_4_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
//TFile *f_mod0_r5_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod0_range_5_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));


TFile *f_mod2_r0_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_0_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile *f_mod2_r1_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_1_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile *f_mod2_r2_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_2_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile *f_mod2_r3_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_3_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
TFile *f_mod2_r4_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_4_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));
//TFile *f_mod2_r5_MC=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_5_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_mc.root",mc.c_str(),nhits_mc));

TFile *f_mod2_r0_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_0_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_mod2_r1_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_1_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_mod2_r2_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_2_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_mod2_r3_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_3_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
TFile *f_mod2_r4_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_4_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_data.root",rd.c_str(),nhits_rd));
//TFile *f_mod2_r5_RD=TFile::Open(Form("comparison_RDMC/%s_h_dist_mod2_range_5_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit.root",rd.c_str(),nhits_rd));
*/


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

/*TH1D *h_mod0_r0_MC=(TH1D*)f_mod0_r0_MC->Get("h_dist_mod0_range0");
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
*/


/*h_theta_e_RD->SetBins(8,0.,0.032);
for(int m=0; m<NBINS; m++){h_theta_e_MC->SetBins(8,0.,0.032);}
*/
//h_theta_mu_RD->SetBins(9,0.0002,0.003);
//h_theta_mu_MC->SetBins(9,0.0002,0.003);

   Double_t cross_section = 0.;
   Double_t entries = 0.;
   Double_t w2_fiducial = 0.;
   Double_t n_mu_ot = 0.;
   Double_t d_tar = 0.;
   Double_t w2_elastic = 0.;
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
    cross_section = 1342.0701;
    entries = 6.91701e+06;
    n_mu_ot=7.9729075e+08;
    d_tar=3.;}
if(mc=="default_bb35b5de"){

    cross_section = 1342.0701;
    entries = 6.91701e+06;
    n_mu_ot=7.9729075e+08;
    d_tar=3.;
}
if(mc=="default_bb35b5de_ideal" ){
    cross_section = 1346.0839;
    entries = 974826;
    n_mu_ot=7.9729075e+08;
    d_tar=3.;}
if(mc=="default_bb35b5de_windowOffset0" ){
    cross_section =1346.0839;
    entries =444736;
    n_mu_ot=7.9729075e+08;
    d_tar=3.;}
if(mc=="my_modifica"){
    cross_section = 1342.0701;
    entries = 6.91701e+06;
    n_mu_ot=7.9729075e+08;
    d_tar=3.;
}


//   Double_t entries =651153;


//   Double_t cross_section = 12.9221;
//   Double_t entries = 4.43846e+07;

//for(int n=0; n < h_opening_MC->GetNbinsX(); n++){cout << n << ") MC: pre " << 1./sqrt(h_opening_MC->GetBinContent(n)) << " binError " << 1./h_opening_MC->GetBinError(n) << endl;}

h_d0_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section*0.973/entries);
h_d1_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section*0.973/entries);
h_d2_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section*0.973/entries);
h_d3_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section*0.973/entries);
h_d4_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section*0.973/entries);
h_d5_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section*0.973/entries);

/*h_mod0_r0_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section*0.973/entries);
h_mod0_r1_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section*0.973/entries);
h_mod0_r2_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section*0.973/entries);

h_mod2_r0_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section*0.973/entries);
h_mod2_r1_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section*0.973/entries);
h_mod2_r2_MC->Scale(n_mu_ot* 5.5*1E+23*d_tar*1E-30 * cross_section*0.973/entries);
*/

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
di.SaveAs(Form("comparison_RDMC/%s_%s_h_dist_basic_chi150_flatArea_the15_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_Bend_%dhit_%dhit.pdf",rd.c_str(),mc.c_str(),nhits_rd,nhits_mc));


}

