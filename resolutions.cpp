#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TSystemDirectory.h"
#include <TStyle.h>

using namespace std;

void resolutions(){

TChain * cbmsim = new TChain("cbmsim");

cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_0-5mrad_1M_2hitFirstModules_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_5-10mrad_1M_2hitFirstModules_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_10-15mrad_1M_2hitFirstModules_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_15-20mrad_1M_2hitFirstModules_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_20-25mrad_1M_2hitFirstModules_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_25-32mrad_1M_2hitFirstModules_NOoutchi2_1M.root");


/*cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_0-5mrad_1M_1hitFirstModules_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_5-10mrad_1M_1hitFirstModules_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_10-15mrad_1M_1hitFirstModules_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_15-20mrad_1M_1hitFirstModules_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_20-25mrad_1M_1hitFirstModules_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_25-32mrad_1M_1hitFirstModules_NOoutchi2_1M.root");
*/

/*cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_0-5mrad_0hit_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_5-10mrad_0hit_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_10-15mrad_0hit_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_15-20mrad_0hit_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_20-25mrad_0hit_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_25-32mrad_0hit_NOoutchi2_1M.root");
*/

/*
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_0-5mrad_1M_1hitFirstModules_NOoutchi2_reassign.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_5-10mrad_1M_1hitFirstModules_NOoutchi2_reassign.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_10-15mrad_1M_1hitFirstModules_NOoutchi2_reassign.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_15-20mrad_1M_1hitFirstModules_NOoutchi2_reassign.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_20-25mrad_1M_1hitFirstModules_NOoutchi2_reassign.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_25-32mrad_1M_1hitFirstModules_NOoutchi2_reassign.root");
*/

/*cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_0-5mrad_1M_2hitFirstModules_NOoutchi2_reassign.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_5-10mrad_1M_2hitFirstModules_NOoutchi2_reassign.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_10-15mrad_1M_2hitFirstModules_NOoutchi2_reassign.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_15-20mrad_1M_2hitFirstModules_NOoutchi2_reassign.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_20-25mrad_1M_2hitFirstModules_NOoutchi2_reassign.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_25-32mrad_1M_2hitFirstModules_NOoutchi2_reassign.root");
*/

        TClonesArray *MCTrack = 0;
        TClonesArray *TrackerStripDigis = 0;
        TClonesArray *TrackerPoints = 0;
        TClonesArray *TrackerStubs = 0;
        MuE::Event *MesmerEvent = 0;
        MUonERecoOutput *ReconstructionOutput = 0;

        cbmsim->SetBranchAddress("MCTrack", &MCTrack);
        cbmsim->SetBranchAddress("TrackerPoints", &TrackerPoints);
        cbmsim->SetBranchAddress("TrackerStripDigis", &TrackerStripDigis);
        cbmsim->SetBranchAddress("TrackerStubs", &TrackerStubs);
        cbmsim->SetBranchAddress("MesmerEvent", &MesmerEvent);
        cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);


   const Int_t NBINS_mu = 14;
   Double_t edges_mu[NBINS_mu + 1] = {0.,0.0001,0.0002,0.0003,0.0004,0.0005,0.0006,0.0007,0.0008,0.0009,0.001,0.002,0.003,0.004,0.005};
   const Int_t NBINS = 14;
   Double_t edges_el[NBINS + 1] = {0.,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.010, 0.015, 0.020, 0.025, 0.032};

TH1D::SetDefaultSumw2(kTRUE);

        std::vector<TH1D*> h_res_vrtx_el(14);
        std::vector<TH1D*> h_res_vrtx_mu(14);
        std::vector<TH1D*> h_res_el_pre(14);
        std::vector<TH1D*> h_res_mu_pre(14);
        std::vector<TH1D*> h_res_vrtx_el_clones(14);
        std::vector<TH1D*> h_res_vrtx_mu_clones(14);
        std::vector<TH1D*> h_res_el_pre_clones(14);
        std::vector<TH1D*> h_res_mu_pre_clones(14);

h_res_vrtx_el.at(0)=new TH1D("res0_vrtx_el", "(the_rec_vrtx-the_true), theta_el[0,1]mrad",100,-0.001,0.001);
h_res_vrtx_el.at(1)=new TH1D("res1_vrtx_el", "(the_rec_vrtx-the_true), theta_el[1,2]mrad",100,-0.001,0.001);
h_res_vrtx_el.at(2)=new TH1D("res2_vrtx_el", "(the_rec_vrtx-the_true), theta_el[2,3]mrad",100,-0.001,0.001);
h_res_vrtx_el.at(3)=new TH1D("res3_vrtx_el", "(the_rec_vrtx-the_true), theta_el[3,4]mrad",100,-0.001,0.001);
h_res_vrtx_el.at(4)=new TH1D("res4_vrtx_el", "(the_rec_vrtx-the_true), theta_el[4,5]mrad",200,-0.002,0.002);
h_res_vrtx_el.at(5)=new TH1D("res5_vrtx_el", "(the_rec_vrtx-the_true), theta_el[5,6]mrad",200,-0.002,0.002);
h_res_vrtx_el.at(6)=new TH1D("res6_vrtx_el", "(the_rec_vrtx-the_true), theta_el[6,7]mrad",200,-0.002,0.002);
h_res_vrtx_el.at(7)=new TH1D("res7_vrtx_el", "(the_rec_vrtx-the_true), theta_el[7,8]mrad",200,-0.002,0.002);
h_res_vrtx_el.at(8)=new TH1D("res8_vrtx_el", "(the_rec_vrtx-the_true), theta_el[8,9]mrad",200,-0.002,0.002);
h_res_vrtx_el.at(9)=new TH1D("res9_vrtx_el", "(the_rec_vrtx-the_true), theta_el[9,10]mrad",200,-0.002,0.002);
h_res_vrtx_el.at(10)=new TH1D("res10_vrtx_el", "(the_rec_vrtx-the_true), theta_el[10,15]mrad",200,-0.01,0.01);
h_res_vrtx_el.at(11)=new TH1D("res11_vrtx_el", "(the_rec_vrtx-the_true), theta_el[15,20]mrad",200,-0.02,0.02);
h_res_vrtx_el.at(12)=new TH1D("res12_vrtx_el", "(the_rec_vrtx-the_true), theta_el[20,25]mrad",200,-0.02,0.02);
h_res_vrtx_el.at(13)=new TH1D("res13_vrtx_el", "(the_rec_vrtx-the_true), theta_el[25,32]mrad",200,-0.02,0.02);


h_res_vrtx_mu.at(0)=new TH1D("res0_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_mu[0.,0.1]mrad",120,-0.0004,0.0004);
h_res_vrtx_mu.at(1)=new TH1D("res1_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_mu[0.1,0.2]mrad",120,-0.0004,0.0004);
h_res_vrtx_mu.at(2)=new TH1D("res2_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_mu[0.2,0.3]mrad",120,-0.0004,0.0004);
h_res_vrtx_mu.at(3)=new TH1D("res3_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_mu[0.3,0.4]mrad",120,-0.0004,0.0004);
h_res_vrtx_mu.at(4)=new TH1D("res4_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_mu[0.4,0.5]mrad",120,-0.0004,0.0004);
h_res_vrtx_mu.at(5)=new TH1D("res5_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_mu[0.5,0.6]mrad",120,-0.0004,0.0004);
h_res_vrtx_mu.at(6)=new TH1D("res6_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_mu[0.6,0.7]mrad",120,-0.0004,0.0004);
h_res_vrtx_mu.at(7)=new TH1D("res7_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_mu[0.7,0.8]mrad",120,-0.0004,0.0004);
h_res_vrtx_mu.at(8)=new TH1D("res8_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_mu[0.8,0.9]mrad",120,-0.0004,0.0004);
h_res_vrtx_mu.at(9)=new TH1D("res9_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_mu[0.9,1]mrad",120,-0.0004,0.0004);
h_res_vrtx_mu.at(10)=new TH1D("res10_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_mu[0.3,0.4]mrad",120,-0.0004,0.0004);
h_res_vrtx_mu.at(11)=new TH1D("res11_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_mu[0.4,0.5]mrad",120,-0.0004,0.0004);
h_res_vrtx_mu.at(12)=new TH1D("res12_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_mu[3,4]mrad",120,-0.0004,0.0004);
h_res_vrtx_mu.at(13)=new TH1D("res13_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_mu[4,5]mrad",120,-0.0004,0.0004);


h_res_el_pre.at(0)=new TH1D("res0_el_pre", "(the_rec-the_true), theta_el[0,1]mrad",100,-0.001,0.001);
h_res_el_pre.at(1)=new TH1D("res1_el_pre", "(the_rec-the_true), theta_el[1,2]mrad",100,-0.001,0.001);
h_res_el_pre.at(2)=new TH1D("res2_el_pre", "(the_rec-the_true), theta_el[2,3]mrad",100,-0.001,0.001);
h_res_el_pre.at(3)=new TH1D("res3_el_pre", "(the_rec-the_true), theta_el[3,4]mrad",100,-0.001,0.001);
h_res_el_pre.at(4)=new TH1D("res4_el_pre", "(the_rec-the_true), theta_el[4,5]mrad",200,-0.002,0.002);
h_res_el_pre.at(5)=new TH1D("res5_el_pre", "(the_rec-the_true), theta_el[5,6]mrad",200,-0.002,0.002);
h_res_el_pre.at(6)=new TH1D("res6_el_pre", "(the_rec-the_true), theta_el[6,7]mrad",200,-0.002,0.002);
h_res_el_pre.at(7)=new TH1D("res7_el_pre", "(the_rec-the_true), theta_el[7,8]mrad",200,-0.002,0.002);
h_res_el_pre.at(8)=new TH1D("res8_el_pre", "(the_rec-the_true), theta_el[8,9]mrad",200,-0.002,0.002);
h_res_el_pre.at(9)=new TH1D("res9_el_pre", "(the_rec-the_true), theta_el[9,10]mrad",200,-0.002,0.002);
h_res_el_pre.at(10)=new TH1D("res10_el_pre", "(the_rec-the_true), theta_el[10,15]mrad",200,-0.01,0.01);
h_res_el_pre.at(11)=new TH1D("res11_el_pre", "(the_rec-the_true), theta_el[15,20]mrad",200,-0.02,0.02);
h_res_el_pre.at(12)=new TH1D("res12_el_pre", "(the_rec-the_true), theta_el[20,25]mrad",200,-0.02,0.02);
h_res_el_pre.at(13)=new TH1D("res13_el_pre", "(the_rec-the_true), theta_el[25,32]mrad",200,-0.02,0.02);


h_res_mu_pre.at(0)=new TH1D("res0_mu_pre", "(thmu_rec-thmu_true), theta_mu[0.,0.1]mrad",120,-0.0004,0.0004);
h_res_mu_pre.at(1)=new TH1D("res1_mu_pre", "(thmu_rec-thmu_true), theta_mu[0.1,0.2]mrad",120,-0.0004,0.0004);
h_res_mu_pre.at(2)=new TH1D("res2_mu_pre", "(thmu_rec-thmu_true), theta_mu[0.2,0.3]mrad",120,-0.0004,0.0004);
h_res_mu_pre.at(3)=new TH1D("res3_mu_pre", "(thmu_rec-thmu_true), theta_mu[0.3,0.4]mrad",120,-0.0004,0.0004);
h_res_mu_pre.at(4)=new TH1D("res4_mu_pre", "(thmu_rec-thmu_true), theta_mu[0.4,0.5]mrad",120,-0.0004,0.0004);
h_res_mu_pre.at(5)=new TH1D("res5_mu_pre", "(thmu_rec-thmu_true), theta_mu[0.5,0.6]mrad",120,-0.0004,0.0004);
h_res_mu_pre.at(6)=new TH1D("res6_mu_pre", "(thmu_rec-thmu_true), theta_mu[0.6,0.7]mrad",120,-0.0004,0.0004);
h_res_mu_pre.at(7)=new TH1D("res7_mu_pre", "(thmu_rec-thmu_true), theta_mu[0.7,0.8]mrad",120,-0.0004,0.0004);
h_res_mu_pre.at(8)=new TH1D("res8_mu_pre", "(thmu_rec-thmu_true), theta_mu[0.8,0.9]mrad",120,-0.0004,0.0004);
h_res_mu_pre.at(9)=new TH1D("res9_mu_pre", "(thmu_rec-thmu_true), theta_mu[0.9,1]mrad",120,-0.0004,0.0004);
h_res_mu_pre.at(10)=new TH1D("res10_mu_pre", "(thmu_rec-thmu_true), theta_mu[1,2]mrad",120,-0.0004,0.0004);
h_res_mu_pre.at(11)=new TH1D("res11_mu_pre", "(thmu_rec-thmu_true), theta_mu[2,3]mrad",120,-0.0004,0.0004);
h_res_mu_pre.at(12)=new TH1D("res12_mu_pre", "(thmu_rec-thmu_true), theta_mu[3,4]mrad",120,-0.0004,0.0004);
h_res_mu_pre.at(13)=new TH1D("res13_mu_pre", "(thmu_rec-thmu_true), theta_mu[4,5]mrad",120,-0.0004,0.0004);


h_res_vrtx_el_clones.at(0)=new TH1D("res0_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[0,1]mrad with clones",100,-0.001,0.001);
h_res_vrtx_el_clones.at(1)=new TH1D("res1_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[1,2]mrad with clones",100,-0.001,0.001);
h_res_vrtx_el_clones.at(2)=new TH1D("res2_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[2,3]mrad with clones",100,-0.001,0.001);
h_res_vrtx_el_clones.at(3)=new TH1D("res3_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[3,4]mrad with clones",100,-0.001,0.001);
h_res_vrtx_el_clones.at(4)=new TH1D("res4_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[4,5]mrad with clones",200,-0.002,0.002);
h_res_vrtx_el_clones.at(5)=new TH1D("res5_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[5,6]mrad with clones",200,-0.002,0.002);
h_res_vrtx_el_clones.at(6)=new TH1D("res6_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[6,7]mrad with clones",200,-0.002,0.002);
h_res_vrtx_el_clones.at(7)=new TH1D("res7_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[7,8]mrad with clones",200,-0.002,0.002);
h_res_vrtx_el_clones.at(8)=new TH1D("res8_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[8,9]mrad with clones",200,-0.002,0.002);
h_res_vrtx_el_clones.at(9)=new TH1D("res9_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[9,10]mrad with clones",200,-0.002,0.002);
h_res_vrtx_el_clones.at(10)=new TH1D("res10_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[10,15]mrad with clones",200,-0.01,0.01);
h_res_vrtx_el_clones.at(11)=new TH1D("res11_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[15,20]mrad with clones",200,-0.02,0.02);
h_res_vrtx_el_clones.at(12)=new TH1D("res12_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[20,25]mrad with clones",200,-0.02,0.02);
h_res_vrtx_el_clones.at(13)=new TH1D("res13_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[25,32]mrad with clones",200,-0.02,0.02);

h_res_vrtx_mu_clones.at(0)=new TH1D("res0_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_mu[0.,0.1]mrad with clones",120,-0.0004,0.0004);
h_res_vrtx_mu_clones.at(1)=new TH1D("res1_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_mu[0.1,0.2]mrad with clones",120,-0.0004,0.0004);
h_res_vrtx_mu_clones.at(2)=new TH1D("res2_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_mu[0.2,0.3]mrad with clones",120,-0.0004,0.0004);
h_res_vrtx_mu_clones.at(3)=new TH1D("res3_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_mu[0.3,0.4]mrad with clones",120,-0.0004,0.0004);
h_res_vrtx_mu_clones.at(4)=new TH1D("res4_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_mu[0.4,0.5]mrad with clones",120,-0.0004,0.0004);
h_res_vrtx_mu_clones.at(5)=new TH1D("res5_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_mu[0.5,0.6]mrad with clones",120,-0.0004,0.0004);
h_res_vrtx_mu_clones.at(6)=new TH1D("res6_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_mu[0.6,0.7]mrad with clones",120,-0.0004,0.0004);
h_res_vrtx_mu_clones.at(7)=new TH1D("res7_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_mu[0.7,0.8]mrad with clones",120,-0.0004,0.0004);
h_res_vrtx_mu_clones.at(8)=new TH1D("res8_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_mu[0.8,0.9]mrad with clones",120,-0.0004,0.0004);
h_res_vrtx_mu_clones.at(9)=new TH1D("res9_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_mu[0.9,1]mrad with clones",120,-0.0004,0.0004);
h_res_vrtx_mu_clones.at(10)=new TH1D("res10_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_mu[1,2]mrad with clones",120,-0.0004,0.0004);
h_res_vrtx_mu_clones.at(11)=new TH1D("res11_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_mu[2,3]mrad with clones",120,-0.0004,0.0004);
h_res_vrtx_mu_clones.at(12)=new TH1D("res12_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_mu[3,4]mrad with clones",120,-0.0004,0.0004);
h_res_vrtx_mu_clones.at(13)=new TH1D("res13_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_mu[4,5]mrad with clones",120,-0.0004,0.0004);

h_res_el_pre_clones.at(0)=new TH1D("res0_el_pre_clones", "(the_rec-the_true), theta_el[0,1]mrad with clones",100,-0.001,0.001);
h_res_el_pre_clones.at(1)=new TH1D("res1_el_pre_clones", "(the_rec-the_true), theta_el[1,2]mrad with clones",100,-0.001,0.001);
h_res_el_pre_clones.at(2)=new TH1D("res2_el_pre_clones", "(the_rec-the_true), theta_el[2,3]mrad with clones",100,-0.001,0.001);
h_res_el_pre_clones.at(3)=new TH1D("res3_el_pre_clones", "(the_rec-the_true), theta_el[3,4]mrad with clones",100,-0.001,0.001);
h_res_el_pre_clones.at(4)=new TH1D("res4_el_pre_clones", "(the_rec-the_true), theta_el[4,5]mrad with clones",200,-0.002,0.002);
h_res_el_pre_clones.at(5)=new TH1D("res5_el_pre_clones", "(the_rec-the_true), theta_el[5,6]mrad with clones",200,-0.002,0.002);
h_res_el_pre_clones.at(6)=new TH1D("res6_el_pre_clones", "(the_rec-the_true), theta_el[6,7]mrad with clones",200,-0.002,0.002);
h_res_el_pre_clones.at(7)=new TH1D("res7_el_pre_clones", "(the_rec-the_true), theta_el[7,8]mrad with clones",200,-0.002,0.002);
h_res_el_pre_clones.at(8)=new TH1D("res8_el_pre_clones", "(the_rec-the_true), theta_el[8,9]mrad with clones",200,-0.002,0.002);
h_res_el_pre_clones.at(9)=new TH1D("res9_el_pre_clones", "(the_rec-the_true), theta_el[9,10]mrad with clones",200,-0.002,0.002);
h_res_el_pre_clones.at(10)=new TH1D("res10_el_pre_clones", "(the_rec-the_true), theta_el[10,15]mrad with clones",200,-0.01,0.01);
h_res_el_pre_clones.at(11)=new TH1D("res11_el_pre_clones", "(the_rec-the_true), theta_el[15,20]mrad with clones",200,-0.02,0.02);
h_res_el_pre_clones.at(12)=new TH1D("res12_el_pre_clones", "(the_rec-the_true), theta_el[20,25]mrad with clones",200,-0.02,0.02);
h_res_el_pre_clones.at(13)=new TH1D("res13_el_pre_clones", "(the_rec-the_true), theta_el[25,32]mrad with clones",200,-0.02,0.02);

h_res_mu_pre_clones.at(0)=new TH1D("res0_mu_pre_clones", "(thmu_rec-thmu_true), theta_mu[0.,0.1]mrad with clones",120,-0.0004,0.0004);
h_res_mu_pre_clones.at(1)=new TH1D("res1_mu_pre_clones", "(thmu_rec-thmu_true), theta_mu[0.1,0.2]mrad with clones",120,-0.0004,0.0004);
h_res_mu_pre_clones.at(2)=new TH1D("res2_mu_pre_clones", "(thmu_rec-thmu_true), theta_mu[0.2,0.3]mrad with clones",120,-0.0004,0.0004);
h_res_mu_pre_clones.at(3)=new TH1D("res3_mu_pre_clones", "(thmu_rec-thmu_true), theta_mu[0.3,0.4]mrad with clones",120,-0.0004,0.0004);
h_res_mu_pre_clones.at(4)=new TH1D("res4_mu_pre_clones", "(thmu_rec-thmu_true), theta_mu[0.4,0.5]mrad with clones",120,-0.0004,0.0004);
h_res_mu_pre_clones.at(5)=new TH1D("res5_mu_pre_clones", "(thmu_rec-thmu_true), theta_mu[0.5,0.6]mrad with clones",120,-0.0004,0.0004);
h_res_mu_pre_clones.at(6)=new TH1D("res6_mu_pre_clones", "(thmu_rec-thmu_true), theta_mu[0.6,0.7]mrad with clones",120,-0.0004,0.0004);
h_res_mu_pre_clones.at(7)=new TH1D("res7_mu_pre_clones", "(thmu_rec-thmu_true), theta_mu[0.7,0.8]mrad with clones",120,-0.0004,0.0004);
h_res_mu_pre_clones.at(8)=new TH1D("res8_mu_pre_clones", "(thmu_rec-thmu_true), theta_mu[0.8,0.9]mrad with clones",120,-0.0004,0.0004);
h_res_mu_pre_clones.at(9)=new TH1D("res9_mu_pre_clones", "(thmu_rec-thmu_true), theta_mu[0.9,1]mrad with clones",120,-0.0004,0.0004);
h_res_mu_pre_clones.at(10)=new TH1D("res10_mu_pre_clones", "(thmu_rec-thmu_true), theta_mu[1,2]mrad with clones",120,-0.0004,0.0004);
h_res_mu_pre_clones.at(11)=new TH1D("res11_mu_pre_clones", "(thmu_rec-thmu_true), theta_mu[2,3]mrad with clones",120,-0.0004,0.0004);
h_res_mu_pre_clones.at(12)=new TH1D("res12_mu_pre_clones", "(thmu_rec-thmu_true), theta_mu[3,4]mrad with clones",120,-0.0004,0.0004);
h_res_mu_pre_clones.at(13)=new TH1D("res13_mu_pre_clones", "(thmu_rec-thmu_true), theta_mu[4,5]mrad with clones",120,-0.0004,0.0004);



TH1D *h_res_vrtx_muin=new TH1D("h_res_vrtx_muin", "(thmu_in_rec-thmu_in_true) ",120,-0.0004,0.0004);

TH1D *h_phi=new TH1D("phi","phi mu in",100,-180.,180.);

//wnorm 1M/1M NLO 300outlier/nooutlier NLO 2hitshared
	double r_wnorm[6]={20.786765103274643,33.313091221576336,42.396733790329876,50.584815206143652,61.828110400735824,106.88513370134392};




for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
                cbmsim->GetEntry(i);
                if(i%10000 == 0) cout<<"Entry "<<i<<endl;

int index=99;
int index_mu=99;

int code_mu_in=-99;
int code_e=-99;
int code_mu=-99;
        TVector3 p_muin_MC;
        TVector3 p_mu_MC;
        TVector3 p_e_MC;
        Double_t the_gen, thmu_gen,theX_gen,theY_gen,thmuX_gen,thmuY_gen;

           int hit_modXmu=0; int hit_modXe=0;
           int hit_modYmu=0; int hit_modYe=0;
           int stereo_mu=0; int stereo_e=0;
     	   int hit_modXmuin=0;
     int hit_modYmuin=0;
     int stereo_muin=0;

        for(int n = 0; n < MCTrack->GetEntries(); n++) {
         const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
         if(MCTr->interactionID()==0 and MCTr->pdgCode()==-13) {code_mu_in=n; p_muin_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); p_muin_MC.Unit();
        for(int s=0; s<TrackerPoints->GetEntries(); s++)
                         {const MUonETrackerPoint *TrackerPt = static_cast<const MUonETrackerPoint*>(TrackerPoints->At(s));
         if(TrackerPt->trackPDGCode()==-13 and TrackerPt->trackID()==n and TrackerPt->stationID()==0){
                                                                                                 if(TrackerPt->moduleID()==0 or TrackerPt->moduleID()==4) hit_modXmuin++;
                                                                                                 else if(TrackerPt->moduleID()==1 or TrackerPt->moduleID()==5) hit_modYmuin++;
                                                                                                 else if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_muin++;}
                }
	}
         if(MCTr->interactionID()==45 and MCTr->pdgCode()==11) {code_e=n; p_e_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); p_e_MC.Unit(); the_gen=p_muin_MC.Angle(p_e_MC);
                                                                                 theX_gen=MCTr->ax();
                                                                                 theY_gen=MCTr->ay();

        for(int s=0; s<TrackerPoints->GetEntries(); s++)
                         {const MUonETrackerPoint *TrackerPt = static_cast<const MUonETrackerPoint*>(TrackerPoints->At(s));
         if(TrackerPt->trackPDGCode()==11 and TrackerPt->trackID()==code_e and TrackerPt->stationID()==1){
                                                                                                 if(TrackerPt->moduleID()==0 or TrackerPt->moduleID()==4) hit_modXe++;
                                                                                                 else if(TrackerPt->moduleID()==1 or TrackerPt->moduleID()==5) hit_modYe++;
                                                                                                 else if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_e++;}
        	}
	 }
         if(MCTr->interactionID()==45 and MCTr->pdgCode()==-13) {code_mu=n; p_mu_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); p_mu_MC.Unit(); thmu_gen=p_muin_MC.Angle(p_mu_MC);
                                                                                 thmuX_gen=MCTr->ax();
                                                                                 thmuY_gen=MCTr->ay();
        for(int s=0; s<TrackerPoints->GetEntries(); s++)
                         {const MUonETrackerPoint *TrackerPt = static_cast<const MUonETrackerPoint*>(TrackerPoints->At(s));
	 if(TrackerPt->trackPDGCode()==-13 and TrackerPt->trackID()==code_mu and TrackerPt->stationID()==1){
                                                                                                 if(TrackerPt->moduleID()==0 or TrackerPt->moduleID()==4) hit_modXmu++;
                                                                                                 else if(TrackerPt->moduleID()==1 or TrackerPt->moduleID()==5) hit_modYmu++;
                                                                                                 else if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_mu++;}
		}
	 }
         i++;
	}

	if(code_mu_in!=-99 and code_e!=-99 and code_mu!=-99 and hit_modXmuin==4 and hit_modYmuin==4 and stereo_muin>1 and hit_modXe==4 and hit_modYe==4 and stereo_e>1 and hit_modXmu==4 and hit_modYmu==4 and stereo_mu>1){


double wnorm=0.;

//cout << "1) the_gen " << the_gen << " thmu_gen " << thmu_gen << " index " << index << " index_mu " << index_mu << endl;
if(the_gen>=0. and the_gen<0.001){index=0; wnorm=r_wnorm[0];}
if(the_gen>=0.001 and the_gen<0.002){index=1; wnorm=r_wnorm[0];}
if(the_gen>=0.002 and the_gen<0.003){index=2; wnorm=r_wnorm[0];}
if(the_gen>=0.003 and the_gen<0.004){index=3; wnorm=r_wnorm[0];}
if(the_gen>=0.004 and the_gen<0.005){index=4; wnorm=r_wnorm[0];}
if(the_gen>=0.005 and the_gen<0.006){index=5; wnorm=r_wnorm[1];}
if(the_gen>=0.006 and the_gen<0.007){index=6; wnorm=r_wnorm[1];}
if(the_gen>=0.007 and the_gen<0.008){index=7; wnorm=r_wnorm[1];}
if(the_gen>=0.008 and the_gen<0.009){index=8; wnorm=r_wnorm[1];}
if(the_gen>=0.009 and the_gen<0.010){index=9; wnorm=r_wnorm[1];}
if(the_gen>=0.010 and the_gen<0.015){index=10; wnorm=r_wnorm[2];}
if(the_gen>=0.015 and the_gen<0.020){index=11; wnorm=r_wnorm[3];}
if(the_gen>=0.020 and the_gen<0.025){index=12; wnorm=r_wnorm[4];}
if(the_gen>=0.025 and the_gen<=0.033){index=13; wnorm=r_wnorm[5];}



if(thmu_gen>=0. and thmu_gen<0.0001)index_mu=0;
if(thmu_gen>=0.0001 and thmu_gen<0.0002)index_mu=1;
if(thmu_gen>=0.0002 and thmu_gen<0.0003)index_mu=2;
if(thmu_gen>=0.0003 and thmu_gen<0.0004)index_mu=3;
if(thmu_gen>=0.0004 and thmu_gen<0.0005)index_mu=4;
if(thmu_gen>=0.0005 and thmu_gen<0.0006)index_mu=5;
if(thmu_gen>=0.0006 and thmu_gen<0.0007)index_mu=6;
if(thmu_gen>=0.0007 and thmu_gen<0.0008)index_mu=7;
if(thmu_gen>=0.0008 and thmu_gen<0.0009)index_mu=8;
if(thmu_gen>=0.0009 and thmu_gen<0.001)index_mu=9;
if(thmu_gen>=0.001 and thmu_gen<0.002)index_mu=10;
if(thmu_gen>=0.002 and thmu_gen<0.003)index_mu=11;
if(thmu_gen>=0.003 and thmu_gen<0.004)index_mu=12;
if(thmu_gen>=0.004 and thmu_gen<=0.005)index_mu=13;


Int_t yes2=0;
Double_t th_inx,th_iny,x0_in,y0_in;
Double_t chi2_muin=999.;
Int_t stubs_muin=0;
Int_t e=0;
Int_t mu=0;
TVector3 p_muin,p_e,p_mu;
Double_t the_rec,theX_rec,theY_rec,thmu_rec,thmuX_rec,thmuY_rec,the_rec_vrtx,thmu_rec_vrtx;

MUonERecoOutputVertex vrtx = ReconstructionOutput->bestVertex();

Double_t chi=vrtx.chi2perDegreeOfFreedom();


vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
int other=0;
int mu_in=0;
vector<double> quality_e; quality_e.reserve(5);
vector<double> quality_mu; quality_mu.reserve(5);

         for (auto&& track : tracks) {

        if(code_mu_in==track.linkedTrackID() and track.sector()==0){
	mu_in++;
         th_inx=track.xSlope();
         th_iny=track.ySlope();
         x0_in=track.x0();
         y0_in=track.y0();
         chi2_muin=track.chi2perDegreeOfFreedom();
         stubs_muin=track.hits().size();
        p_muin.SetXYZ(th_inx,th_iny,1.0);
        p_muin=p_muin.Unit();
                        }
        else if(track.processIDofLinkedTrack()==45 and track.sector()==1)
                {
                 if(code_e==track.linkedTrackID()) {yes2++; e++; theX_rec=track.xSlope(); theY_rec=track.ySlope();p_e.SetXYZ(theX_rec,theY_rec,1.0); p_e=p_e.Unit();
							if(track.fractionOfHitsSharedWithLinkedTrack()>=0.65)the_rec=p_e.Angle(p_muin);
                                                    quality_e.push_back(track.fractionOfHitsSharedWithLinkedTrack());}

                 if(code_mu==track.linkedTrackID()) {yes2++; mu++;thmuX_rec=track.xSlope(); thmuY_rec=track.ySlope();p_mu.SetXYZ(thmuX_rec,thmuY_rec,1.0); p_mu=p_mu.Unit();
							if(track.fractionOfHitsSharedWithLinkedTrack()>=0.65)thmu_rec=p_mu.Angle(p_muin);
                                                    quality_mu.push_back(track.fractionOfHitsSharedWithLinkedTrack());}
                }
        else if(track.processIDofLinkedTrack()!=45 and track.sector()==1){other++;}

         }//for


if(stubs_muin>=5 and chi2_muin<5){//and thmu_gen>0.0002){

                        h_res_vrtx_muin->Fill(p_muin.Theta()-p_muin_MC.Theta(),MesmerEvent->wgt_LO*wnorm);

if(find_if(quality_e.begin(),quality_e.end(),[](double i){return i>=0.65;})!=end(quality_e) and find_if(quality_mu.begin(),quality_mu.end(),[](double i){return i>=0.65;})!=end(quality_mu)){



	if(chi!=0 and e==1 and mu==1){
                        h_res_el_pre.at(index)->Fill(the_rec-the_gen,MesmerEvent->wgt_LO*wnorm);
                        h_res_mu_pre.at(index_mu)->Fill(thmu_rec-thmu_gen,MesmerEvent->wgt_LO*wnorm);
			}

	if(chi!=0 and e>=1 and mu>=1){
                        h_res_el_pre_clones.at(index)->Fill(the_rec-the_gen,MesmerEvent->wgt_LO*wnorm);
                        h_res_mu_pre_clones.at(index_mu)->Fill(thmu_rec-thmu_gen,MesmerEvent->wgt_LO*wnorm);}


}


if(vrtx.outgoingMuon().fractionOfHitsSharedWithLinkedTrack()>=0.65 and vrtx.outgoingElectron().fractionOfHitsSharedWithLinkedTrack()>=0.65){


	if(chi!=0 and e==1 and mu==1){
			if(vrtx.outgoingMuon().linkedTrackID()==code_mu)h_res_vrtx_mu.at(index_mu)->Fill(vrtx.muonTheta()-thmu_gen,MesmerEvent->wgt_LO*wnorm);
			if(vrtx.outgoingMuon().linkedTrackID()==code_e)h_res_vrtx_el.at(index)->Fill(vrtx.muonTheta()-the_gen,MesmerEvent->wgt_LO*wnorm);
                        if(vrtx.outgoingElectron().linkedTrackID()==code_mu)h_res_vrtx_mu.at(index_mu)->Fill(vrtx.electronTheta()-thmu_gen,MesmerEvent->wgt_LO*wnorm);
                        if(vrtx.outgoingElectron().linkedTrackID()==code_e)h_res_vrtx_el.at(index)->Fill(vrtx.electronTheta()-the_gen,MesmerEvent->wgt_LO*wnorm);}

	if(chi!=0 and e>=1 and mu>=1){
                        if(vrtx.outgoingMuon().linkedTrackID()==code_mu)h_res_vrtx_mu_clones.at(index_mu)->Fill(vrtx.muonTheta()-thmu_gen,MesmerEvent->wgt_LO*wnorm);
                        if(vrtx.outgoingMuon().linkedTrackID()==code_e)h_res_vrtx_el_clones.at(index)->Fill(vrtx.muonTheta()-the_gen,MesmerEvent->wgt_LO*wnorm);
                        if(vrtx.outgoingElectron().linkedTrackID()==code_mu)h_res_vrtx_mu_clones.at(index_mu)->Fill(vrtx.electronTheta()-thmu_gen,MesmerEvent->wgt_LO*wnorm);
                        if(vrtx.outgoingElectron().linkedTrackID()==code_e)h_res_vrtx_el_clones.at(index)->Fill(vrtx.electronTheta()-the_gen,MesmerEvent->wgt_LO*wnorm);}
			}//if quality vrtx


                }//if mu_in
        }//if generated
}//for


TCanvas n0("n0","n0",2100,3500);
n0.Divide(3,5);
for(int m=0; m<NBINS; m++){
n0.cd(m+1);
h_res_vrtx_el.at(m)->SetMinimum(1.);
h_res_vrtx_el.at(m)->Draw("hist");
h_res_el_pre.at(m)->SetLineColor(kRed+1);
h_res_el_pre.at(m)->Draw("hist same");
h_res_el_pre.at(m)->SaveAs(Form("quality_tracks/root/res_prevrtx_el_%d_2hitFirstModules_LO.root",static_cast<char>(m)));
h_res_vrtx_el.at(m)->SaveAs(Form("quality_tracks/root/res_vrtx_el_%d_2hitFirstModules_LO.root",static_cast<char>(m)));
//gPad->SetLogy();
}
n0.SaveAs("quality_tracks/h_res_pre_el_2hitFirstModules_LO.pdf");


TCanvas n1("n1","n1",2100,3500);
n1.Divide(3,5);
for(int m=0; m<NBINS; m++){
n1.cd(m+1);
h_res_vrtx_el_clones.at(m)->Draw("hist");
h_res_el_pre_clones.at(m)->SetLineColor(kRed+1);
h_res_el_pre_clones.at(m)->Draw("hist same");
h_res_el_pre_clones.at(m)->SaveAs(Form("quality_tracks/root/res_prevrtx_clones_el_%d_2hitFirstModules_LO.root",static_cast<char>(m)));
h_res_vrtx_el_clones.at(m)->SaveAs(Form("quality_tracks/root/res_vrtx_clones_el_%d_2hitFirstModules_LO.root",static_cast<char>(m)));
}
n1.SaveAs("quality_tracks/h_res_vrtx_el_clones_2hitFirstModules_LO.pdf");


TCanvas n2("n2","n2",2100,3500);
n2.Divide(3,5);
for(int m=0; m<NBINS_mu; m++){
n2.cd(m+1);
h_res_vrtx_mu.at(m)->SetMinimum(1.);
h_res_vrtx_mu.at(m)->Draw("hist");
h_res_mu_pre.at(m)->SetLineColor(kRed+1);
h_res_mu_pre.at(m)->Draw("hist same");
h_res_mu_pre.at(m)->SaveAs(Form("quality_tracks/root/res_prevrtx_AngleMu_%d_2hitFirstModules_LO.root",static_cast<char>(m)));
h_res_vrtx_mu.at(m)->SaveAs(Form("quality_tracks/root/res_vrtx_AngleMu_%d_2hitFirstModules_LO.root",static_cast<char>(m)));
//gPad->SetLogy();
}
n2.SaveAs("quality_tracks/h_res_pre_mu_2hitFirstModules_LO.pdf");

TCanvas n3("n3","n3",2100,3500);
n3.Divide(3,5);
for(int m=0; m<NBINS_mu; m++){
n3.cd(m+1);
h_res_vrtx_mu_clones.at(m)->Draw("hist");
h_res_mu_pre_clones.at(m)->SetLineColor(kRed+1);
h_res_mu_pre_clones.at(m)->Draw("hist same");
h_res_mu_pre_clones.at(m)->SaveAs(Form("quality_tracks/root/res_prevrtx_clones_AngleMu_%d_2hitFirstModules_LO.root",static_cast<char>(m)));
h_res_vrtx_mu_clones.at(m)->SaveAs(Form("quality_tracks/root/res_vrtx_clones_AngleMu_%d_2hitFirstModules_LO.root",static_cast<char>(m)));
}
n3.SaveAs("quality_tracks/h_res_vrtx_mu_clones_2hitFirstModules_LO.pdf");

TCanvas in("in","in",700,700);
h_res_vrtx_muin->Draw();
in.SaveAs("quality_tracks/h_res_pmuin_2hitFirstModules_LO.pdf");




TCanvas d1("d1","d1",700,700);
h_phi->Draw("hist");
d1.SaveAs("phi_180.pdf");




}

