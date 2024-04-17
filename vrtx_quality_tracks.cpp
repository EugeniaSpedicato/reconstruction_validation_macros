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

void vrtx_quality_tracks(){

TChain * cbmsim = new TChain("cbmsim");

//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/theta_0-32mrad_1M_2hit_NOoutchi2_1M.root");
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/theta_0-32mrad_1M_1hit_NOoutchi2.root");

/*
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_0-5mrad_1M_1hit_NOoutchi2.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_5-10mrad_1M_1hit_NOoutchi2.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_10-15mrad_1M_1hit_NOoutchi2.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_15-20mrad_1M_1hit_NOoutchi2.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_20-25mrad_1M_1hit_NOoutchi2.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_25-32mrad_1M_1hit_NOoutchi2.root");
*/

/*
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_0-5mrad_1M_2hit_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_5-10mrad_1M_2hit_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_10-15mrad_1M_2hit_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_15-20mrad_1M_2hit_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_20-25mrad_1M_2hit_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_25-32mrad_1M_2hit_NOoutchi2_1M.root");
*/


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


                auto pos_on_track = [](double q, double m, double z){return (q + m*z);};
   const Int_t NBINS = 14;
   const Int_t NBINS_mu = 14;
   Double_t edges_el[NBINS + 1] = {0.,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.010, 0.015, 0.020, 0.025, 0.032};
   Double_t edges_mu[NBINS_mu + 1] = {0.,0.0001,0.0002,0.0003,0.0004,0.0005,0.0006,0.0007,0.0008,0.0009,0.001,0.002,0.003,0.004,0.005};

TH1D::SetDefaultSumw2(kTRUE);

        std::vector<TH1D*> h_res_vrtx_el(6);
        std::vector<TH1D*> h_res_vrtx_mu(6);
        std::vector<TH1D*> h_res_el_pre(6);
        std::vector<TH1D*> h_res_mu_pre(6);
        std::vector<TH1D*> h_res_vrtx_el_clones(6);
        std::vector<TH1D*> h_res_vrtx_mu_clones(6);
        std::vector<TH1D*> h_res_el_pre_clones(6);
        std::vector<TH1D*> h_res_mu_pre_clones(6);

h_res_vrtx_el.at(0)=new TH1D("res0_vrtx_el", "(the_rec_vrtx-the_true), theta_el[0,5]mrad",150,-0.03,0.03);
h_res_vrtx_el.at(1)=new TH1D("res1_vrtx_el", "(the_rec_vrtx-the_true), theta_el[5,10]mrad",150,-0.03,0.03);
h_res_vrtx_el.at(2)=new TH1D("res2_vrtx_el", "(the_rec_vrtx-the_true), theta_el[10,15]mrad",500,-0.1,0.1);
h_res_vrtx_el.at(3)=new TH1D("res3_vrtx_el", "(the_rec_vrtx-the_true), theta_el[15,20]mrad",500,-0.1,0.1);
h_res_vrtx_el.at(4)=new TH1D("res4_vrtx_el", "(the_rec_vrtx-the_true), theta_el[20,25]mrad",500,-0.1,0.1);
h_res_vrtx_el.at(5)=new TH1D("res5_vrtx_el", "(the_rec_vrtx-the_true), theta_el[25,32]mrad",500,-0.1,0.1);

h_res_vrtx_mu.at(0)=new TH1D("res0_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_el[0,5]mrad",600,-0.003,0.003);
h_res_vrtx_mu.at(1)=new TH1D("res1_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_el[5,10]mrad",120,-0.0006,0.0006);
h_res_vrtx_mu.at(2)=new TH1D("res2_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_el[10,15]mrad",120,-0.0006,0.0006);
h_res_vrtx_mu.at(3)=new TH1D("res3_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_el[15,20]mrad",120,-0.0006,0.0006);
h_res_vrtx_mu.at(4)=new TH1D("res4_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_el[20,25]mrad",120,-0.0006,0.0006);
h_res_vrtx_mu.at(5)=new TH1D("res5_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_el[25,32]mrad",120,-0.0006,0.0006);


h_res_el_pre.at(0)=new TH1D("res0_el_pre", "(the_rec-the_true), theta_el[0,5]mrad",150,-0.03,0.03);
h_res_el_pre.at(1)=new TH1D("res1_el_pre", "(the_rec-the_true), theta_el[5,10]mrad",150,-0.03,0.03);
h_res_el_pre.at(2)=new TH1D("res2_el_pre", "(the_rec-the_true), theta_el[10,15]mrad",500,-0.1,0.1);
h_res_el_pre.at(3)=new TH1D("res3_el_pre", "(the_rec-the_true), theta_el[15,20]mrad",500,-0.1,0.1);
h_res_el_pre.at(4)=new TH1D("res4_el_pre", "(the_rec-the_true), theta_el[20,25]mrad",500,-0.1,0.1);
h_res_el_pre.at(5)=new TH1D("res5_el_pre", "(the_rec-the_true), theta_el[25,32]mrad",500,-0.1,0.1);

h_res_mu_pre.at(0)=new TH1D("res0_mu_pre", "(thmu_rec-thmu_true), theta_el[0,5]mrad",600,-0.003,0.003);
h_res_mu_pre.at(1)=new TH1D("res1_mu_pre", "(thmu_rec-thmu_true), theta_el[5,10]mrad",120,-0.0006,0.0006);
h_res_mu_pre.at(2)=new TH1D("res2_mu_pre", "(thmu_rec-thmu_true), theta_el[10,15]mrad",120,-0.0006,0.0006);
h_res_mu_pre.at(3)=new TH1D("res3_mu_pre", "(thmu_rec-thmu_true), theta_el[15,20]mrad",120,-0.0006,0.0006);
h_res_mu_pre.at(4)=new TH1D("res4_mu_pre", "(thmu_rec-thmu_true), theta_el[20,25]mrad",120,-0.0006,0.0006);
h_res_mu_pre.at(5)=new TH1D("res5_mu_pre", "(thmu_rec-thmu_true), theta_el[25,32]mrad",120,-0.0006,0.0006);



h_res_vrtx_el_clones.at(0)=new TH1D("res0_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[0,5]mrad with clones",150,-0.03,0.03);
h_res_vrtx_el_clones.at(1)=new TH1D("res1_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[5,10]mrad with clones",150,-0.03,0.03);
h_res_vrtx_el_clones.at(2)=new TH1D("res2_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[10,15]mrad with clones",500,-0.1,0.1);
h_res_vrtx_el_clones.at(3)=new TH1D("res3_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[15,20]mrad with clones",500,-0.1,0.1);
h_res_vrtx_el_clones.at(4)=new TH1D("res4_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[20,25]mrad with clones",500,-0.1,0.1);
h_res_vrtx_el_clones.at(5)=new TH1D("res5_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[25,32]mrad with clones",500,-0.1,0.1);

h_res_vrtx_mu_clones.at(0)=new TH1D("res0_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_el[0,5]mrad with clones",200,-0.001,0.001);
h_res_vrtx_mu_clones.at(1)=new TH1D("res1_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_el[5,10]mrad with clones",120,-0.0006,0.0006);
h_res_vrtx_mu_clones.at(2)=new TH1D("res2_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_el[10,15]mrad with clones",120,-0.0006,0.0006);
h_res_vrtx_mu_clones.at(3)=new TH1D("res3_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_el[15,20]mrad with clones",120,-0.0006,0.0006);
h_res_vrtx_mu_clones.at(4)=new TH1D("res4_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_el[20,25]mrad with clones",120,-0.0006,0.0006);
h_res_vrtx_mu_clones.at(5)=new TH1D("res5_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_el[25,32]mrad with clones",120,-0.0006,0.0006);


h_res_el_pre_clones.at(0)=new TH1D("res0_el_pre_clones", "(the_rec-the_true), theta_el[0,5]mrad with clones",150,-0.03,0.03);
h_res_el_pre_clones.at(1)=new TH1D("res1_el_pre_clones", "(the_rec-the_true), theta_el[5,10]mrad with clones",150,-0.03,0.03);
h_res_el_pre_clones.at(2)=new TH1D("res2_el_pre_clones", "(the_rec-the_true), theta_el[10,15]mrad with clones",500,-0.1,0.1);
h_res_el_pre_clones.at(3)=new TH1D("res3_el_pre_clones", "(the_rec-the_true), theta_el[15,20]mrad with clones",500,-0.1,0.1);
h_res_el_pre_clones.at(4)=new TH1D("res4_el_pre_clones", "(the_rec-the_true), theta_el[20,25]mrad with clones",500,-0.1,0.1);
h_res_el_pre_clones.at(5)=new TH1D("res5_el_pre_clones", "(the_rec-the_true), theta_el[25,32]mrad with clones",500,-0.1,0.1);

h_res_mu_pre_clones.at(0)=new TH1D("res0_mu_pre_clones", "(thmu_rec-thmu_true), theta_el[0,5]mrad with clones",200,-0.001,0.001);
h_res_mu_pre_clones.at(1)=new TH1D("res1_mu_pre_clones", "(thmu_rec-thmu_true), theta_el[5,10]mrad with clones",120,-0.0006,0.0006);
h_res_mu_pre_clones.at(2)=new TH1D("res2_mu_pre_clones", "(thmu_rec-thmu_true), theta_el[10,15]mrad with clones",120,-0.0006,0.0006);
h_res_mu_pre_clones.at(3)=new TH1D("res3_mu_pre_clones", "(thmu_rec-thmu_true), theta_el[15,20]mrad with clones",120,-0.0006,0.0006);
h_res_mu_pre_clones.at(4)=new TH1D("res4_mu_pre_clones", "(thmu_rec-thmu_true), theta_el[20,25]mrad with clones",120,-0.0006,0.0006);
h_res_mu_pre_clones.at(5)=new TH1D("res5_mu_pre_clones", "(thmu_rec-thmu_true), theta_el[25,32]mrad with clones",120,-0.0006,0.0006);




TH1D *h_res_vrtx_muin=new TH1D("h_res_vrtx_muin", "(thmu_in_rec-thmu_in_true) ",120,-0.0006,0.0006);

TH2D *h_2d=new TH2D("h_2d","Electron reco vrtx angle Vs muon reco vrtx angle",1000,0.,0.03,200,0.,0.005);
TH2D *h_2d_gen=new TH2D("h_2d_gen","Electron true angle Vs muon true angle",1000,0.,0.03,200,0.,0.005);

TH2D *h_2d_wrong_gen=new TH2D("h_2d_wrong_gen","Electron gen angle Vs muon gen angle wrong assignment",1000,0.,0.03,200,0.,0.005);
TH2D *h_2d_wrong=new TH2D("h_2d_wrong","Electron reco vrtx angle Vs muon reco vrtx angle wrong assignment",1000,0.,0.03,200,0.,0.005);


TH1D *h_process_el=new TH1D("h_process_el","Process of particle to which mismatched outgoing elastic electron is associated",50,0,50);
TH1D *h_process_mu=new TH1D("h_process_mu","Process of particle to which mismatched outgoing elastic muon is associated",50,0,50);

TH1D *chi_wrong=new TH1D("chi_wrong","vrtx chi2 when wrong assignment",300,0,100.);
TH1D *chi_correct=new TH1D("chi_correct","vrtx chi2 when correct assignment",300,0,100.);
TH1D *chi_swapped=new TH1D("chi_swapped","vrtx chi2 when swapped assignment",300,0,100.);

TH1D *h_z_wrong=new TH1D("h_z_wrong","zPositionFit for wrong events",800,800.,1000.);
TH1D *h_z_correct=new TH1D("h_z_correct","zPositionFit for correct events",800,800.,1000.);
TH1D *h_z_swapped=new TH1D("h_z_swapped","zPositionFit for swapped events",800,800.,1000.);

TH1D *theta_e_clone=new TH1D("theta_e_clone", "Reco elastic event: Electron scattering reco angles from MESMER with clones",NBINS,edges_el);
TH1D *theta_mu_clone=new TH1D("theta_mu_clone", "Reco elastic event: Muon scattering reco angles from MESMER with clones",NBINS_mu,edges_mu);

TH1D *theta_e_single_clone=new TH1D("theta_single_e_clone", "Reco single electron: Electron scattering reco angles from MESMER with clones",NBINS,edges_el);
TH1D *theta_mu_single_clone=new TH1D("theta_single_mu_clone", "Reco single muon: Muon scattering reco angles from MESMER with clones",NBINS_mu,edges_mu);

TH1D *theta_e_clone_quality=new TH1D("theta_e_clone_quality", "Reco elastic event: Electron scattering reco angles from MESMER with clones quality check",NBINS,edges_el);
TH1D *theta_mu_clone_quality=new TH1D("theta_mu_clone_quality", "Reco elastic event: Muon scattering reco angles from MESMER with clones quality check",NBINS_mu,edges_mu);

TH1D *theta_e_single_clone_quality=new TH1D("theta_single_e_clone_quality", "Reco single electron: Electron scattering reco angles from MESMER with clones quality check",NBINS,edges_el);
TH1D *theta_mu_single_clone_quality=new TH1D("theta_single_mu_clone_quality", "Reco single muon: Muon scattering reco angles from MESMER with clones quality check",NBINS_mu,edges_mu);

TH1D *theta_e_gen=new TH1D("theta_e_gen", "Electron scattering generated angles from MESMER",NBINS,edges_el);
TH1D *theta_mu_gen=new TH1D("theta_mu_gen", "Muon scattering generated angles from MESMER",NBINS_mu,edges_mu);


//wnorm 1M/1M NLO 300outlier/nooutlier NLO 2hitshared
	double r_wnorm[6]={20.786765103274643,33.313091221576336,42.396733790329876,50.584815206143652,61.828110400735824,106.88513370134392};

//1000chi2 1M single range
//      double r_wnorm[6]={169.10075957722117,169.10075957722117,169.10075957722117,169.10075957722117,169.10075957722117,169.10075957722117};


double sum_wgt=0.;
double sumW2[6]={0.};

double n_correct=0.;
double n_swapped=0.;
double n_wrong=0.;
double tot=0.;
double all=0.;
double n_novrtx=0.;
double n_yesvrtx=0.;
double n_yes_bestvrtx=0.;


double no_clones=0.;
double clones=0.;
double no_clones_q=0.;
double clones_q=0.;


for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
                cbmsim->GetEntry(i);
                if(i%100 == 0) cout<<"Entry "<<i<<endl;

int index=7;

Double_t code_mu_in=-99;
Double_t code_e=-99;
Double_t code_mu=-99;
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

if(the_gen>=0 and the_gen<0.005){index=0; wnorm=r_wnorm[0];}
if(the_gen>=0.005 and the_gen<0.010){index=1; wnorm=r_wnorm[1];}
if(the_gen>=0.010 and the_gen<0.015){index=2; wnorm=r_wnorm[2];}
if(the_gen>=0.015 and the_gen<0.020){index=3; wnorm=r_wnorm[3];}
if(the_gen>=0.020 and the_gen<0.025){index=4; wnorm=r_wnorm[4];}
if(the_gen>=0.025 and the_gen<0.032){index=5; wnorm=r_wnorm[5];}

all+=MesmerEvent->wgt_full*wnorm;


double opening_angle=p_mu_MC.Angle(p_e_MC);


Int_t yes_mu=0;
Int_t yes_e=0;
Int_t yes2=0;
Double_t th_inx,th_iny,x0_in,y0_in;
Double_t chi2_muin=999.;
Int_t stubs_muin=0;
Int_t e=0;
Int_t mu=0;
TVector3 p_muin,p_e,p_mu;
Double_t the_rec,theX_rec,theY_rec,thmu_rec,thmuX_rec,thmuY_rec,the_rec_vrtx,thmu_rec_vrtx;
MUonERecoOutputVertex vrtx = ReconstructionOutput->bestVertex();

vector<MUonERecoOutputVertex> vrtx_all = ReconstructionOutput->reconstructedVertices();
Double_t chi=vrtx.chi2perDegreeOfFreedom();

if(vrtx_all.size()==0){n_novrtx+=MesmerEvent->wgt_full*wnorm;}
else{n_yesvrtx+=MesmerEvent->wgt_full*wnorm; if(chi!=0)n_yes_bestvrtx+=MesmerEvent->wgt_full*wnorm;}

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
        if(track.processIDofLinkedTrack()==45 and track.sector()==1)
                {
                 if(code_e==track.linkedTrackID()) {yes2++; e++; theX_rec=track.xSlope(); theY_rec=track.ySlope();p_e.SetXYZ(theX_rec,theY_rec,1.0); p_e=p_e.Unit();
							if(track.fractionOfHitsSharedWithLinkedTrack()>=0.65)the_rec=p_e.Angle(p_muin);
                                                    quality_e.push_back(track.fractionOfHitsSharedWithLinkedTrack());}

                 if(code_mu==track.linkedTrackID()) {yes2++; mu++;thmuX_rec=track.xSlope(); thmuY_rec=track.ySlope();p_mu.SetXYZ(thmuX_rec,thmuY_rec,1.0); p_mu=p_mu.Unit();
							if(track.fractionOfHitsSharedWithLinkedTrack()>=0.65)thmu_rec=p_mu.Angle(p_muin);
                                                    quality_mu.push_back(track.fractionOfHitsSharedWithLinkedTrack());}
                }
        if(track.processIDofLinkedTrack()!=45 and track.sector()==1){other++;}

         }//for


//if(stubs_muin==6 and abs(posxIN)<=1.5 and abs(posyIN)<=1.5 and chi2_muin<2 and thmu_gen>0.0002){
if(stubs_muin>=5 and chi2_muin<2 and thmu_gen>0.0002){

tot+=MesmerEvent->wgt_full*wnorm;
theta_mu_gen->Fill(thmu_gen,MesmerEvent->wgt_full*wnorm);
theta_e_gen->Fill(the_gen,MesmerEvent->wgt_full*wnorm);


if(find_if(quality_e.begin(),quality_e.end(),[](double i){return i>=0.65;})!=end(quality_e) and find_if(quality_mu.begin(),quality_mu.end(),[](double i){return i>=0.65;})!=end(quality_mu)){

 if(e>=1){theta_e_single_clone->Fill(the_gen,MesmerEvent->wgt_full*wnorm);}
 if(mu>=1){theta_mu_single_clone->Fill(thmu_gen,MesmerEvent->wgt_full*wnorm);}

 if(e==1 and mu==1 and chi!=0) no_clones+=MesmerEvent->wgt_full*wnorm;
 if(((e>=1 and mu>1) or (e>1 and mu>=1)) and chi!=0){clones+=MesmerEvent->wgt_full*wnorm;}

 if(e>=1 and mu>=1 and chi!=0){
                   	theta_mu_clone->Fill(thmu_gen,MesmerEvent->wgt_full*wnorm);
                        theta_e_clone->Fill(the_gen,MesmerEvent->wgt_full*wnorm);}

//quality checks
 if(chi!=0 and mu>=1 and vrtx.outgoingMuon().fractionOfHitsSharedWithLinkedTrack()>=0.65 and ((vrtx.outgoingMuon().linkedTrackID()==code_e and vrtx.outgoingElectron().linkedTrackID()==code_mu) or (vrtx.outgoingMuon().linkedTrackID()==code_mu and vrtx.outgoingElectron().linkedTrackID()==code_e)) ){
		theta_mu_single_clone->Fill(thmu_gen,MesmerEvent->wgt_full*wnorm);}
 else if(chi!=0 and e>=1 and vrtx.outgoingElectron().fractionOfHitsSharedWithLinkedTrack()>=0.65 and ((vrtx.outgoingMuon().linkedTrackID()==code_e and vrtx.outgoingElectron().linkedTrackID()==code_mu) or (vrtx.outgoingMuon().linkedTrackID()==code_mu and vrtx.outgoingElectron().linkedTrackID()==code_e)) ){
		theta_e_single_clone->Fill(the_gen,MesmerEvent->wgt_full*wnorm);}
 if( (vrtx.outgoingMuon().linkedTrackID()!=code_e and vrtx.outgoingMuon().linkedTrackID()!=code_mu) or (vrtx.outgoingElectron().linkedTrackID()!=code_mu and vrtx.outgoingElectron().linkedTrackID()!=code_e) ){
        n_wrong+=MesmerEvent->wgt_full*wnorm;h_2d_wrong_gen->Fill(the_gen,thmu_gen,MesmerEvent->wgt_full*wnorm);h_2d_wrong->Fill(the_rec_vrtx,thmu_rec_vrtx,MesmerEvent->wgt_full*wnorm);
        chi_wrong->Fill(chi,MesmerEvent->wgt_full*wnorm);h_z_wrong->Fill(vrtx.zPositionFit(),MesmerEvent->wgt_full*wnorm);
		for(int n = 0; n < MCTrack->GetEntries(); n++) {
                 const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
                if(n==vrtx.outgoingMuon().linkedTrackID()) h_process_mu->Fill(MCTr->interactionID(),MesmerEvent->wgt_full*wnorm);//cout << "Muon is linked to pdgCode " << MCTr->pdgCode() << " from interaction " << MCTr->interactionID() << endl;
                if(n==vrtx.outgoingElectron().linkedTrackID()) h_process_el->Fill(MCTr->interactionID(),MesmerEvent->wgt_full*wnorm);//cout << "Electron is linked to pdgCode " << MCTr->pdgCode() << " from interaction " << MCTr->interactionID() << endl;
                }
	continue;}

if(vrtx.outgoingMuon().fractionOfHitsSharedWithLinkedTrack()>=0.65 and vrtx.outgoingElectron().fractionOfHitsSharedWithLinkedTrack()>=0.65){

   if((vrtx.outgoingMuon().linkedTrackID()==code_e and vrtx.outgoingElectron().linkedTrackID()==code_mu) or (vrtx.outgoingMuon().linkedTrackID()==code_mu and vrtx.outgoingElectron().linkedTrackID()==code_e))n_correct+=MesmerEvent->wgt_full*wnorm;

/*
		        chi_correct->Fill(chi,MesmerEvent->wgt_full*wnorm);
			h_z_correct->Fill(vrtx.zPositionFit(),MesmerEvent->wgt_full*wnorm);
			h_2d->Fill(vrtx.electronTheta(),vrtx.muonTheta(),MesmerEvent->wgt_full*wnorm);
                        h_2d_gen->Fill(the_gen,thmu_gen,MesmerEvent->wgt_full*wnorm);
*/
                        h_res_vrtx_muin->Fill(p_muin.Theta()-p_muin_MC.Theta(),MesmerEvent->wgt_full*wnorm);


if(chi!=0 and e==1 and mu==1){
no_clones_q+=MesmerEvent->wgt_full*wnorm;
			h_res_el_pre.at(index)->Fill(the_rec-the_gen,MesmerEvent->wgt_full*wnorm);
			h_res_mu_pre.at(index)->Fill(thmu_rec-thmu_gen,MesmerEvent->wgt_full*wnorm);
			if(vrtx.outgoingMuon().linkedTrackID()==code_mu)h_res_vrtx_mu.at(index)->Fill(vrtx.muonTheta()-thmu_gen,MesmerEvent->wgt_full*wnorm);
			if(vrtx.outgoingMuon().linkedTrackID()==code_e)h_res_vrtx_mu.at(index)->Fill(vrtx.electronTheta()-the_gen,MesmerEvent->wgt_full*wnorm);
                        if(vrtx.outgoingElectron().linkedTrackID()==code_mu)h_res_vrtx_el.at(index)->Fill(vrtx.muonTheta()-thmu_gen,MesmerEvent->wgt_full*wnorm);
                        if(vrtx.outgoingElectron().linkedTrackID()==code_e)h_res_vrtx_el.at(index)->Fill(vrtx.electronTheta()-the_gen,MesmerEvent->wgt_full*wnorm);

}

if(chi!=0 and ((e>=1 and mu>1) or (e>1 and mu>=1))){clones_q+=MesmerEvent->wgt_full*wnorm;}

if(chi!=0 and e>=1 and mu>=1){
			theta_mu_clone_quality->Fill(thmu_gen,MesmerEvent->wgt_full*wnorm);
			theta_e_clone_quality->Fill(the_gen,MesmerEvent->wgt_full*wnorm);
//per il pre un po di modifiche con i cloni. Devi scegliere tracce migliori e fare angolo
                        h_res_el_pre_clones.at(index)->Fill(the_rec-the_gen,MesmerEvent->wgt_full*wnorm);
                        h_res_mu_pre_clones.at(index)->Fill(thmu_rec-thmu_gen,MesmerEvent->wgt_full*wnorm);
//----------------------------------------------
                        if(vrtx.outgoingMuon().linkedTrackID()==code_mu)h_res_vrtx_mu_clones.at(index)->Fill(vrtx.muonTheta()-thmu_gen,MesmerEvent->wgt_full*wnorm);
                        if(vrtx.outgoingMuon().linkedTrackID()==code_e)h_res_vrtx_mu_clones.at(index)->Fill(vrtx.muonTheta()-the_gen,MesmerEvent->wgt_full*wnorm);
                        if(vrtx.outgoingElectron().linkedTrackID()==code_mu)h_res_vrtx_el_clones.at(index)->Fill(vrtx.electronTheta()-thmu_gen,MesmerEvent->wgt_full*wnorm);
                        if(vrtx.outgoingElectron().linkedTrackID()==code_e)h_res_vrtx_el_clones.at(index)->Fill(vrtx.electronTheta()-the_gen,MesmerEvent->wgt_full*wnorm);

				}

			}//if quality vrtx

		 }// if quality reco
                }//if mu_in

        }//if generated

}//for


cout << "fractions elastic events without quality without clones " << no_clones/tot*100 << "%" << endl;
cout << "fractions elastic events without quality with clones " << clones/tot*100 << "%" << endl;

cout << "fractions elastic events with quality without clones " << no_clones_q/tot*100 << "%" << endl;
cout << "fractions elastic events with quality with clones " << clones_q/tot*100 << "%" << endl;



cout<< "Number of reco events without vertex " << n_novrtx/all*100 << "%"<< endl;
cout<< "Number of reco events with vertex " << n_yesvrtx/all*100 << "%"<< endl;
cout<< "Number of reco events with best vertex " <<  n_yes_bestvrtx/all*100 << "%"<< endl;

cout<< " Number of vrtx created with good incoming muon over all number of reconstructible events " << tot/all*100 << "%"<< endl;
cout<< " Number of vrtx's track id assignment correct " << n_correct << " when vrtx is created -> " << n_correct/tot*100 << "%"<< endl;
cout<< " Number of vrtx's track id assignment swapped " << n_swapped << " when vrtx is created -> " << n_swapped/tot*100 << "%"<< endl;
cout<< " Number of vrtx's track id assignment which are not correct and are not just swapped " << n_wrong << " when vrtx is created -> " << n_wrong/tot*100 << "%"<< endl;


TCanvas n0("n0","n0",1000,1000);
n0.Divide(2,3);
for(int m=0; m<6; m++){
n0.cd(m+1);
h_res_vrtx_el.at(m)->SetMinimum(1.);
h_res_vrtx_el.at(m)->Draw("hist");
h_res_el_pre.at(m)->SetLineColor(kRed+1);
h_res_el_pre.at(m)->Draw("hist same");
h_res_el_pre.at(m)->SaveAs(Form("quality_tracks/root/res_prevrtx_el_%d_2hitFirstModules.root",static_cast<char>(m)));
h_res_vrtx_el.at(m)->SaveAs(Form("quality_tracks/root/res_vrtx_el_%d_2hitFirstModules.root",static_cast<char>(m)));
//gPad->SetLogy();
}
n0.SaveAs("quality_tracks/h_res_vrtx_el_1hitFirstModules.pdf");


TCanvas n1("n1","n1",1000,1000);
n1.Divide(2,3);
for(int m=0; m<6; m++){
n1.cd(m+1);
h_res_vrtx_el_clones.at(m)->Draw("hist");
h_res_el_pre_clones.at(m)->SetLineColor(kRed+1);
h_res_el_pre_clones.at(m)->Draw("hist same");
h_res_el_pre_clones.at(m)->SaveAs(Form("quality_tracks/root/res_prevrtx_clones_el_%d_2hitFirstModules.root",static_cast<char>(m)));
h_res_vrtx_el_clones.at(m)->SaveAs(Form("quality_tracks/root/res_vrtx_clones_el_%d_2hitFirstModules.root",static_cast<char>(m)));
}
n1.SaveAs("quality_tracks/h_res_vrtx_el_clones_1hitFirstModules.pdf");


TCanvas n2("n2","n2",1000,1000);
n2.Divide(2,3);
for(int m=0; m<6; m++){
n2.cd(m+1);
h_res_vrtx_mu.at(m)->SetMinimum(1.);
h_res_vrtx_mu.at(m)->Draw("hist");
h_res_mu_pre.at(m)->SetLineColor(kRed+1);
h_res_mu_pre.at(m)->Draw("hist same");
h_res_mu_pre.at(m)->SaveAs(Form("quality_tracks/root/res_prevrtx_mu_%d_2hitFirstModules.root",static_cast<char>(m)));
h_res_vrtx_mu.at(m)->SaveAs(Form("quality_tracks/root/res_vrtx_mu_%d_2hitFirstModules.root",static_cast<char>(m)));
//gPad->SetLogy();
}
n2.SaveAs("quality_tracks/h_res_vrtx_mu_1hitFirstModules.pdf");

TCanvas n3("n3","n3",1000,1000);
n3.Divide(2,3);
for(int m=0; m<6; m++){
n3.cd(m+1);
h_res_vrtx_mu_clones.at(m)->Draw("hist");
h_res_mu_pre_clones.at(m)->SetLineColor(kRed+1);
h_res_mu_pre_clones.at(m)->Draw("hist same");
h_res_mu_pre_clones.at(m)->SaveAs(Form("quality_tracks/root/res_prevrtx_clones_mu_%d_2hitFirstModules.root",static_cast<char>(m)));
h_res_vrtx_mu_clones.at(m)->SaveAs(Form("quality_tracks/root/res_vrtx_clones_mu_%d_2hitFirstModules.root",static_cast<char>(m)));
}
n3.SaveAs("quality_tracks/h_res_vrtx_mu_clones_1hitFirstModules.pdf");


TCanvas n4("n4","n4",1000,1000);
n4.Divide(1,5);
n4.cd(1);
h_res_vrtx_muin->Draw("hist");
n4.cd(2);
h_2d->Draw("COLZ");
n4.cd(3);
h_2d_gen->Draw("COLZ");
n4.cd(4);
h_2d_wrong_gen->Draw("COLZ");
n4.cd(5);
h_2d_wrong->Draw("COLZ");
n4.SaveAs("quality_tracks/muin_2d_1hitFirstModules.pdf");


TCanvas n5("n5","n5",1000,1000);
n5.Divide(1,2);
n5.cd(1);
h_process_mu->Draw("hist");
gPad->SetLogy();
n5.cd(2);
h_process_el->Draw("hist");
n5.SaveAs("quality_tracks/h_process_1hitFirstModules.pdf");


TCanvas n6("n6","n6",700,700);
n6.Divide(1,2);
n6.cd(1);
chi_correct->SetMinimum(1.);
chi_correct->Draw("hist");
chi_swapped->SetLineColor(kGreen+1);
chi_swapped->Draw("hist same");
chi_wrong->SetLineColor(kRed+1);
chi_wrong->Draw("hist same");
gPad->SetLogy();

n6.cd(2);
h_z_correct->SetMinimum(1.);
h_z_correct->Draw("hist");
h_z_swapped->SetLineColor(kGreen+1);
h_z_swapped->Draw("hist same");
h_z_wrong->SetLineColor(kRed+1);
h_z_wrong->Draw("hist same");
gPad->SetLogy();
n6.SaveAs("quality_tracks/vrtx_chi2_h_z_1hitFirstModules.pdf");



TH1D * h0_original = (TH1D*) theta_e_clone->Clone();
TH1D * h0gen = (TH1D*) theta_e_gen->Clone();
h0_original->Divide(h0_original,h0gen,1,1,"B");
TH1D * h1_original = (TH1D*) theta_mu_clone->Clone();
TH1D * h1gen = (TH1D*) theta_mu_gen->Clone();
h1_original->Divide(h1_original,h1gen,1,1,"B");

TH1D * h0 = (TH1D*) theta_e_clone_quality->Clone();
h0->Divide(h0,h0gen,1,1,"B");
TH1D * h1 = (TH1D*) theta_mu_clone_quality->Clone();
h1->Divide(h1,h1gen,1,1,"B");

TH1D * h0_quality = (TH1D*) theta_e_clone_quality->Clone();
TH1D * h0_no = (TH1D*) theta_e_clone->Clone();
h0_quality->Divide(h0_quality,h0_no,1,1,"B");
TH1D * h1_quality = (TH1D*) theta_mu_clone_quality->Clone();
TH1D * h1_no = (TH1D*) theta_mu_clone->Clone();
h1_quality->Divide(h1_quality,h1_no,1,1,"B");


TCanvas n7("n7","n7",700,700);
n7.Divide(2,2);
n7.cd(1);
TGaxis::SetMaxDigits(3);
h0_original->SetTitle("Reco eff. elastic event as a funtion of el angle");
h0_original->SetLineColor(kOrange+2);
h0_original->SetMinimum(0.6);
h0_original->Draw("E");
h0->Draw("E same");
h0->SaveAs("proposal/quality_el_gen_2hitFirstModules.root");
gStyle->SetOptStat(0);
n7.cd(2);
TGaxis::SetMaxDigits(3);
h1_original->SetTitle("Reco eff. elastic event as a funtion of mu angle");
h1_original->SetLineColor(kOrange+2);
h1_original->SetMinimum(0.6);
h1_original->Draw("E");
h1->Draw("E same");
h1->SaveAs("proposal/quality_mu_gen_2hitFirstModules.root");
gStyle->SetOptStat(0);
n7.cd(3);
h0_quality->SetTitle("Ratio between sample with quality cut/without quality track VS el angle");
TGaxis::SetMaxDigits(3);
h0_quality->SetMinimum(0.6);
h0_quality->Draw("E");
h0_quality->SaveAs("proposal/quality_el_clone_2hitFirstModules.root");
gStyle->SetOptStat(0);
n7.cd(4);
h1_quality->SetTitle("Ratio between sample with quality cut/without quality track VS el angle");
TGaxis::SetMaxDigits(3);
h1_quality->SetMinimum(0.6);
h1_quality->Draw("E");
h1_quality->SaveAs("proposal/quality_mu_clone_2hitFirstModules.root");
gStyle->SetOptStat(0);
n7.SaveAs("prova_1hitFirstModules.pdf");

}

