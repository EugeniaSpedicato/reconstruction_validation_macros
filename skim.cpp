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

#include "Station_hits.h"
#include "releff.h"

using namespace std;

void skim(){

int NMODULES = 12;
int iskim = 41;

TChain * cbmsim = new TChain("cbmsim");
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/minbias_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/minbias_1M_2.root");
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_0-32mrad_1M_1hit.root");

        TClonesArray *MCTrack = 0;
   TClonesArray *TrackerStubs = 0;
   MUonERecoOutput *ReconstructionOutput = 0;

   cbmsim->SetBranchAddress("MCTrack", &MCTrack);
   cbmsim->SetBranchAddress("TrackerStubs", &TrackerStubs);
   cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);


 auto pos_on_track = [](double q, double m, double z){return (q + m*z);};

string n_mult="ALLmu";
string hname = "h_skim/histos_"+n_mult+"_minbias.root";
  TFile *hfile = new TFile(hname.c_str(),"RECREATE");

  TFile *o_file;
  TTree *o_tree;

if (iskim > 0) {
    //
    // OUTPUT
    //
    // max size of skimmed output files: 2 GB
    TTree::SetMaxTreeSize(2000000000LL);
    //  TTree::SetMaxTreeSize(10000000LL);
    o_file = TFile::Open(("h_skim/output_"+n_mult+"_minbias.root").c_str(), "recreate");
    if (!o_file) {cout << "Error opening file" << std::endl; exit(-1);}
    o_tree = new TTree("cbmsim", "");

        TClonesArray *MCTrack2 = 0;
   TClonesArray *TrackerStubs2 = 0;
   MUonERecoOutput *ReconstructionOutput2 = 0;

   o_tree->SetBranchAddress("MCTrack", &MCTrack2);
   o_tree->SetBranchAddress("TrackerStubs", &TrackerStubs2);
   o_tree->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput2);
  }



const Double_t nmax = 70;
  Int_t nb = nmax;

TH1D *h_nstubsPerModule = new TH1D("h_nstubsPerModule","nstubs on modules",NMODULES,0,NMODULES);
  TH1D *h_nstubsPerModule_S0 = new TH1D("h_nstubsPerModule_S0","nstubs on modules - Station 0",NMODULES,0,NMODULES);
  TH1D *h_nstubsPerModule_S1 = new TH1D("h_nstubsPerModule_S1","nstubs on modules - Station 1",NMODULES,0,NMODULES);


  TH2D *h2_nStubs = new TH2D("h2_nStubs","N stubs on S1 vs S0, with ALLmu",nb,0,nmax,nb,0,nmax);
  TH1D *h_nStubs_0 = new TH1D("h_nStubs_0","Tot N stubs First Station - S0, with ALLmu",nb,0,nmax);
  TH1D *h_nStubs_1 = new TH1D("h_nStubs_1","Tot N stubs Second Station - S1, with ALLmu",nb,0,nmax);

  TH1D *h_nTotStubs = new TH1D("h_nTotStubs","Total N of Stubs with ALLmu",nb,0,nmax);

  TH1D *h_nTotStubs_after = new TH1D("h_nTotStubs_after","Total N of Stubs after skim with ALLmu",nb,0,nmax);

  TH1D *h_nTotStubs_zero_0 = new TH1D("h_nTotStubs_zero_0","Total N of Stubs - events with zero hits in S_0",nb,0,nmax);
  TH1D *h_nTotStubs_zero_1 = new TH1D("h_nTotStubs_zero_1","Total N of Stubs - events with zero hits in S_1",nb,0,nmax);
  TH1D *h_nTotStubs_noise = new TH1D("h_nTotStubs_noise","Total N of Stubs noise events",nb,0,nmax);
  TH1D *h_nTotStubs_missing_0 = new TH1D("h_nTotStubs_missing_0","Total N of Stubs muons missing S_0",nb,0,nmax);
  TH1D *h_nTotStubs_missing_1 = new TH1D("h_nTotStubs_missing_1","Total N of Stubs muons missing S_1",nb,0,nmax);
  TH1D *h_nTotStubs_not_trackable = new TH1D("h_nTotStubs_not_trackable","Total N of Stubs - not trackable events",nb,0,nmax);
  TH1D *h_nTotStubs_passing_golden = new TH1D("h_nTotStubs_passing_golden","Total N of Stubs for best passing mu",nb,0,nmax);
  TH1D *h_nTotStubs_passing_1mu = new TH1D("h_nTotStubs_passing_1mu","Total N of Stubs for passing mu",nb,0,nmax);
  TH1D *h_nTotStubs_passing_2mu = new TH1D("h_nTotStubs_passing_2mu","Total N of Stubs for passing 2mu",nb,0,nmax);
  TH1D *h_nTotStubs_passing_3mu = new TH1D("h_nTotStubs_passing_3mu","Total N of Stubs for passing 3mu",nb,0,nmax);

  TH1D *h_nTotStubs_golden = new TH1D("h_nTotStubs_golden","Total N of Stubs for Golden signal selection; Number of Stubs; Events",nb,0,nmax);  
  TH1D *h_nTotStubs_umsel = new TH1D("h_nTotStubs_umsel","Total N of Stubs for Umberto's selection",nb,0,nmax);
  TH1D *h_nTotStubs_single_cand = new TH1D("h_nTotStubs_single_cand","Total N of Stubs for candidate single mu signal",nb,0,nmax);
  TH1D *h_nTotStubs_single_clean = new TH1D("h_nTotStubs_single_clean","Total N of Stubs for clean single mu signal",nb,0,nmax);
  TH1D *h_nTotStubs_pileup_any = new TH1D("h_nTotStubs_pileup_any","Total N of Stubs for signal+(any)mu pileup",nb,0,nmax);
  TH1D *h_nTotStubs_pileup_2mu = new TH1D("h_nTotStubs_pileup_2mu","Total N of Stubs for 2mu pileup",nb,0,nmax);
  TH1D *h_nTotStubs_pileup_3mu = new TH1D("h_nTotStubs_pileup_3mu","Total N of Stubs for 3mu pileup",nb,0,nmax);
  TH1D *h_nTotStubs_pileup_4mu = new TH1D("h_nTotStubs_pileup_4mu","Total N of Stubs for 4mu pileup",nb,0,nmax);
  TH1D *h_nTotStubs_pileup_many = new TH1D("h_nTotStubs_pileup_many","Total N of Stubs for >=4mu pileup",nb,0,nmax);
  TH1D *h_nTotStubs_pileup_skim = new TH1D("h_nTotStubs_pileup_skim","Total N of Stubs for skimmed pileup mu",nb,0,nmax);

  TH1D *process_st0 = new TH1D("process_st0","Interaction ID of reco tracks in station 0",50,0,50);
  TH1D *process_st1 = new TH1D("process_st1","Interaction ID of reco tracks in station 1",50,0,50);

Long64_t N_input_events = 0;
  Long64_t N_input_stubs_0 = 0;
  Long64_t N_input_stubs_1 = 0;
  Long64_t N_input_stubs = 0;    // total S0+S1
  Long64_t N_stubs_skim = 0;

  Long64_t N_zero_0 = 0;
  Long64_t N_noise_0 = 0;
  Long64_t N_passing_golden_0 = 0;
  Long64_t N_trackable_0 = 0;
  Long64_t N_single_cand_0 = 0;
  Long64_t N_single_clean_0 = 0;
  Long64_t N_pileup_any_0 = 0;
  Long64_t N_pileup_2mu_0 = 0;
  Long64_t N_pileup_3mu_0 = 0;
  Long64_t N_pileup_4mu_0 = 0;
  Long64_t N_pileup_many_0 = 0;
  
  Long64_t N_zero_1 = 0;
  Long64_t N_noise_1 = 0;
  Long64_t N_passing_golden_1 = 0;
  Long64_t N_passing_cand_1 = 0;
  Long64_t N_passing_clean_1 = 0;
  Long64_t N_2tracks_1 = 0;
  Long64_t N_golden_1 = 0;
  
  Long64_t N_signal_plus_pileup_1 = 0;
  Long64_t N_signal_plus_tripileup_1 = 0;
  Long64_t N_signal_plus_fourpileup_1 = 0;
  Long64_t N_pileup_2mu_1 = 0;
  Long64_t N_pileup_3mu_1 = 0;
  Long64_t N_pileup_4mu_1 = 0;
 
  Long64_t N_zero_events = 0;
  Long64_t N_noise_events = 0;
  Long64_t N_missing_0_events = 0;
  Long64_t N_missing_1_events = 0;
  Long64_t N_not_trackable_events = 0;
  Long64_t N_golden_events = 0;
  Long64_t N_umsel_events = 0;
  Long64_t N_passing_golden_events = 0;
  Long64_t N_passing_golden_S0_zero_S1_events = 0;
  Long64_t N_passing_golden_S1_zero_S0_events = 0;
  Long64_t N_passing_golden_S0_noise_S1_events = 0;
  Long64_t N_passing_golden_S1_noise_S0_events = 0;
  Long64_t N_passing_golden_S0_trackable_S1_events = 0;
  Long64_t N_passing_golden_S1_trackable_S0_events = 0;
  Long64_t N_passing_1mu_events = 0;
  Long64_t N_passing_2mu_events = 0;
  Long64_t N_passing_3mu_events = 0;
  Long64_t N_passing_4mu_events = 0;
  Long64_t N_single_cand_events = 0;
  Long64_t N_single_clean_events = 0;
  Long64_t N_pileup_any_events = 0;
  Long64_t N_pileup_2mu_events = 0;
  Long64_t N_pileup_3mu_events = 0;
  Long64_t N_pileup_4mu_events = 0;
  Long64_t N_pileup_many_events = 0;
  Long64_t N_pileup_skim_events = 0;
  Long64_t N_loose_events = 0;




for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {

 cbmsim->GetEntry(i);
 if(i%100 == 0) cout<<"Entry "<<i<<endl;

int mu_gen=0;
int e_out=0;
int mu_out=0;
double xz=10.;
double yz=10.;
double xz_e=10.;
double yz_e=10.;
double xz_mu=10.;
double yz_mu=10.;


           int hit_modXmu=0; int hit_modXe=0;
           int hit_modYmu=0; int hit_modYe=0;
           int stereo_mu=0; int stereo_e=0;

        for(int n = 0; n < MCTrack->GetEntries(); n++) {
         const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
         if(MCTr->interactionID()==0 and MCTr->pdgCode()==13){mu_gen++; xz=pos_on_track(MCTr->bx(),MCTr->ax(),911.2); yz=pos_on_track(MCTr->by(),MCTr->ay(),911.2);}

         if(MCTr->interactionID()==45 and MCTr->pdgCode()==11) { e_out=n; xz_e=pos_on_track(MCTr->bx(),MCTr->ax(),911.2+93.7693); yz_e=pos_on_track(MCTr->by(),MCTr->ay(),911.2+93.7693);}
         if(MCTr->interactionID()==45 and MCTr->pdgCode()==-13) {mu_out=n; xz_mu=pos_on_track(MCTr->bx(),MCTr->ax(),911.2+93.7693); yz_mu=pos_on_track(MCTr->by(),MCTr->ay(),911.2+93.7693);}
        }


//cout << "expected number of incoming muons: " << mu_gen << endl;

//if( abs(xz_e)<4. and abs(yz_e)<4. and abs(xz_mu)<4. and abs(yz_mu)<4.){

if(mu_gen>=0){
N_input_events++;

N_input_stubs += TrackerStubs->GetEntries();

h_nTotStubs->Fill(TrackerStubs->GetEntries());

std::array<int,12> nstubs = {0,0,0,0,0,0,0,0,0,0,0,0};


  int nStubs_0 = 0;
  int nStubs_1 = 0;

  int nTotStubs = TrackerStubs->GetEntries();


   for(int t=0; t<TrackerStubs->GetEntries(); t++)
   {const MUonETrackerStub *stubs = static_cast<const MUonETrackerStub*>(TrackerStubs->At(t));

 	 int imod = stubs->moduleID()+stubs->stationID()*6;
	 nstubs.at(imod)++;
	 h_nstubsPerModule->Fill(imod);
	 if (imod <6) {
	  nStubs_0++;
	  h_nstubsPerModule_S0->Fill(imod);}
 	 else {
	  nStubs_1++;
	  h_nstubsPerModule_S1->Fill(imod);
	    }
//cout << "stubs->moduleID() " << stubs->moduleID() << " stubs->stationID() " << stubs->stationID() << " -> imod " << imod << endl;
    }// end loop on all stubs in the two stations

//for(int s=0; s<12; s++){ cout <<" mod " << s << ": " << nstubs.at(s) << endl;}

    N_input_stubs_0 += nStubs_0;
    N_input_stubs_1 += nStubs_1;

    h_nStubs_0->Fill(nStubs_0);
    h_nStubs_1->Fill(nStubs_1);
    h2_nStubs->Fill(nStubs_0, nStubs_1);


    //////////////////////////////////////////////////////////////
    // DEFINE the structures specifying the stations' hit patterns
    //////////////////////////////////////////////////////////////
    Station_hits S_0 = fill_station_hits(nstubs,0);
    Station_hits S_1 = fill_station_hits(nstubs,6);
    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////

    bool ok_skim = false;

    // conditions on station 0
    bool is_zero_0 = nStubs_0 == 0;
    if (is_zero_0) {
 N_zero_0++;
 h_nTotStubs_zero_0->Fill(nTotStubs);
    }

    bool is_noise_0 = (S_0.fired_xy_modules + S_0.fired_uv_modules)<=1;
    if (is_noise_0) N_noise_0++;


    bool is_passing_golden_0 = nStubs_0 == 6 && S_0.fired_xy_modules == 4 && S_0.fired_uv_modules == 2;
    if (is_passing_golden_0) {
 N_passing_golden_0++;
    }


    bool is_trackable_0 = S_0.fired_xy_modules == 4 && S_0.fired_uv_modules >0;
    if (is_trackable_0) N_trackable_0++;
    bool is_single_cand_0 = is_trackable_0 && (S_0.multifired_xy_modules + S_0.multifired_uv_modules)<=1;
    if (is_single_cand_0) N_single_cand_0++;
    bool is_single_clean_0 = is_single_cand_0 && nStubs_0 <8;
    if (is_single_clean_0) N_single_clean_0++;

    bool is_pileup_any_0 = S_0.multifired_xy_modules == 4 && S_0.nhits_uv >=2;
    if (is_pileup_any_0) N_pileup_any_0++;

    bool is_pileup_2mu_0 = S_0.multifired_xy_modules == 4 && S_0.nhits_uv >=2 &&
 (S_0.threefired_xy_modules + S_0.threefired_uv_modules)<=1 &&
 (S_0.nhits_xy+S_0.nhits_uv)<=13 ;
    if (is_pileup_2mu_0) N_pileup_2mu_0++;

    bool is_pileup_3mu_0 = S_0.threefired_xy_modules == 4 && S_0.nhits_uv >=3 &&
 (S_0.fourfired_xy_modules + S_0.fourfired_uv_modules)<=1 &&
 (S_0.nhits_xy+S_0.nhits_uv)<=19 ;
    if (is_pileup_3mu_0) N_pileup_3mu_0++;

    bool is_pileup_4mu_0 = S_0.fourfired_xy_modules == 4 && S_0.nhits_uv >=4 &&
 (S_0.fivefired_xy_modules + S_0.fivefired_uv_modules)<=1 &&
 (S_0.nhits_xy+S_0.nhits_uv)<=25 ;
    if (is_pileup_4mu_0) N_pileup_4mu_0++;

    // (4 or more muons)
    bool is_pileup_many_0 = S_0.fourfired_xy_modules == 4 && S_0.nhits_uv >=4;  // (nhits >=20)
    if (is_pileup_many_0) N_pileup_many_0++;

    // conditions on station 1
    //
    bool is_zero_1 = nStubs_1 == 0;
    if (is_zero_1) {
 N_zero_1++;
 h_nTotStubs_zero_1->Fill(nTotStubs);
    }

    bool is_noise_1 = (S_1.fired_xy_modules + S_1.fired_uv_modules)<=1;
    if (is_noise_1) N_noise_1++;

    bool is_passing_golden_1 = nStubs_1 == 6 && S_1.fired_xy_modules == 4 && S_1.fired_uv_modules == 2;
    if (is_passing_golden_1) {
 N_passing_golden_1++;
    }


    bool is_trackable_1 = S_1.fired_xy_modules == 4 && S_1.fired_uv_modules >0;
    bool is_passing_cand_1 = is_trackable_1 && (S_1.multifired_xy_modules + S_1.multifired_uv_modules)<=1;
    if (is_passing_cand_1) N_passing_cand_1++;
    bool is_passing_clean_1 = is_passing_cand_1 && nStubs_1 <8;
    if (is_passing_clean_1) N_passing_clean_1++;

    bool is_2nd_pattern_1 = (S_1.multifired_xy_modules>1 && S_1.nhits_uv >1) ||
   (S_1.multifired_xy_modules>2 && S_1.nhits_uv >0);
    bool is_2tracks_1 = is_trackable_1 && is_2nd_pattern_1;
    if (is_2tracks_1) N_2tracks_1++;

    bool is_golden_1 = nStubs_1 == 12 && S_1.multifired_xy_modules == 4 && S_1.multifired_uv_modules == 2;
    if (is_golden_1) N_golden_1++;

    bool is_pileup_2mu_1 = S_1.multifired_xy_modules == 4 && S_1.nhits_uv >=2 &&
 (S_1.threefired_xy_modules + S_1.threefired_uv_modules)<=1 &&
 (S_1.nhits_xy+S_1.nhits_uv)<=13 ;
    if (is_pileup_2mu_1) N_pileup_2mu_1++;

    bool is_pileup_3mu_1 = S_1.threefired_xy_modules == 4 && S_1.nhits_uv >=3 &&
 (S_1.fourfired_xy_modules + S_1.fourfired_uv_modules)<=1 &&
 (S_1.nhits_xy+S_1.nhits_uv)<=19 ;
    if (is_pileup_3mu_1) N_pileup_3mu_1++;

    bool is_pileup_4mu_1 = S_1.fourfired_xy_modules == 4 && S_1.nhits_uv >=4 &&
 (S_1.fivefired_xy_modules + S_1.fivefired_uv_modules)<=1 &&
 (S_1.nhits_xy+S_1.nhits_uv)<=25 ;
    if (is_pileup_4mu_1) N_pileup_4mu_1++;

    bool is_signal_plus_pileup_1 =
 S_1.multifired_xy_modules == 4 &&
 ((S_1.threefired_xy_modules >=2 && S_1.nhits_uv >= 3)
 ||
  (S_1.threefired_xy_modules >=3 && S_1.nhits_uv >= 2));
    if (is_signal_plus_pileup_1) N_signal_plus_pileup_1++;

    bool is_signal_plus_tripileup_1 =
 S_1.threefired_xy_modules == 4 &&
 ((S_1.fourfired_xy_modules >=2 && S_1.nhits_uv >= 4)
  ||
 (S_1.fourfired_xy_modules >=3 && S_1.nhits_uv >= 3));
    if (is_signal_plus_tripileup_1) N_signal_plus_tripileup_1++;

    bool is_signal_plus_fourpileup_1 =
 S_1.fourfired_xy_modules == 4 &&
 ((S_1.fivefired_xy_modules >=2 && S_1.nhits_uv >= 5)
  ||
 (S_1.fivefired_xy_modules >=3 && S_1.nhits_uv >= 4));
    if (is_signal_plus_fourpileup_1) N_signal_plus_fourpileup_1++;


    // condition on S0+S1
    //
    bool is_zero_event = is_zero_0 && is_zero_1;
    if (is_zero_event) N_zero_events++;

    bool is_noise_event = is_noise_0 && is_noise_1;
    if (is_noise_event) {
 N_noise_events++;
 h_nTotStubs_noise->Fill(nTotStubs);
    }

    bool is_missing_0_event = is_noise_0 &&  is_trackable_1;
    if (is_missing_0_event) {
 N_missing_0_events++;
 h_nTotStubs_missing_0->Fill(nTotStubs);
    }

    bool is_missing_1_event = is_trackable_0 &&  is_noise_1;
    if (is_missing_1_event) {
 N_missing_1_events++;
 h_nTotStubs_missing_1->Fill(nTotStubs);
    }

    bool is_not_trackable_event = !is_trackable_0 || !is_trackable_1;
    if (is_not_trackable_event) {
 N_not_trackable_events++;
 h_nTotStubs_not_trackable->Fill(nTotStubs);
    }

    bool is_passing_golden_event = is_passing_golden_0 && is_passing_golden_1;
    if (is_passing_golden_event) {
 N_passing_golden_events++;
 h_nTotStubs_passing_golden->Fill(nTotStubs);
    }

    bool is_passing_golden_S0_zero_S1_event = is_passing_golden_0 && is_zero_1;
    if (is_passing_golden_S0_zero_S1_event) N_passing_golden_S0_zero_S1_events++;
    bool is_passing_golden_S1_zero_S0_event = is_passing_golden_1 && is_zero_0;
    if (is_passing_golden_S1_zero_S0_event) N_passing_golden_S1_zero_S0_events++;

    bool is_passing_golden_S0_noise_S1_event = is_passing_golden_0 && is_noise_1;
    if (is_passing_golden_S0_noise_S1_event) N_passing_golden_S0_noise_S1_events++;
    bool is_passing_golden_S1_noise_S0_event = is_passing_golden_1 && is_noise_0;
    if (is_passing_golden_S1_noise_S0_event) N_passing_golden_S1_noise_S0_events++;

    bool is_passing_golden_S0_trackable_S1_event = is_passing_golden_0 && is_trackable_1;
    if (is_passing_golden_S0_trackable_S1_event) {
 N_passing_golden_S0_trackable_S1_events++;
    }
    bool is_passing_golden_S1_trackable_S0_event = is_passing_golden_1 && is_trackable_0;
    if (is_passing_golden_S1_trackable_S0_event) {
 N_passing_golden_S1_trackable_S0_events++;
    }


    bool is_passing_1mu_event = is_single_clean_0 && is_passing_clean_1;
    if (is_passing_1mu_event) {
 N_passing_1mu_events++;
 h_nTotStubs_passing_1mu->Fill(nTotStubs);
    }
    bool is_passing_2mu_event = is_pileup_2mu_0 && is_pileup_2mu_1;
    if (is_passing_2mu_event) {
 N_passing_2mu_events++;
 h_nTotStubs_passing_2mu->Fill(nTotStubs);
    }
    bool is_passing_3mu_event = is_pileup_3mu_0 && is_pileup_3mu_1;
    if (is_passing_3mu_event) {
 N_passing_3mu_events++;
 h_nTotStubs_passing_3mu->Fill(nTotStubs);
    }
    bool is_passing_4mu_event = is_pileup_4mu_0 && is_pileup_4mu_1;
    if (is_passing_4mu_event) N_passing_4mu_events++;

    bool is_single_cand_event = is_single_cand_0 && is_2tracks_1;
    if (is_single_cand_event) {
 N_single_cand_events++;
 h_nTotStubs_single_cand->Fill(nTotStubs);
    }

    bool is_single_clean_event = is_single_clean_0 && is_2tracks_1;
    if (is_single_clean_event) {
 N_single_clean_events++;
 h_nTotStubs_single_clean->Fill(nTotStubs);
 if (iskim == 1) ok_skim = true;
    }

    bool is_umsel_event = nStubs_0 >=5 && nStubs_1 >=5 && (nStubs_1 - nStubs_0) >= 5;
    if (is_umsel_event) {
 N_umsel_events++;
 h_nTotStubs_umsel->Fill(nTotStubs);
    }

    bool is_golden_event = is_passing_golden_0 && is_golden_1;
    if (is_golden_event) {
 N_golden_events++;
 h_nTotStubs_golden->Fill(nTotStubs);
    }

    bool is_pileup_any_event = is_pileup_any_0 && is_signal_plus_pileup_1;
    if (is_pileup_any_event) {
 N_pileup_any_events++;
 h_nTotStubs_pileup_any->Fill(nTotStubs);
    }    
    bool is_pileup_2mu_event = is_pileup_2mu_0 && is_signal_plus_pileup_1;
    if (is_pileup_2mu_event) {
 N_pileup_2mu_events++;
 h_nTotStubs_pileup_2mu->Fill(nTotStubs);
    }
    bool is_pileup_3mu_event = is_pileup_3mu_0 && is_signal_plus_tripileup_1;
    if (is_pileup_3mu_event) {
 N_pileup_3mu_events++;
 h_nTotStubs_pileup_3mu->Fill(nTotStubs);
    }
    bool is_pileup_4mu_event = is_pileup_4mu_0 && is_signal_plus_fourpileup_1;
    if (is_pileup_4mu_event) {
 N_pileup_4mu_events++;
 h_nTotStubs_pileup_4mu->Fill(nTotStubs);
    }
    bool is_pileup_many_event = is_pileup_many_0 && is_signal_plus_fourpileup_1;
    if (is_pileup_many_event) {
 N_pileup_many_events++;
 h_nTotStubs_pileup_many->Fill(nTotStubs);
    }

    bool is_pileup_skim_event = is_pileup_2mu_event || is_pileup_3mu_event || is_pileup_many_event;
    if (is_pileup_skim_event) {
 N_pileup_skim_events++;
 h_nTotStubs_pileup_skim->Fill(nTotStubs);
 if (iskim == 40) ok_skim = true;
    }

    bool is_loose_event = is_single_clean_event || is_pileup_2mu_event || is_pileup_3mu_event || is_pileup_many_event;
    if (is_loose_event) {
 N_loose_events++;
 if (iskim == 41) ok_skim = true;
    }

 if (ok_skim) {
/*   for(int t=0; t<TrackerStubs->GetEntries(); t++)
   {const MUonETrackerStub *stubs = static_cast<const MUonETrackerStub*>(TrackerStubs->At(t));
if(stubs->stationID()==1){
        const vector<MUonETrackDeposit> dep= stubs->seedingClusterLinkedTrackDeposits();
	for(int d=0; d<dep.size(); d++) {cout << d << ") seedingClusterLinkedTrackDeposits trackID" << dep.at(d).trackID() <<endl;
                for(int n = 0; n < MCTrack->GetEntries(); n++) {
         const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
         if(n==dep.at(d).trackID()){cout << "MCTr->interactionID() " << MCTr->interactionID() << endl;}
                }
         }
	}
    }*/

h_nTotStubs_after->Fill(TrackerStubs->GetEntries());

vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
         for (auto&& track : tracks) {
//if(e_out==track.linkedTrackID()){process_st0->Fill(track.processIDofLinkedTrack());}
//if(mu_out==track.linkedTrackID()){process_st1->Fill(track.processIDofLinkedTrack());}
if(track.sector()==0){process_st0->Fill(track.processIDofLinkedTrack());}
if(track.sector()==1){process_st1->Fill(track.processIDofLinkedTrack());}
}

N_stubs_skim += nTotStubs;
o_tree->Fill();

}

		}//if mu_in
	}//end for loop

TCanvas s("s","s",700,700);
s.Divide(1,2);
s.cd(1);
process_st0->Draw("hist");
s.cd(2);
process_st1->Draw("hist");
s.SaveAs("process_ALLmu.pdf");

TCanvas a("a","a",700,700);
h_nTotStubs->Draw("hist");
h_nTotStubs->SetTitle(" ");
h_nTotStubs_after->SetLineColor(kRed);
h_nTotStubs_after->Draw("hist same");
gPad->SetLogy();
gStyle->SetOptStat(0);
a.SaveAs("nstubs_ALLmu.pdf");

  cout <<"\n"<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
  cout<<"\n Station 0:"<<endl;
  cout<<"Number of events without any hit in S_0 = "<<endl;
  cout<<"\t" << N_zero_0 <<endl; 
  cout<<"Number of Empty or Noise events in S_0 = "<<endl;
  cout<<"\t" << N_noise_0 <<endl; 
  cout<<"Number of golden events in Station 0 (all 6 modules with one hit) = " <<endl;
  cout<<"\t" << N_passing_golden_0 <<endl;
  cout<<"Number of trackable events (all 4 XY modules and at least 1 UV) = "<<endl;
  cout<<"\t "<<N_trackable_0<<endl;
  cout<<"Number of candidate single mu events in S0 = "<<endl;
  cout<<"\t "<<N_single_cand_0 <<endl;
  cout<<"Number of clean single mu events in S_0 = "<<endl;
  cout<<"\t "<<N_single_clean_0 <<endl;
  cout<<"Number of (any) pileup events in S_0 = "<<endl;
  cout<<"\t "<<N_pileup_any_0 <<endl;
  cout<<"Number of 2mu pileup events in S_0 = "<<endl;
  cout<<"\t "<<N_pileup_2mu_0 <<endl;  
  cout<<"Number of 3mu pileup events in S_0 = "<<endl;
  cout<<"\t "<<N_pileup_3mu_0 <<endl;  
   cout<<"Number of >=4mu pileup events in S_0 = "<<endl;
  cout<<"\t "<<N_pileup_many_0 <<endl;  
 
  cout<<"\n Station 1:"<<endl;
  cout<<"Number of events without any hit in S_1 = "<<endl;
  cout<<"\t" << N_zero_1 <<endl; 
  cout<<"Number of Empty or Noise events in S_1 = "<<endl;
  cout<<"\t" << N_noise_1 <<endl; 
  cout<<"Number of passing events with one hit in each module in S_1 = " <<endl;
  cout<<"\t" << N_passing_golden_1 <<endl;
  cout<<"Number of candidate single mu events in S_1 = "<<endl;
  cout<<"\t "<<N_passing_cand_1 <<endl;
  cout<<"Number of clean single mu events in S_1 = "<<endl;
  cout<<"\t "<<N_passing_clean_1 <<endl;
  cout<<"Number of interaction candidates (2 tracks) in S_1 = "<<endl;
  cout<<"\t "<<N_2tracks_1<<endl;
  cout<<"Number of golden events in S_1 (2 hits in all modules in S_1)= " <<endl;
  cout<<"\t" << N_golden_1 <<endl;
  cout<<"Number of signal+pileup candidates in S_1 = "<<endl;
  cout<<"\t "<<N_signal_plus_pileup_1<<endl;
  
  cout<<"\n EVENTS:   "<<endl;
  cout<<"Number of events without any hit in both S_0 and S_1 = "<<endl;
  cout<<"\t" << N_zero_events <<endl; 
  cout<<"Number of Noise events (empty or <=1 module in both stations) = " << endl;
  cout<<"\t "<< N_noise_events <<endl;
  cout<<"Number of events missing S_0 (empty or <=1 module) and trackable in S_1 = " << endl;
  cout<<"\t "<< N_missing_0_events <<endl;
  cout<<"Number of events missing S_1 (empty or <=1 module) and trackable in S_0 = " << endl;
  cout<<"\t "<< N_missing_1_events <<endl;
  cout<<"Number of events not trackable (in S_0 or S_1) = "<<endl;
  cout<<"\t "<< N_not_trackable_events <<endl;
  cout<<"Number of clean 1mu passing events = " <<endl;
  cout<<"\t "<< N_passing_1mu_events <<endl;
  cout<<"Number of golden 1mu passing events = " << endl;
  cout<<"\t "<< N_passing_golden_events <<endl;

  cout<<"Number of passing 1mu golden S_0 and trackable in S_1 = " << endl;
  cout<<"\t "<< N_passing_golden_S0_trackable_S1_events <<endl;
  cout<<"\tNumber of passing 1mu golden S_0 and zero hits in S_1 = " << endl;
  cout<<"\t\t "<< N_passing_golden_S0_zero_S1_events <<endl;
  cout<<"\tNumber of passing 1mu golden S_0 and <=1 module hits in S_1 = " << endl;
  cout<<"\t\t "<< N_passing_golden_S0_noise_S1_events <<endl;
  //
  cout<<"Number of passing 1mu golden S_1 and trackable in S_0 = " << endl;
  cout<<"\t "<< N_passing_golden_S1_trackable_S0_events <<endl;
  cout<<"\tNumber of passing 1mu golden S_1 and zero hits in S_0 = " << endl;
  cout<<"\t\t "<< N_passing_golden_S1_zero_S0_events <<endl;
  cout<<"\tNumber of passing 1mu golden S_1 and <=1 module hits in S_0 = " << endl;
  cout<<"\t\t "<< N_passing_golden_S1_noise_S0_events <<endl;

  cout<<"Number of candidate Single mu interaction events = " << endl;
  cout<<"\t "<< N_single_cand_events <<endl;
  cout<<"Number of clean Single mu interaction events = " << endl;
  cout<<"\t "<< N_single_clean_events <<endl;
  cout<<"Number of golden Single mu interaction events = " << endl;
  cout<<"\t "<< N_golden_events <<endl;
  cout<<"Number of events passing Umbertos selection:  N0>=5 && N1>=5 && (N1-N0)>=5  ="<<endl;
  cout<<"\t "<< N_umsel_events <<endl;
  cout<<"Number of any pileup interaction events = " << endl;
  cout<<"\t "<< N_pileup_any_events <<endl;
  cout<<"Number of 2mu pileup interaction events = " << endl;
  cout<<"\t "<< N_pileup_2mu_events <<endl;
  cout<<"Number of 3mu pileup interaction events = " << endl;
  cout<<"\t "<< N_pileup_3mu_events <<endl;
  cout<<"Number of >=4mu pileup interaction events = " << endl;
  cout<<"\t "<< N_pileup_many_events <<endl;
  cout<<"Total number of pileup (2,3,>=4) interaction events = " << endl;
  cout<<"\t "<< N_pileup_skim_events <<endl;
  cout<<"Total number of Single + Pileup interaction events = " << endl;
  cout<<"\t "<< N_loose_events <<endl;
  cout <<      "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;

  cout<<"N INPUT Events    = "<< N_input_events <<endl;
  cout<<"N input stubs S_0 = "<< N_input_stubs_0 <<endl;  
  cout<<"N input stubs S_1 = "<< N_input_stubs_1 <<endl;  
  cout<<"N input stubs tot = "<< N_input_stubs <<endl;  
  cout<<endl;
  cout<<"Fraction missing S_0             = "<< double(N_missing_0_events)/double(N_input_events) <<endl;
  cout<<"Fraction missing S_1             = "<< double(N_missing_1_events)/double(N_input_events) <<endl;
  cout<<"Fraction not fully trackable     = "<< double(N_not_trackable_events)/double(N_input_events) <<endl;
  cout<<"Fraction Passing Golden mu       = "<< double(N_passing_golden_events)/double(N_input_events) <<endl;
//line 689 Giovanni's skimming

  pair<double,double> ineff_zero_S1 = binomial_eff(N_passing_golden_S0_zero_S1_events, N_passing_golden_0);
  pair<double,double> ineff_noise_S1 = binomial_eff(N_passing_golden_S0_noise_S1_events, N_passing_golden_0);
  pair<double,double> eff_S1 = binomial_eff(N_passing_golden_S0_trackable_S1_events, N_passing_golden_0);

  pair<double,double> ineff_zero_S0 = binomial_eff(N_passing_golden_S1_zero_S0_events, N_passing_golden_1);
  pair<double,double> ineff_noise_S0 = binomial_eff(N_passing_golden_S1_noise_S0_events, N_passing_golden_1);
  pair<double,double> eff_S0 = binomial_eff(N_passing_golden_S1_trackable_S0_events, N_passing_golden_1);

  cout<<"\t estimated Efficiency reconstruction S_1      = "<< eff_S1.first << " +- " << eff_S1.second <<endl; 
  cout<<"\t\t estimated inefficiency (zero hits in S_1)  = "<< ineff_zero_S1.first << " +- " << ineff_zero_S1.second <<endl; 
  cout<<"\t\t estimated inefficiency (<=1 module in S_1) = "<< ineff_noise_S1.first << " +- " << ineff_noise_S1.second <<endl; 
  cout<<"\t estimated Efficiency reconstruction S_0    = "<< eff_S0.first << " +- " << eff_S0.second<<endl;
  cout<<"\t\t estimated inefficiency (zero hits in S_0)  = "<< ineff_zero_S0.first << " +- " << ineff_zero_S0.second <<endl; 
  cout<<"\t\t estimated inefficiency (<=1 module in S_0) = "<< ineff_noise_S0.first << " +- " << ineff_noise_S0.second <<endl; 

  ofstream ofeff1("relative_efficiency_S1_goodS0_"+n_mult+".txt");
  ofeff1<<"RELATIVE EFFICIENCY Station S_1 (in events with Golden Passing Mu in S_0)"<<endl;
  ofeff1<< N_passing_golden_0 <<endl;
  ofeff1<< N_passing_golden_S0_trackable_S1_events << endl;
  ofeff1<< eff_S1.first << " +- " << eff_S1.second <<endl;
  ofeff1<< N_passing_golden_S0_noise_S1_events << endl;
  ofeff1<< ineff_noise_S1.first << " +- " << ineff_noise_S1.second <<endl;
  ofeff1<< N_passing_golden_S0_zero_S1_events << endl;
  ofeff1<< ineff_zero_S1.first << " +- " << ineff_zero_S1.second <<endl;
  ofeff1.close();
  
  ofstream ofeff0("relative_efficiency_S0_goodS1_"+n_mult+".txt");
  ofeff0<<"RELATIVE EFFICIENCY Station S_0 (in events with Golden Passing Mu in S_1)"<<endl;
  ofeff0<< N_passing_golden_1 <<endl;
  ofeff0<< N_passing_golden_S1_trackable_S0_events << endl;
  ofeff0<< eff_S0.first << " +- " << eff_S0.second <<endl;
  ofeff0<< N_passing_golden_S1_noise_S0_events << endl;
  ofeff0<< ineff_noise_S0.first << " +- " << ineff_noise_S0.second <<endl;
  ofeff0<< N_passing_golden_S1_zero_S0_events << endl;
  ofeff0<< ineff_zero_S0.first << " +- " << ineff_zero_S0.second <<endl;
  ofeff0.close();

  cout<<"Fraction Passing 1mu             = "<< double(N_passing_1mu_events)/double(N_input_events) <<endl;
  cout<<"Fraction Single Mu Int.          = "<< double(N_single_clean_events)/double(N_input_events) <<endl;
  cout<<"Fraction Golden sel              = "<< double(N_golden_events)/double(N_input_events) <<endl;
  cout<<"Fraction Umberto's sel           = "<< double(N_umsel_events)/double(N_input_events) <<endl;
  cout<<"Fraction 2 Mu PU Int.            = "<< double(N_pileup_2mu_events)/double(N_input_events) <<endl;
  cout<<"Fraction 3 Mu PU Int.            = "<< double(N_pileup_3mu_events)/double(N_input_events) <<endl;
  cout<<"Fraction >=4 Mu PU Int.          = "<< double(N_pileup_many_events)/double(N_input_events) <<endl;
  cout<<"Total Fraction PU Mu Int.        = "<< double(N_pileup_skim_events)/double(N_input_events) <<endl;
  cout<<"Total Fraction Single+PU Mu Int. = "<< double(N_loose_events)/double(N_input_events) <<endl;

hfile->Write();

  if (iskim > 0) {
    o_file->Write();
  }


}//end of file
