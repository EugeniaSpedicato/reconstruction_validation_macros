#include <array>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"

#include "Station_hits.h"
#include "releff.h"

using namespace std;


void fifth_new_skim_eff97(){
int NMODULES = 12;
int iskim = 61;

///////
  string filelist = "runs/";
  string inputdir = "inputs/";
  string n_mult="fifth_ALLPUmu_eff97";
  string outdir = "outputs/";

    filelist = filelist  + "files.txt";

TChain * cbmsim = new TChain("cbmsim");
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/minbias_1M_new_085PUmean.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/gen_digi/commit_2f4e96f4_MCsignal_SIM-DIGI.root");

  //////////////////////////////////////////////////////////////////////
  auto n_entries = cbmsim->GetEntries();
  cout << "Total number of input events  = " << n_entries << endl;

   TClonesArray *MCTrack = 0;
   TClonesArray *TrackerStubs = 0;
   //MUonERecoOutput *ReconstructionOutput = 0;

   cbmsim->SetBranchAddress("MCTrack", &MCTrack);
   cbmsim->SetBranchAddress("TrackerStubs", &TrackerStubs);
   //cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);

  //////////////////////////////////////////////////////////////////////
   TClonesArray *o_MCTrack = 0;
   TClonesArray *o_TrackerStubs = 0;
   //MUonERecoOutput *o_ReconstructionOutput = 0;

TFile *o_file;
  TTree *o_tree;

  Long64_t o_nevt = 0;
if (iskim > 0) {
    //
    // OUTPUT
    //
    // max size of skimmed output files: 2 GB
    TTree::SetMaxTreeSize(2000000000LL);
    //  TTree::SetMaxTreeSize(10000000LL);
    o_file = TFile::Open((outdir+n_mult+".root").c_str(), "recreate");
    if (!o_file) {cout << "Error opening file" << std::endl; exit(-1);}
    o_tree = new TTree("cbmsim", "");
    o_tree->Branch("MCTrack", &o_MCTrack);
    o_tree->Branch("TrackerStubs", &o_TrackerStubs);
    //o_tree->Branch("ReconstructionOutput", &o_ReconstructionOutput);
        }

  // histograms: find X range (file number - SuperId - time)
  //
  // assume maximum 70 stubs per event in the two stations
  string hname = outdir+"/histos_"+n_mult+".root";
  TFile *hfile = new TFile(hname.c_str(),"RECREATE");
    const Double_t nmax = 70;
  Int_t nb = nmax;

  // stub counts per module and station 
  TH1D *h_nstubsPerModule = new TH1D("h_nstubsPerModule","nstubs on modules",NMODULES,0,NMODULES);
  TH1D *h_nstubsPerModule_S0 = new TH1D("h_nstubsPerModule_S0","nstubs on modules - Station 0",NMODULES,0,NMODULES);
  TH1D *h_nstubsPerModule_S1 = new TH1D("h_nstubsPerModule_S1","nstubs on modules - Station 1",NMODULES,0,NMODULES);
  // for different selections
  TH1D *h_nstubsPerModule_presel = new TH1D("h_nstubsPerModule_presel","nstubs on modules (preselected events)",NMODULES,0,NMODULES);
  TH1D *h_nstubsPerModule_trackable_0_onehit_1 = new TH1D("h_nstubsPerModule_trackable_0_onehit_1","nstubs on modules (S0_trackable_S1_onehit events))",NMODULES,0,NMODULES);
  TH1D *h_nstubsPerModule_trackable_1_onehit_0 = new TH1D("h_nstubsPerModule_trackable_1_onehit_0","nstubs on modules (S1_trackable_S0_onehit events))",NMODULES,0,NMODULES);
  TH1D *h_nstubsPerModule_trackable = new TH1D("h_nstubsPerModule_trackable","nstubs on modules (trackable events)",NMODULES,0,NMODULES);
  TH1D *h_nstubsPerModule_fired12 = new TH1D("h_nstubsPerModule_fired12","nstubs on modules (fired12 events)",NMODULES,0,NMODULES);
  TH1D *h_nstubsPerModule_single_clean = new TH1D("h_nstubsPerModule_single_clean","nstubs on modules (clean single mu interaction events)",NMODULES,0,NMODULES);
  TH1D *h_nstubsPerModule_pileup234_skim = new TH1D("h_nstubsPerModule_pileup234_skim","nstubs on modules (pileup234_skim events)",NMODULES,0,NMODULES);
  TH1D *h_nstubsPerModule_fired12_plus_single_clean = new TH1D("h_nstubsPerModule_fired12_plus_single_clean","nstubs on modules (fired12_plus_single_clean events)",NMODULES,0,NMODULES);

  // stub multiplicities for different selections
  TH1D *h_nTotStubs = new TH1D("h_nTotStubs","Total N of Stubs; Number of Stubs; Events",nb,0,nmax);
  TH2D *h2_nStubs = new TH2D("h2_nStubs","N stubs on S1 vs S0",nb,0,nmax,nb,0,nmax);
  TH1D *h_nStubs_0 = new TH1D("h_nStubs_0","Tot N stubs First Station - S0",nb,0,nmax);
  TH1D *h_nStubs_1 = new TH1D("h_nStubs_1","Tot N stubs Second Station - S1",nb,0,nmax);

  TH1D *h_nTotStubs_presel = new TH1D("h_nTotStubs_presel","Total N of Stubs (preselected events); Number of Stubs; Events",nb,0,nmax);
  TH2D *h2_nStubs_presel = new TH2D("h2_nStubs_presel","N stubs on S1 vs S0 (preselected events)",nb,0,nmax,nb,0,nmax);
  TH1D *h_nStubs_0_presel = new TH1D("h_nStubs_0_presel","Tot N stubs First Station - S0 (preselected events)",nb,0,nmax);
  TH1D *h_nStubs_1_presel = new TH1D("h_nStubs_1_presel","Tot N stubs Second Station - S1 (preselected events)",nb,0,nmax);

  TH1D *h_nTotStubs_S0_trackable_S1_onehit = new TH1D("h_nTotStubs_S0_trackable_S1_onehit","Total N of Stubs (S0_trackable_S1_onehit events",nb,0,nmax);
  TH2D *h2_nStubs_S0_trackable_S1_onehit = new TH2D("h2_nStubs_S0_trackable_S1_onehit","N stubs on S1 vs S0 (S0_trackable_S1_onehit events)",nb,0,nmax,nb,0,nmax);
  TH1D *h_nStubs_0_S0_trackable_S1_onehit = new TH1D("h_nStubs_0_S0_trackable_S1_onehit","Tot N stubs First Station (S0_trackable_S1_onehit events)",nb,0,nmax);
  TH1D *h_nStubs_1_S0_trackable_S1_onehit = new TH1D("h_nStubs_1_S0_trackable_S1_onehit","Tot N stubs Second Station (S0_trackable_S1_onehit events)",nb,0,nmax);

  TH1D *h_nTotStubs_S1_trackable_S0_onehit = new TH1D("h_nTotStubs_S1_trackable_S0_onehit","Total N of Stubs (S1_trackable_S0_onehit events",nb,0,nmax);
  TH2D *h2_nStubs_S1_trackable_S0_onehit = new TH2D("h2_nStubs_S1_trackable_S0_onehit","N stubs on S1 vs S0 (S1_trackable_S0_onehit events)",nb,0,nmax,nb,0,nmax);
  TH1D *h_nStubs_0_S1_trackable_S0_onehit = new TH1D("h_nStubs_0_S1_trackable_S0_onehit","Tot N stubs First Station (S1_trackable_S0_onehit events)",nb,0,nmax);
  TH1D *h_nStubs_1_S1_trackable_S0_onehit = new TH1D("h_nStubs_1_S1_trackable_S0_onehit","Tot N stubs Second Station (S1_trackable_S0_onehit events)",nb,0,nmax);

  TH1D *h_nTotStubs_trackable = new TH1D("h_nTotStubs_trackable","Total N of Stubs (trackable events); Number of Stubs; Events",nb,0,nmax);
  TH2D *h2_nStubs_trackable = new TH2D("h2_nStubs_trackable","N stubs on S1 vs S0 (trackable events)",nb,0,nmax,nb,0,nmax);
  TH1D *h_nStubs_0_trackable = new TH1D("h_nStubs_0_trackable","Tot N stubs First Station - S0 (trackable events)",nb,0,nmax);
  TH1D *h_nStubs_1_trackable = new TH1D("h_nStubs_1_trackable","Tot N stubs Second Station - S1 (trackable events)",nb,0,nmax);
  
  TH1D *h_nTotStubs_zero_0 = new TH1D("h_nTotStubs_zero_0","Total N of Stubs - events with zero hits in S_0",nb,0,nmax);
  TH1D *h_nTotStubs_zero_1 = new TH1D("h_nTotStubs_zero_1","Total N of Stubs - events with zero hits in S_1",nb,0,nmax);
  TH1D *h_nTotStubs_noise = new TH1D("h_nTotStubs_noise","Total N of Stubs noise events",nb,0,nmax);
  TH1D *h_nTotStubs_missing_0 = new TH1D("h_nTotStubs_missing_0","Total N of Stubs muons missing S_0",nb,0,nmax);
  TH1D *h_nTotStubs_missing_1 = new TH1D("h_nTotStubs_missing_1","Total N of Stubs muons missing S_1",nb,0,nmax);
  TH1D *h_nTotStubs_not_trackable = new TH1D("h_nTotStubs_not_trackable","Total N of Stubs - not trackable events",nb,0,nmax);
  TH1D *h_nTotStubs_presel_not_trackable = new TH1D("h_nTotStubs_presel_not_trackable","Total N of Stubs (preselected events) - not trackable events",nb,0,nmax);
  
  TH1D *h_nTotStubs_passing_1mu = new TH1D("h_nTotStubs_passing_1mu","Total N of Stubs for passing mu",nb,0,nmax);
  TH1D *h_nTotStubs_passing_2mu = new TH1D("h_nTotStubs_passing_2mu","Total N of Stubs for passing 2mu",nb,0,nmax);
  TH1D *h_nTotStubs_passing_3mu = new TH1D("h_nTotStubs_passing_3mu","Total N of Stubs for passing 3mu",nb,0,nmax);
  TH1D *h_nTotStubs_passing_4mu = new TH1D("h_nTotStubs_passing_4mu","Total N of Stubs for passing 4mu",nb,0,nmax);

  TH1D *h_nTotStubs_golden = new TH1D("h_nTotStubs_golden","Total N of Stubs for Golden signal selection; Number of Stubs; Events",nb,0,nmax);  
  TH1D *h_nTotStubs_umsel = new TH1D("h_nTotStubs_umsel","Total N of Stubs for Umberto's selection",nb,0,nmax);
  TH1D *h_nTotStubs_single_cand = new TH1D("h_nTotStubs_single_cand","Total N of Stubs for candidate single mu signal",nb,0,nmax);
  TH1D *h_nTotStubs_single_clean = new TH1D("h_nTotStubs_single_clean","Total N of Stubs for clean single mu signal",nb,0,nmax);
  TH2D *h2_nStubs_single_clean = new TH2D("h2_nStubs_single_clean","N stubs on S1 vs S0 (single_clean events)",nb,0,nmax,nb,0,nmax);
  TH1D *h_nStubs_0_single_clean = new TH1D("h_nStubs_0_single_clean","Tot N stubs First Station - S0 (single_clean events)",nb,0,nmax);
  TH1D *h_nStubs_1_single_clean = new TH1D("h_nStubs_1_single_clean","Tot N stubs Second Station - S1 (single_clean events)",nb,0,nmax);

  TH1D *h_nTotStubs_pileup_any = new TH1D("h_nTotStubs_pileup_any","Total N of Stubs for signal+(any)mu pileup",nb,0,nmax);
  TH1D *h_nTotStubs_pileup_2mu = new TH1D("h_nTotStubs_pileup_2mu","Total N of Stubs for 2mu pileup",nb,0,nmax);
  TH1D *h_nTotStubs_pileup_3mu = new TH1D("h_nTotStubs_pileup_3mu","Total N of Stubs for 3mu pileup",nb,0,nmax);
  TH1D *h_nTotStubs_pileup_4mu = new TH1D("h_nTotStubs_pileup_4mu","Total N of Stubs for 4mu pileup",nb,0,nmax);
  TH1D *h_nTotStubs_pileup_many = new TH1D("h_nTotStubs_pileup_many","Total N of Stubs for >=4mu pileup",nb,0,nmax);
  TH1D *h_nTotStubs_pileup_skim = new TH1D("h_nTotStubs_pileup_skim","Total N of Stubs (pileup_skim events); Number of Stubs; Events",nb,0,nmax);
  TH2D *h2_nStubs_pileup_skim = new TH2D("h2_nStubs_pileup_skim","N stubs on S1 vs S0 (pileup_skim events)",nb,0,nmax,nb,0,nmax);
  TH1D *h_nStubs_0_pileup_skim = new TH1D("h_nStubs_0_pileup_skim","Tot N stubs First Station - S0 (pileup_skim events)",nb,0,nmax);
  TH1D *h_nStubs_1_pileup_skim = new TH1D("h_nStubs_1_pileup_skim","Tot N stubs Second Station - S1 (pileup_skim events)",nb,0,nmax);

  TH1D *h_nTotStubs_pileup234_skim = new TH1D("h_nTotStubs_pileup234_skim","Total N of Stubs (pileup234_skim events); Number of Stubs; Events",nb,0,nmax);
  TH2D *h2_nStubs_pileup234_skim = new TH2D("h2_nStubs_pileup234_skim","N stubs on S1 vs S0 (pileup234_skim events)",nb,0,nmax,nb,0,nmax);
  TH1D *h_nStubs_0_pileup234_skim = new TH1D("h_nStubs_0_pileup234_skim","Tot N stubs First Station - S0 (pileup234_skim events)",nb,0,nmax);
  TH1D *h_nStubs_1_pileup234_skim = new TH1D("h_nStubs_1_pileup234_skim","Tot N stubs Second Station - S1 (pileup234_skim events)",nb,0,nmax);

  TH1D *h_nTotStubs_fired12 = new TH1D("h_nTotStubs_fired12","Total N of Stubs (fired12 events); Number of Stubs; Events",nb,0,nmax);
  TH2D *h2_nStubs_fired12 = new TH2D("h2_nStubs_fired12","N stubs on S1 vs S0 (fired12 events)",nb,0,nmax,nb,0,nmax);
  TH1D *h_nStubs_0_fired12 = new TH1D("h_nStubs_0_fired12","Tot N stubs First Station - S0 (fired12 events)",nb,0,nmax);
  TH1D *h_nStubs_1_fired12 = new TH1D("h_nStubs_1_fired12","Tot N stubs Second Station - S1 (fired12 events)",nb,0,nmax);

  TH1D *h_nTotStubs_fired12_plus_single_clean = new TH1D("h_nTotStubs_fired12_plus_single_clean","Total N of Stubs (fired12_plus_single_clean events); Number of Stubs; Events",nb,0,nmax);
  TH2D *h2_nStubs_fired12_plus_single_clean = new TH2D("h2_nStubs_fired12_plus_single_clean","N stubs on S1 vs S0 (fired12_plus_single_clean events)",nb,0,nmax,nb,0,nmax);
  TH1D *h_nStubs_0_fired12_plus_single_clean = new TH1D("h_nStubs_0_fired12_plus_single_clean","Tot N stubs First Station - S0 (fired12_plus_single_clean events)",nb,0,nmax);
  TH1D *h_nStubs_1_fired12_plus_single_clean = new TH1D("h_nStubs_1_fired12_plus_single_clean","Tot N stubs Second Station - S1 (fired12_plus_single_clean events)",nb,0,nmax);


array<TH1D*,12> stub_per_mod_trackable;
for(int m=0; m<12; m++){
string name="stub_per_mod_trackable"+to_string(m);
string title="Number of stubs when N_stubs<15 per module "+to_string(m);
stub_per_mod_trackable.at(m)=new TH1D(name.c_str(),title.c_str(),20,0,20);
}

array<TH1D*,12> stub_per_mod_single_clean;
for(int m=0; m<12; m++){
string name="stub_per_mod_single_clean"+to_string(m);
string title="Number of stubs when _single_clean per module "+to_string(m);
stub_per_mod_single_clean.at(m)=new TH1D(name.c_str(),title.c_str(),20,0,20);
}

array<TH1D*,12> stub_per_mod_pileup234;
for(int m=0; m<12; m++){
string name="stub_per_mod_pileup234"+to_string(m);
string title="Number of stubs when _pileup234 per module "+to_string(m);
stub_per_mod_pileup234.at(m)=new TH1D(name.c_str(),title.c_str(),20,0,20);
}



  Long64_t N_input_events = 0;
  Long64_t N_input_stubs_0 = 0;
  Long64_t N_input_stubs_1 = 0;
  Long64_t N_input_stubs = 0;    // total S0+S1
  Long64_t N_stubs_skim = 0;

  Long64_t N_presel_events = 0;
  Long64_t N_presel_stubs_0 = 0;
  Long64_t N_presel_stubs_1 = 0;
  Long64_t N_presel_stubs = 0;    // total S0+S1

  //  Long64_t N_trackable_0_onehit_1 = 0;
  Long64_t N_trackable_0_onehit_1_stubs_0 = 0;
  Long64_t N_trackable_0_onehit_1_stubs_1 = 0;
  Long64_t N_trackable_0_onehit_1_stubs = 0;  // total S0+S1

  //  Long64_t N_trackable_1_onehit_0 = 0;
  Long64_t N_trackable_1_onehit_0_stubs_0 = 0;
  Long64_t N_trackable_1_onehit_0_stubs_1 = 0;
  Long64_t N_trackable_1_onehit_0_stubs = 0;  // total S0+S1

  //  Long64_t N_single_clean_events = 0;
  Long64_t N_single_clean_stubs_0 = 0;
  Long64_t N_single_clean_stubs_1 = 0;
  Long64_t N_single_clean_stubs = 0;    // total S0+S1

  //  Long64_t N_pileup_skim_events = 0;
  Long64_t N_pileup_skim_stubs_0 = 0;
  Long64_t N_pileup_skim_stubs_1 = 0;
  Long64_t N_pileup_skim_stubs = 0;    // total S0+S1

  //  Long64_t N_pileup234_skim_events = 0;
  Long64_t N_pileup234_skim_stubs_0 = 0;
  Long64_t N_pileup234_skim_stubs_1 = 0;
  Long64_t N_pileup234_skim_stubs = 0;    // total S0+S1
  
  Long64_t N_zero_0 = 0;
  Long64_t N_noise_0 = 0;
  //
  Long64_t N_onehit_0 = 0;
  Long64_t N_twohit_0 = 0;
  Long64_t N_threehit_0 = 0;
  Long64_t N_fourhit_0 = 0;
  Long64_t N_trackable_0 = 0;
  Long64_t N_trackable_0_onehit_1 = 0;
  Long64_t N_trackable_0_twohit_1 = 0;
  Long64_t N_trackable_0_threehit_1 = 0;
  Long64_t N_trackable_0_fourhit_1 = 0;
  Long64_t N_trackable_1_onehit_0 = 0;
  Long64_t N_trackable_events = 0;
  //
  Long64_t N_passing_1mu_golden_0 = 0;
  Long64_t N_passing_2mu_golden_0 = 0;
  Long64_t N_passing_3mu_golden_0 = 0;
  Long64_t N_passing_4mu_golden_0 = 0;
  Long64_t N_single_cand_0 = 0;
  Long64_t N_single_clean_0 = 0;
  Long64_t N_pileup_any_0 = 0;
  Long64_t N_pileup_2mu_0 = 0;
  Long64_t N_pileup_3mu_0 = 0;
  Long64_t N_pileup_4mu_0 = 0;
  Long64_t N_pileup_many_0 = 0;
  
  Long64_t N_zero_1 = 0;
  Long64_t N_noise_1 = 0;
  Long64_t N_passing_1mu_golden_1 = 0;
  Long64_t N_passing_2mu_golden_1 = 0;
  Long64_t N_passing_3mu_golden_1 = 0;
  Long64_t N_passing_4mu_golden_1 = 0;
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
  Long64_t N_presel_not_trackable_events = 0;
  Long64_t N_not_trackable_events = 0;
  Long64_t N_passing_golden_S0_zero_S1_events = 0;
  Long64_t N_passing_golden_S1_zero_S0_events = 0;

  Long64_t N_passing_1mu_events = 0;
  Long64_t N_passing_2mu_events = 0;
  Long64_t N_passing_3mu_events = 0;
  Long64_t N_passing_4mu_events = 0;

  Long64_t N_passing_1mu_golden_events = 0;
  Long64_t N_passing_2mu_golden_events = 0;
  Long64_t N_passing_3mu_golden_events = 0;
  Long64_t N_passing_4mu_golden_events = 0;

  Long64_t N_passing_golden_S0_trackable_S1_events = 0;
  Long64_t N_passing_golden_S1_trackable_S0_events = 0;
  Long64_t N_presel_passing_golden_S0_events = 0;
  Long64_t N_presel_passing_golden_S1_events = 0;
  //
  Long64_t N_presel_passing_2mu_golden_S0_events = 0;
  Long64_t N_presel_passing_2mu_golden_S1_events = 0;
  Long64_t N_presel_passing_2mu_golden_S0_trackable_S1_events = 0;
  Long64_t N_presel_passing_2mu_golden_S1_trackable_S0_events = 0;
  //
  Long64_t N_presel_passing_3mu_golden_S0_events = 0;
  Long64_t N_presel_passing_3mu_golden_S1_events = 0;
  Long64_t N_presel_passing_3mu_golden_S0_trackable_S1_events = 0;
  Long64_t N_presel_passing_3mu_golden_S1_trackable_S0_events = 0;
  //
  Long64_t N_presel_passing_4mu_golden_S0_events = 0;
  Long64_t N_presel_passing_4mu_golden_S1_events = 0;
  Long64_t N_presel_passing_4mu_golden_S0_trackable_S1_events = 0;
  Long64_t N_presel_passing_4mu_golden_S1_trackable_S0_events = 0;

  Long64_t N_golden_events = 0;
  Long64_t N_umsel_events = 0;

  Long64_t N_single_cand_events = 0;
  Long64_t N_single_clean_events = 0;
  Long64_t N_pileup_any_events = 0;
  Long64_t N_pileup_2mu_events = 0;
  Long64_t N_pileup_3mu_events = 0;
  Long64_t N_pileup_4mu_events = 0;
  Long64_t N_pileup_many_events = 0;
  Long64_t N_pileup_skim_events = 0;
  Long64_t N_pileup234_skim_events = 0;
  Long64_t N_loose_events = 0;
  Long64_t N_loose234_events = 0;

  Long64_t N_fired6_0 = 0;
  Long64_t N_fired6_1 = 0;
  Long64_t N_fired12_events = 0;
  Long64_t N_fired12_plus_single_clean_events = 0;


  for(Long64_t i = 0; i < n_entries; i++) {
  //  for(Long64_t i = 0; i < 100001; i++) {

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


        for(int n = 0; n < MCTrack->GetEntries(); n++) {
         const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
         if(MCTr->interactionID()==0 and MCTr->pdgCode()==-13){mu_gen++;}

         if(MCTr->interactionID()==45 and MCTr->pdgCode()==11) { e_out=n;}
         if(MCTr->interactionID()==45 and MCTr->pdgCode()==-13) {mu_out=n;}
        }


if(mu_gen>0){

    N_input_events++;
    std::array<int,12> nstubs = {0,0,0,0,0,0,0,0,0,0,0,0};
    int nTotStubs = 0;  
    int nStubs_0 = 0;
    int nStubs_1 = 0;


  for(int t=0; t<TrackerStubs->GetEntries(); t++)
   {const MUonETrackerStub *stubs = static_cast<const MUonETrackerStub*>(TrackerStubs->At(t));

        int imod = stubs->moduleID()+stubs->stationID()*6;
      int a=1;

//      if(imod == 0 or imod == 4 or imod==6 or imod==10) a=-1;
      float local=a*(stubs->seedClusterCenterStrip() + 0.5 + 0.5 * stubs->bend()) * 9.144 / 1016 - 0.5 * 9.144;

//simulating modules efficiency of 97%
        bool hit= (rand()%100) <= 97;

         if (imod <6 and abs(local)<3 ) {
        nstubs.at(imod)+=static_cast<int>(hit);
        nTotStubs+=static_cast<int>(hit);
          nStubs_0+=static_cast<int>(hit);
         if(static_cast<int>(hit)==1)h_nstubsPerModule->Fill(imod);
          if(static_cast<int>(hit)==1)h_nstubsPerModule_S0->Fill(imod);}
         else if(imod>=6){
        nstubs.at(imod)+=static_cast<int>(hit);
        nTotStubs+=static_cast<int>(hit);
          nStubs_1+=static_cast<int>(hit);
         if(static_cast<int>(hit)==1)h_nstubsPerModule->Fill(imod);
          if(static_cast<int>(hit)==1)h_nstubsPerModule_S1->Fill(imod);
            }
	}

    N_input_stubs += nTotStubs;
    h_nTotStubs->Fill(nTotStubs);

    // end loop on all stubs in the two stations
    //////////////////////// ///////////////////////////////

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

    int fired_modules_0 = S_0.fired_xy_modules + S_0.fired_uv_modules;
    bool is_fired6_0 = fired_modules_0 == 6;
    if (is_fired6_0) N_fired6_0++;
    
    if (fired_modules_0 >= 1) {
      N_onehit_0++;
    }
    if (fired_modules_0 >= 2) {
      N_twohit_0++;
    }
    if (fired_modules_0 >= 3) {
      N_threehit_0++;
    }
    if (fired_modules_0 >= 4) {
      N_fourhit_0++;
    }
            
    bool is_trackable_0 = S_0.fired_xy_modules == 4 && S_0.fired_uv_modules >0;
    if (is_trackable_0) {
      N_trackable_0++;
    }

    bool is_passing_golden_0 = nStubs_0 == 6 && S_0.fired_xy_modules == 4 && S_0.fired_uv_modules == 2;
    if (is_passing_golden_0) {
      N_passing_1mu_golden_0++;
    }

    bool is_passing_2mu_golden_0 = nStubs_0 == 12 && S_0.multifired_xy_modules == 4 && S_0.multifired_uv_modules == 2;
    if (is_passing_2mu_golden_0) {
      N_passing_2mu_golden_0++;
    }
    
    bool is_passing_3mu_golden_0 = nStubs_0 == 18 && S_0.threefired_xy_modules == 4 && S_0.threefired_uv_modules == 2;
    if (is_passing_3mu_golden_0) {
      N_passing_3mu_golden_0++;
    }
    bool is_passing_4mu_golden_0 = nStubs_0 == 24 && S_0.fourfired_xy_modules == 4 && S_0.fourfired_uv_modules == 2;
    if (is_passing_4mu_golden_0) {
      N_passing_4mu_golden_0++;
    }

    bool is_passing_2mu_trackable_0 = S_0.multifired_xy_modules == 4 && S_0.nhits_uv >=2;
    bool is_passing_3mu_trackable_0 = S_0.threefired_xy_modules == 4 && S_0.nhits_uv >=3;
    bool is_passing_4mu_trackable_0 = S_0.fourfired_xy_modules == 4 && S_0.nhits_uv >=4;

    bool is_single_cand_0 = is_trackable_0 && (S_0.multifired_xy_modules + S_0.multifired_uv_modules)<=1;
    if (is_single_cand_0) N_single_cand_0++;
    bool is_single_clean_0 = is_single_cand_0 && nStubs_0 <8;
    if (is_single_clean_0) N_single_clean_0++;

    bool is_fired6_plus_single_clean_0 =
	 is_fired6_0 && S_0.multifired_xy_modules == 4 && S_0.multifired_uv_modules >=1 && nStubs_0 <=13;   
    
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

    int fired_modules_1 = S_1.fired_xy_modules + S_1.fired_uv_modules;
    bool is_fired6_1 = fired_modules_1 == 6;
    if (is_fired6_1) N_fired6_1++;
    
    bool is_passing_golden_1 = nStubs_1 == 6 && S_1.fired_xy_modules == 4 && S_1.fired_uv_modules == 2;
    if (is_passing_golden_1) {
      N_passing_1mu_golden_1++;
    }
    bool is_passing_2mu_golden_1 = nStubs_1 == 12 && S_1.multifired_xy_modules == 4 && S_1.multifired_uv_modules == 2;
    if (is_passing_2mu_golden_1) {
      N_passing_2mu_golden_1++;
    }
    bool is_passing_3mu_golden_1 = nStubs_1 == 18 && S_1.threefired_xy_modules == 4 && S_1.threefired_uv_modules == 2;
    if (is_passing_3mu_golden_1) {
      N_passing_3mu_golden_1++;
    }
    bool is_passing_4mu_golden_1 = nStubs_1 == 24 && S_1.fourfired_xy_modules == 4 && S_1.fourfired_uv_modules == 2;
    if (is_passing_4mu_golden_1) {
      N_passing_4mu_golden_1++;
    }

    bool is_passing_2mu_trackable_1 = S_1.multifired_xy_modules == 4 && S_1.nhits_uv >=2;
    bool is_passing_3mu_trackable_1 = S_1.threefired_xy_modules == 4 && S_1.nhits_uv >=3;
    bool is_passing_4mu_trackable_1 = S_1.fourfired_xy_modules == 4 && S_1.nhits_uv >=4;

    bool is_trackable_1 = S_1.fired_xy_modules == 4 && S_1.fired_uv_modules >0;
    bool is_passing_cand_1 = is_trackable_1 && (S_1.multifired_xy_modules + S_1.multifired_uv_modules)<=1;
    if (is_passing_cand_1) N_passing_cand_1++;
    bool is_passing_clean_1 = is_passing_cand_1 && nStubs_1 <8;
    if (is_passing_clean_1) N_passing_clean_1++;
   
    bool is_2nd_pattern_1 = (S_1.multifired_xy_modules>1 && S_1.nhits_uv >1) ||
                            (S_1.multifired_xy_modules>2 && S_1.nhits_uv >0);   
    bool is_2tracks_1 = is_trackable_1 && is_2nd_pattern_1;
    if (is_2tracks_1) N_2tracks_1++;

    bool is_fired6_plus_trackable_1 = is_fired6_1 && S_1.multifired_xy_modules == 4 && S_1.nhits_uv >= 3;

    bool is_2nd_pattern_overfired6_1 = (S_1.threefired_xy_modules >=2 && S_1.nhits_uv >=4) ||
                                       (S_1.threefired_xy_modules >=3 && S_1.nhits_uv >=3);

    bool is_fired6_2tracks_1 = is_fired6_plus_trackable_1 && is_2nd_pattern_overfired6_1;
      
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

    // at least one hit on both stations
    bool is_presel_event = !is_zero_0 && !is_zero_1;
    if (is_presel_event) {
      N_presel_events++;
      N_presel_stubs += nTotStubs;
      N_presel_stubs_0 += nStubs_0;
      N_presel_stubs_1 += nStubs_1;

      h_nTotStubs_presel->Fill(nTotStubs);
      h_nStubs_0_presel->Fill(nStubs_0);
      h_nStubs_1_presel->Fill(nStubs_1);
      h2_nStubs_presel->Fill(nStubs_0, nStubs_1);
      

      for (int k=0; k<NMODULES; k++) {
	h_nstubsPerModule_presel->Fill(k, nstubs[k]);
      }
    }

    // S0 trackable and increasing requests in S1 
    if (is_trackable_0 && fired_modules_1 >= 1) {
      N_trackable_0_onehit_1++;
      N_trackable_0_onehit_1_stubs += nTotStubs;
      N_trackable_0_onehit_1_stubs_0 += nStubs_0;
      N_trackable_0_onehit_1_stubs_1 += nStubs_1;
      
      h_nTotStubs_S0_trackable_S1_onehit->Fill(nTotStubs);
      h_nStubs_0_S0_trackable_S1_onehit->Fill(nStubs_0);
      h_nStubs_1_S0_trackable_S1_onehit->Fill(nStubs_1);
      h2_nStubs_S0_trackable_S1_onehit->Fill(nStubs_0,nStubs_1);


      for (int k=0; k<NMODULES; k++) {
	h_nstubsPerModule_trackable_0_onehit_1->Fill(k, nstubs[k]);
      }
    }
    if (is_trackable_0 && fired_modules_1 >= 2) {
      N_trackable_0_twohit_1++;
    }
    if (is_trackable_0 && fired_modules_1 >= 3) {
      N_trackable_0_threehit_1++;
    }
    if (is_trackable_0 && fired_modules_1 >= 4) {
      N_trackable_0_fourhit_1++;
    }

    // S1 trackable and at least one hit in S0 
    if (is_trackable_1 && fired_modules_0 >= 1) {
      N_trackable_1_onehit_0++;
      N_trackable_1_onehit_0_stubs += nTotStubs;
      N_trackable_1_onehit_0_stubs_0 += nStubs_0;
      N_trackable_1_onehit_0_stubs_1 += nStubs_1;
      
      h_nTotStubs_S1_trackable_S0_onehit->Fill(nTotStubs);
      h_nStubs_0_S1_trackable_S0_onehit->Fill(nStubs_0);
      h_nStubs_1_S1_trackable_S0_onehit->Fill(nStubs_1);
      h2_nStubs_S1_trackable_S0_onehit->Fill(nStubs_0,nStubs_1);


      for (int k=0; k<NMODULES; k++) {
	h_nstubsPerModule_trackable_1_onehit_0->Fill(k, nstubs[k]);
      }
    }
      
    // at least a muon trackable in both stations
    bool is_trackable_event = is_trackable_0 && is_trackable_1;
    if (is_trackable_event) {
      N_trackable_events++;
      h_nTotStubs_trackable->Fill(nTotStubs);

      h_nStubs_0_trackable->Fill(nStubs_0);
      h_nStubs_1_trackable->Fill(nStubs_1);
      h2_nStubs_trackable->Fill(nStubs_0, nStubs_1);

      for (int k=0; k<NMODULES; k++) {
	h_nstubsPerModule_trackable->Fill(k, nstubs[k]);
        stub_per_mod_trackable.at(k)->Fill(nstubs.at(k));
      }
    }
    
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

    bool is_presel_not_trackable_event = is_presel_event && is_not_trackable_event;
    if (is_presel_not_trackable_event) {
      N_presel_not_trackable_events++;
      h_nTotStubs_presel_not_trackable->Fill(nTotStubs);
    }

    bool is_fired12_event = is_fired6_0 && is_fired6_1;
    if (is_fired12_event) {
      N_fired12_events++;
      h_nTotStubs_fired12->Fill(nTotStubs);
      h_nStubs_0_fired12->Fill(nStubs_0);
      h_nStubs_1_fired12->Fill(nStubs_1);
      h2_nStubs_fired12->Fill(nStubs_0, nStubs_1);

      for (int k=0; k<NMODULES; k++) {
	h_nstubsPerModule_fired12->Fill(k, nstubs[k]);
      }
    }
    
    bool is_passing_golden_event = is_passing_golden_0 && is_passing_golden_1;
    if (is_passing_golden_event) {
      N_passing_1mu_golden_events++;
    }

    bool is_presel_passing_golden_S0_event = is_presel_event && is_passing_golden_0;
    if (is_presel_passing_golden_S0_event) {
      N_presel_passing_golden_S0_events++;
    }
    
    bool is_presel_passing_golden_S1_event = is_presel_event && is_passing_golden_1;
    if (is_presel_passing_golden_S1_event) {
      N_presel_passing_golden_S1_events++;
    }

    //*************** efficiency 2 mu ***********************
    bool is_presel_passing_2mu_golden_S0_event = is_presel_event && is_passing_2mu_golden_0;
    if (is_presel_passing_2mu_golden_S0_event) {
      N_presel_passing_2mu_golden_S0_events++;
    }
    bool is_presel_passing_2mu_golden_S1_event = is_presel_event && is_passing_2mu_golden_1;
    if (is_presel_passing_2mu_golden_S1_event) {
      N_presel_passing_2mu_golden_S1_events++;
    }
    bool is_presel_passing_2mu_golden_S0_trackable_S1_event = is_presel_passing_2mu_golden_S0_event && is_passing_2mu_trackable_1;
    if (is_presel_passing_2mu_golden_S0_trackable_S1_event) {
      N_presel_passing_2mu_golden_S0_trackable_S1_events++;
    }
    bool is_presel_passing_2mu_golden_S1_trackable_S0_event = is_presel_passing_2mu_golden_S1_event && is_passing_2mu_trackable_0;
    if (is_presel_passing_2mu_golden_S1_trackable_S0_event) {
      N_presel_passing_2mu_golden_S1_trackable_S0_events++;
    }

    bool is_passing_2mu_golden_event = is_passing_2mu_golden_0 && is_passing_2mu_golden_1;
    if (is_passing_2mu_golden_event) {
      N_passing_2mu_golden_events++;
    }
    
    //*************** efficiency 3 mu ***********************
    bool is_presel_passing_3mu_golden_S0_event = is_presel_event && is_passing_3mu_golden_0;
    if (is_presel_passing_3mu_golden_S0_event) {
      N_presel_passing_3mu_golden_S0_events++;
    }
    bool is_presel_passing_3mu_golden_S1_event = is_presel_event && is_passing_3mu_golden_1;
    if (is_presel_passing_3mu_golden_S1_event) {
      N_presel_passing_3mu_golden_S1_events++;
    }
    bool is_presel_passing_3mu_golden_S0_trackable_S1_event = is_presel_passing_3mu_golden_S0_event && is_passing_3mu_trackable_1;
    if (is_presel_passing_3mu_golden_S0_trackable_S1_event) {
      N_presel_passing_3mu_golden_S0_trackable_S1_events++;
    }
    bool is_presel_passing_3mu_golden_S1_trackable_S0_event = is_presel_passing_3mu_golden_S1_event && is_passing_3mu_trackable_0;
    if (is_presel_passing_3mu_golden_S1_trackable_S0_event) {
      N_presel_passing_3mu_golden_S1_trackable_S0_events++;
    }

    bool is_passing_3mu_golden_event = is_passing_3mu_golden_0 && is_passing_3mu_golden_1;
    if (is_passing_3mu_golden_event) {
      N_passing_3mu_golden_events++;
    }
    
    //*************** efficiency 4 mu ***********************
    bool is_presel_passing_4mu_golden_S0_event = is_presel_event && is_passing_4mu_golden_0;
    if (is_presel_passing_4mu_golden_S0_event) {
      N_presel_passing_4mu_golden_S0_events++;
    }
    bool is_presel_passing_4mu_golden_S1_event = is_presel_event && is_passing_4mu_golden_1;
    if (is_presel_passing_4mu_golden_S1_event) {
      N_presel_passing_4mu_golden_S1_events++;
    }
    bool is_presel_passing_4mu_golden_S0_trackable_S1_event = is_presel_passing_4mu_golden_S0_event && is_passing_4mu_trackable_1;
    if (is_presel_passing_4mu_golden_S0_trackable_S1_event) {
      N_presel_passing_4mu_golden_S0_trackable_S1_events++;
    }
    bool is_presel_passing_4mu_golden_S1_trackable_S0_event = is_presel_passing_4mu_golden_S1_event && is_passing_4mu_trackable_0;
    if (is_presel_passing_4mu_golden_S1_trackable_S0_event) {
      N_presel_passing_4mu_golden_S1_trackable_S0_events++;
    }

    bool is_passing_4mu_golden_event = is_passing_4mu_golden_0 && is_passing_4mu_golden_1;
    if (is_passing_4mu_golden_event) {
      N_passing_4mu_golden_events++;
    }
    //***************
    
    bool is_passing_golden_S0_zero_S1_event = is_passing_golden_0 && is_zero_1;
    if (is_passing_golden_S0_zero_S1_event) N_passing_golden_S0_zero_S1_events++;
    bool is_passing_golden_S1_zero_S0_event = is_passing_golden_1 && is_zero_0;
    if (is_passing_golden_S1_zero_S0_event) N_passing_golden_S1_zero_S0_events++;
  
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
    if (is_passing_4mu_event) {
      N_passing_4mu_events++;
      h_nTotStubs_passing_4mu->Fill(nTotStubs);
    }
    
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

      N_single_clean_stubs += nTotStubs;
      N_single_clean_stubs_0 += nStubs_0;
      N_single_clean_stubs_1 += nStubs_1;


      h_nStubs_0_single_clean->Fill(nStubs_0);
      h_nStubs_1_single_clean->Fill(nStubs_1);
      h2_nStubs_single_clean->Fill(nStubs_0, nStubs_1);
      

      for (int k=0; k<NMODULES; k++) {
	h_nstubsPerModule_single_clean->Fill(k, nstubs[k]);
	stub_per_mod_single_clean.at(k)->Fill(nstubs.at(k));
      }
    }

    bool is_fired12_plus_single_clean_event = is_fired6_plus_single_clean_0 && is_fired6_2tracks_1;
    if (is_fired12_plus_single_clean_event) {
      N_fired12_plus_single_clean_events++;
      h_nTotStubs_fired12_plus_single_clean->Fill(nTotStubs);

      h_nStubs_0_fired12_plus_single_clean->Fill(nStubs_0);
      h_nStubs_1_fired12_plus_single_clean->Fill(nStubs_1);
      h2_nStubs_fired12_plus_single_clean->Fill(nStubs_0, nStubs_1);

      for (int k=0; k<NMODULES; k++) {
	h_nstubsPerModule_fired12_plus_single_clean->Fill(k, nstubs[k]);
      }
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

      N_pileup_skim_stubs += nTotStubs;
      N_pileup_skim_stubs_0 += nStubs_0;
      N_pileup_skim_stubs_1 += nStubs_1;

      h_nStubs_0_pileup_skim->Fill(nStubs_0);
      h_nStubs_1_pileup_skim->Fill(nStubs_1);
      h2_nStubs_pileup_skim->Fill(nStubs_0, nStubs_1);
      
    }

    bool is_pileup234_skim_event = is_pileup_2mu_event || is_pileup_3mu_event || is_pileup_4mu_event;
    if (is_pileup234_skim_event) {
      N_pileup234_skim_events++;
      h_nTotStubs_pileup234_skim->Fill(nTotStubs);     
      if (iskim == 60) ok_skim = true;

      N_pileup234_skim_stubs += nTotStubs;
      N_pileup234_skim_stubs_0 += nStubs_0;
      N_pileup234_skim_stubs_1 += nStubs_1;

      h_nStubs_0_pileup234_skim->Fill(nStubs_0);
      h_nStubs_1_pileup234_skim->Fill(nStubs_1);
      h2_nStubs_pileup234_skim->Fill(nStubs_0, nStubs_1);



      for (int k=0; k<NMODULES; k++) {
	h_nstubsPerModule_pileup234_skim->Fill(k, nstubs[k]);
	stub_per_mod_pileup234.at(k)->Fill(nstubs.at(k));
      }
    }

    bool is_loose_event = is_single_clean_event || is_pileup_2mu_event || is_pileup_3mu_event || is_pileup_many_event;
    if (is_loose_event) {
      N_loose_events++;
      if (iskim == 41) ok_skim = true;
    }

    bool is_loose234_event = is_single_clean_event || is_pileup_2mu_event || is_pileup_3mu_event || is_pileup_4mu_event;
    if (is_loose234_event) {
      N_loose234_events++;
      if (iskim == 61) ok_skim = true;
    }

    if (ok_skim) {
      N_stubs_skim += nTotStubs;
    }
   }//mu_in
  } //  for(Long64_t i = 0; i < n_entries; i++)



  cout <<"\n"<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
  cout<<"\n Station 0:"<<endl;
  cout<<"Number of events without any hit in S_0 = "<<endl;
  cout<<"\t" << N_zero_0 <<endl; 
  cout<<"Number of Empty or Noise events in S_0 = "<<endl;
  cout<<"\t" << N_noise_0 <<endl;
  cout<<"Number of events with >=1 fired S_0 module = "<<endl;
  cout<<"\t" << N_onehit_0 <<endl;
  cout<<"Number of events with >=2 fired S_0 modules = "<<endl;
  cout<<"\t" << N_twohit_0 <<endl;
  cout<<"Number of events with >=3 fired S_0 modules = "<<endl;
  cout<<"\t" << N_threehit_0 <<endl;
  cout<<"Number of events with >=4 fired S_0 modules = "<<endl;
  cout<<"\t" << N_fourhit_0 <<endl;
  cout<<"Number of events with all 6 fired S_0 modules = "<<endl;
  cout<<"\t" << N_fired6_0 <<endl;
  cout<<"Number of trackable events (all 4 XY modules and at least 1 UV) = "<<endl;
  cout<<"\t "<<N_trackable_0<<endl;
  cout<<"Number of candidate single mu events in S0 = "<<endl;
  cout<<"\t "<<N_single_cand_0 <<endl;
  cout<<"Number of clean single mu events in S_0 = "<<endl;
  cout<<"\t "<<N_single_clean_0 <<endl;
  cout<<"Number of golden 1mu passing events in Station 0 (all 6 modules with one hit) = " <<endl;
  cout<<"\t" << N_passing_1mu_golden_0 <<endl;
  cout<<"Number of golden 2mu passing events in Station 0 (all 6 modules with two hits) = " <<endl;
  cout<<"\t" << N_passing_2mu_golden_0 <<endl;
  cout<<"Number of golden 3mu passing events in Station 0 (all 6 modules with three hits) = " <<endl;
  cout<<"\t" << N_passing_3mu_golden_0 <<endl;
  cout<<"Number of golden 4mu passing events in Station 0 (all 6 modules with four hits) = " <<endl;
  cout<<"\t" << N_passing_4mu_golden_0 <<endl;
  cout<<"Number of (any) pileup events in S_0 = "<<endl;
  cout<<"\t "<<N_pileup_any_0 <<endl;
  cout<<"Number of 2mu pileup events in S_0 = "<<endl;
  cout<<"\t "<<N_pileup_2mu_0 <<endl;  
  cout<<"Number of 3mu pileup events in S_0 = "<<endl;
  cout<<"\t "<<N_pileup_3mu_0 <<endl;  
  cout<<"Number of 4mu pileup events in S_0 = "<<endl;
  cout<<"\t "<<N_pileup_4mu_0 <<endl;  
  cout<<"Number of >=4mu pileup events in S_0 = "<<endl;
  cout<<"\t "<<N_pileup_many_0 <<endl;  
 
  cout<<"\n Station 1:"<<endl;
  cout<<"Number of events without any hit in S_1 = "<<endl;
  cout<<"\t" << N_zero_1 <<endl; 
  cout<<"Number of Empty or Noise events in S_1 = "<<endl;
  cout<<"\t" << N_noise_1 <<endl; 
  cout<<"Number of events with all 6 fired S_1 modules = "<<endl;
  cout<<"\t" << N_fired6_1 <<endl;
  cout<<"Number of golden 1mu passing events in Station 1 (all 6 modules with one hit) = " <<endl;
  cout<<"\t" << N_passing_1mu_golden_1 <<endl;
  cout<<"Number of golden 2mu passing events in Station 1 (all 6 modules with two hits) = " <<endl;
  cout<<"\t" << N_passing_2mu_golden_1 <<endl;
  cout<<"Number of golden 3mu passing events in Station 1 (all 6 modules with three hits) = " <<endl;
  cout<<"\t" << N_passing_3mu_golden_1 <<endl;
  cout<<"Number of golden 4mu passing events in Station 1 (all 6 modules with four hits) = " <<endl;
  cout<<"\t" << N_passing_4mu_golden_1 <<endl;
  cout<<"Number of candidate single mu events in S_1 = "<<endl;
  cout<<"\t "<<N_passing_cand_1 <<endl;
  cout<<"Number of clean single mu events in S_1 = "<<endl;
  cout<<"\t "<<N_passing_clean_1 <<endl;
  cout<<"Number of interaction candidates (2 tracks) in S_1 = "<<endl;
  cout<<"\t "<< N_2tracks_1<<endl;
  cout<<"Number of golden events in S_1 (2 hits in all modules in S_1)= " <<endl;
  cout<<"\t" << N_golden_1 <<endl;
  cout<<"Number of signal+pileup candidates in S_1 = "<<endl;
  cout<<"\t "<<N_signal_plus_pileup_1<<endl;
  
  cout<<"\n EVENTS:   "<<endl;
  cout<<"Number all input events = " << endl;
  cout<<"\t" << N_input_events <<endl;
  cout<<"Number of events without any hit in both S_0 and S_1 = "<<endl;
  cout<<"\t" << N_zero_events <<endl; 
  cout<<"Number of Noise events (empty or <=1 module in both stations) = " << endl;
  cout<<"\t "<< N_noise_events <<endl;
  cout<<"Number of events missing S_0 (empty or <=1 module) and trackable in S_1 = " << endl;
  cout<<"\t "<< N_missing_0_events <<endl;
  cout<<"Number of events missing S_1 (empty or <=1 module) and trackable in S_0 = " << endl;
  cout<<"\t "<< N_missing_1_events <<endl;
  cout<<"Number of Preselected events = "<<endl;
  cout<<"\t "<< N_presel_events <<endl;
  //
  cout<<"Number of trackable S_0 events with >=1 fired S_1 module = "<<endl;
  cout<<"\t" << N_trackable_0_onehit_1 <<endl;
  cout<<"Number of trackable S_0 events with >=2 fired S_1 module = "<<endl;
  cout<<"\t" << N_trackable_0_twohit_1 <<endl;
  cout<<"Number of trackable S_0 events with >=3 fired S_1 module = "<<endl;
  cout<<"\t" << N_trackable_0_threehit_1 <<endl;
  cout<<"Number of trackable S_0 events with >=4 fired S_1 module = "<<endl;
  cout<<"\t" << N_trackable_0_fourhit_1 <<endl;
  cout<<"Number of trackable events = "<<endl; // both stations trackable
  cout<<"\t "<< N_trackable_events <<endl;

  cout<<"Number of events with ALL 12 fired modules = "<<endl;
  cout<<"\t" << N_fired12_events <<endl;
  
  cout<<"Number of clean 1mu passing events = " <<endl;
  cout<<"\t "<< N_passing_1mu_events <<endl;

  cout<<"Number of passing 1mu golden S_0 in preselected events = "<< endl;
  cout<<"\t "<< N_presel_passing_golden_S0_events <<endl;
  cout<<"Number of passing 1mu golden S_0 and trackable in S_1 = " << endl;
  cout<<"\t "<< N_passing_golden_S0_trackable_S1_events <<endl;
  cout<<"\tNumber of passing 1mu golden S_0 and zero hits in S_1 = " << endl;
  cout<<"\t\t "<< N_passing_golden_S0_zero_S1_events <<endl;
  //
  cout<<"Number of passing 1mu golden S_1 in preselected events = "<< endl;
  cout<<"\t "<< N_presel_passing_golden_S1_events <<endl;
  cout<<"Number of passing 1mu golden S_1 and trackable in S_0 = " << endl;
  cout<<"\t "<< N_passing_golden_S1_trackable_S0_events <<endl;
  cout<<"\tNumber of passing 1mu golden S_1 and zero hits in S_0 = " << endl;
  cout<<"\t\t "<< N_passing_golden_S1_zero_S0_events <<endl;

  cout<<"Number of golden 1mu passing events (all 12 modules with one hit) = " <<endl;
  cout<<"\t" << N_passing_1mu_golden_events <<endl;
  cout<<"Number of golden 2mu passing events (all 12 modules with two hits) = " <<endl;
  cout<<"\t" << N_passing_2mu_golden_events <<endl;
  cout<<"Number of golden 3mu passing events (all 12 modules with three hits) = " <<endl;
  cout<<"\t" << N_passing_3mu_golden_events <<endl;
  cout<<"Number of golden 4mu passing events (all 12 modules with four hits) = " <<endl;
  cout<<"\t" << N_passing_4mu_golden_events <<endl;
  
  cout<<"Number of candidate Single mu interaction events = " << endl;
  cout<<"\t "<< N_single_cand_events <<endl;
  cout<<"Number of clean Single mu interaction events = " << endl;
  cout<<"\t "<< N_single_clean_events <<endl;
  cout<<"Number of golden Single mu interaction events = " << endl;
  cout<<"\t "<< N_golden_events <<endl;
  cout<<"Number of any pileup interaction events = " << endl;
  cout<<"\t "<< N_pileup_any_events <<endl;
  cout<<"Number of 2mu pileup interaction events = " << endl;
  cout<<"\t "<< N_pileup_2mu_events <<endl;
  cout<<"Number of 3mu pileup interaction events = " << endl;
  cout<<"\t "<< N_pileup_3mu_events <<endl;
  cout<<"Number of 4mu pileup interaction events = " << endl;
  cout<<"\t "<< N_pileup_4mu_events <<endl;
  cout<<"Number of >=4mu pileup interaction events = " << endl;
  cout<<"\t "<< N_pileup_many_events <<endl;
  cout<<"Total number of pileup (2,3,4) interaction events = " << endl;
  cout<<"\t "<< N_pileup234_skim_events <<endl;
  cout<<"Total number of pileup (2,3,>=4) interaction events = " << endl;
  cout<<"\t "<< N_pileup_skim_events <<endl;
  cout<<"Total number of Single + Pileup (2,3,4) interaction events = " << endl;
  cout<<"\t "<< N_loose234_events <<endl;
  cout<<"Total number of Single + Pileup (2,3,>=4) interaction events = " << endl;
  cout<<"\t "<< N_loose_events <<endl;
  cout<<"Number of events with 1mu (12stubs) overlapping Single Clean events = "<<endl;
  cout<<"\t" << N_fired12_plus_single_clean_events <<endl;
  cout <<      "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;

  cout<<"N INPUT Events    = "<< N_input_events <<endl;
  cout<<"N input stubs S_0 = "<< N_input_stubs_0 <<endl;  
  cout<<"N input stubs S_1 = "<< N_input_stubs_1 <<endl;  
  cout<<"N input stubs tot = "<< N_input_stubs <<endl;  
  cout<<endl;
  cout<<"Fraction missing S_0             = "<< double(N_missing_0_events)/double(N_input_events) <<endl;
  cout<<"Fraction missing S_1             = "<< double(N_missing_1_events)/double(N_input_events) <<endl;
  cout<<"Fraction not trackable           = "<< double(N_not_trackable_events)/double(N_input_events) <<endl;

  cout<<endl<<"N PRESELECTED Events     = "<< N_presel_events <<endl;
  cout<<"N stubs S_0 (presel.ev.) = "<< N_presel_stubs_0 <<endl;  
  cout<<"N stubs S_1 (presel.ev.) = "<< N_presel_stubs_1 <<endl;  
  cout<<"N stubs tot (presel.ev.) = "<< N_presel_stubs <<endl;  

  cout<<endl<<"N PRESELECTED S_0 Trackable Events     = "<< N_trackable_0_onehit_1 <<endl;
  cout<<"N stubs S_0 (presel.S0-trackable ev.) = "<< N_trackable_0_onehit_1_stubs_0 <<endl;  
  cout<<"N stubs S_1 (presel.S0-trackable ev.) = "<< N_trackable_0_onehit_1_stubs_1 <<endl;  
  cout<<"N stubs tot (presel.S0-trackable ev.) = "<< N_trackable_0_onehit_1_stubs <<endl;  

  cout<<endl<<"N TRACKABLE Events       = "<< N_trackable_events <<endl;
  
  cout<<"Fraction PRESELECTED / all       = "<< double(N_presel_events)/double(N_input_events) <<endl;
  cout<<"Fraction PRESEL. S_0-Trackable / PRESEL = "<< double(N_trackable_0_onehit_1)/double(N_presel_events) <<endl;
  cout<<"Fraction TRACKABLE / PRESEL      = "<< double(N_trackable_events)/double(N_presel_events) <<endl;

  ////////////////////////////////////////////////////

  cout<<endl;
  cout<<"***** RATES w.r.t. all MERGED events ***************************************************" << endl;
  cout<<"Fraction Passing 1mu             = "<< double(N_passing_1mu_events)/double(N_input_events) <<endl;
  cout<<"Fraction Passing GOLDEN 1mu      = "<< double(N_passing_1mu_golden_events)/double(N_input_events) <<endl;
  cout<<"Fraction Passing GOLDEN 2mu      = "<< double(N_passing_2mu_golden_events)/double(N_input_events) <<endl;
  cout<<"Fraction Passing GOLDEN 3mu      = "<< double(N_passing_3mu_golden_events)/double(N_input_events) <<endl;
  cout<<"Fraction Passing GOLDEN 4mu      = "<< double(N_passing_4mu_golden_events)/double(N_input_events) <<endl;
  cout<<"Fraction Single Mu Int.          = "<< double(N_single_clean_events)/double(N_input_events) <<endl;
  cout<<"Fraction Golden sel              = "<< double(N_golden_events)/double(N_input_events) <<endl;
  cout<<"Fraction 2 Mu PU Int.            = "<< double(N_pileup_2mu_events)/double(N_input_events) <<endl;
  cout<<"Fraction 3 Mu PU Int.            = "<< double(N_pileup_3mu_events)/double(N_input_events) <<endl;
  cout<<"Fraction 4 Mu PU Int.            = "<< double(N_pileup_4mu_events)/double(N_input_events) <<endl;
  cout<<"Fraction >=4 Mu PU Int.          = "<< double(N_pileup_many_events)/double(N_input_events) <<endl;
  cout<<"Total Fraction PU 2+3+4 Mu Int.  = "<< double(N_pileup234_skim_events)/double(N_input_events) <<endl;
  cout<<"Total Fraction PU Mu Int.        = "<< double(N_pileup_skim_events)/double(N_input_events) <<endl;
  cout<<"Total Fraction Single+PU 2+3+4   = "<< double(N_loose234_events)/double(N_input_events) <<endl;
  cout<<"Total Fraction Single+PU 2+3+4+  = "<< double(N_loose_events)/double(N_input_events) <<endl;
  cout<<endl;
  cout<<"===========================================================================================" <<endl;
  cout<<"Fraction of PRESELECTED events (with at least one stub in S0 and one stub in S1): "
      << double(N_presel_events)/double(N_input_events) << endl;
  cout<<"===========================================================================================" <<endl;
  cout<<endl;
  cout<<"***** RATES w.r.t. PRESELECTED events ***************************************************" << endl;
  cout<<"Fraction Passing 1mu             = "<< double(N_passing_1mu_events)/double(N_presel_events) <<endl;
  cout<<"Fraction Passing GOLDEN 1mu      = "<< double(N_passing_1mu_golden_events)/double(N_presel_events) <<endl;
  cout<<"Fraction Passing GOLDEN 2mu      = "<< double(N_passing_2mu_golden_events)/double(N_presel_events) <<endl;
  cout<<"Fraction Passing GOLDEN 3mu      = "<< double(N_passing_3mu_golden_events)/double(N_presel_events) <<endl;
  cout<<"Fraction Passing GOLDEN 4mu      = "<< double(N_passing_4mu_golden_events)/double(N_presel_events) <<endl;
  cout<<"Fraction Single Mu Int.          = "<< double(N_single_clean_events)/double(N_presel_events) <<endl;
  cout<<"Fraction Golden sel              = "<< double(N_golden_events)/double(N_presel_events) <<endl;
  cout<<"Fraction 2 Mu PU Int.            = "<< double(N_pileup_2mu_events)/double(N_presel_events) <<endl;
  cout<<"Fraction 3 Mu PU Int.            = "<< double(N_pileup_3mu_events)/double(N_presel_events) <<endl;
  cout<<"Fraction 4 Mu PU Int.            = "<< double(N_pileup_4mu_events)/double(N_presel_events) <<endl;
  cout<<"Fraction >=4 Mu PU Int.          = "<< double(N_pileup_many_events)/double(N_presel_events) <<endl;
  cout<<"Total Fraction PU 2+3+4 Mu Int.  = "<< double(N_pileup234_skim_events)/double(N_presel_events) <<endl;
  cout<<"Total Fraction PU Mu Int.        = "<< double(N_pileup_skim_events)/double(N_presel_events) <<endl;
  cout<<"Total Fraction Single+PU 2+3+4   = "<< double(N_loose234_events)/double(N_presel_events) <<endl;
  cout<<"Total Fraction Single+PU 2+3+4+  = "<< double(N_loose_events)/double(N_presel_events) <<endl;
  cout<<endl;
  cout<<"===========================================================================================" <<endl;
  cout<<"Fraction of S0-Trackable events (PRESELECTED, with at least one stub in S1) w.r.t. all MERGED: "
      << double(N_trackable_0_onehit_1)/double(N_input_events) << endl;
  cout<<"Fraction of S0-Trackable events (PRESELECTED, with at least one stub in S1) w.r.t. PRESELECTED: "
      << double(N_trackable_0_onehit_1)/double(N_presel_events) << endl;
  cout<<"===========================================================================================" <<endl;
  cout<<endl;
  cout<<"***** RATES w.r.t. S0-Trackable PRESELECTED events **************************************" << endl;
  cout<<"Fraction Passing 1mu             = "<< double(N_passing_1mu_events)/double(N_trackable_0_onehit_1) <<endl;
  cout<<"Fraction Passing GOLDEN 1mu      = "<< double(N_passing_1mu_golden_events)/double(N_trackable_0_onehit_1) <<endl;
  cout<<"Fraction Passing GOLDEN 2mu      = "<< double(N_passing_2mu_golden_events)/double(N_trackable_0_onehit_1) <<endl;
  cout<<"Fraction Passing GOLDEN 3mu      = "<< double(N_passing_3mu_golden_events)/double(N_trackable_0_onehit_1) <<endl;
  cout<<"Fraction Passing GOLDEN 4mu      = "<< double(N_passing_4mu_golden_events)/double(N_trackable_0_onehit_1) <<endl;
  cout<<"Fraction Single Mu Int.          = "<< double(N_single_clean_events)/double(N_trackable_0_onehit_1) <<endl;
  cout<<"Fraction Golden sel              = "<< double(N_golden_events)/double(N_trackable_0_onehit_1) <<endl;
  cout<<"Fraction 2 Mu PU Int.            = "<< double(N_pileup_2mu_events)/double(N_trackable_0_onehit_1) <<endl;
  cout<<"Fraction 3 Mu PU Int.            = "<< double(N_pileup_3mu_events)/double(N_trackable_0_onehit_1) <<endl;
  cout<<"Fraction 4 Mu PU Int.            = "<< double(N_pileup_4mu_events)/double(N_trackable_0_onehit_1) <<endl;
  cout<<"Fraction >=4 Mu PU Int.          = "<< double(N_pileup_many_events)/double(N_trackable_0_onehit_1) <<endl;
  cout<<"Total Fraction PU 2+3+4 Mu Int.  = "<< double(N_pileup234_skim_events)/double(N_trackable_0_onehit_1) <<endl;
  cout<<"Total Fraction PU Mu Int.        = "<< double(N_pileup_skim_events)/double(N_trackable_0_onehit_1) <<endl;
  cout<<"Total Fraction Single+PU 2+3+4   = "<< double(N_loose234_events)/double(N_trackable_0_onehit_1) <<endl;
  cout<<"Total Fraction Single+PU 2+3+4+  = "<< double(N_loose_events)/double(N_trackable_0_onehit_1) <<endl;
  cout<<endl;
  cout<<"===========================================================================================" <<endl;
  cout<<"Fraction of TRACKABLE events w.r.t. all Merged events: "
      << double(N_trackable_events)/double(N_input_events) << endl;
  cout<<"Fraction of TRACKABLE events w.r.t. PRESELECTED events: "
      << double(N_trackable_events)/double(N_presel_events) << endl;
  cout<<"Fraction of TRACKABLE events w.r.t. PRESELECTED && S0-Trackable events: "
      << double(N_trackable_events)/double(N_trackable_0_onehit_1) << endl;
  cout<<"===========================================================================================" <<endl;
  cout<<endl;
  cout<<"***** RATES w.r.t. TRACKABLE events ***************************************************" << endl;
  cout<<"Fraction Passing 1mu             = "<< double(N_passing_1mu_events)/double(N_trackable_events) <<endl;
  cout<<"Fraction Passing GOLDEN 1mu      = "<< double(N_passing_1mu_golden_events)/double(N_trackable_events) <<endl;
  cout<<"Fraction Passing GOLDEN 2mu      = "<< double(N_passing_2mu_golden_events)/double(N_trackable_events) <<endl;
  cout<<"Fraction Passing GOLDEN 3mu      = "<< double(N_passing_3mu_golden_events)/double(N_trackable_events) <<endl;
  cout<<"Fraction Passing GOLDEN 4mu      = "<< double(N_passing_4mu_golden_events)/double(N_trackable_events) <<endl;
  cout<<"Fraction Single Mu Int.          = "<< double(N_single_clean_events)/double(N_trackable_events) <<endl;
  cout<<"Fraction Golden sel              = "<< double(N_golden_events)/double(N_trackable_events) <<endl;
  cout<<"Fraction 2 Mu PU Int.            = "<< double(N_pileup_2mu_events)/double(N_trackable_events) <<endl;
  cout<<"Fraction 3 Mu PU Int.            = "<< double(N_pileup_3mu_events)/double(N_trackable_events) <<endl;
  cout<<"Fraction 4 Mu PU Int.            = "<< double(N_pileup_4mu_events)/double(N_trackable_events) <<endl;
  cout<<"Fraction >=4 Mu PU Int.          = "<< double(N_pileup_many_events)/double(N_trackable_events) <<endl;
  cout<<"Total Fraction PU 2+3+4 Mu Int.  = "<< double(N_pileup234_skim_events)/double(N_trackable_events) <<endl;
  cout<<"Total Fraction PU Mu Int.        = "<< double(N_pileup_skim_events)/double(N_trackable_events) <<endl;
  cout<<"Total Fraction Single+PU 2+3+4   = "<< double(N_loose234_events)/double(N_trackable_events) <<endl;
  cout<<"Total Fraction Single+PU 2+3+4+  = "<< double(N_loose_events)/double(N_trackable_events) <<endl;

  cout<<"===========================================================================================" <<endl;
  cout<<"Fraction of 1mu PASSING GOLDEN events w.r.t. all Merged events: "
      << double(N_passing_1mu_golden_events)/double(N_input_events) << endl;
  cout<<"Fraction of 1mu PASSING GOLDEN events w.r.t. PRESELECTED events: "
      << double(N_passing_1mu_golden_events)/double(N_presel_events) << endl;
  cout<<"Fraction of 1mu PASSING GOLDEN events w.r.t. PRESELECTED && S0-Trackable events: "
      << double(N_passing_1mu_golden_events)/double(N_trackable_0_onehit_1) << endl;
  cout<<"Fraction of 1mu PASSING GOLDEN events w.r.t. TRACKABLE events: "
      << double(N_passing_1mu_golden_events)/double(N_trackable_events) << endl;

  cout<<"===========================================================================================" <<endl;
  cout<<endl;
  cout<<"***** RATES w.r.t. 1mu PASSING GOLDEN events **********************************************" << endl;
  cout<<"Fraction Passing 1mu             = "<< double(N_passing_1mu_events)/double(N_passing_1mu_golden_events) <<endl;
  cout<<"Fraction Passing GOLDEN 1mu      = "<< double(N_passing_1mu_golden_events)/double(N_passing_1mu_golden_events) <<endl;
  cout<<"Fraction Passing GOLDEN 2mu      = "<< double(N_passing_2mu_golden_events)/double(N_passing_1mu_golden_events) <<endl;
  cout<<"Fraction Passing GOLDEN 3mu      = "<< double(N_passing_3mu_golden_events)/double(N_passing_1mu_golden_events) <<endl;
  cout<<"Fraction Passing GOLDEN 4mu      = "<< double(N_passing_4mu_golden_events)/double(N_passing_1mu_golden_events) <<endl;
  cout<<"Fraction Single Mu Int.          = "<< double(N_single_clean_events)/double(N_passing_1mu_golden_events) <<endl;
  cout<<"Fraction Golden sel              = "<< double(N_golden_events)/double(N_passing_1mu_golden_events) <<endl;
  cout<<"Fraction 2 Mu PU Int.            = "<< double(N_pileup_2mu_events)/double(N_passing_1mu_golden_events) <<endl;
  cout<<"Fraction 3 Mu PU Int.            = "<< double(N_pileup_3mu_events)/double(N_passing_1mu_golden_events) <<endl;
  cout<<"Fraction 4 Mu PU Int.            = "<< double(N_pileup_4mu_events)/double(N_passing_1mu_golden_events) <<endl;
  cout<<"Fraction >=4 Mu PU Int.          = "<< double(N_pileup_many_events)/double(N_passing_1mu_golden_events) <<endl;
  cout<<"Total Fraction PU 2+3+4 Mu Int.  = "<< double(N_pileup234_skim_events)/double(N_passing_1mu_golden_events) <<endl;
  cout<<"Total Fraction PU Mu Int.        = "<< double(N_pileup_skim_events)/double(N_passing_1mu_golden_events) <<endl;
  cout<<"Total Fraction Single+PU 2+3+4   = "<< double(N_loose234_events)/double(N_passing_1mu_golden_events) <<endl;
  cout<<"Total Fraction Single+PU 2+3+4+  = "<< double(N_loose_events)/double(N_passing_1mu_golden_events) <<endl;

  cout<<"===========================================================================================" <<endl;
  cout<<"Fraction of events with ALL 12 modules fired w.r.t. all Merged events: "
      << double(N_fired12_events)/double(N_input_events) << endl;
  cout<<"Fraction of events with ALL 12 modules fired w.r.t. PRESELECTED events: "
      << double(N_fired12_events)/double(N_presel_events) << endl;
  cout<<"Fraction of events with ALL 12 modules fired w.r.t. PRESELECTED && S0-Trackable events: "
      << double(N_fired12_events)/double(N_trackable_0_onehit_1) << endl;
  cout<<"Fraction of events with ALL 12 modules fired w.r.t. TRACKABLE events: "
      << double(N_fired12_events)/double(N_trackable_events) << endl;

  cout<<"===========================================================================================" <<endl;
  cout<<endl;
  cout<<"***** RATES w.r.t. events with ALL 12 modules fired ***************************************" << endl;
  cout<<"Fraction Passing 1mu             = "<< double(N_passing_1mu_events)/double(N_fired12_events) <<endl;
  cout<<"Fraction Passing GOLDEN 1mu      = "<< double(N_passing_1mu_golden_events)/double(N_fired12_events) <<endl;
  cout<<"Fraction Passing GOLDEN 2mu      = "<< double(N_passing_2mu_golden_events)/double(N_fired12_events) <<endl;
  cout<<"Fraction Passing GOLDEN 3mu      = "<< double(N_passing_3mu_golden_events)/double(N_fired12_events) <<endl;
  cout<<"Fraction Passing GOLDEN 4mu      = "<< double(N_passing_4mu_golden_events)/double(N_fired12_events) <<endl;
  cout<<"Fraction Single Mu Int.          = "<< double(N_single_clean_events)/double(N_fired12_events) <<endl;
  cout<<"Fraction Golden sel              = "<< double(N_golden_events)/double(N_fired12_events) <<endl;
  cout<<"Fraction 2 Mu PU Int.            = "<< double(N_pileup_2mu_events)/double(N_fired12_events) <<endl;
  cout<<"Fraction 3 Mu PU Int.            = "<< double(N_pileup_3mu_events)/double(N_fired12_events) <<endl;
  cout<<"Fraction 4 Mu PU Int.            = "<< double(N_pileup_4mu_events)/double(N_fired12_events) <<endl;
  cout<<"Fraction >=4 Mu PU Int.          = "<< double(N_pileup_many_events)/double(N_fired12_events) <<endl;
  cout<<"Total Fraction PU 2+3+4 Mu Int.  = "<< double(N_pileup234_skim_events)/double(N_fired12_events) <<endl;
  cout<<"Total Fraction PU Mu Int.        = "<< double(N_pileup_skim_events)/double(N_fired12_events) <<endl;
  cout<<"Total Fraction Single+PU 2+3+4   = "<< double(N_loose234_events)/double(N_fired12_events) <<endl;
  cout<<"Total Fraction Single+PU 2+3+4+  = "<< double(N_loose_events)/double(N_fired12_events) <<endl;
  cout<<"Fraction Single Mu Int.+FIRED12  = "<< double(N_fired12_plus_single_clean_events)/double(N_fired12_events) <<endl;

  // reset stat errors on the per module histos
  h_nstubsPerModule->Sumw2(false);
  h_nstubsPerModule_presel->Sumw2(false);
  h_nstubsPerModule_trackable_0_onehit_1->Sumw2(false);
  h_nstubsPerModule_trackable_1_onehit_0->Sumw2(false);
  h_nstubsPerModule_trackable->Sumw2(false);
  h_nstubsPerModule_fired12->Sumw2(false);
  h_nstubsPerModule_single_clean->Sumw2(false);
  h_nstubsPerModule_pileup234_skim->Sumw2(false);
  h_nstubsPerModule_fired12_plus_single_clean->Sumw2(false);
  
  hfile->Write();

  if (iskim > 0) {
    cout<<"\n Number of output skimmed events          = "<<o_nevt<<endl;
    cout<<" Number of stubs in output skimmed events = "<<N_stubs_skim<<endl;
  }

  return 0;
}













