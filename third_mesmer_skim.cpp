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


void third_mesmer_skim(){

int NMODULES = 12;
int iskim = 41;

///////
  string filelist = "runs/";
  string inputdir = "inputs/";
  string n_mult="third_1PUmu";
  string outdir = "outputs_mesmer/";

    filelist = filelist  + "files.txt";


TChain * cbmsim = new TChain("cbmsim");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/theta_0-2mrad_1M_1hit_NOoutchi2.root");
// /mnt/raid10/DATA/espedica/fairmu/theta_0-32mrad_1M_1hit_NOoutchi2.root");
// /mnt/raid10/DATA/espedica/fairmu/theta_0-2mrad_1M_1hit_NOoutchi2.root");

 auto pos_on_track = [](double q, double m, double z){return (q + m*z);};

  //////////////////////////////////////////////////////////////////////
  auto n_entries = cbmsim->GetEntries();
  cout << "Total number of input events  = " << n_entries << endl;

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

  //////////////////////////////////////////////////////////////////////
   TClonesArray *o_MCTrack = 0;
   TClonesArray *o_TrackerStripDigis = 0;
   TClonesArray *o_TrackerPoints = 0;
   TClonesArray *o_TrackerStubs = 0;
   MuE::Event *o_MesmerEvent = 0;
   MUonERecoOutput *o_ReconstructionOutput = 0;

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
/*    o_file = TFile::Open((outdir+n_mult+".root").c_str(), "recreate");
    if (!o_file) {cout << "Error opening file" << std::endl; exit(-1);}
    o_tree = new TTree("cbmsim", "");

    o_tree->Branch("MCTrack", &o_MCTrack);
    o_tree->SetBranchAddress("TrackerPoints", &o_TrackerPoints);
    o_tree->SetBranchAddress("TrackerStripDigis", &o_TrackerStripDigis);
    o_tree->SetBranchAddress("MesmerEvent", &o_MesmerEvent);
    o_tree->Branch("TrackerStubs", &o_TrackerStubs);
    o_tree->Branch("ReconstructionOutput", &o_ReconstructionOutput);*/
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
  
  //  TH1D *h_nTotStubs_passing_golden = new TH1D("h_nTotStubs_passing_golden","Total N of Stubs for best passing mu",nb,0,nmax);
  TH1D *h_nTotStubs_passing_1mu = new TH1D("h_nTotStubs_passing_1mu","Total N of Stubs for passing mu",nb,0,nmax);
  TH1D *h_nTotStubs_passing_2mu = new TH1D("h_nTotStubs_passing_2mu","Total N of Stubs for passing 2mu",nb,0,nmax);
  TH1D *h_nTotStubs_passing_3mu = new TH1D("h_nTotStubs_passing_3mu","Total N of Stubs for passing 3mu",nb,0,nmax);

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

  TH1D *h_nTotStubs_pileup23_skim = new TH1D("h_nTotStubs_pileup23_skim","Total N of Stubs (pileup23_skim events); Number of Stubs; Events",nb,0,nmax);
  TH2D *h2_nStubs_pileup23_skim = new TH2D("h2_nStubs_pileup23_skim","N stubs on S1 vs S0 (pileup23_skim events)",nb,0,nmax,nb,0,nmax);
  TH1D *h_nStubs_0_pileup23_skim = new TH1D("h_nStubs_0_pileup23_skim","Tot N stubs First Station - S0 (pileup23_skim events)",nb,0,nmax);
  TH1D *h_nStubs_1_pileup23_skim = new TH1D("h_nStubs_1_pileup23_skim","Tot N stubs Second Station - S1 (pileup23_skim events)",nb,0,nmax);

  TH1D *h_nTotStubs_pileup234_skim = new TH1D("h_nTotStubs_pileup234_skim","Total N of Stubs (pileup234_skim events); Number of Stubs; Events",nb,0,nmax);
  TH2D *h2_nStubs_pileup234_skim = new TH2D("h2_nStubs_pileup234_skim","N stubs on S1 vs S0 (pileup234_skim events)",nb,0,nmax,nb,0,nmax);
  TH1D *h_nStubs_0_pileup234_skim = new TH1D("h_nStubs_0_pileup234_skim","Tot N stubs First Station - S0 (pileup234_skim events)",nb,0,nmax);
  TH1D *h_nStubs_1_pileup234_skim = new TH1D("h_nStubs_1_pileup234_skim","Tot N stubs Second Station - S1 (pileup234_skim events)",nb,0,nmax);

  TH1D *h_nTotStubs_after_1PUmu = new TH1D("h_nTotStubs_after","Total N of Stubs after skim with 1PUmu",nb,0,nmax);
  TH1D *h_nTotStubs_after = new TH1D("h_nTotStubs_after","Total N of Stubs after skim with 1PUmu",nb,0,nmax);

array<TH1D*,12> stub_per_mod;
for(int m=0; m<12; m++){
string name="stub_per_mod"+to_string(m);
string title="Number of stubs when N_stubs<15 per module "+to_string(m);
stub_per_mod.at(m)=new TH1D(name.c_str(),title.c_str(),20,0,20);
}

array<TH1D*,12> stub_per_mod_after;
for(int m=0; m<12; m++){
string name="stub_after_per_mod"+to_string(m);
string title="Number of stubs after skim per module "+to_string(m);
stub_per_mod_after.at(m)=new TH1D(name.c_str(),title.c_str(),20,0,20);
}

  TH1D *process_st0 = new TH1D("process_st0","Interaction ID of reco tracks in station 0",50,0,50);
  TH1D *process_st1 = new TH1D("process_st1","Interaction ID of reco tracks in station 1",50,0,50);

   const Int_t NBINS = 6;
   Double_t edges[NBINS + 1] = {0.0, 0.005, 0.010, 0.015, 0.020, 0.025, 0.032};
TH1D *theta_e_gen=new TH1D("theta_e_gen", "Electron scattering generated angles from MESMER",64,0.,0.032);
TH1D *theta_mu_gen=new TH1D("theta_mu_gen", "Muon scattering generated angles from MESMER",20,0.,0.005);
TH1D *theta_e_gen_after=new TH1D("theta_e_gen_after", "Electron scattering generated angles from MESMER after skim",64,0.,0.032);
TH1D *theta_mu_gen_after=new TH1D("theta_mu_gen_after", "Muon scattering generated angles from MESMER after skim",20,0.,0.005);

TH2D *theta_2D=new TH2D("theta_2D", "Muon electron scattering angles from MESMER",50,0.,0.002,250,0.,0.005);
TH2D *theta_2D_after=new TH2D("theta_2D", "Muon electron scattering angles from MESMER when N_stubs<15",50,0.,0.002,250,0.,0.005);

  // FOR STATION EFFICIENCIES
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

  //  Long64_t N_single_clean_events = 0;
  Long64_t N_single_clean_stubs_0 = 0;
  Long64_t N_single_clean_stubs_1 = 0;
  Long64_t N_single_clean_stubs = 0;    // total S0+S1

  //  Long64_t N_pileup_skim_events = 0;
  Long64_t N_pileup_skim_stubs_0 = 0;
  Long64_t N_pileup_skim_stubs_1 = 0;
  Long64_t N_pileup_skim_stubs = 0;    // total S0+S1

  //  Long64_t N_pileup23_skim_events = 0;
  Long64_t N_pileup23_skim_stubs_0 = 0;
  Long64_t N_pileup23_skim_stubs_1 = 0;
  Long64_t N_pileup23_skim_stubs = 0;    // total S0+S1

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
  Long64_t N_golden_events = 0;
  Long64_t N_umsel_events = 0;
  Long64_t N_passing_1mu_golden_events = 0;
  Long64_t N_passing_2mu_golden_events = 0;
  Long64_t N_passing_3mu_golden_events = 0;
  Long64_t N_passing_4mu_golden_events = 0;
  Long64_t N_passing_golden_S0_zero_S1_events = 0;
  Long64_t N_passing_golden_S1_zero_S0_events = 0;
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
   //
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
  Long64_t N_pileup23_skim_events = 0;
  Long64_t N_pileup234_skim_events = 0;
  Long64_t N_loose_events = 0;
  Long64_t N_loose23_events = 0;
  Long64_t N_loose234_events = 0;

  Long64_t N_fired6_0 = 0;
  Long64_t N_fired6_1 = 0;
  Long64_t N_fired12_events = 0;
  // i passing 1,2,3,4 mu golden sono un sottinsieme di questi
  // i golden interactions 1+2 sono un sottinsieme di questi
  // single_clean ?
  // pileup 2,3,4 mu ?
  Long64_t N_fired12_plus_single_clean_events = 0;


Double_t n_tot=0.;
Double_t n_filtered=0.;


for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {

 cbmsim->GetEntry(i);
 if(i%100 == 0) cout<<"Entry "<<i<<endl;

int mu_gen=0;
int e_out=-99;
int mu_out=-99;
double xz=10.;
double yz=10.;
double xz_e=10.;
double yz_e=10.;
double xz_mu=10.;
double yz_mu=10.;
double the_gen=0.;
double thmu_gen=0.;
TVector3 p_muin_MC,p_e_MC,p_mu_MC;

           int hit_modXmu=0; int hit_modXe=0;
           int hit_modYmu=0; int hit_modYe=0;
           int stereo_mu=0; int stereo_e=0;
           int hit_modXmuin=0;
           int hit_modYmuin=0;
           int stereo_muin=0; 

        for(int n = 0; n < MCTrack->GetEntries(); n++) {
         const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
         if(MCTr->interactionID()==0 and MCTr->pdgCode()==13){mu_gen++; p_muin_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); p_muin_MC.Unit();
        for(int s=0; s<TrackerPoints->GetEntries(); s++)
                         {const MUonETrackerPoint *TrackerPt = static_cast<const MUonETrackerPoint*>(TrackerPoints->At(s));
         if(TrackerPt->trackPDGCode()==13 and TrackerPt->trackID()==n and TrackerPt->stationID()==0){
                                                                                                 if(TrackerPt->moduleID()==0 or TrackerPt->moduleID()==4) hit_modXmuin++;
                                                                                                 if(TrackerPt->moduleID()==1 or TrackerPt->moduleID()==5) hit_modYmuin++;
                                                                                                 if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_muin++;}
                }
        }
         if(MCTr->interactionID()==45 and MCTr->pdgCode()==11) { e_out=n; p_e_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); p_e_MC.Unit(); the_gen=p_muin_MC.Angle(p_e_MC);

        for(int s=0; s<TrackerPoints->GetEntries(); s++)
                         {const MUonETrackerPoint *TrackerPt = static_cast<const MUonETrackerPoint*>(TrackerPoints->At(s));
         if(TrackerPt->trackPDGCode()==11 and TrackerPt->trackID()==e_out and TrackerPt->stationID()==1){
                                                                                                 if(TrackerPt->moduleID()==0 or TrackerPt->moduleID()==4) hit_modXe++;
                                                                                                 if(TrackerPt->moduleID()==1 or TrackerPt->moduleID()==5) hit_modYe++;
                                                                                                 if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_e++;}
                }
         }
         if(MCTr->interactionID()==45 and MCTr->pdgCode()==-13) {mu_out=n; p_mu_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); p_mu_MC.Unit(); thmu_gen=p_muin_MC.Angle(p_mu_MC);
        for(int s=0; s<TrackerPoints->GetEntries(); s++)
                         {const MUonETrackerPoint *TrackerPt = static_cast<const MUonETrackerPoint*>(TrackerPoints->At(s));
         if(TrackerPt->trackPDGCode()==-13 and TrackerPt->trackID()==mu_out and TrackerPt->stationID()==1){
                                                                                                 if(TrackerPt->moduleID()==0 or TrackerPt->moduleID()==4) hit_modXmu++;
                                                                                                 if(TrackerPt->moduleID()==1 or TrackerPt->moduleID()==5) hit_modYmu++;
                                                                                                 if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_mu++;}
                }
         }
	}

if(mu_gen==1 and e_out!=-99 and mu_out!=-99 and hit_modXmuin==4 and hit_modYmuin==4 and stereo_muin>1 and hit_modXe==4 and hit_modYe==4 and stereo_e>1 and hit_modXmu==4 and hit_modYmu==4 and stereo_mu>1 and thmu_gen>0.0002 ){

    N_input_events++;


    int nTotStubs = TrackerStubs->GetEntries();
    N_input_stubs += nTotStubs;
    h_nTotStubs->Fill(nTotStubs,MesmerEvent->wgt_full);

    std::array<int,12> nstubs = {0,0,0,0,0,0,0,0,0,0,0,0};

    int nStubs_0 = 0;
    int nStubs_1 = 0;

theta_e_gen->Fill(the_gen,MesmerEvent->wgt_full);
theta_mu_gen->Fill(thmu_gen,MesmerEvent->wgt_full);
n_tot+=MesmerEvent->wgt_full;

    ///////////////////////////////////////////////////////
    // loop on all stubs in the two stations
   for(int t=0; t<TrackerStubs->GetEntries(); t++)
   {const MUonETrackerStub *stubs = static_cast<const MUonETrackerStub*>(TrackerStubs->At(t));

         int imod = stubs->moduleID()+stubs->stationID()*6;
        nstubs.at(imod)++;
         h_nstubsPerModule->Fill(imod,MesmerEvent->wgt_full);
         if (imod <6) {
          nStubs_0++;
          h_nstubsPerModule_S0->Fill(imod,MesmerEvent->wgt_full);}
         else {
          nStubs_1++;
          h_nstubsPerModule_S1->Fill(imod,MesmerEvent->wgt_full);
            }
        }

 theta_2D->Fill(the_gen,thmu_gen,MesmerEvent->wgt_full);
if(TrackerStubs->GetEntries()<=14) {theta_2D_after->Fill(the_gen,thmu_gen,MesmerEvent->wgt_full);
for(int s=0; s<12; s++){ stub_per_mod.at(s)->Fill(nstubs.at(s),MesmerEvent->wgt_full);}
}
/*        for(int s=0; s<TrackerPoints->GetEntries(); s++){
        const MUonETrackerPoint *TrackerPt = static_cast<const MUonETrackerPoint*>(TrackerPoints->At(s));
if(TrackerPt->trackPDGCode()==11 and TrackerPt->trackID()==e_out and TrackerPt->stationID()==1){cout << "e) totalEnergyDeposit " << TrackerPt->totalEnergyDeposit() << endl;}
if(TrackerPt->trackPDGCode()==-13 and TrackerPt->trackID()==mu_out and TrackerPt->stationID()==1){cout << "mu) totalEnergyDeposit " << TrackerPt->totalEnergyDeposit() << endl;}
        }

   for(int t=0; t<TrackerStubs->GetEntries(); t++){
	const MUonETrackerStub *stubs = static_cast<const MUonETrackerStub*>(TrackerStubs->At(t));
	vector<MUonETrackDeposit> SeedingDeposits = stubs->SeedingClusterLinkedTrackDeposits();
        vector<MUonETrackDeposit> CorrelationDeposits = stubs->CorrelationClusterLinkedTrackDeposits();

		for(int j=0; j<SeedingDeposits.size(); j++){cout <<"Seed deposit " << SeedingDeposits.at(j) << endl;}
                for(int j=0; j<CorrelationDeposits.size(); j++){cout <<"Correlation deposit " << CorrelationDeposits.at(j) << endl;}
	}
}*/

    // end loop on all stubs in the two stations
    //////////////////////// ///////////////////////////////



//for(int s=0; s<12; s++){ stub_per_mod.at(s)->Fill(nstubs.at(s),MesmerEvent->wgt_full);}

    N_input_stubs_0 += nStubs_0;
    N_input_stubs_1 += nStubs_1;

    h_nStubs_0->Fill(nStubs_0,MesmerEvent->wgt_full);
    h_nStubs_1->Fill(nStubs_1,MesmerEvent->wgt_full);
    h2_nStubs->Fill(nStubs_0, nStubs_1,MesmerEvent->wgt_full);

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
      h_nTotStubs_zero_0->Fill(nTotStubs,MesmerEvent->wgt_full);
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
    //    if (is_fired6_plus_single_clean_0) N_fired6_plus_single_clean_0++;
    
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

    // (5 or more muons)
    //    bool is_pileup_5plus_0 = S_0.fivefired_xy_modules == 4 && S_0.nhits_uv >=5;
    //    if (is_pileup_5plus_0) N_pileup_5plus_0++;
   
    // conditions on station 1
    //
    bool is_zero_1 = nStubs_1 == 0;
    if (is_zero_1) {
      N_zero_1++;
      h_nTotStubs_zero_1->Fill(nTotStubs,MesmerEvent->wgt_full);
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
    //    if (is_fired6_2tracks_1) N_fired6_2tracks_1++;
      
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

      h_nTotStubs_presel->Fill(nTotStubs,MesmerEvent->wgt_full);
      h_nStubs_0_presel->Fill(nStubs_0,MesmerEvent->wgt_full);
      h_nStubs_1_presel->Fill(nStubs_1,MesmerEvent->wgt_full);
      h2_nStubs_presel->Fill(nStubs_0, nStubs_1,MesmerEvent->wgt_full);
      
    }

    // S0 trackable and increasing requests in S1 
    if (is_trackable_0 && fired_modules_1 >= 1) {
      N_trackable_0_onehit_1++;
      N_trackable_0_onehit_1_stubs += nTotStubs;
      N_trackable_0_onehit_1_stubs_0 += nStubs_0;
      N_trackable_0_onehit_1_stubs_1 += nStubs_1;
      
      h_nTotStubs_S0_trackable_S1_onehit->Fill(nTotStubs,MesmerEvent->wgt_full);
      h_nStubs_0_S0_trackable_S1_onehit->Fill(nStubs_0,MesmerEvent->wgt_full);
      h_nStubs_1_S0_trackable_S1_onehit->Fill(nStubs_1,MesmerEvent->wgt_full);
      h2_nStubs_S0_trackable_S1_onehit->Fill(nStubs_0,nStubs_1,MesmerEvent->wgt_full);

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
      
    // at least a muon trackable in both stations
    bool is_trackable_event = is_trackable_0 && is_trackable_1;
    if (is_trackable_event) {
      N_trackable_events++;
      h_nTotStubs_trackable->Fill(nTotStubs,MesmerEvent->wgt_full);

      h_nStubs_0_trackable->Fill(nStubs_0,MesmerEvent->wgt_full);
      h_nStubs_1_trackable->Fill(nStubs_1,MesmerEvent->wgt_full);
      h2_nStubs_trackable->Fill(nStubs_0, nStubs_1,MesmerEvent->wgt_full);
    }
    
    bool is_noise_event = is_noise_0 && is_noise_1;
    if (is_noise_event) {
      N_noise_events++;
      h_nTotStubs_noise->Fill(nTotStubs,MesmerEvent->wgt_full);
    }
    
    bool is_missing_0_event = is_noise_0 &&  is_trackable_1;
    if (is_missing_0_event) {
      N_missing_0_events++;
      h_nTotStubs_missing_0->Fill(nTotStubs,MesmerEvent->wgt_full);
    }

    bool is_missing_1_event = is_trackable_0 &&  is_noise_1;
    if (is_missing_1_event) {
      N_missing_1_events++;
      h_nTotStubs_missing_1->Fill(nTotStubs,MesmerEvent->wgt_full);
    }

    bool is_not_trackable_event = !is_trackable_0 || !is_trackable_1;
    if (is_not_trackable_event) {
      N_not_trackable_events++;
      h_nTotStubs_not_trackable->Fill(nTotStubs,MesmerEvent->wgt_full); 
    }

    bool is_presel_not_trackable_event = is_presel_event && is_not_trackable_event;
    if (is_presel_not_trackable_event) {
      N_presel_not_trackable_events++;
      h_nTotStubs_presel_not_trackable->Fill(nTotStubs,MesmerEvent->wgt_full);
    }

    bool is_fired12_event = is_fired6_0 && is_fired6_1;
    if (is_fired12_event) N_fired12_events++;
    
    bool is_passing_golden_event = is_passing_golden_0 && is_passing_golden_1;
    if (is_passing_golden_event) {
      N_passing_1mu_golden_events++;
      //      h_nTotStubs_passing_golden->Fill(nTotStubs);
    }

    bool is_presel_passing_golden_S0_event = is_presel_event && is_passing_golden_0;
    if (is_presel_passing_golden_S0_event) {
      N_presel_passing_golden_S0_events++;
    }
    
    bool is_presel_passing_golden_S1_event = is_presel_event && is_passing_golden_1;
    if (is_presel_passing_golden_S1_event) {
      N_presel_passing_golden_S1_events++;
    }

    // *************** efficiency 2 mu ***********************
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
    
    // *************** efficiency 3 mu ***********************
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
    
    // *************** efficiency 4 mu ***********************
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
    // ***************
    
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
      h_nTotStubs_passing_1mu->Fill(nTotStubs,MesmerEvent->wgt_full);
    }
    bool is_passing_2mu_event = is_pileup_2mu_0 && is_pileup_2mu_1;
    if (is_passing_2mu_event) {
      N_passing_2mu_events++;
      h_nTotStubs_passing_2mu->Fill(nTotStubs,MesmerEvent->wgt_full);
    }
    bool is_passing_3mu_event = is_pileup_3mu_0 && is_pileup_3mu_1;
    if (is_passing_3mu_event) {
      N_passing_3mu_events++;
      h_nTotStubs_passing_3mu->Fill(nTotStubs,MesmerEvent->wgt_full);
    }
    bool is_passing_4mu_event = is_pileup_4mu_0 && is_pileup_4mu_1;
    if (is_passing_4mu_event) N_passing_4mu_events++;

    bool is_single_cand_event = is_single_cand_0 && is_2tracks_1;
    if (is_single_cand_event) {
      N_single_cand_events++;
      h_nTotStubs_single_cand->Fill(nTotStubs,MesmerEvent->wgt_full);
    }

    bool is_single_clean_event = is_single_clean_0 && is_2tracks_1;
    if (is_single_clean_event) {
      N_single_clean_events++;
      h_nTotStubs_single_clean->Fill(nTotStubs,MesmerEvent->wgt_full);
      if (iskim == 1) ok_skim = true;

      N_single_clean_stubs += nTotStubs;
      N_single_clean_stubs_0 += nStubs_0;
      N_single_clean_stubs_1 += nStubs_1;

      h_nStubs_0_single_clean->Fill(nStubs_0,MesmerEvent->wgt_full);
      h_nStubs_1_single_clean->Fill(nStubs_1,MesmerEvent->wgt_full);
      h2_nStubs_single_clean->Fill(nStubs_0, nStubs_1,MesmerEvent->wgt_full);
      
    }

    bool is_fired12_plus_single_clean_event = is_fired6_plus_single_clean_0 && is_fired6_2tracks_1;
    if (is_fired12_plus_single_clean_event) {
      N_fired12_plus_single_clean_events++;
      //      h_nTotStubs_single_clean->Fill(nTotStubs);
      //      N_single_clean_stubs += nTotStubs;
      //      N_single_clean_stubs_0 += nStubs_0;
      //      N_single_clean_stubs_1 += nStubs_1;
      //
      //      h_nStubs_0_single_clean->Fill(nStubs_0);
      //      h_nStubs_1_single_clean->Fill(nStubs_1);
      //      h2_nStubs_single_clean->Fill(nStubs_0, nStubs_1);
    }
    
    bool is_umsel_event = nStubs_0 >=5 && nStubs_1 >=5 && (nStubs_1 - nStubs_0) >= 5;
    if (is_umsel_event) {
      N_umsel_events++;
      h_nTotStubs_umsel->Fill(nTotStubs,MesmerEvent->wgt_full);
    }
    
    bool is_golden_event = is_passing_golden_0 && is_golden_1;
    if (is_golden_event) {
      N_golden_events++;
      h_nTotStubs_golden->Fill(nTotStubs,MesmerEvent->wgt_full);
    }
    
    bool is_pileup_any_event = is_pileup_any_0 && is_signal_plus_pileup_1;
    if (is_pileup_any_event) {
      N_pileup_any_events++;
      h_nTotStubs_pileup_any->Fill(nTotStubs,MesmerEvent->wgt_full);
    }    
    bool is_pileup_2mu_event = is_pileup_2mu_0 && is_signal_plus_pileup_1;
    if (is_pileup_2mu_event) {
      N_pileup_2mu_events++;
      h_nTotStubs_pileup_2mu->Fill(nTotStubs,MesmerEvent->wgt_full);
    }
    bool is_pileup_3mu_event = is_pileup_3mu_0 && is_signal_plus_tripileup_1;
    if (is_pileup_3mu_event) {
      N_pileup_3mu_events++;
      h_nTotStubs_pileup_3mu->Fill(nTotStubs,MesmerEvent->wgt_full);
    }
    bool is_pileup_4mu_event = is_pileup_4mu_0 && is_signal_plus_fourpileup_1;
    if (is_pileup_4mu_event) {
      N_pileup_4mu_events++;
      h_nTotStubs_pileup_4mu->Fill(nTotStubs,MesmerEvent->wgt_full);
    }
    bool is_pileup_many_event = is_pileup_many_0 && is_signal_plus_fourpileup_1;
    if (is_pileup_many_event) {
      N_pileup_many_events++;
      h_nTotStubs_pileup_many->Fill(nTotStubs),MesmerEvent->wgt_full;     
    }

    bool is_pileup_skim_event = is_pileup_2mu_event || is_pileup_3mu_event || is_pileup_many_event;
    if (is_pileup_skim_event) {
      N_pileup_skim_events++;
      h_nTotStubs_pileup_skim->Fill(nTotStubs,MesmerEvent->wgt_full);     
      if (iskim == 40) ok_skim = true;

      N_pileup_skim_stubs += nTotStubs;
      N_pileup_skim_stubs_0 += nStubs_0;
      N_pileup_skim_stubs_1 += nStubs_1;

      h_nStubs_0_pileup_skim->Fill(nStubs_0,MesmerEvent->wgt_full);
      h_nStubs_1_pileup_skim->Fill(nStubs_1,MesmerEvent->wgt_full);
      h2_nStubs_pileup_skim->Fill(nStubs_0, nStubs_1,MesmerEvent->wgt_full);
      
    }

    bool is_pileup23_skim_event = is_pileup_2mu_event || is_pileup_3mu_event;
    if (is_pileup23_skim_event) {
      N_pileup23_skim_events++;
      h_nTotStubs_pileup23_skim->Fill(nTotStubs,MesmerEvent->wgt_full);     
      if (iskim == 50) ok_skim = true;

      N_pileup23_skim_stubs += nTotStubs;
      N_pileup23_skim_stubs_0 += nStubs_0;
      N_pileup23_skim_stubs_1 += nStubs_1;

      h_nStubs_0_pileup23_skim->Fill(nStubs_0,MesmerEvent->wgt_full);
      h_nStubs_1_pileup23_skim->Fill(nStubs_1,MesmerEvent->wgt_full);
      h2_nStubs_pileup23_skim->Fill(nStubs_0, nStubs_1,MesmerEvent->wgt_full);
      
    }

    bool is_pileup234_skim_event = is_pileup_2mu_event || is_pileup_3mu_event || is_pileup_4mu_event;
    if (is_pileup234_skim_event) {
      N_pileup234_skim_events++;
      h_nTotStubs_pileup234_skim->Fill(nTotStubs,MesmerEvent->wgt_full);     
      if (iskim == 60) ok_skim = true;

      N_pileup234_skim_stubs += nTotStubs;
      N_pileup234_skim_stubs_0 += nStubs_0;
      N_pileup234_skim_stubs_1 += nStubs_1;

      h_nStubs_0_pileup234_skim->Fill(nStubs_0,MesmerEvent->wgt_full);
      h_nStubs_1_pileup234_skim->Fill(nStubs_1,MesmerEvent->wgt_full);
      h2_nStubs_pileup234_skim->Fill(nStubs_0, nStubs_1,MesmerEvent->wgt_full);
      
    }
    
    bool is_loose_event = is_single_clean_event || is_pileup_2mu_event || is_pileup_3mu_event || is_pileup_many_event;
    if (is_loose_event) {
      N_loose_events++;
      if (iskim == 41) ok_skim = true;
    }

    bool is_loose23_event = is_single_clean_event || is_pileup_2mu_event || is_pileup_3mu_event;
    if (is_loose23_event) {
      N_loose23_events++;
      if (iskim == 51) ok_skim = true;
    }

    bool is_loose234_event = is_single_clean_event || is_pileup_2mu_event || is_pileup_3mu_event || is_pileup_4mu_event;
    if (is_loose234_event) {
      N_loose234_events++;
      if (iskim == 61) ok_skim = true;
    }

    if (ok_skim) {
n_filtered+=MesmerEvent->wgt_full;

h_nTotStubs_after->Fill(TrackerStubs->GetEntries(),MesmerEvent->wgt_full);

theta_e_gen_after->Fill(the_gen,MesmerEvent->wgt_full);
theta_mu_gen_after->Fill(thmu_gen,MesmerEvent->wgt_full);
for(int s=0; s<12; s++){ stub_per_mod_after.at(s)->Fill(nstubs.at(s),MesmerEvent->wgt_full);}

vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
         for (auto&& track : tracks) {
//if(e_out==track.linkedTrackID()){process_st0->Fill(track.processIDofLinkedTrack());}
//if(mu_out==track.linkedTrackID()){process_st1->Fill(track.processIDofLinkedTrack());}
if(track.sector()==0){process_st0->Fill(track.processIDofLinkedTrack(),MesmerEvent->wgt_full);}
if(track.sector()==1){process_st1->Fill(track.processIDofLinkedTrack(),MesmerEvent->wgt_full);}
}


      /*  o_MCTrack = MCTrack;
        o_TrackerStubs = TrackerStubs;
        o_ReconstructionOutput = ReconstructionOutput;

      o_tree->Fill();*/
      o_nevt++;

      N_stubs_skim += nTotStubs;
    }
   }// if mu in
  } //  for(Long64_t i = 0; i < n_entries; i++)

/*
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
  cout<<"Number of golden 1mu passing events = " << endl;
  cout<<"\t "<< N_passing_1mu_golden_events <<endl;

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
  cout<<"Total number of pileup (2,3) interaction events = " << endl;
  cout<<"\t "<< N_pileup23_skim_events <<endl;
  cout<<"Total number of pileup (2,3,4) interaction events = " << endl;
  cout<<"\t "<< N_pileup234_skim_events <<endl;
  cout<<"Total number of pileup (2,3,>=4) interaction events = " << endl;
  cout<<"\t "<< N_pileup_skim_events <<endl;
  cout<<"Total number of Single + Pileup (2,3) interaction events = " << endl;
  cout<<"\t "<< N_loose23_events <<endl;
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
  cout<<"Fraction Passing Golden mu       = "<< double(N_passing_1mu_golden_events)/double(N_input_events) <<endl;

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

  // Average Station efficiencies
  pair<double,double> ineff_zero_S1 = binomial_eff(N_passing_golden_S0_zero_S1_events, N_passing_1mu_golden_0);
  pair<double,double> eff_S1 = binomial_eff(N_passing_golden_S0_trackable_S1_events, N_passing_1mu_golden_0);
  pair<double,double> eff_presel_S1 = binomial_eff(N_passing_golden_S0_trackable_S1_events, N_presel_passing_golden_S0_events);

  pair<double,double> ineff_zero_S0 = binomial_eff(N_passing_golden_S1_zero_S0_events, N_passing_1mu_golden_1);
  pair<double,double> eff_S0 = binomial_eff(N_passing_golden_S1_trackable_S0_events, N_passing_1mu_golden_1);
  pair<double,double> eff_presel_S0 = binomial_eff(N_passing_golden_S1_trackable_S0_events, N_presel_passing_golden_S1_events);

  cout<<endl
      <<"estimated Efficiency reconstruction S_1 (after presel) = "<< eff_presel_S1.first << " +- " << eff_presel_S1.second <<endl; 
  cout<<"estimated Efficiency reconstruction S_1 (all events)   = "<< eff_S1.first << " +- " << eff_S1.second <<endl; 
  cout<<"\t estimated inefficiency (zero hits in S_1)  = "<< ineff_zero_S1.first << " +- " << ineff_zero_S1.second <<endl; 
  cout<<"estimated Efficiency reconstruction S_0 (after presel) = "<< eff_presel_S0.first << " +- " << eff_presel_S0.second<<endl;
  cout<<"estimated Efficiency reconstruction S_0 (all events)   = "<< eff_S0.first << " +- " << eff_S0.second<<endl;
  cout<<"\t estimated inefficiency (zero hits in S_0)  = "<< ineff_zero_S0.first << " +- " << ineff_zero_S0.second <<endl; 


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
  cout<<"Total Fraction PU 2+3 Mu Int.    = "<< double(N_pileup23_skim_events)/double(N_input_events) <<endl;
  cout<<"Total Fraction PU 2+3+4 Mu Int.  = "<< double(N_pileup234_skim_events)/double(N_input_events) <<endl;
  cout<<"Total Fraction PU Mu Int.        = "<< double(N_pileup_skim_events)/double(N_input_events) <<endl;
  cout<<"Total Fraction Single+PU 2+3     = "<< double(N_loose23_events)/double(N_input_events) <<endl;
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
  cout<<"Total Fraction PU 2+3 Mu Int.    = "<< double(N_pileup23_skim_events)/double(N_presel_events) <<endl;
  cout<<"Total Fraction PU 2+3+4 Mu Int.  = "<< double(N_pileup234_skim_events)/double(N_presel_events) <<endl;
  cout<<"Total Fraction PU Mu Int.        = "<< double(N_pileup_skim_events)/double(N_presel_events) <<endl;
  cout<<"Total Fraction Single+PU 2+3     = "<< double(N_loose23_events)/double(N_presel_events) <<endl;
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
  cout<<"Total Fraction PU 2+3 Mu Int.    = "<< double(N_pileup23_skim_events)/double(N_trackable_0_onehit_1) <<endl;
  cout<<"Total Fraction PU 2+3+4 Mu Int.  = "<< double(N_pileup234_skim_events)/double(N_trackable_0_onehit_1) <<endl;
  cout<<"Total Fraction PU Mu Int.        = "<< double(N_pileup_skim_events)/double(N_trackable_0_onehit_1) <<endl;
  cout<<"Total Fraction Single+PU 2+3     = "<< double(N_loose23_events)/double(N_trackable_0_onehit_1) <<endl;
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
  cout<<"Total Fraction PU 2+3 Mu Int.    = "<< double(N_pileup23_skim_events)/double(N_trackable_events) <<endl;
  cout<<"Total Fraction PU 2+3+4 Mu Int.  = "<< double(N_pileup234_skim_events)/double(N_trackable_events) <<endl;
  cout<<"Total Fraction PU Mu Int.        = "<< double(N_pileup_skim_events)/double(N_trackable_events) <<endl;
  cout<<"Total Fraction Single+PU 2+3     = "<< double(N_loose23_events)/double(N_trackable_events) <<endl;
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
  cout<<"Total Fraction PU 2+3 Mu Int.    = "<< double(N_pileup23_skim_events)/double(N_passing_1mu_golden_events) <<endl;
  cout<<"Total Fraction PU 2+3+4 Mu Int.  = "<< double(N_pileup234_skim_events)/double(N_passing_1mu_golden_events) <<endl;
  cout<<"Total Fraction PU Mu Int.        = "<< double(N_pileup_skim_events)/double(N_passing_1mu_golden_events) <<endl;
  cout<<"Total Fraction Single+PU 2+3     = "<< double(N_loose23_events)/double(N_passing_1mu_golden_events) <<endl;
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
  cout<<"Fraction Single Mu Int.+FIRED12  = "<< double(N_fired12_plus_single_clean_events)/double(N_fired12_events) <<endl;
  cout<<"Fraction Golden sel              = "<< double(N_golden_events)/double(N_fired12_events) <<endl;
  cout<<"Fraction 2 Mu PU Int.            = "<< double(N_pileup_2mu_events)/double(N_fired12_events) <<endl;
  cout<<"Fraction 3 Mu PU Int.            = "<< double(N_pileup_3mu_events)/double(N_fired12_events) <<endl;
  cout<<"Fraction 4 Mu PU Int.            = "<< double(N_pileup_4mu_events)/double(N_fired12_events) <<endl;
  cout<<"Fraction >=4 Mu PU Int.          = "<< double(N_pileup_many_events)/double(N_fired12_events) <<endl;
  cout<<"Total Fraction PU 2+3 Mu Int.    = "<< double(N_pileup23_skim_events)/double(N_fired12_events) <<endl;
  cout<<"Total Fraction PU 2+3+4 Mu Int.  = "<< double(N_pileup234_skim_events)/double(N_fired12_events) <<endl;
  cout<<"Total Fraction PU Mu Int.        = "<< double(N_pileup_skim_events)/double(N_fired12_events) <<endl;
  cout<<"Total Fraction Single+PU 2+3     = "<< double(N_loose23_events)/double(N_fired12_events) <<endl;
  cout<<"Total Fraction Single+PU 2+3+4   = "<< double(N_loose234_events)/double(N_fired12_events) <<endl;
  cout<<"Total Fraction Single+PU 2+3+4+  = "<< double(N_loose_events)/double(N_fired12_events) <<endl;
  */
  hfile->Write();

  if (iskim > 0) {
    cout<<"\n Number of output skimmed events          = "<<o_nevt<<endl;
    cout<<" Number of stubs in output skimmed events = "<<N_stubs_skim<<endl;
  }



TCanvas s("s","s",700,700);
s.Divide(1,2);
s.cd(1);
process_st0->Draw("hist");
s.cd(2);
process_st1->Draw("hist");
s.SaveAs("outputs_mesmer/process_1PUmu_02cut_0_2_mrad.pdf");

TH1D * h2 = (TH1D*) theta_e_gen_after->Clone();
TH1D * h2gen = (TH1D*) theta_e_gen->Clone();
h2->Divide(h2,h2gen,1,1,"B");
TH1D * h3 = (TH1D*) theta_mu_gen_after->Clone();
TH1D * h3gen = (TH1D*) theta_mu_gen->Clone();
h3->Divide(h3,h3gen,1,1,"B");

Int_t nxTh = theta_2D->GetNbinsX();
Int_t nyTh = theta_2D->GetNbinsY();
for (Int_t i=1; i<nxTh+1; i++) {
for (Int_t j=1; j<nyTh+1; j++) {
if (theta_2D->GetBinContent(i,j)<1) theta_2D->SetBinContent(i,j,0);}}

Int_t nxTh_after = theta_2D_after->GetNbinsX();
Int_t nyTh_after = theta_2D_after->GetNbinsY();
for (Int_t i=1; i<nxTh_after+1; i++) {
for (Int_t j=1; j<nyTh_after+1; j++) {
if (theta_2D_after->GetBinContent(i,j)<1) theta_2D_after->SetBinContent(i,j,0);}}

TCanvas a("a","a",700,700);
a.Divide(1,3);
a.cd(1);
h_nTotStubs->Draw("hist");
h_nTotStubs->SetTitle(" ");
h_nTotStubs->SetMinimum(1.);
h_nTotStubs_after->SetLineColor(kRed);
h_nTotStubs_after->Draw("hist same");
gPad->SetLogy();
gStyle->SetOptStat(0);
a.cd(2);
theta_e_gen->Draw("hist");
theta_e_gen->SetTitle(" ");
theta_e_gen->SetMinimum(1.);
theta_e_gen_after->SetLineColor(kRed);
theta_e_gen_after->Draw("hist same");
gStyle->SetOptStat(0);
a.cd(3);
theta_mu_gen->SetTitle(" ");
theta_mu_gen->SetMinimum(1.);
theta_mu_gen->Draw("hist");
theta_mu_gen_after->SetLineColor(kRed);
theta_mu_gen_after->Draw("hist same");
gPad->SetLogy();
gStyle->SetOptStat(0);
a.SaveAs("outputs_mesmer/nstubs_1PUmu_02cut_0_2_mrad.pdf");

TCanvas b("b","b",700,700);
b.Divide(1,2);
b.cd(1);
theta_2D->Draw();
b.cd(2);
theta_2D_after->Draw();
b.SaveAs("outputs_mesmer/2D_plots_02cut_0_2_mrad.pdf");


cout << "electron " << endl;
for(int b=1; b<64; b++){
cout << b << " ) " << h2->GetBinContent(b) << " +- " << h2->GetBinError(b)<< " ---> " << h2->GetBinError(b)/h2->GetBinContent(b)*100 << "%" <<endl ;
}
cout << "muon " << endl;
for(int b=1; b<25; b++){
cout << b << " ) " << h3->GetBinContent(b) << " +- " << h3->GetBinError(b)<< " ----> " << h3->GetBinError(b)/h3->GetBinContent(b)*100 << "%" << endl;
}

TCanvas a1("a1","a1",700,700);
a1.Divide(1,2);
a1.cd(1);
h2->SetMinimum(0.845);
h2->SetMaximum(1.02);
h2->Draw("E");
gStyle->SetOptStat(0);
a1.cd(2);
h3->SetMinimum(0.945);;
h3->SetMaximum(1.02);
h3->Draw("E");
gStyle->SetOptStat(0);
a1.SaveAs("outputs_mesmer/eff_skim_02cut_0_2_mrad.pdf");
h2->SaveAs("outputs_mesmer/eff_skim_e_02cut_0_2_mrad.root");
h3->SaveAs("outputs_mesmer/eff_skim_mu_02cut_0_2_mrad.root");

TCanvas a2("a2","a2",700,700);
a2.Divide(3,4);
for(int m=0; m<12; m++){
a2.cd(m+1);
stub_per_mod.at(m)->SetMinimum(1.0);
stub_per_mod.at(m)->SetLineWidth(3);
stub_per_mod.at(m)->Draw("hist");
//stub_per_mod_after.at(m)->SetLineColor(kRed);
//stub_per_mod_after.at(m)->Draw("hist same");
gPad->SetLogy();
}
a2.SaveAs("outputs_mesmer/stub_mod_02cut_0_2_mrad.pdf");



}













