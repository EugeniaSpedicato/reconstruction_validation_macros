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


int counting_single_int(int nhits,string run, string type, string save){

  int nthreads = 7;


  ROOT::EnableImplicitMT(nthreads);


TChain * cbmsim = new TChain("cbmsim");
TChain * cbmsim_g = new TChain("cbmsim");

if(run=="run8" and nhits==0 and type=="single_mu_int_0"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run8/single_muon_interaction_0/reco/single_muon_interaction_0_dataProcessor.root");
}
else if(run=="run8" and nhits==0 and type=="single_mu_int_1"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run8/single_muon_interaction_1/reco/single_muon_interaction_1_dataProcessor.root");
}
else if(run=="run11" and nhits==0 and type=="unbiased"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run11/unbiased/unbiased_500files_dataProcessor.root");
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run11/unbiased/unbiased_500files_dataProcessor_1.root");
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run11/unbiased/unbiased_500files_dataProcessor_2.root");
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run11/unbiased/unbiased_500files_dataProcessor_3.root");
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run11/unbiased/unbiased_500files_dataProcessor_4.root");
}
else if(run=="run11" and nhits==0 and type=="unbiased_1"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run11/unbiased/unbiased_2files.root");
}
else if(run=="run16" and nhits==0 and type=="single_mu_int_0"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run16/single_muon_interaction_0/4files_commit62871d8b.root");
}
else if(run=="run16" and nhits==0 and type=="single_mu_int_1"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run16/single_muon_interaction_1/4files_commit62871d8b.root");
}
else if(run=="run11" and nhits==0 and type=="single_mu_int_0"){
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run11/single_muon_interaction_0/4files_commit62871d8b.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run11/single_muon_interaction_0/allRun_0hit_noMuPid.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run11/single_muon_interaction_0/allRun_0hit_noMuPid_1.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run11/single_muon_interaction_0/allRun_0hit_noMuPid_2.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run11/single_muon_interaction_0/allRun_0hit_noMuPid_3.root");
}
else if(run=="run11" and nhits==0 and type=="single_mu_int_1"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run11/single_muon_interaction_1/4files_commit62871d8b.root");
}
else if(run=="run20" and nhits==0 and type=="single_mu_int_1"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run20/single_muon_interaction_1/2files_commit62871d8b_pt1.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run20/single_muon_interaction_1/2files_commit62871d8b_pt2.root");
}
else if(run=="run20" and nhits==0 and type=="single_mu_int_0"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run20/single_muon_interaction_0/2files_commit62871d8b_pt1.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run20/single_muon_interaction_0/2files_commit62871d8b_pt2.root");
}
else if(run=="run20_newAl" and nhits==0 and type=="single_mu_int_0"){
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run20/single_muon_interaction_0/2files_commit62871d8b_pt1_newAl.root");
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run20/single_muon_interaction_0/2files_commit62871d8b_pt2_newAl.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run20/single_muon_interaction_0/4files_commit62871d8b_newAl2.root");
}
else if(run=="run20_newAl" and nhits==0 and type=="single_mu_int_1"){
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run20/single_muon_interaction_1/2files_commit62871d8b_pt1_newAl.root");
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run20/single_muon_interaction_1/2files_commit62871d8b_pt2_newAl.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run20/single_muon_interaction_1/4files_commit62871d8b_newAl2.root");
}

else if(run=="run17" and nhits==0 and type=="single_mu_int_0"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run17/single_muon_interaction_0/1files_commit62871d8b.root");
}
else if(run=="run17" and nhits==0 and type=="single_mu_int_1"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run17/single_muon_interaction_1/1files_commit62871d8b.root");
}



else if(run=="run17_0800" and nhits==0 and type=="single_mu_int_0"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run17/single_muon_interaction_0/1files_commit62871d8b_0800.root");
}
else if(run=="run17_0800" and nhits==0 and type=="single_mu_int_1"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run17/single_muon_interaction_1/1files_commit62871d8b_0800.root");
}



else if(run=="run18" and nhits==0 and type=="single_mu_int_0"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run18/single_muon_interaction_0/1files_commit62871d8b.root");
}
else if(run=="run18" and nhits==0 and type=="single_mu_int_1"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run18/single_muon_interaction_1/1files_commit62871d8b.root");
}

else if(run=="run18_0041" and nhits==0 and type=="single_mu_int_0"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run18/single_muon_interaction_0/1files_commit62871d8b_0041.root");
}
else if(run=="run18_0041" and nhits==0 and type=="single_mu_int_1"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run18/single_muon_interaction_1/1files_commit62871d8b_0041.root");
}

else if(run=="run18_0041_newAl" and nhits==0 and type=="single_mu_int_0"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run18/single_muon_interaction_0/1files_commit62871d8b_0041_newAl.root");
}
else if(run=="run18_0041_newAl" and nhits==0 and type=="single_mu_int_1"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run18/single_muon_interaction_1/1files_commit62871d8b_0041_newAl.root");
}



else if(run=="run19" and nhits==0 and type=="single_mu_int_0"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run19/single_muon_interaction_0/1files_commit62871d8b.root");
}
else if(run=="run19" and nhits==0 and type=="single_mu_int_1"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run19/single_muon_interaction_1/1files_commit62871d8b.root");
}
else if(run=="run19_pt1" and nhits==0 and type=="single_mu_int_0"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run19/single_muon_interaction_0/1files_commit62871d8b_pt1.root");
}
else if(run=="run19_pt1" and nhits==0 and type=="single_mu_int_1"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run19/single_muon_interaction_1/1files_commit62871d8b_pt1.root");
}


else if(run=="run19_1015" and nhits==0 and type=="single_mu_int_0"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run19/single_muon_interaction_0/1files_commit62871d8b_1015.root");
}
else if(run=="run19_1015" and nhits==0 and type=="single_mu_int_1"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run19/single_muon_interaction_1/1files_commit62871d8b_1015.root");
}


else if(run=="run19_1020" and nhits==0 and type=="single_mu_int_0"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run19/single_muon_interaction_0/1files_commit62871d8b_1020.root");
}
else if(run=="run19_1020" and nhits==0 and type=="single_mu_int_1"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run19/single_muon_interaction_1/1files_commit62871d8b_1020.root");
}


else if(run=="run19_1517" and nhits==0 and type=="single_mu_int_0"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run19/single_muon_interaction_0/1files_commit62871d8b_1020.root");
}
else if(run=="run19_1517" and nhits==0 and type=="single_mu_int_1"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run19/single_muon_interaction_1/1files_commit62871d8b_1020.root");
}



else if(run=="run20_1130" and nhits==0 and type=="single_mu_int_0"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run20/single_muon_interaction_0/1files_commit62871d8b_1130.root");
}
else if(run=="run20_1130" and nhits==0 and type=="single_mu_int_1"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run20/single_muon_interaction_1/1files_commit62871d8b_1130.root");
}



else if(run=="run20_1400" and nhits==0 and type=="single_mu_int_0"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run20/single_muon_interaction_0/1files_commit62871d8b_1400.root");
}
else if(run=="run20_1400" and nhits==0 and type=="single_mu_int_1"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run20/single_muon_interaction_1/1files_commit62871d8b_1400.root");
}


else if(run=="run20_1230" and nhits==0 and type=="single_mu_int_0"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run20/single_muon_interaction_0/1files_commit62871d8b_1230.root");
}
else if(run=="run20_1230" and nhits==0 and type=="single_mu_int_1"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run20/single_muon_interaction_1/1files_commit62871d8b_1230.root");
}

else if(run=="run20_1447" and nhits==0 and type=="single_mu_int_0"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run20/single_muon_interaction_0/1files_commit62871d8b_1447.root");
}
else if(run=="run20_1447" and nhits==0 and type=="single_mu_int_1"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run20/single_muon_interaction_1/1files_commit62871d8b_1447.root");
}


else if(run=="run21" and nhits==0 and type=="single_mu_int_0"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run21/single_muon_interaction_0/1files_commit62871d8b.root");
}
else if(run=="run21" and nhits==0 and type=="single_mu_int_1"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run21/single_muon_interaction_1/1files_commit62871d8b.root");
}

else if(run=="run23" and nhits==0 and type=="single_mu_int_0"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run23/single_muon_interaction_0/1files_commit62871d8b.root");
}
else if(run=="run23" and nhits==0 and type=="single_mu_int_1"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run23/single_muon_interaction_1/1files_commit62871d8b.root");
}

else if(run=="run24" and nhits==0 and type=="single_mu_int_0"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run24/single_muon_interaction_0/1files_commit62871d8b.root");
}

else if(run=="run50" and nhits==0 and type=="single_mu_int_0"){
cbmsim->Add("/home/espedica/clement_fair_install/instFairRoot/share/MUonE/macros/run50_reconstruction.root");
}

else if(run=="run48_3" and nhits==0 and type=="single_mu_int_0"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run48/single_muon_interaction_0/3files_noMuPid_alignment48.root");
}
else if(run=="run48_3" and nhits==0 and type=="single_mu_int_1"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run48/single_muon_interaction_1/3files_noMuPid_alignment48.root");
}
else if(run=="run11_mark" and nhits==0 and type=="single_mu_int_0"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run11/single_muon_interaction_0/2files_0hit_noMuPid.root");
}
else if(run=="run11_mark" and nhits==0 and type=="single_mu_int_1"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run11/single_muon_interaction_1/2files_0hit_noMuPid.root");
}
else if(run=="run29_0" and nhits==0 and type=="single_mu_int_0"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run29/single_muon_interaction_0/muedaq03-1753800115_0hit.root");
}
else if(run=="run29_1" and nhits==0 and type=="single_mu_int_1"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run29/single_muon_interaction_1/muedaq03-1753800115_0hit.root");
}




cout << "cbmsim->GetEntries() " << cbmsim->GetEntries() <<endl;

ROOT::TTreeProcessorMT tp1(*cbmsim,nthreads);

      MUonERecoOutputAnalysis *ReconstructionOutput = 0;



TH1::SetDefaultSumw2(kTRUE);
   const Int_t NBINS = 6;
   Double_t edges[NBINS + 1] = {0.0, 0.005, 0.010, 0.015, 0.020, 0.025, 0.032};

   const Int_t NBINS2 = 12;
   Double_t edges2[NBINS2 + 1] = {0.,0.0002,0.0004,0.0006,0.0008,0.001,0.00125,0.0015,0.00175,0.002,0.003,0.004,0.005};


   ROOT::TThreadedObject<TH1D> nstub_mod_s0("h_nstub_mod_s0","nstubs per mod station0",6,0,6);
   ROOT::TThreadedObject<TH1D> nstub_mod_s1("h_nstub_mod_s1","nstubs per mod station1",6,0,6);
   ROOT::TThreadedObject<TH1D> nstub_mod_s2("h_nstub_mod_s2","nstubs per mod station2",6,0,6);

   ROOT::TThreadedObject<TH1D> h_nstubs_s0("h_nstubs_s0","nstubs station0",30,0,30);
   ROOT::TThreadedObject<TH1D> h_nstubs_s1("h_nstubs_s1","nstubs station1",30,0,30);
   ROOT::TThreadedObject<TH1D> h_nstubs_s2("h_nstubs_s2","nstubs station2",30,0,30);

   ROOT::TThreadedObject<TH1D> h_ntracks_s0("h_ntracks_s0","ntracks station0",30,0,30);
   ROOT::TThreadedObject<TH1D> h_ntracks_s1("h_ntracks_s1","ntracks station1",30,0,30);
   ROOT::TThreadedObject<TH1D> h_ntracks_s2("h_ntracks_s2","ntracks station2",30,0,30);

   ROOT::TThreadedObject<TH1D> h_nvrtx_1("h_nvrtx_s2","nvrtxs station01",10,0,10);
   ROOT::TThreadedObject<TH1D> h_nvrtx_2("h_nvrtx_s2","nvrtxs station12",10,-5,5);

   ROOT::TThreadedObject<TH1D> vrtx_chi2_1("vrtx_chi2_1","vrtx_chi2 station01",80,0,20);
   ROOT::TThreadedObject<TH1D> vrtx_chi2_2("vrtx_chi2_2","vrtx_chi2 station12",80,0,20);

   ROOT::TThreadedObject<TH1D> count_tr_1("count_chi1","count_chi1",2,0,2);
   ROOT::TThreadedObject<TH1D> count_tr_2("count_chi2","count_chi2",2,0,2);

   ROOT::TThreadedObject<TH1D> good_vrtx_1("good_vrtx_1","good_vrtx_1",2,0,2);
   ROOT::TThreadedObject<TH1D> good_vrtx_2("good_vrtx_2","good_vrtx_2",2,0,2);

   ROOT::TThreadedObject<TH1D> z_pos_1("z_pos_1","z_pos_1",1000,400,900);
   ROOT::TThreadedObject<TH1D> z_pos_2("z_pos_2","z_pos_2",1000,400,900);
   ROOT::TThreadedObject<TH1D> z_pos_BestAll("z_pos_BestAll","z_pos_BestAll",1000,400,1000);
   ROOT::TThreadedObject<TH1D> z_pos_all("z_pos_all","z_pos_all",1000,400,1000);

 auto myFunction = [&](TTreeReader &myReader) {


     TTreeReaderValue<std::vector<MUonERecoOutputTrackAnalysis>> RVtracks(myReader, "ReconstructedTracks");
     TTreeReaderValue<MUonERecoOutputVertexAnalysis> vrtx(myReader, "BestVertex");
     TTreeReaderValue<std::vector<MUonERecoOutputVertexAnalysis>> RVall_vrtx(myReader, "ReconstructedVertices");
     TTreeReaderValue<std::vector<MUonERecoOutputHitAnalysis>> RVstubs(myReader, "ReconstructedHits");
//     TTreeReaderValue<std::vector<MUonETrackerStub>> tr_stubs(myReader,"TrackerStubs");

     TTreeReaderValue<Bool_t> trig_smi0(myReader, "TriggerSingleMuonInteraction0");
     TTreeReaderValue<Bool_t> trig_smi1(myReader, "TriggerSingleMuonInteraction1");

count_tr_1->Fill(-99); h_nvrtx_1->Fill(-99); vrtx_chi2_1->Fill(-99);
count_tr_2->Fill(-99); h_nvrtx_2->Fill(-99); vrtx_chi2_2->Fill(-99);
good_vrtx_1->Fill(-99);
good_vrtx_2->Fill(-99);
z_pos_1->Fill(-99);
z_pos_2->Fill(-99);
z_pos_BestAll->Fill(-99);

     while (myReader.Next()) {



 double chi=vrtx->chi2perDegreeOfFreedom();

/*if(type!="unbiased" or type!="unbiased_1"){
 if(chi!=0 and vrtx->stationIndex()==1){ h_nvrtx_1->Fill(1); vrtx_chi2_1->Fill(chi); if(chi<20) good_vrtx_1->Fill(1);}
 if(chi!=0 and vrtx->stationIndex()==2){ h_nvrtx_2->Fill(1); vrtx_chi2_2->Fill(chi); if(chi<20) good_vrtx_2->Fill(1);}
}
else{*/
if(*trig_smi0==1)count_tr_1->Fill(1);
if(*trig_smi1==1)count_tr_2->Fill(1);
// if(chi!=0 and vrtx->stationIndex()==1 and *trig_smi0==1){ h_nvrtx_1->Fill(1); vrtx_chi2_1->Fill(chi); if(chi<20 and vrtx->zPositionFit()>659.6 and vrtx->zPositionFit()<669.6) {good_vrtx_1->Fill(1); z_pos_1->Fill(vrtx->zPositionFit());} }
// if(chi!=0 and vrtx->stationIndex()==2 and *trig_smi1==1){ h_nvrtx_2->Fill(1); vrtx_chi2_2->Fill(chi); if(chi<20 and vrtx->zPositionFit()>775.7 and vrtx->zPositionFit()<785.7) {good_vrtx_2->Fill(1); z_pos_2->Fill(vrtx->zPositionFit());} }

 if(chi!=0) z_pos_BestAll->Fill(vrtx->zPositionFit());

 if(chi!=0 and vrtx->stationIndex()==1 and *trig_smi0==1){ h_nvrtx_1->Fill(1); vrtx_chi2_1->Fill(chi); if(chi<20) {good_vrtx_1->Fill(1); z_pos_1->Fill(vrtx->zPositionFit());} }
 if(chi!=0 and vrtx->stationIndex()==2 and *trig_smi1==1){ h_nvrtx_2->Fill(1); vrtx_chi2_2->Fill(chi); if(chi<20) {good_vrtx_2->Fill(1); z_pos_2->Fill(vrtx->zPositionFit());} }

//}

auto all_vrtx= *RVall_vrtx;
	for(auto&& av : all_vrtx){
		if(av.chi2perDegreeOfFreedom()<20) z_pos_all->Fill(av.zPositionFit());
	}

int sec0=0;
int sec1=0;
int sec2=0;


         auto tracks = *RVtracks;
         for (auto&& track : tracks) {
        if(track.sector()==0) sec0++;
        if(track.sector()==1) sec1++;
        if(track.sector()==2) sec2++;
        }

h_ntracks_s0->Fill(sec0);
h_ntracks_s1->Fill(sec1);
h_ntracks_s2->Fill(sec2);

int stub0=0;int stub1=0;int stub2=0;

         auto hits = *RVstubs;
         for (auto&& hit : hits) {
        if(hit.stationID()==0){stub0++;nstub_mod_s0->Fill(hit.moduleID());}
        if(hit.stationID()==1){stub1++;nstub_mod_s1->Fill(hit.moduleID());}
        if(hit.stationID()==2){stub2++;nstub_mod_s2->Fill(hit.moduleID());}
        }

h_nstubs_s0->Fill(stub0);
h_nstubs_s1->Fill(stub1);
h_nstubs_s2->Fill(stub2);


 } //end of general while
}; //end of ySlopefunction

  tp1.Process(myFunction);


  auto nstub_mod_s0M=nstub_mod_s0.Merge();
  auto nstub_mod_s1M=nstub_mod_s1.Merge();
  auto nstub_mod_s2M=nstub_mod_s2.Merge();
  auto h_nstubs_s0M=h_nstubs_s0.Merge();
  auto h_nstubs_s1M=h_nstubs_s1.Merge();
  auto h_nstubs_s2M=h_nstubs_s2.Merge();
  auto h_ntracks_s0M=h_ntracks_s0.Merge();
  auto h_ntracks_s1M=h_ntracks_s1.Merge();
  auto h_ntracks_s2M=h_ntracks_s2.Merge();
  auto h_nvrtx_1M=h_nvrtx_1.Merge();
  auto h_nvrtx_2M=h_nvrtx_2.Merge();
  auto vrtx_chi2_1M=vrtx_chi2_1.Merge();
  auto vrtx_chi2_2M=vrtx_chi2_2.Merge();
  auto count_tr_1M=count_tr_1.Merge();
  auto count_tr_2M=count_tr_2.Merge();
  auto good_vrtx_1M=good_vrtx_1.Merge();
  auto good_vrtx_2M=good_vrtx_2.Merge();

  auto z_pos_1M=z_pos_1.Merge();
  auto z_pos_2M=z_pos_2.Merge();
  auto z_pos_BestAllM=z_pos_BestAll.Merge();
  auto z_pos_allM=z_pos_all.Merge();

    std::ofstream out("report.txt"); // apre (o crea) il file in scrittura
    if (!out) {
        std::cerr << "Errore nell'aprire il file!" << std::endl;
        return 1;
    }

    out << "Entries " << cbmsim->GetEntries() << std::endl;
    out << "Events with trigger 0 " << count_tr_1M->Integral() << std::endl;
    out << "Events with trigger 1 " << count_tr_2M->Integral() << std::endl;
    out << "Events with reco vrtx in st01 " << h_nvrtx_1M->Integral()
        << " with chi<20 " << good_vrtx_1->Integral() << std::endl;
    out << "Events with reco vrtx in st12 " << h_nvrtx_2M->Integral()
        << " with chi<20 " << good_vrtx_2->Integral() << std::endl;

    out.close();
/*
cout <<"Entries " << cbmsim->GetEntries() << endl;
cout <<"Events with trigger 0 " << count_tr_1M->Integral() << endl;
cout <<"Events with trigger 1 " << count_tr_2M->Integral() << endl;
cout <<"Events with reco vrtx in st01 " << h_nvrtx_1M->Integral() << " with chi<20 " << good_vrtx_1->Integral() << endl;
cout <<"Events with reco vrtx in st12 " << h_nvrtx_2M->Integral() << " with chi<20 " << good_vrtx_2->Integral() << endl;
*/
if(save=="save"){
                TFile *outfile = new TFile(Form("/home/espedica/clement_fair_install/instFairRoot/share/MUonE/macros/counting/%s_%s_hits%i.root",run.c_str(),type.c_str(),nhits), "RECREATE");

nstub_mod_s0M->Write();
nstub_mod_s1M->Write();
nstub_mod_s2M->Write();
h_nstubs_s0M->Write();
h_nstubs_s1M->Write();
h_nstubs_s2M->Write();
h_ntracks_s0M->Write();
h_ntracks_s1M->Write();
h_ntracks_s2M->Write();
vrtx_chi2_1M->Write();
vrtx_chi2_2M->Write();
z_pos_1M->Write();
z_pos_2M->Write();
z_pos_BestAllM->Write();
z_pos_allM->Write();
                outfile->Close();
}
return 0;
}
