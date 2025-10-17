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
#include "ROOT/TThreadedObject.hxx"


using namespace std;


int extra_hit(int nhits,string run, string type, string save){

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
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run11/single_muon_interaction_0/4files_commit62871d8b.root");
}
else if(run=="run11" and nhits==0 and type=="single_mu_int_1"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run11/single_muon_interaction_1/4files_commit62871d8b.root");
}
else if(run=="run11_scint" and nhits==0 and type=="single_mu_int_0"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run11/single_muon_interaction_0/4files_commit62871d8b_scintillator.root");
}
else if(run=="run16_scint" and nhits==0 and type=="single_mu_int_0"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run16/single_muon_interaction_0/4files_commit62871d8b_scintillator.root");
}
else if(run=="run11_scint" and nhits==0 and type=="single_mu_int_1"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run11/single_muon_interaction_1/1files_commit62871d8b_scintillator.root");
}
else if(run=="run16_scint" and nhits==0 and type=="single_mu_int_1"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run16/single_muon_interaction_1/2files_commit62871d8b_scintillator.root");
}


cout << "cbmsim->GetEntries() " << cbmsim->GetEntries() <<endl;

ROOT::TTreeProcessorMT tp1(*cbmsim,nthreads);

      MUonERecoOutputAnalysis *ReconstructionOutput = 0;



TH1::SetDefaultSumw2(kTRUE);
   const Int_t NBINS = 6;
   Double_t edges[NBINS + 1] = {0.0, 0.005, 0.010, 0.015, 0.020, 0.025, 0.032};

   const Int_t NBINS2 = 12;
   Double_t edges2[NBINS2 + 1] = {0.,0.0002,0.0004,0.0006,0.0008,0.001,0.00125,0.0015,0.00175,0.002,0.003,0.004,0.005};



   ROOT::TThreadedObject<TH1D> h_mod_extra("h_mod_extra","module_extra_hit",7,0,7);

   ROOT::TThreadedObject<TH1D> h_dist_extra_0("h_dist_extra_0","Distance extra - track position in mod 0",90,-0.3,0.3);
   ROOT::TThreadedObject<TH1D> h_pos_extra_0("h_pos_extra_0","Position extra in mod 0",1046, -5, 5);

   ROOT::TThreadedObject<TH1D> h_dist_extra_1("h_dist_extra_1","Distance extra - track position in mod 1",90,-0.3,0.3);
   ROOT::TThreadedObject<TH1D> h_pos_extra_1("h_pos_extra_1","Position extra in mod 1",1046, -5, 5);

   ROOT::TThreadedObject<TH1D> h_dist_extra_2("h_dist_extra_2","Distance extra - track position in mod 2",90,-0.3,0.3);
   ROOT::TThreadedObject<TH1D> h_pos_extra_2("h_pos_extra_2","Position extra in mod 2",1046, -5, 5);

   ROOT::TThreadedObject<TH1D> h_dist_extra_3("h_dist_extra_3","Distance extra - track position in mod 3",90,-0.3,0.3);
   ROOT::TThreadedObject<TH1D> h_pos_extra_3("h_pos_extra_3","Position extra in mod 3",1046, -5, 5);

   ROOT::TThreadedObject<TH1D> h_dist_extra_4("h_dist_extra_4","Distance extra - track position in mod 4",90,-0.3,0.3);
   ROOT::TThreadedObject<TH1D> h_pos_extra_4("h_pos_extra_4","Position extra in mod 4",1046, -5, 5);

   ROOT::TThreadedObject<TH1D> h_dist_extra_5("h_dist_extra_5","Distance extra - track position in mod 5",90,-0.3,0.3);
   ROOT::TThreadedObject<TH1D> h_pos_extra_5("h_pos_extra_5","Position extra in mod 5",1046, -5, 5);

 auto myFunction = [&](TTreeReader &myReader) {


     TTreeReaderValue<std::vector<MUonERecoOutputTrackAnalysis>> RVtracks(myReader, "ReconstructedTracks");
     TTreeReaderValue<MUonERecoOutputVertexAnalysis> vrtx(myReader, "BestVertex");
     TTreeReaderValue<std::vector<MUonERecoOutputVertexAnalysis>> RVall_vrtx(myReader, "ReconstructedVertices");
     TTreeReaderValue<std::vector<MUonERecoOutputHitAnalysis>> RVstubs(myReader, "ReconstructedHits");
//     TTreeReaderValue<std::vector<MUonETrackerStub>> tr_stubs(myReader,"TrackerStubs");

     TTreeReaderValue<Bool_t> trig_smi0(myReader, "TriggerSingleMuonInteraction0");
     TTreeReaderValue<Bool_t> trig_smi1(myReader, "TriggerSingleMuonInteraction1");

     while (myReader.Next()) {



 double chi=vrtx->chi2perDegreeOfFreedom();


int sec0=0;
int sec1=0;
int sec2=0;


         auto tracks = *RVtracks;
         for (auto&& track : tracks) {
        if(track.sector()==0) sec0++;
        if(track.sector()==1) sec1++;
        if(track.sector()==2) sec2++;
        }

int stub0=0;int stub1=0;int stub2=0;

         auto hits = *RVstubs;
         for (auto&& hit : hits) {
        if(hit.stationID()==0){stub0++;}
        if(hit.stationID()==1){stub1++;}
        if(hit.stationID()==2){stub2++;}
        }

MUonERecoOutputHitAnalysis extra_hit;
if(sec0==1 and stub0==7){
         for (auto&& track : tracks) {
		if(track.sector()==0 and track.hits().size()==6){
	         for (auto&& hit : hits) {

		         if(hit.stationID()==0){ if(std::find_if(track.hits().begin(),track.hits().end(), [&](auto p){return p.position()==hit.position();})==track.hits().end() ) extra_hit=hit;}

		for (auto&& thit : track.hits()) {
			if(thit.moduleID()==extra_hit.moduleID()){
			if(thit.moduleID()==0){h_dist_extra_0->Fill(extra_hit.position()-thit.position()); h_pos_extra_0->Fill(extra_hit.position());}
			if(thit.moduleID()==1){h_dist_extra_1->Fill(extra_hit.position()-thit.position()); h_pos_extra_1->Fill(extra_hit.position());}
			if(thit.moduleID()==2){h_dist_extra_2->Fill(extra_hit.position()-thit.position()); h_pos_extra_2->Fill(extra_hit.position());}
			if(thit.moduleID()==3){h_dist_extra_3->Fill(extra_hit.position()-thit.position()); h_pos_extra_3->Fill(extra_hit.position());}
			if(thit.moduleID()==4){h_dist_extra_4->Fill(extra_hit.position()-thit.position()); h_pos_extra_4->Fill(extra_hit.position());}
			if(thit.moduleID()==5){h_dist_extra_5->Fill(extra_hit.position()-thit.position()); h_pos_extra_5->Fill(extra_hit.position());}
			}
		 }
		}
	}
 }
}

h_mod_extra->Fill(extra_hit.moduleID());




 } //end of general while
}; //end of ySlopefunction

  tp1.Process(myFunction);




if(save=="save"){
                TFile *outfile = new TFile(Form("/home/espedica/clement_fair_install/instFairRoot/share/MUonE/macros/counting/extrahit_%s_%s_hits%i.root",run.c_str(),type.c_str(),nhits), "RECREATE");

std::shared_ptr<TH1D> h_mod_extraM = h_mod_extra.Merge();
h_mod_extraM->Write();

std::shared_ptr<TH1D> h_dist_extra_0M = h_dist_extra_0.Merge();
h_dist_extra_0M->Write();
std::shared_ptr<TH1D> h_pos_extra_0M = h_pos_extra_0.Merge();
h_pos_extra_0M->Write();

std::shared_ptr<TH1D> h_dist_extra_1M = h_dist_extra_1.Merge();
h_dist_extra_1M->Write();
std::shared_ptr<TH1D> h_pos_extra_1M = h_pos_extra_1.Merge();
h_pos_extra_1M->Write();

std::shared_ptr<TH1D> h_dist_extra_2M = h_dist_extra_2.Merge();
h_dist_extra_2M->Write();
std::shared_ptr<TH1D> h_pos_extra_2M = h_pos_extra_2.Merge();
h_pos_extra_2M->Write();

std::shared_ptr<TH1D> h_dist_extra_3M = h_dist_extra_3.Merge();
h_dist_extra_3M->Write();
std::shared_ptr<TH1D> h_pos_extra_3M = h_pos_extra_3.Merge();
h_pos_extra_3M->Write();

std::shared_ptr<TH1D> h_dist_extra_4M = h_dist_extra_4.Merge();
h_dist_extra_4M->Write();
std::shared_ptr<TH1D> h_pos_extra_4M = h_pos_extra_4.Merge();
h_pos_extra_4M->Write();

std::shared_ptr<TH1D> h_dist_extra_5M = h_dist_extra_5.Merge();
h_dist_extra_5M->Write();
std::shared_ptr<TH1D> h_pos_extra_5M = h_pos_extra_5.Merge();
h_pos_extra_5M->Write();

TCanvas a("a","a",1000,1500);
a.Divide(2,3);
a.cd(1);
h_dist_extra_0M->Draw("hist");
a.cd(2);
h_dist_extra_1M->Draw("hist");
a.cd(3);
h_dist_extra_2M->Draw("hist");
a.cd(4);
h_dist_extra_3M->Draw("hist");
a.cd(5);
h_dist_extra_4M->Draw("hist");
a.cd(6);
h_dist_extra_5M->Draw("hist");
//a.SaveAs(Form("counting/dist_extra_%s_%s_hits%i.pdf",run.c_str(),type.c_str(),nhits));
a.Write();

TCanvas b("b","b",1000,1500);
b.Divide(2,3);
b.cd(1);
h_pos_extra_0M->Draw("hist");
b.cd(2);
h_pos_extra_1M->Draw("hist");
b.cd(3);
h_pos_extra_2M->Draw("hist");
b.cd(4);
h_pos_extra_3M->Draw("hist");
b.cd(5);
h_pos_extra_4M->Draw("hist");
b.cd(6);
h_pos_extra_5M->Draw("hist");
//b.SaveAs(Form("counting/pos_extra_%s_%s_hits%i.pdf",run.c_str(),type.c_str(),nhits));
b.Write();
                outfile->Close();
}
return 0;
}
