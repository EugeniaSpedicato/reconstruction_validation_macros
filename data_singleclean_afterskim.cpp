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

void data_singleclean_afterskim(){

TChain * cbmsim = new TChain("cbmsim");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_run6_singlemu_1hitshared_nochi2_ALL.root");
// /mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_run6_singlemu_1hitshared_nochi2.root");

        TClonesArray *MCTrack = 0;
        MUonERecoOutput *ReconstructionOutput = 0;

        cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);

TTreeReader myReader(cbmsim);

  TH1D *h_reco_tracks = new TH1D("h_reco_tracks","number of reco tracks", 20,0,20);
  TH1D *h_reco_tracks_no = new TH1D("h_reco_tracks_no","number of reco tracks no vrtx", 20,0,20);

  TH1D *h_reco_tracks_0 = new TH1D("h_reco_tracks_0","number of reco tracks station0", 20,0,20);
  TH1D *h_reco_tracks_1 = new TH1D("h_reco_tracks_1","number of reco tracks station1", 20,0,20);

  TH1D *h_reco_tracks_0_no = new TH1D("h_reco_tracks_0_no","number of reco tracks station0 no vrtx", 20,0,20);
  TH1D *h_reco_tracks_1_no = new TH1D("h_reco_tracks_1_no","number of reco tracks station1 no vrtx", 20,0,20);

  TH1D *h_th1=new TH1D("h_th1","Angle single particle when not correctly reconstructed",200,0.,0.001);
  TH1D *h_linkedID=new TH1D("h_linkedID","Linked track ID when not correctly reconstructed",20,0,20);

  TH2D *h_reco_tracks_2d = new TH2D("h_reco_tracks_2d","number of reco tracks station0 VS station1", 20,0,20, 20,0,20);
  TH2D *h_reco_tracks_2d_no = new TH2D("h_reco_tracks_2d_no","number of reco tracks station0 VS station1 no vrtx", 20,0,20, 20,0,20);

  TH1D *h_nStubs_0_1PUmu_no = new TH1D("h_nStubs_0_1PUmu_no","number stubs station 0 no vrtx",30,0,30);
  TH1D *h_nStubs_1_1PUmu_no = new TH1D("h_nStubs_1_1PUmu_no","number stubs station 1 no vrtx",30,0,30);

  TH1D *h_nStubs_0_1PUmu = new TH1D("h_nStubs_0_1PUmu","number stubs station 0",30,0,30);
  TH1D *h_nStubs_1_1PUmu = new TH1D("h_nStubs_1_1PUmu","number stubs station 1",30,0,30);

  TH1D *h_Xpos = new TH1D("h_Xpos","Stubs'X position station0 when number stubs<7",100,-10.,10.);
  TH1D *h_Xpos_7 = new TH1D("h_Xpos_7","Stubs'X position station0 when number stubs==7",100,-10.,10.);
  TH1D *h_Ypos = new TH1D("h_Ypos","Stubs'Y position station0 when number stubs<7",100,-10.,10.);
  TH1D *h_Ypos_7 = new TH1D("h_Ypos_7","Stubs'Y position station0 when number stubs==7",100,-10.,10.);

for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
		cbmsim->GetEntry(i);
		if(i%1000 == 0) cout<<"Entry "<<i<<endl;


int nStubs_0 = 0;
int nStubs_1 = 0;


vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();

int n_0=0;
int n_1=0;


for(int t=0; t<tracks.size(); t++){
if(tracks.at(t).sector()==0)n_0++;
if(tracks.at(t).sector()==1)n_1++;
}


std::vector<MUonERecoOutputHit> stubs=ReconstructionOutput->reconstructedHits();

  for (auto & hit : stubs) {
      if (hit.stationID()==0) nStubs_0++;
      if (hit.stationID()==1) nStubs_1++;
}
  for (auto & hit : stubs) {
      if (hit.stationID()==0){
                                if(hit.moduleID()==0 || hit.moduleID()==4){if(nStubs_0==7){h_Xpos_7->Fill(hit.positionPerpendicular());}else{h_Xpos->Fill(hit.positionPerpendicular());}}
                                if(hit.moduleID()==1 || hit.moduleID()==5){if(nStubs_0==7){h_Ypos_7->Fill(hit.positionPerpendicular());}else{h_Ypos->Fill(hit.positionPerpendicular());}}
				}
		}
TVector3 p_muin,p1;
double th1=-99.;
int link1=-99;
for(int t=0; t<tracks.size(); t++){
if(n_0==1 and n_1==1 and ReconstructionOutput->reconstructedVertices().size()==0){
        if(tracks.at(t).sector()==0){p_muin.SetXYZ(tracks.at(t).xSlope(),tracks.at(t).ySlope(),1.0); p_muin=p_muin.Unit();}
        else{p1.SetXYZ(tracks.at(t).xSlope(),tracks.at(t).ySlope(),1.0); p1=p1.Unit();th1=p1.Angle(p_muin); link1=tracks.at(t).linkedTrackID();}
  }
}

if(ReconstructionOutput->reconstructedVertices().size()!=0){h_reco_tracks_0->Fill(n_0);h_reco_tracks_1->Fill(n_1);h_reco_tracks->Fill(tracks.size());
                                                                h_nStubs_0_1PUmu->Fill(nStubs_0);h_nStubs_1_1PUmu->Fill(nStubs_1);h_reco_tracks_2d->Fill(n_0,n_1);}
else{h_reco_tracks_0_no->Fill(n_0);h_reco_tracks_1_no->Fill(n_1);h_th1->Fill(th1); h_linkedID->Fill(link1);h_reco_tracks_no->Fill(tracks.size());
                                                                h_nStubs_0_1PUmu_no->Fill(nStubs_0);h_nStubs_1_1PUmu_no->Fill(nStubs_1);h_reco_tracks_2d_no->Fill(n_0,n_1);}

} //end of general for


TCanvas a("a","a",1000,700);
a.Divide(3,2);
a.cd(1);
h_reco_tracks_0_no->SetLineColor(kRed);
h_reco_tracks_0_no->Draw("hist");
h_reco_tracks_0->Draw("hist same");
a.cd(2);
h_reco_tracks_1_no->SetLineColor(kRed);
h_reco_tracks_1_no->Draw("hist");
h_reco_tracks_1->Draw("hist same");
a.cd(3);
h_linkedID->Draw("hist");
gPad->SetLogy();
a.cd(4);
h_th1->Draw("hist");
a.cd(5);
h_reco_tracks_no->SetLineColor(kRed);
h_reco_tracks_no->Draw("hist");
h_reco_tracks->Draw("hist same");
a.SaveAs("realdata_singlemuclean/h_reco_tracks_1PUmu.pdf");

TCanvas b("b","b",1000,700);
b.Divide(2,2);
b.cd(1);
h_reco_tracks_2d->SetXTitle("n_reco_tracks station0");
h_reco_tracks_2d->SetYTitle("n_reco_tracks station1");
h_reco_tracks_2d->Draw("COLZ");
b.cd(2);
h_reco_tracks_2d_no->SetXTitle("n_reco_tracks station0");
h_reco_tracks_2d_no->SetYTitle("n_reco_tracks station1");
h_reco_tracks_2d_no->Draw("COLZ");
b.cd(3);
h_nStubs_0_1PUmu_no->SetLineColor(kRed);
h_nStubs_0_1PUmu_no->Draw("hist");
h_nStubs_0_1PUmu->Draw("hist same");
b.cd(4);
h_nStubs_1_1PUmu_no->SetLineColor(kRed);
h_nStubs_1_1PUmu_no->Draw("hist");
h_nStubs_1_1PUmu->Draw("hist same");
b.SaveAs("realdata_singlemuclean/stubs_tracks_1PUmu.pdf");

TCanvas p("p","p");
p.Divide(1,2);
p.cd(1);
h_Xpos->Draw("hist");
h_Xpos_7->SetLineColor(kRed);
h_Xpos_7->Draw("hist same");
p.cd(2);
h_Ypos->Draw("hist");
h_Ypos_7->SetLineColor(kRed);
h_Ypos_7->Draw("hist same");
p.SaveAs("realdata_singlemuclean/position_st0.pdf");
}
