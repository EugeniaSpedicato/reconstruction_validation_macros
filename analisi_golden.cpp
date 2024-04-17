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


for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
		cbmsim->GetEntry(i);
		if(i%1000 == 0) cout<<"Entry "<<i<<endl;


int nStubs_0 = 0;
int nStubs_1 = 0;


vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();

int n_0=0;
int n_1=0;

//number of reco tracks in sector0 and sector1 (sector==station)
for(int t=0; t<tracks.size(); t++){
if(tracks.at(t).sector()==0)n_0++;
if(tracks.at(t).sector()==1)n_1++;
}


std::vector<MUonERecoOutputHit> stubs=ReconstructionOutput->reconstructedHits();
// number of stubs per sector
  for (auto & hit : stubs) {
      if (hit.stationID()==0) nStubs_0++;
      if (hit.stationID()==1) nStubs_1++;
	}

TVector3 p_muin,p1;
double th1=-99.;
int link1=-99;


//golden mu sector0
if(n_0==1 and nStubs_0==6){

	for(int t=0; t<tracks.size(); t++){
	//direction mu_in
	if(tracks.at(t).sector()==0){p_muin.SetXYZ(tracks.at(t).xSlope(),tracks.at(t).ySlope(),1.0); p_muin=p_muin.Unit(); double theta_in=p_muin.Angle()}
	}

  //what happens in station1
  cout << "Number of tracks Station 1 when golden muon in Station 0: " << n_1 << " with number of stubs: " << nStubs_1 << endl;
   }



 } //end of general for

}//end of file
