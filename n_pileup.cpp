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

using namespace std;


void n_pileup(){

TChain * cbmsim = new TChain("cbmsim");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/minbias_1M_2.root");
//1M_2.root");//10k_new.root");

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

double n_1=0.;
double n_2=0.;
double n_3=0.;
double n_4=0.;


for(Long64_t i = 0; i < 500000; i++) {

 cbmsim->GetEntry(i);
 if(i%100 == 0) cout<<"Entry "<<i<<endl;
int mu_gen=0;
        for(int n = 0; n < MCTrack->GetEntries(); n++) {
         const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
         if(MCTr->interactionID()==0 and MCTr->pdgCode()==13){mu_gen++;}
	}
if(mu_gen==1)n_1++;
if(mu_gen==2)n_2++;
if(mu_gen==3)n_3++;
if(mu_gen==4)n_4++;

	}//end for

cout << "Pileup n=1: " << n_1 << " events on " << cbmsim->GetEntries() << " entries -> " << n_1/cbmsim->GetEntries() * 100 << "%" <<endl;
cout << "Pileup n=2: " << n_2 << " events on " << cbmsim->GetEntries() << " entries -> " << n_2/cbmsim->GetEntries() * 100 << "%" <<endl;
cout << "Pileup n=3: " << n_3 << " events on " << cbmsim->GetEntries() << " entries -> " << n_3/cbmsim->GetEntries() * 100 << "%" <<endl;
cout << "Pileup n=4: " << n_4 << " events on " << cbmsim->GetEntries() << " entries -> " << n_4/cbmsim->GetEntries() * 100 << "%" <<endl;

}//end file
