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

void analysis_example2_mu(){

//input file name
 	TFile *inputfile = new TFile("example2_muonfilter.root");//TRMesmer_box_100k.root");

//TTree name
TTree* cbmsim = (TTree*) inputfile->Get("cbmsim");

//branches present in the TTrree
TClonesArray *MCTrack = 0;
TClonesArray *SignalTracks = 0;
TClonesArray *TrackerStripDigis = 0;
TClonesArray *TrackerPoints = 0;
TClonesArray *TrackerStubs = 0;
MuE::Event *MesmerEvent = 0;
MUonERecoOutput *ReconstructionOutput = 0;

cbmsim->SetBranchAddress("MCTrack", &MCTrack);
cbmsim->SetBranchAddress("SignalTracks", &SignalTracks);
cbmsim->SetBranchAddress("TrackerPoints", &TrackerPoints);
cbmsim->SetBranchAddress("TrackerStripDigis", &TrackerStripDigis);
cbmsim->SetBranchAddress("TrackerStubs", &TrackerStubs);
cbmsim->SetBranchAddress("MesmerEvent", &MesmerEvent);
cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);


//Loop on the events of the TTree
for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
		cbmsim->GetEntry(i);
		if(i%1000 == 0) cout<<"Entry "<<i<<endl;

//Let's look inside the MC tracks container, are there the electron and muon from elastic interaction? to see the ROOT ID for interactions go to https://root.cern/doc/v610/TMCProcess_8h.html, elastic interaction is interactionID==45 while pair production is interactionID==5. While to see the PDG particle ID https://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf (electron is 11 and muon is 13)

TVector3 pmuin,pe,pmu;
int e=0; int mu=0;
int code_e = -99; int code_mu=-99;

	for(int n = 0; n < MCTrack->GetEntries(); n++) {

	 const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));

 if(MCTr->interactionID()==45 and MCTr->pdgCode()==11) {
 e=1;
 code_e=n;}  //important for MC truth => association with reconstructed track}

 if(MCTr->interactionID()==45 and MCTr->pdgCode()==13) {
 mu=1;
 code_mu=n;} //important for MC truth => association with reconstructed track}
	}

//if electron and muon have been generated, is there any reconstructed track and vertex?
if(e==1 and mu==1){
vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();


//loop over reconstructed tracks

for(int j=0; j<tracks.size();j++)
{
 //are there any particle in the statione behind calo?
 if(tracks.at(j).sector()==2){
	cout << "There is a particle in the muon filter with linkID processID " << tracks.at(j).processIDofLinkedTrack() << endl;
	for(int n = 0; n < MCTrack->GetEntries(); n++) {
         const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
	 if(n==tracks.at(j).linkedTrackID()) {cout << "pdg Code is " << MCTr->pdgCode() << " at Z " << MCTr->startZ() << endl;}
	 }
	}

}

 }//close if(e==1 and mu==1)
}//close for loop

}//end of file
