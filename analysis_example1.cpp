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

void analysis_example1(){

//input file name
 	TFile *inputfile = new TFile("example1.root");

//TTree name
TTree* cbmsim = (TTree*) inputfile->Get("cbmsim");

//branches present in the TTrree
TClonesArray *MCTrack = 0;
TClonesArray *SignalTracks = 0;
TClonesArray *TrackerStripDigis = 0;
TClonesArray *TrackerPoints = 0;
TClonesArray *CalorimeterPoints = 0;
TClonesArray *TrackerStubs = 0;
MUonECalorimeterDigiDeposit *CalorimeterDigiDeposit = 0;
MuE::Event *MesmerEvent = 0;
MUonERecoOutput *ReconstructionOutput = 0;

cbmsim->SetBranchAddress("MCTrack", &MCTrack);
cbmsim->SetBranchAddress("SignalTracks", &SignalTracks);
cbmsim->SetBranchAddress("TrackerPoints", &TrackerPoints);
cbmsim->SetBranchAddress("CalorimeterPoints", &CalorimeterPoints);
cbmsim->SetBranchAddress("TrackerStripDigis", &TrackerStripDigis);
cbmsim->SetBranchAddress("TrackerStubs", &TrackerStubs);
cbmsim->SetBranchAddress("CalorimeterDigiDeposit", &CalorimeterDigiDeposit);
cbmsim->SetBranchAddress("MesmerEvent", &MesmerEvent);
cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);


TH1D *theta_e=new TH1D("theta_e", "Reconstructed electron scattering angles",100,0.,0.032);
TH1D *theta_mu=new TH1D("theta_mu", "Reconstructed muon scattering angles",100,0.,0.005);


//Loop on the events of the TTree
for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
		cbmsim->GetEntry(i);
		if(i%1000 == 0) cout<<"Entry "<<i<<endl;

//Let's look inside the MC tracks container, are there the electron and muon from elastic interaction? to see the ROOT ID for interactions go to https://root.cern/doc/v610/TMCProcess_8h.html, elastic interaction is interactionID==45 while pair production is interactionID==5. While to see the PDG particle ID https://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf (electron is 11 and muon is 13)

int e=0; int mu=0;
int code_e = -99; int code_mu=-99;

	for(int n = 0; n < MCTrack->GetEntries(); n++) {

	 const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));

if(MCTr->interactionID()==0 and MCTr->pdgCode()==13){
	   cout << "I'm the incoming muon " << endl; }

 if(MCTr->interactionID()==45 and MCTr->pdgCode()==11) {
 e=1; 
 code_e=n;  //important for MC truth => association with reconstructed track
 cout << "I'm the MC electron" << endl; }

 if(MCTr->interactionID()==45 and MCTr->pdgCode()==13) {
 mu=1; 
 code_mu=n; //important for MC truth => association with reconstructed track
 cout << "I'm the MC positive muon" << endl; }
	}

//if electron and muon have been generated, is there any reconstructed track and vertex?
if(e==1 and mu==1){

//To know what is inside the reconstruction output (name of the objects, variable, attributes etc.) go to https://gitlab.cern.ch/mgoncerz/fairmuone/-/blob/master/MUonE/MUonEReconstruction/MUonERecoOutput.h?ref_type=heads

vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
vector<MUonERecoOutputVertex> vrtx = ReconstructionOutput->reconstructedVertices();


//How to surf inside ECAL output
for(int c=0; c<25; c++){
double column = cbmsim->GetLeaf("Crystals.crystalColumnID")->GetValue(c);
double row = cbmsim->GetLeaf("Crystals.crystalRowID")->GetValue(c);
double value = cbmsim->GetLeaf("Crystals.EnergyDeposit")->GetValue(c);
cout << "Energy dep. " << value << " in crystal (row,column): -> (" <<row << ", " << column <<")"<< endl;
}


//loop over vertex container, but you can choose the best vertex just defining the one below, so you don't need to loop over all the vrtx reconstructed. Best vrtx is the one with best chi2.
MUonERecoOutputVertex best_vrtx = ReconstructionOutput->bestVertex();

for(int j=0; j<vrtx.size();j++)
{
if(vrtx.at(j).stationIndex()==1)
{
	cout << "There is a reconstructed vrtx with chi2 = " << vrtx.at(j).chi2perDegreeOfFreedom() << endl;
	}
}

cout << "The best vrtx has a chi2 = " << best_vrtx.chi2perDegreeOfFreedom() << endl;


std::vector<MUonERecoOutputHit> stubs=ReconstructionOutput->ReconstructedHits();

//reconstructedHits();
//vector<MUonETrackDeposit> hits_



//loop over reconstructed tracks

cout << "There are " << tracks.size() << " reconstructed tracks" << endl;

TVector3 p_muin,p_mu,p_e;
double thmu_rec=0; double the_rec=0;
for(int j=0; j<tracks.size();j++)
{
//you can ask the processID of the tracksa to see from which interaction comes. If ==0 it is the incoming muon. The sector is the station where the particle is reconstructed
if(tracks.at(j).processIDofLinkedTrack()==0 and tracks.at(j).sector()==0){
		cout << "This is the reconstructed incoming muon in station " << tracks.at(j).sector()  << endl;
	        double th_inx=tracks.at(j).xSlope();
        	double th_iny=tracks.at(j).ySlope();
        	p_muin.SetXYZ(th_inx,th_iny,1.0);
	        p_muin=p_muin.Unit();
	}
//are there any particle from elastic scattering?


if(tracks.at(j).processIDofLinkedTrack()==45 and tracks.at(j).sector()==1){
		cout << "This is one of the elastic particle reconstructed in station " << tracks.at(j).sector()  << endl;
//Is it the muon or the electron? You can use linkedTrackID() which associate the reconstructed track to one of the MCTracks container. The variable represent the number of that particle in the container. So if the outgoing muon was the second particle inside MCTracks container, then linkedTrackID==2. So if you want to use MC truth, you have to save those number when you loop on MCTracks.
//We also look at the reco scattering angle. NB it needs to be wrt the incoming muon direction
if(code_mu==tracks.at(j).linkedTrackID()){cout << "I'm the reconstructed muon" << endl;
					double thmuX_rec=tracks.at(j).xSlope(); double thmuY_rec=tracks.at(j).ySlope();
					p_mu.SetXYZ(thmuX_rec,thmuY_rec,1.0);
					p_mu=p_mu.Unit();
					thmu_rec=p_mu.Angle(p_muin);}
if(code_e==tracks.at(j).linkedTrackID()){cout << "I'm the reconstructed electron" << endl;
					double theX_rec=tracks.at(j).xSlope(); double theY_rec=tracks.at(j).ySlope();
					p_e.SetXYZ(theX_rec,theY_rec,1.0); p_e=p_e.Unit();
					the_rec=p_e.Angle(p_muin);}
	}
}


theta_mu->Fill(thmu_rec,MesmerEvent->wgt_full);
theta_e->Fill(the_rec,MesmerEvent->wgt_full);

 }//close if(e==1 and mu==1)
}//close for loop

TCanvas c("c","c",700,700);
c.Divide(1,2);
c.cd(1);
theta_mu->Draw("hist");
c.cd(2);
theta_e->Draw("hist");
c.SaveAs("scattering_angles.pdf");
}//end of file
