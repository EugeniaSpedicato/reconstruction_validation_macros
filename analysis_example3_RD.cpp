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
#include <yaml-cpp/yaml.h>
#include "geometry.h"

using namespace std;


void analysis_example3_RD(){

 	TFile *inputfile = new TFile("example3_RD.root");

       TTree* cbmsim = (TTree*) inputfile->Get("cbmsim");

        MUonERecoOutput *ReconstructionOutput = 0;


   cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);

    TH1D *h_mult0=new TH1D("h_mult0","Multiplicity of track in station0",10,0,10);
    TH1D *h_mult1=new TH1D("h_mult1","Multiplicity of track in station1",10,0,10);



for(Long64_t i = 0; i < 1000000; i++) {//cbmsim->GetEntries(); i++) {
		cbmsim->GetEntry(i);
		if(i%1000 == 0) cout<<"Entry "<<i<<endl;

vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();


    int sec0=0; int sec1=0;


    for(int j=0; j<tracks.size();j++)
    {
	std::vector<MUonERecoOutputHit> hits_=tracks.at(j).hits();
        if(tracks.at(j).sector()==0)sec0++;
        if(tracks.at(j).sector()==1)sec1++;
	}

h_mult0->Fill(sec0);
h_mult1->Fill(sec1);


}


    TCanvas n("n","n",700,700);
    h_mult0->Draw();
    h_mult1->SetLineColor(kRed);
    h_mult1->Draw("same");
    n.SaveAs("info.pdf");

}


