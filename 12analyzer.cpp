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


void RealDataAnalyzer(){

 	TFile *inputfile = new TFile("/mnt/raid10/DATA/espedica/fairmu/dataReconstruction_entire.root");
//dataReconstruction_merged_3232_3233_small4M_outlier.root");
//dataReconstruction_merged_3238_3239_0hit_1.root");
//dataReconstruction_merged_3232_3233_small4M_outlier.root");
        TTree* cbmsim = (TTree*) inputfile->Get("cbmsim");

        MUonERecoOutput *ReconstructionOutput = 0;


   cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);


        std::vector<TH1D*> residual1(12);
for(int m=0; m<12; m++){
	string name="residual1_"+to_string(m);
                string title="Residual station1 of module "+to_string(m);
if(m==2 or m==3 or m==8 or m==9) residual1.at(m)=new TH1D(name.c_str(),title.c_str(),600,-0.06,0.06);
else residual1.at(m)=new TH1D(name.c_str(),title.c_str(),600,-0.06,0.06);
			}


        std::vector<TH1D*> localX1(12);
for(int m=0; m<12; m++){
	        string name="local1_"+to_string(m);
                string title="LocalX station1 of module "+to_string(m);
		localX1.at(m)=new TH1D(name.c_str(),title.c_str(),100,-5.,5.);}

TH1D *h_nstubs_1=new TH1D("h_nstubs_1", "number of stubs for a track with 12 modules", 13,0,13);


// X or Y position on the track at a given Z
                auto pos_on_track = [](double q, double m, double z){return (q + m*z);};

for(Long64_t i = 0; i < 1000000; i++) {//cbmsim->GetEntries(); i++) {
		cbmsim->GetEntry(i);
		if(i%1000 == 0) cout<<"Entry "<<i<<endl;

vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
vector<MUonERecoOutputVertex> vrtx = ReconstructionOutput->reconstructedVertices();
//vector<MUonERecoOutputAdaptiveFitterVertex> vrtx = ReconstructionOutput->adaptiveFitterVertices();


    int sec0=0;

double pos1[12];

    for(int j=0; j<tracks.size();j++)
    {if(tracks.at(j).sector()==0)sec0++;}

    for(int j=0; j<tracks.size();j++)
    {

     	if(tracks.at(j).sector()==0 and sec0==1){

std::vector<MUonERecoOutputHit> hits_=tracks.at(j).hits();
h_nstubs_1->Fill(hits_.size());
  for(int h=0;h<hits_.size();h++)
  {
        if(hits_.at(h).stationID()==0)residual1.at(hits_.at(h).moduleID())->Fill(hits_.at(h).perpendicularResiduum());
   	if(hits_.at(h).stationID()==1)residual1.at(hits_.at(h).moduleID()+6)->Fill(hits_.at(h).perpendicularResiduum());

  }
 }
}

}

TCanvas n3("n3","n3",700,700);
n3.Divide(2,6);
for(int m=0; m<12; m++){
n3.cd(m+1);
residual1[m]->Draw();
gPad->SetLogy();}
n3.SaveAs("res1.pdf");

TCanvas n4("n4","n4",700,700);
h_nstubs_1->Draw();
gPad->SetLogy();
n4.SaveAs("nstubs_12.pdf");

}
