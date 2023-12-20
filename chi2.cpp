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
#include "TMatrix.h"
#include "TSystemDirectory.h"
#include <TStyle.h>
#include <yaml-cpp/yaml.h>
#include "geometry.h"

using namespace std;


void RealDataAnalyzer(){

 	TFile *inputfile = new TFile("/mnt/raid10/DATA/espedica/fairmu/dataReconstruction_3234-3235_separate_alin.root");
///mnt/raid10/DATA/espedica/fairmu/dataReconstruction_3234-3235_try_shift.root");
//TRMesmer_noInt_10k.root");
//dataReconstruction_3234-3235_new_straight.root");
//new_straight.root");
//separate_alin.root");
//dataReconstruction_merged_3238_3239_0hit_1.root");
//dataReconstruction_merged_3232_3233_small4M_outlier.root");
        TTree* cbmsim = (TTree*) inputfile->Get("cbmsim");

        MUonERecoOutput *ReconstructionOutput = 0;


   cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);


        std::vector<TH1D*> h_chi2(2);
h_chi2.at(0)=new TH1D("h_chi20","chi2 mu_in station 0 ",200,0.,100.);
h_chi2.at(1)=new TH1D("h_chi21","chi2 mu_out station 1 ",200,0.,100.);


// X or Y position on the track at a given Z
                auto pos_on_track = [](double q, double m, double z){return (q + m*z);};



for(Long64_t i = 0; i < 1000000; i++) {//cbmsim->GetEntries(); i++) {
		cbmsim->GetEntry(i);
		if(i%1000 == 0) cout<<"Entry "<<i<<endl;

vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
vector<MUonERecoOutputVertex> vrtx = ReconstructionOutput->reconstructedVertices();

    int sec0=0; int sec1=0;

        TVector3 in;
        TVector3 out1;
        TVector3 out2;


double x0=99;double y0=99;double x1=99;double y1=99;
double mx0=99;double my0=99;double mx1=99;double my1=99;double qx0=99;double qy0=99;double qx1=99;double qy1=99;
double chi2_0,chi2_1;

    for(int j=0; j<tracks.size();j++)
    {
std::vector<MUonERecoOutputHit> hits_=tracks.at(j).hits();
        if(tracks.at(j).sector()==0)
					{sec0++;
					 qx0=tracks.at(j).x0();
					 qy0=tracks.at(j).y0();
                                         mx0=tracks.at(j).xSlope();
                                         my0=tracks.at(j).ySlope();
					 chi2_0=tracks.at(j).chi2();}
        if(tracks.at(j).sector()==1)
					 {sec1++;
                                         qx1=tracks.at(j).x0();
                                         qy1=tracks.at(j).y0();
                                         mx1=tracks.at(j).xSlope();
                                         my1=tracks.at(j).ySlope();
                                         chi2_1=tracks.at(j).chi2();
}
	}

if(sec0==1 and sec1==1){

	h_chi2.at(0)->Fill(chi2_0);
        h_chi2.at(1)->Fill(chi2_1);

	}



}


TCanvas n1("n1","n1",1000,1000);
n1.Divide(1,2);
n1.cd(1);
h_chi2.at(0)->Draw();
n1.cd(2);
h_chi2.at(1)->Draw();
n1.SaveAs("chi2_6.pdf");
h_chi2.at(0)->SaveAs("chi2_muin_6.root");
h_chi2.at(1)->SaveAs("chi2_muout_6.root");


}


