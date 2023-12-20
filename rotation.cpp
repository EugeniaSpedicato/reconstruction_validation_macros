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

	std::vector<TH1D*> h_dif(2);
h_dif.at(0)=new TH1D("h_div0","deltaX station 0 - station 1 at Z=921",100,-0.5,0.5);//100,0.,0.2);
h_dif.at(1)=new TH1D("h_div1","deltaY station 0 - station 1 at Z=921",100,-0.5,0.5);//100,0.2,0.5);

        std::vector<TH1D*> h_imp(4);
h_imp.at(0)=new TH1D("h_imp0","X profile station 0 ",800,-5.,5.);
h_imp.at(1)=new TH1D("h_imp1","Y profile station 0 ",800,-5.,5.);
h_imp.at(2)=new TH1D("h_imp2","X profile station 1 ",800,-5.,5.);
h_imp.at(3)=new TH1D("h_imp3","Y profile station 1 ",800,-5.,5.);


        std::vector<TH1D*> h_m(4);
h_m.at(0)=new TH1D("h_m0","Delta mx1-mx0",200,-0.006,0.006);//-0.005,0.);
h_m.at(1)=new TH1D("h_m1","Delta my1-my0",200,-0.006,0.006);//-0.006,0.);
h_m.at(2)=new TH1D("h_m2","Delta mx1-mx0 after rotation",200,-0.006,0.006);//-0.005,0.);
h_m.at(3)=new TH1D("h_m3","Delta my1-my0 after rotation",200,-0.006,0.006);//-0.006,0.);

double z_mod[6]={18.0218,21.8693,55.3635,56.6205,89.9218,93.7693};

        std::vector<TH1D*> localX1(6);
for(int m=0; m<6; m++){
	        string name="local1_"+to_string(m);
                string title="Residual station1 of module "+to_string(m);
		localX1.at(m)=new TH1D(name.c_str(),title.c_str(),600,-0.06,0.06);}

        std::vector<TH1D*> localX2(6);
for(int m=0; m<6; m++){
                string name="local2_"+to_string(m);
                string title="Residual of stub in station2 and expected position from track station1 of module "+to_string(m);
                localX2.at(m)=new TH1D(name.c_str(),title.c_str(),400,-0.4,0.4);}//500,-5.,5.);}


// X or Y position on the track at a given Z
                auto pos_on_track = [](double q, double m, double z){return (q + m*z);};



                //interface to parse the alignment parameters from FairMUonE (yaml file)
                YAML::Node alignmentFile = YAML::LoadFile("/home/espedica/fair_install/instFairRoot/share/MUonE/common/alignment/TR_separate.yaml");//TR_alin.yaml");

                for(int station_index = 0; station_index < alignmentFile.size(); station_index++) {
                        auto const& station = alignmentFile[station_index];
                        for(int module_index = 0; module_index < station.size(); module_index++) {
                                int linkID = module_index + 6*station_index;
                                auto const& module = station[module_index];
                                XOFFSET[linkID] = module["xOffset"].as<double>();
                                YOFFSET[linkID] = module["yOffset"].as<double>();
                                ZOFFSET[linkID] = module["zOffset"].as<double>();
                                ALPHAOFFSET[linkID] = module["angleOffset"].as<double>()*TMath::DegToRad();
                                TILTOFFSET[linkID]  = module["tiltOffset"].as<double>()*TMath::DegToRad();

                                COSALPHA[linkID] = TMath::Cos(ALPHA[linkID] + ALPHAOFFSET[linkID]);
                                SINALPHA[linkID] = TMath::Sin(ALPHA[linkID] + ALPHAOFFSET[linkID]);

                                COSTILT[linkID] = TMath::Cos(TILT[linkID] + TILTOFFSET[linkID]);
                                SINTILT[linkID] = TMath::Sin(TILT[linkID] + TILTOFFSET[linkID]);
                        }
                }







for(Long64_t i = 0; i < 1000000; i++) {//cbmsim->GetEntries(); i++) {
		cbmsim->GetEntry(i);
		if(i%1000 == 0) cout<<"Entry "<<i<<endl;

vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
vector<MUonERecoOutputVertex> vrtx = ReconstructionOutput->reconstructedVertices();
//vector<MUonERecoOutputAdaptiveFitterVertex> vrtx = ReconstructionOutput->adaptiveFitterVertices();


    int sec0=0; int sec1=0;

        TVector3 in;
        TVector3 out1;
        TVector3 out2;


double x0=99;double y0=99;double x1=99;double y1=99;
double mx0=99;double my0=99;double mx1=99;double my1=99;double qx0=99;double qy0=99;double qx1=99;double qy1=99;
double exp_pos[6];
double pos1[6];
double pos2[6];

    for(int j=0; j<tracks.size();j++)
    {
std::vector<MUonERecoOutputHit> hits_=tracks.at(j).hits();
        if(tracks.at(j).sector()==0)
					{sec0++;
					 qx0=tracks.at(j).x0();
					 qy0=tracks.at(j).y0();
                                         mx0=tracks.at(j).xSlope();
                                         my0=tracks.at(j).ySlope();}
        if(tracks.at(j).sector()==1)
					 {sec1++;
                                         qx1=tracks.at(j).x0();
                                         qy1=tracks.at(j).y0();
                                         mx1=tracks.at(j).xSlope();
                                         my1=tracks.at(j).ySlope();
}
	}

if(sec0==1 and sec1==1){
TVector3 k_in(-0.002133,-0.00393,1.0);
TVector3 k_out(mx1,my1,1.0);

Double_t pz=k_in.Z();
Double_t py=k_in.Y();
Double_t px=k_in.X();

Double_t ptz=sqrt(px*px+pz*pz);


// costruzione della matrice di rotazione
Double_t psi=atan2(px,pz);
Double_t phi=atan2(py,ptz);

TMatrixD R(3,3);
R[0][0]=cos(psi);
R[0][1]=0;
R[0][2]=sin(psi);
    R[1][0]=-sin(phi)*sin(psi);
    R[1][1]=cos(phi);
    R[1][2]=sin(phi)*cos(psi);
        R[2][0]=-cos(phi)*sin(psi);
        R[2][1]=-sin(phi);
        R[2][2]=cos(phi)*cos(psi);

TMatrixD pO(3,1);
    pO[0][0]=k_out.X();
    pO[1][0]=k_out.Y();
    pO[2][0]=k_out.Z();


TMatrixD pN(R, TMatrixD::kMult,pO);
//cout << "[0][0] " << pN[0][0] << " VS " << mx0 << endl;
//cout << "[1][0] " << pN[1][0] << " VS " << my0 << endl;
//cout << "[2][0] " << pN[2][0] << " VS 1.0" << endl;

h_m[0]->Fill(mx1-mx0);
h_m[1]->Fill(my1-my0);
h_m[2]->Fill(pN[0][0]-mx0);
h_m[3]->Fill(pN[1][0]-my0);

                                         x0=pos_on_track(qx0,mx0,912.7);
                                         y0=pos_on_track(qy0,my0,912.7);
                                         x1=pos_on_track(qx1,pN[0][0],912.7);
                                         y1=pos_on_track(qy1,pN[1][0],912.7);

h_dif[0]->Fill(x1-x0);
h_dif[1]->Fill(y1-y0);


}



}


TCanvas n1("n1","n1",1000,1000);
n1.Divide(2,2);
n1.cd(1);
h_m[0]->Draw();
n1.cd(2);
h_m[1]->Draw();
n1.cd(3);
h_m[2]->Draw();
n1.cd(4);
h_m[3]->Draw();
n1.SaveAs("beam_fmu.pdf");


TCanvas n2("n2","n2",1000,1000);
n2.Divide(1,2);
n2.cd(1);
h_dif[0]->Draw();
n2.cd(2);
h_dif[1]->Draw();
n2.SaveAs("imp.pdf");

}


