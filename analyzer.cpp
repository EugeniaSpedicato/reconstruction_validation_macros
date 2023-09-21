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

 	TFile *inputfile = new TFile("/mnt/raid10/DATA/espedica/fairmu/dataReconstruction_3234-3235_filtered.root");//dataReconstruction_target.root");
//dataReconstruction_merged_3232_3233_small4M_outlier.root");
//dataReconstruction_merged_3238_3239_0hit_1.root");
//dataReconstruction_merged_3232_3233_small4M_outlier.root");
        TTree* cbmsim = (TTree*) inputfile->Get("cbmsim");

        MUonERecoOutput *ReconstructionOutput = 0;


   cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);

    TH1D *h_mult0=new TH1D("h_mult0","Multiplicity of track in station0",10,0,10);
    TH1D *h_mult1=new TH1D("h_mult1","Multiplicity of track in station1",10,0,10);

    TH2D *h_2d_1=new TH2D("h_2d_0","XY impact point station 0 (black) and station 1 (red)" ,800,-5.,5.,800,-5.,5.);
    TH2D *h_2d_2=new TH2D("h_2d_1","XY impact point station 0 (black) and station 1 (red)" ,800,-5.,5.,800,-5.,5.);

	std::vector<TH1D*> h_dif(2);
h_dif.at(0)=new TH1D("h_div0","deltaX station 0 - station 1 at Z=921",100,-0.03,0.03);
h_dif.at(1)=new TH1D("h_div1","deltaY station 0 - station 1 at Z=921",100,-0.03,0.03);

        std::vector<TH1D*> h_imp(4);
h_imp.at(0)=new TH1D("h_imp0","X profile station 0 ",800,-5.,5.);
h_imp.at(1)=new TH1D("h_imp1","Y profile station 0 ",800,-5.,5.);
h_imp.at(2)=new TH1D("h_imp2","X profile station 1 ",800,-5.,5.);
h_imp.at(3)=new TH1D("h_imp3","Y profile station 1 ",800,-5.,5.);

    TH1D *h_Z=new TH1D("h_Z","Z of the vrtx found",800,800.,1000.);

        std::vector<TH1D*> residual1(6);
for(int m=0; m<6; m++){
	string name="residual1_"+to_string(m);
                string title="Residual station1 of module "+to_string(m);
if(m==2 or m==3) residual1.at(m)=new TH1D(name.c_str(),title.c_str(),600,-0.06,0.06);//100,-0.1,0.1);
else residual1.at(m)=new TH1D(name.c_str(),title.c_str(),600,-0.06,0.06);//500,-0.005,0.005);
			}

        std::vector<TH1D*> residual2(6);
for(int m=0; m<6; m++){
        string name="residual2_"+to_string(m);
                string title="Residual station2 of module "+to_string(m);
if(m==2 or m==3) residual2.at(m)=new TH1D(name.c_str(),title.c_str(),600,-0.06,0.06);//100,-0.1,0.1);
else residual2.at(m)=new TH1D(name.c_str(),title.c_str(),600,-0.06,0.06);//500,-0.005,0.005);
                        }

        std::vector<TH1D*> localX1(6);
for(int m=0; m<6; m++){
	        string name="local1_"+to_string(m);
                string title="LocalX station1 of module "+to_string(m);
		localX1.at(m)=new TH1D(name.c_str(),title.c_str(),100,-5.,5.);}

        std::vector<TH1D*> localX2(6);
for(int m=0; m<6; m++){
                string name="local2_"+to_string(m);
                string title="LocalX station2 of module "+to_string(m);
                localX2.at(m)=new TH1D(name.c_str(),title.c_str(),100,-5.,5.);}


TH1D *h_nstubs_1=new TH1D("h_nstubs_1", "number of stubs for a track in station1", 13,0,13);
TH1D *h_nstubs_2=new TH1D("h_nstubs_2", "number of stubs for a track in station1", 13,0,13);


// X or Y position on the track at a given Z
                auto pos_on_track = [](double q, double m, double z){return (q + m*z);};

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
double exp_pos[6];
double pos1[6];
double pos2[6];

    for(int j=0; j<tracks.size();j++)
    {

        if(tracks.at(j).sector()==0)
					{sec0++;
					 x0=pos_on_track(tracks.at(j).x0(),tracks.at(j).xSlope(),912.7);
					 y0=pos_on_track(tracks.at(j).y0(),tracks.at(j).ySlope(),912.7);}
        if(tracks.at(j).sector()==1)
					 {sec1++;
					 x1=pos_on_track(tracks.at(j).x0(),tracks.at(j).xSlope(),912.7);
					 y1=pos_on_track(tracks.at(j).y0(),tracks.at(j).ySlope(),912.7);}
	}

    for(int j=0; j<tracks.size();j++)
    {

     	if(tracks.at(j).sector()==0 and sec0==1){

std::vector<MUonERecoOutputHit> hits_=tracks.at(j).hits();
h_nstubs_1->Fill(hits_.size());
  for(int h=0;h<hits_.size();h++)
  {
        if(sec1==1 and hits_.at(h).stationID()==0) residual1.at(hits_.at(h).moduleID())->Fill(hits_.at(h).perpendicularResiduum());
        localX1.at(hits_.at(h).moduleID())->Fill(hits_.at(h).positionPerpendicular());


/*pos1[hits_.at(h).moduleID()]=hits_.at(h).position();
   if(hits_.at(h).moduleID()%2==0)exp_posX1[hits_.at(h).moduleID()]=pos_on_track(tracks.at(j).x0(),tracks.at(j).xSlope(),hits_.at(h).z());
   if(hits_.at(h).moduleID()%2!=0)exp_posY1[hits_.at(h).moduleID()]=pos_on_track(tracks.at(j).y0(),tracks.at(j).ySlope(),hits_.at(h).z());
*/
  }
 }

   	if(tracks.at(j).sector()==1 and sec1==1){
std::vector<MUonERecoOutputHit> hits_=tracks.at(j).hits();
h_nstubs_2->Fill(hits_.size());
  for(int h=0;h<hits_.size();h++)
  {
	if(sec0==1 and hits_.at(h).stationID()==1) residual2.at(hits_.at(h).moduleID())->Fill(hits_.at(h).perpendicularResiduum());
	localX2.at(hits_.at(h).moduleID())->Fill(hits_.at(h).positionPerpendicular());
  }
 }

}

h_mult0->Fill(sec0);
h_mult1->Fill(sec1);


if(x0!=99 and x1!=99 and sec0==1 and sec1==1 and y0!=99 and y1!=99 and sec0==1 and sec1==1){
h_imp[0]->Fill(x0);h_imp[1]->Fill(y0);h_dif[0]->Fill(x0-x1);
h_imp[2]->Fill(x1);h_imp[3]->Fill(y1);h_dif[1]->Fill(y0-y1);
h_2d_1->Fill(x0,y0);
h_2d_2->Fill(x1,y1);}



}


    TCanvas n("n","n",700,700);
    n.Divide(1,2);
    n.cd(1);
    h_mult0->Draw();
    h_mult1->SetLineColor(kRed);
    h_mult1->Draw("same");
    n.cd(2);
    h_2d_1->Draw();
    h_2d_2->SetMarkerColor(kRed);
    h_2d_2->Draw("same");
    n.SaveAs("info.pdf");

TCanvas n2("n2","n2",1000,1000);
n2.Divide(2,2);
n2.cd(1);
h_dif[0]->Draw();
n2.cd(2);
h_dif[1]->Draw();
n2.cd(3);
h_imp[0]->Draw();
h_imp[2]->SetLineColor(kRed);
h_imp[2]->Draw("same");
n2.cd(4);
h_imp[1]->Draw();
h_imp[3]->SetLineColor(kRed);
h_imp[3]->Draw("same");
n2.SaveAs("beam_fmu.pdf");

TCanvas n3("n3","n3",1000,1000);
n3.Divide(2,3);
for(int m=0; m<6; m++){
n3.cd(m+1);
residual1[m]->Draw();
gPad->SetLogy();}
n3.SaveAs("res1.pdf");


TCanvas n4("n4","n4",1000,1000);
n4.Divide(2,3);
for(int m=0; m<6; m++){
n4.cd(m+1);
residual2[m]->Draw();
gPad->SetLogy();}
n4.SaveAs("res2.pdf");

TCanvas n5("n5","n5",1000,1000);
n5.Divide(2,3);
for(int m=0; m<6; m++){
n5.cd(m+1);
localX1[m]->Draw();
localX2[m]->SetLineColor(kRed);
localX2[m]->Draw("same");
gPad->SetLogy();}
n5.SaveAs("local.pdf");

TCanvas n6("n6","n6",700,700);
h_nstubs_1->Draw();
h_nstubs_2->SetLineColor(kRed);
h_nstubs_2->Draw("same");
gPad->SetLogy();
n6.SaveAs("nstubs_6_6.pdf");

}


