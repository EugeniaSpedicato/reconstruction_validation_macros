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


void RealDataAnalyzer(){

 	TFile *inputfile = new TFile("/mnt/raid10/DATA/espedica/fairmu/dataReconstruction_3234-3235_new12_st0_6.root");
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

    TH1D *h_mult0=new TH1D("h_mult0","Multiplicity of track in station0",10,0,10);
    TH1D *h_mult1=new TH1D("h_mult1","Multiplicity of track in station1",10,0,10);

    TH2D *h_2d_1=new TH2D("h_2d_0","XY impact point station 0 (black) and station 1 (red)" ,800,-5.,5.,800,-5.,5.);
    TH2D *h_2d_2=new TH2D("h_2d_1","XY impact point station 0 (black) and station 1 (red)" ,800,-5.,5.,800,-5.,5.);

	std::vector<TH1D*> h_dif(2);
h_dif.at(0)=new TH1D("h_div0","deltaX station 0 - station 1 at Z=921",100,-0.01,0.01);//100,0.,0.2);
h_dif.at(1)=new TH1D("h_div1","deltaY station 0 - station 1 at Z=921",100,-0.01,0.01);//100,0.2,0.5);

        std::vector<TH1D*> h_imp(4);
h_imp.at(0)=new TH1D("h_imp0","X profile station 0 ",800,-5.,5.);
h_imp.at(1)=new TH1D("h_imp1","Y profile station 0 ",800,-5.,5.);
h_imp.at(2)=new TH1D("h_imp2","X profile station 1 ",800,-5.,5.);
h_imp.at(3)=new TH1D("h_imp3","Y profile station 1 ",800,-5.,5.);


        std::vector<TH1D*> h_m(2);
h_m.at(0)=new TH1D("h_m0","Delta mx1-mx0",200,-0.001,0.001);//-0.005,0.);
h_m.at(1)=new TH1D("h_m1","Delta my1-my0",200,-0.001,0.001);//-0.006,0.);



    TH1D *h_Z=new TH1D("h_Z","Z of the vrtx found",800,800.,1000.);

double z_mod[6]={18.0218,21.8693,55.3635,56.6205,89.9218,93.7693};

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
                string title="Residual station1 of module "+to_string(m);
		localX1.at(m)=new TH1D(name.c_str(),title.c_str(),600,-0.06,0.06);}

        std::vector<TH1D*> localX2(6);
for(int m=0; m<6; m++){
                string name="local2_"+to_string(m);
                string title="Residual of stub in station2 and expected position from track station1 of module "+to_string(m);
                localX2.at(m)=new TH1D(name.c_str(),title.c_str(),400,-0.4,0.4);}//500,-5.,5.);}


TH1D *h_nstubs_1=new TH1D("h_nstubs_1", "number of stubs for a track in station1", 13,0,13);
TH1D *h_nstubs_2=new TH1D("h_nstubs_2", "number of stubs for a track in station1", 13,0,13);


// X or Y position on the track at a given Z
                auto pos_on_track = [](double q, double m, double z){return (q + m*z);};



                //interface to parse the alignment parameters from FairMUonE (yaml file)
                YAML::Node alignmentFile = YAML::LoadFile("/home/espedica/fair_install/instFairRoot/share/MUonE/common/alignment/TR_alin.yaml");//TR_3234-3235_new6.yaml");//TR_alin.yaml");

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
					 x0=pos_on_track(tracks.at(j).x0(),tracks.at(j).xSlope(),912.7);
					 y0=pos_on_track(tracks.at(j).y0(),tracks.at(j).ySlope(),912.7);
                                         mx0=tracks.at(j).xSlope();
                                         my0=tracks.at(j).ySlope();}
        if(tracks.at(j).sector()==1)
					 {sec1++;
                                         qx1=tracks.at(j).x0();
                                         qy1=tracks.at(j).y0();
					 x1=pos_on_track(tracks.at(j).x0(),tracks.at(j).xSlope(),912.7);
					 y1=pos_on_track(tracks.at(j).y0(),tracks.at(j).ySlope(),912.7);
                                         mx1=tracks.at(j).xSlope();
                                         my1=tracks.at(j).ySlope();
  /*for(int h=0;h<hits_.size();h++)
  {
  localX2.at(hits_.at(h).moduleID())->Fill(hits_.at(h).positionPerpendicular());
	}*/

}
	}




    for(int j=0; j<tracks.size();j++)
    {

     	if(tracks.at(j).sector()==0 and sec0==1 and sec1==1){

/*	ROOT::Math::SMatrix<Double_t, 4> cov= tracks.at(j).linearFitCovarianceMatrix();

	for(int r=0; r<4; r++){
		cout << "riga " << r << endl;
		for(int c=0; c<4; c++){
			cout << cov[r][c] << " ";
				}
			cout << endl;
			}

cout << "tracks.at(j).xSlopeError()^2 " << tracks.at(j).xSlopeError()*tracks.at(j).xSlopeError() << endl;
*/
std::vector<MUonERecoOutputHit> hits_=tracks.at(j).hits();
h_nstubs_1->Fill(hits_.size());
  for(int h=0;h<hits_.size();h++)
  {
/*cout << "staz 1: error mod " << hits_.at(h).moduleID() << ": " <<
sqrt(((hits_.at(h).z())*tracks.at(j).xSlopeError())*((hits_.at(h).z())*tracks.at(j).xSlopeError()) + tracks.at(j).x0Error()*tracks.at(j).x0Error()) << endl;
cout << "staz 1: error con correlazini mod " << hits_.at(h).moduleID() << ": "
<< sqrt((hits_.at(h).z()*tracks.at(j).xSlopeError())*(hits_.at(h).z()*tracks.at(j).xSlopeError()) + tracks.at(j).x0Error()*tracks.at(j).x0Error() + 2*cov[0][2]) << endl;
*/
        if(sec1==1 and hits_.at(h).stationID()==0) residual1.at(hits_.at(h).moduleID())->Fill(hits_.at(h).perpendicularResiduum());
        //if(hits_.at(h).moduleID()%2==0)
	localX1.at(hits_.at(h).moduleID())->Fill(-hits_.at(h).positionPerpendicular()+(pos_on_track(qx0,mx0,hits_.at(h).z())-XOFFSET[hits_.at(h).moduleID()])*COSALPHA[hits_.at(h).moduleID()]+(pos_on_track(qy0,my0,hits_.at(h).z())-YOFFSET[hits_.at(h).moduleID()])*SINALPHA[hits_.at(h).moduleID()]);

/*if(hits_.at(h).moduleID()%2==0)cout << hits_.at(h).moduleID() <<
") perp pos: " << hits_.at(h).positionPerpendicular()
<< ", calculated: " << (pos_on_track(qx0,mx0,hits_.at(h).z())-XOFFSET[hits_.at(h).moduleID()])*COSALPHA[hits_.at(h).moduleID()]
+(pos_on_track(qy0,my0,hits_.at(h).z())-YOFFSET[hits_.at(h).moduleID()])*SINALPHA[hits_.at(h).moduleID()]
<< ", z " << hits_.at(h).z() <<endl;
*/
  }
 }

   	if(tracks.at(j).sector()==1 and sec1==1 and sec0==1){


std::vector<MUonERecoOutputHit> hits_=tracks.at(j).hits();
h_nstubs_2->Fill(hits_.size());
  for(int h=0;h<hits_.size();h++)
  {
	if(sec0==1 and hits_.at(h).stationID()==1) residual2.at(hits_.at(h).moduleID())->Fill(hits_.at(h).perpendicularResiduum());
localX2.at(hits_.at(h).moduleID())->Fill(-hits_.at(h).positionPerpendicular()+(pos_on_track(qx0,mx0,hits_.at(h).z())-XOFFSET[hits_.at(h).moduleID()+6])*COSALPHA[hits_.at(h).moduleID()+6]+(pos_on_track(qy0,my0,hits_.at(h).z())-YOFFSET[hits_.at(h).moduleID()+6])*SINALPHA[hits_.at(h).moduleID()+6]);
	//if(hits_.at(h).moduleID()%2!=0)localX2.at(hits_.at(h).moduleID())->Fill(-hits_.at(h).positionPerpendicular()+(pos_on_track(qy0,my0,hits_.at(h).z())-YOFFSET[hits_.at(h).moduleID()+6])*SINALPHA[hits_.at(h).moduleID()+6]);

  }
 }

}

h_mult0->Fill(sec0);
h_mult1->Fill(sec1);


if(x0!=99 and x1!=99 and sec0==1 and sec1==1 and y0!=99 and y1!=99 and sec0==1 and sec1==1){
cout << "event "<< i << endl;
h_imp[0]->Fill(x0);h_imp[1]->Fill(y0);h_dif[0]->Fill(x0-x1);
h_imp[2]->Fill(x1);h_imp[3]->Fill(y1);h_dif[1]->Fill(y0-y1);
h_2d_1->Fill(x0,y0);
h_2d_2->Fill(x1,y1);
h_m[0]->Fill(mx0-mx1);h_m[1]->Fill(my0-my1);
//cout << "--------- event "<< i << endl;
//cout << "x0-x1 " << x0-x1 << ", y0-y1 " << y0-y1 << endl;
}



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
h_m[0]->Draw();
n2.cd(4);
h_m[1]->Draw();
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
n5.Divide(2,4);
n5.cd(1);
localX2[0]->SetLineColor(kRed);
localX2[0]->Draw();
gPad->SetLogy();
n5.cd(2);
localX2[1]->SetLineColor(kRed);
localX2[1]->Draw();
gPad->SetLogy();
n5.cd(3);
localX2[2]->SetLineColor(kRed);
localX2[2]->Draw();
gPad->SetLogy();
n5.cd(4);
localX2[3]->SetLineColor(kRed);
localX2[3]->Draw();
gPad->SetLogy();
n5.cd(5);
localX2[4]->SetLineColor(kRed);
localX2[4]->Draw();
gPad->SetLogy();
n5.cd(6);
localX2[5]->SetLineColor(kRed);
localX2[5]->Draw();
gPad->SetLogy();
/*n5.cd(5);
localX1[0]->Draw();
gPad->SetLogy();
n5.cd(6);
localX1[1]->Draw();
gPad->SetLogy();
n5.cd(7);
localX1[4]->Draw();
gPad->SetLogy();
n5.cd(8);
localX1[5]->Draw();
gPad->SetLogy();
*/
/*for(int m=0; m<6; m++){
n5.cd(m+1);
localX2[m]->SetLineColor(kRed);
localX2[m]->Draw();}*/
//gPad->SetLogy();}
n5.SaveAs("local.pdf");
/*
TCanvas n6("n6","n6",700,700);
h_nstubs_1->Draw();
h_nstubs_2->SetLineColor(kRed);
h_nstubs_2->Draw("same");
gPad->SetLogy();
n6.SaveAs("nstubs_6_6.pdf");
*/
}


