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

void study_data(){

TChain * cbmsim = new TChain("cbmsim");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_2hit_1.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_2hit_2.root");


        TClonesArray *MCTrack = 0;
        MuE::Event *MesmerEvent = 0;
        MUonERecoOutput *ReconstructionOutput = 0;

        cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);

                auto pos_on_track = [](double q, double m, double z){return (q + m*z);};

double signal=0.;

double reco_1=0.; double reco1_1=0.; double more_reco_1=0.; double reco0_1=0.;double reco3_1=0.;
double reco_2=0.; double reco1_2=0.; double more_reco_2=0.; double reco0_2=0.;double reco3_2=0.;
double reco_3=0.; double reco1_3=0.; double more_reco_3=0.; double reco0_3=0.;double reco3_3=0.;
double reco_4=0.; double reco1_4=0.; double more_reco_4=0.; double reco0_4=0.;double reco3_4=0.;
double reco_5=0.; double reco1_5=0.; double more_reco_5=0.; double reco0_5=0.;double reco3_5=0.;
double reco_6=0.; double reco1_6=0.; double more_reco_6=0.; double reco0_6=0.;double reco3_6=0.;

double error1=0.;double error2=0.;double error3=0.;double error4=0.;double error5=0.;double error6=0.;

double reco_v=0.; double more_reco_v=0.; double reco0_v=0.;
int yes2=0; int yes_v=0;
int code_mu=-99; int code_e=-99; int code_mu_in=-99;
int TrackIdreco=-99;
double z_fix=912.7;

   TH1D* h_z_kin_pre=new TH1D("h_z_kin_pre"," Z selected events kin fit pre-selection",450, 890.,940.);
   TH1D* h_z_pos_pre=new TH1D("h_z_pos_pre","Z selected events pos fit pre-selection",450, 890.,940.);
   TH1D* h_z_kin=new TH1D("h_z_kin"," Z selected events kin fit",450, 890.,940.);
   TH1D* h_z_pos=new TH1D("h_z_pos","Z selected events pos fit",450, 890.,940.);
   TH1D* th_in=new TH1D("th_in","Incoming muon theta",100,0.,0.01);


for(Long64_t i = 0; i < 1000000; i++) {//cbmsim->GetEntries(); i++) {
		cbmsim->GetEntry(i);
		if(i%1000 == 0) cout<<"Entry "<<i<<endl;

vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
MUonERecoOutputVertex vrtx = ReconstructionOutput->bestVertex();

 double chi=vrtx.chi2perDegreeOfFreedom();

int yes_mu=0;
int yes_e=0;
double th_inx,th_iny,x0_in,y0_in;
double chi2_muin;
int sec0=0;
int sec1=0;
int stubs_muin=0.;
double th_muin=0.;

    for(int j=0; j<tracks.size();j++)
    {
        if(tracks.at(j).sector()==0) sec0++;
        if(tracks.at(j).sector()==1) sec1++;
        }


for(int j=0; j<tracks.size();j++)
{
        if(tracks.at(j).sector()==0 and sec0==1){
	std::vector<MUonERecoOutputHit> hits_=tracks.at(j).hits();
	stubs_muin=hits_.size();
        th_inx=tracks.at(j).xSlope();
        th_iny=tracks.at(j).ySlope();
        x0_in=tracks.at(j).x0();
        y0_in=tracks.at(j).y0();
        chi2_muin=tracks.at(j).chi2perDegreeOfFreedom();
	TVector3 p_muin(th_inx,th_iny,1.0);
        p_muin=p_muin.Unit();
        th_muin=p_muin.Theta();
                        }
	if(tracks.at(j).sector()==1 and sec1==2)
	{yes2++;}
}

//estrapolo posizione negli ultimi due moduli della seconda stazione
/*double posxIN=pos_on_track(x0_in,th_inx,899.92180);
double posyIN=pos_on_track(y0_in,th_iny,903.76930);*/

// posizione locale negli ultimi due moduli della seconda stazione
double posxIN=99.;
double posyIN=99.;

//posxIN=pos_on_track(x0_in,th_inx,z_fix);
//posyIN=pos_on_track(y0_in,th_iny,z_fix);


std::vector<MUonERecoOutputHit> stubs=ReconstructionOutput->reconstructedHits();

int stub0 = 0;
int stub1 = 0;
for(int s=0; s<stubs.size(); s++){
if(stubs.at(s).stationID()==0){stub0++; if(stubs.at(s).moduleID()==4){posxIN=stubs.at(s).position(); } else if(stubs.at(s).moduleID()==5){posyIN=stubs.at(s).position();}   }
if(stubs.at(s).stationID()==1)stub1++;
}

 if(sec0==1 and stubs_muin==6 and abs(posxIN)<=1.5 and abs(posyIN)<=1.5 and chi2_muin<=2 and stub0==6)th_in->Fill(th_muin);

 if(sec0==1 and stubs_muin==6 and abs(posxIN)<=1.5 and abs(posyIN)<=1.5 and chi2_muin<=2 and stub0==6 and th_muin<0.004){// and stub1<=15){

 signal++;

if(chi!=0){

        h_z_kin_pre->Fill(vrtx.zKinematicFit());
        h_z_pos_pre->Fill(vrtx.zPositionFit());

 MUonERecoOutputTrack mu_in = vrtx.incomingMuon();
 MUonERecoOutputTrack mu_out = vrtx.outgoingMuon();
 MUonERecoOutputTrack e_out = vrtx.outgoingElectron();
        TVector3 p_muin(mu_in.xSlope(),mu_in.ySlope(),1.0);
        TVector3 p_mu(mu_out.xSlope(),mu_out.ySlope(),1.0);
        TVector3 p_e(e_out.xSlope(),e_out.ySlope(),1.0);
	p_e.Unit();p_mu.Unit();p_muin.Unit();

                                                double dotProduct_v = p_mu.Dot(p_e);
                                                TVector3 crossProduct_v = p_mu.Cross(p_e);
                                                double T_v = p_muin.Dot(crossProduct_v);
                                                TVector3 im_v= p_muin.Cross(p_mu);
                                                TVector3 ie_v= p_muin.Cross(p_e);
                                                T_v = T_v>0? 1:-1;
                                                double acoplanarity_v= T_v*(TMath::Pi()- acos( ((im_v).Dot(ie_v))/(im_v.Mag()*ie_v.Mag()) ));

  if(abs(acoplanarity_v)<=1 and chi<20 and vrtx.muonTheta()>=0.0002 and stub1<=15){//vrtx.electronTheta()>=0.0005 and vrtx.electronTheta()<=0.02){//vrtx.electronTheta()<=0.032){
//first
 if(yes2>=2){

	h_z_kin->Fill(vrtx.zKinematicFit());
	h_z_pos->Fill(vrtx.zPositionFit());


}

if(vrtx.electronTheta()>0.0 and vrtx.electronTheta()<=0.005){
if(yes2>=2){reco_1++;}

if(yes2==2 and tracks.size()==3) {reco3_1++;}

}

//second
if(vrtx.electronTheta()>0.005 and vrtx.electronTheta()<=0.01){
if(yes2>=2){reco_2++;}

if(yes2==2 and tracks.size()==3) {reco3_2++;}

}

//third
if(vrtx.electronTheta()>0.01 and vrtx.electronTheta()<=0.015){
if(yes2>=2){reco_3++;}

if(yes2==2 and tracks.size()==3) {reco3_3++;}

}

//fourth
if(vrtx.electronTheta()>0.015 and vrtx.electronTheta()<=0.02){
if(yes2>=2){reco_4++;}

if(yes2==2 and tracks.size()==3) {reco3_4++;}

}

//fifth
if(vrtx.electronTheta()>0.02 and vrtx.electronTheta()<=0.025){
if(yes2>=2){reco_5++;}

if(yes2==2 and tracks.size()==3) {reco3_5++;}

}

//sixth
if(vrtx.electronTheta()>0.025 and vrtx.electronTheta()<=0.032){
if(yes2>=2){reco_6++;}

if(yes2==2 and tracks.size()==3) {reco3_6++;}

}

			}//aco e chi
		}//chi!=0
	}//mu_in
yes2=0;
} //end of general for



TCanvas a("a","a",700,700);
a.Divide(1,2);
a.cd(1);
h_z_pos_pre->SetLineColor(kOrange+10);
h_z_pos_pre->Draw("hist");
a.cd(2);
h_z_pos->SetLineColor(kOrange+10);
h_z_pos->Draw("hist");

a.SaveAs("data_pos.pdf");


}
