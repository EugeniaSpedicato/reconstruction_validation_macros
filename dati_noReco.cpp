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

        TFile *inputfile = new TFile("/mnt/raid10/DATA/espedica/fairmu/dataReconstruction_3234-3235_new_1hit.root");

//dataReconstruction_3234-3235_new12_st0_6_1shared_5M.root");
//minbias_1shared_5M.root");
//dataReconstruction_3234-3235_new12_st0_6_1shared_5M.root");
        TTree* cbmsim = (TTree*) inputfile->Get("cbmsim");

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

TH2D* h_xy=new TH2D("h_xy","Beam Spot Z=target",100,-5.,5.,100,-5.,5.);
TH2D* h_2d=new TH2D("h2D","theta mu vs theta E with all cuts", 300,0.,0.035,50,0.,0.005);

TH1::SetDefaultSumw2(kTRUE);

   const Int_t NBINS = 6;
   Double_t edges[NBINS + 1] = {0.0, 0.005, 0.010, 0.015, 0.020, 0.025, 0.032};
   TH1D* d_eff = new TH1D("d_eff_real", "Efficiency as a function of the electron's angle",NBINS,edges);
   TH1D* theta_e = new TH1D("theta_e", "Electron scattering angles from MESMER",10,0.,0.035);
   TH1D* theta_mu = new TH1D("theta_mu", "Muon scattering angles from MESMER",20,0.,0.005);

for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
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
                        }
	if(tracks.at(j).sector()==1 and sec1==2)
	{yes2++;}
}

double posxIN=pos_on_track(x0_in,th_inx,z_fix);
double posyIN=pos_on_track(y0_in,th_iny,z_fix);

h_xy->Fill(posxIN,posyIN);

 if(sec0==1 and stubs_muin==6 and abs(posxIN)<=1.5 and abs(posyIN)<=1.5 and chi2_muin<2){
 signal++;

if(chi!=0){
 MUonERecoOutputTrack mu_in = vrtx.incomingMuon();
 MUonERecoOutputTrack mu_out = vrtx.outgoingMuon();
 MUonERecoOutputTrack e_out = vrtx.outgoingElectron();
        TVector3 p_muin(mu_in.xSlope(),mu_in.ySlope(),1.0);
        TVector3 p_mu(mu_out.xSlope(),mu_out.ySlope(),1.0);
        TVector3 p_e(e_out.xSlope(),e_out.ySlope(),1.0);

                                                double dotProduct_v = p_mu.Dot(p_e);
                                                TVector3 crossProduct_v = p_mu.Cross(p_e);
                                                double T_v = p_muin.Dot(crossProduct_v);
                                                TVector3 im_v= p_muin.Cross(p_mu);
                                                TVector3 ie_v= p_muin.Cross(p_e);
                                                T_v = T_v>0? 1:-1;
                                                double acoplanarity_v= T_v*(TMath::Pi()- acos( ((im_v).Dot(ie_v))/(im_v.Mag()*ie_v.Mag()) ));

  if(abs(acoplanarity_v)<=1 and chi<20 and vrtx.muonTheta()>0.0002){
//first
 if(yes2>=2){d_eff->Fill(vrtx.electronTheta());

 h_2d->Fill(vrtx.electronTheta(),vrtx.muonTheta());
 theta_mu->Fill(vrtx.muonTheta());
 theta_e->Fill(vrtx.electronTheta());}

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


cout <<endl;
cout << "Theta range: [0,5] mrad" << endl;
cout << "Su " << signal << " golden muons, " << reco_1 << " +- " << sqrt(reco_1) << " sono ricostruiti, con un rapporto del " << reco_1/signal*100 << "%"<< endl;
cout << "Su " << signal << " golden muons, " << reco3_1 << " sono ricostruiti con 3 tracce, con un rapporto del " << (reco3_1/signal)*100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con piu tracce (stesso id), " << more_reco_1 << " sono ricostruiti, con un rapporto del " << more_reco_1/signal*100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con 0 tracce di segnale reco, " << reco0_1 << ", con un rapporto del " << reco0_1/signal*100 << "%"<< endl;
//cout << "Su " << signal << " eventi di con 1 sola traccia di segnale reco, " << reco1_1 << ", con un rapporto del " << reco1_1/signal*100 << "%"<< endl;
cout <<endl;
cout << "Theta range: [5,10] mrad" << endl;
cout << "Su " << signal << " golden muons, " << reco_2 << " +- " << sqrt(reco_2) << " sono ricostruiti, con un rapporto del " << reco_2/signal *100 << "%"<< endl;
cout << "Su " << signal << " golden muons, " << reco3_2 << " sono ricostruiti con 3 tracce, con un rapporto del " << (reco3_2/signal)*100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con piu tracce (stesso id), " << more_reco_2 << " sono ricostruiti, con un rapporto del " <<  more_reco_3/signal *100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con 0 tracce di segnale reco, " << reco0_2 << ", con un rapporto del " << reco0_2/signal *100 << "%"<< endl;
//cout << "Su " << signal << " eventi di con 1 sola traccia di segnale reco, " << reco1_2 << ", con un rapporto del " << reco1_2/signal *100 << "%"<< endl;
cout <<endl;
cout << "Theta range: [10,15] mrad" << endl;
cout << "Su " << signal << " golden muons, " << reco_3 << " +- " << sqrt(reco_3) << " sono ricostruiti, con un rapporto del " << reco_3/signal *100 << "%"<< endl;
cout << "Su " << signal << " golden muons, " << reco3_3 << " sono ricostruiti con 3 tracce, con un rapporto del " << (reco3_3/signal)*100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con piu tracce (stesso id), " << more_reco_3 << " sono ricostruiti, con un rapporto del " <<  more_reco_3/signal *100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con 0 tracce di segnale reco, " << reco0_3 << ", con un rapporto del " << reco0_3/signal *100 << "%"<< endl;
//cout << "Su " << signal << " eventi di con 1 sola traccia di segnale reco, " << reco1_3 << ", con un rapporto del " << reco1_3/signal *100 << "%"<< endl;
cout <<endl;
cout << "Theta range: [15,20] mrad" << endl;
cout << "Su " << signal << " golden muons, " << reco_4 << " +- " << sqrt(reco_4) << " sono ricostruiti, con un rapporto del " << reco_4/signal *100 << "%"<< endl;
cout << "Su " << signal << " golden muons, " << reco3_4 << " sono ricostruiti con 3 tracce, con un rapporto del " << (reco3_4/signal)*100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con piu tracce (stesso id), " << more_reco_4 << " sono ricostruiti, con un rapporto del " <<  more_reco_4/signal *100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con 0 tracce di segnale reco, " << reco0_4 << ", con un rapporto del " << reco0_4/signal *100 << "%"<< endl;
//cout << "Su " << signal << " eventi di con 1 sola traccia di segnale reco, " << reco1_4 << ", con un rapporto del " << reco1_4/signal *100 << "%"<< endl;
cout <<endl;
cout << "Theta range: [20,25] mrad" << endl;
cout << "Su " << signal << " golden muons, " << reco_5 << " +- " << sqrt(reco_5) << " sono ricostruiti, con un rapporto del " << reco_5/signal *100 << "%"<< endl;
cout << "Su " << signal << " golden muons, " << reco3_5 << " sono ricostruiti con 3 tracce, con un rapporto del " << (reco3_5/signal)*100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con piu tracce (stesso id), " << more_reco_5 << " sono ricostruiti, con un rapporto del " <<  more_reco_5/signal *100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con 0 tracce di segnale reco, " << reco0_5 << ", con un rapporto del " << reco0_5/signal *100 << "%"<< endl;
//cout << "Su " << signal << " eventi di con 1 sola traccia di segnale reco, " << reco1_5 << ", con un rapporto del " << reco1_5/signal *100 << "%"<< endl;
cout <<endl;
cout << "Theta range: [25,32] mrad" << endl;
cout << "Su " << signal << " golden muons, " << reco_6 << " +- " << sqrt(reco_6) << " sono ricostruiti, con un rapporto del " << reco_6/signal *100 << "%"<< endl;
cout << "Su " << signal << " golden muons, " << reco3_6 << " sono ricostruiti con 3 tracce, con un rapporto del " << (reco3_6/signal)*100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con piu tracce (stesso id), " << more_reco_6 << " sono ricostruiti, con un rapporto del " <<  more_reco_6/signal *100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con 0 tracce di segnale reco, " << reco0_6<< ", con un rapporto del " << reco0_6/signal *100 << "%"<< endl;
//cout << "Su " << signal << " eventi di con 1 sola traccia di segnale reco, " << reco1_6 << ", con un rapporto del " << reco1_6/signal *100 << "%"<< endl;

TCanvas a("a","a",700,700);
d_eff->Draw("E");
a.SaveAs("d_eff_RD.pdf");
d_eff->SaveAs("d_eff_RD.root");

TCanvas b("b","b",700,700);
h_2d->Draw();
h_2d->SaveAs("2D_RD.root");

TCanvas c("c","c",700,700);
theta_mu->Draw("E");
theta_mu->SaveAs("theta_mu_RD.root");

TCanvas d("d","d",700,700);
theta_e->Draw("E");
theta_e->SaveAs("theta_e_RD.root");


}
