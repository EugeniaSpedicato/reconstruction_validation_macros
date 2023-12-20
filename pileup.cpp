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

        TFile *inputfile = new TFile("/mnt/raid10/DATA/espedica/fairmu/dataReconstruction_3234-3235_new12.root");
//dataReconstruction_3234-3235_new12.root");
//dataReconstruction_3234-3235_int_new_nobend.root");
//dataReconstruction_3234-3235_new12.root");
//TRMesmer_3cm.root");
//dataReconstruction_3234-3235_filtered.root");//dataReconstruction_try_sigma.root");
//TRMesmer_3cm.root");
//dataReconstruction_try.root");
        TTree* cbmsim = (TTree*) inputfile->Get("cbmsim");

        MUonERecoOutput *ReconstructionOutput = 0;


   cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);

TH2D *h_2d=new TH2D("h_2d","Electron VS muone angle with |aco|<0.5 and chi2 cut" ,700,0.,0.07,90,0.,0.003);
TH2D *h_2d_vrtx=new TH2D("h_2d_vrtx","Electron VS muone angle from best vrtx with |aco|<0.5 and chi2 cut",700,0.,0.07,90,0.,0.003);

TH2D *h_2d_chi=new TH2D("h_2d_chi","Electron VS muone angle chi2_vrtx cut" ,700,0.,0.07,90,0.,0.003);
TH2D *h_2d_vrtx_chi=new TH2D("h_2d_vrtx_chi","Electron VS muone angle from best vrtx chi2_vrtx cut",700,0.,0.07,90,0.,0.003);

TH1D* h_aco=new TH1D("aco","Acoplanarity of muone+electron from sig minbias pre vrtx",600,-3.14,3.14);
TH1D* h_aco_v=new TH1D("aco_v","Acoplanarity of muone+electron from sig minbias after vrtx",600,-3.14,3.14);

TH1D* h_aco_chi=new TH1D("acoc","Acoplanarity of muone+electron from sig minbias pre vrtx with chi2<20 cut",600,-3.14,3.14);
TH1D* h_aco_chi_v=new TH1D("aco_vc","Acoplanarity of muone+electron from sig minbias after vrtx with chi2<20 cut",600,-3.14,3.14);

TH1D* h_stubs=new TH1D("h_stubs","Number of stubs in strange events, where 12 are used", 50,10,60);
TH1D* h_stubs0=new TH1D("h_stubs0","Number of stubs in strange events, where 12 are used STATION0", 50,0,50);
TH1D* h_stubs1=new TH1D("h_stubs1","Number of stubs in strange events, where 12 are used STATION1", 50,0,50);



// X or Y position on the track at a given Z
                auto pos_on_track = [](double q, double m, double z){return (q + m*z);};

double all=0.;
double e_c=0.;
double mu_c=0.;
double mu_in_c=0.;
double v=0.;
double mu_st0=0.;
double elastic=0.;

for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
		cbmsim->GetEntry(i);
		if(i%1000 == 0) cout<<"Entry "<<i<<endl;

std::vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
//std::vector<MUonERecoOutputHit> stubs=ReconstructionOutput->ReconstructedHits();
//vector<MUonERecoOutputVertex> vrtx = ReconstructionOutput->reconstructedVertices();
//vector<MUonERecoOutputAdaptiveFitterVertex> vrtx = ReconstructionOutput->adaptiveFitterVertices();

double chi=99.;
double Z_sig=0.;

//cout << " ReconstructedHitsMultiplicity " << ReconstructionOutput->reconstructedHitsMultiplicity() << endl;
MUonERecoOutputVertex vrtx = ReconstructionOutput->bestVertex();
MUonERecoOutputTrack muin=vrtx.incomingMuon();
MUonERecoOutputTrack mu=vrtx.outgoingMuon();
MUonERecoOutputTrack e=vrtx.outgoingElectron();

TVector3 th_muin_v(muin.xSlope(),muin.ySlope(),1.0);
th_muin_v.Unit();

chi=vrtx.chi2perDegreeOfFreedom();

double chi2_muin_vrtx=vrtx.incomingMuon().chi2perDegreeOfFreedom();
double chi2_e_vrtx=vrtx.outgoingElectron().chi2perDegreeOfFreedom();
double chi2_mu_vrtx=vrtx.outgoingMuon().chi2perDegreeOfFreedom();

double x_v=0.;
double y_v=0.;
double z_v=0.;

if(chi!=99 and chi!=0) {//cout << "Best vertex chi " << chi << endl;
                        x_v=vrtx.x(); y_v=vrtx.y(); z_v=vrtx.z();}

    int sec0=0; int sec1=0;

TVector3 thin1;
TVector3 electron;
TVector3 muon;
std::array<TVector3,2> thin2_v;
std::array<double,2> theta;

double th_muin, th_mu, th_el;
double acoplanarity,acoplanarity_v;
double x0[2];
double y0[2];
double mx[2];
double my[2];
int stubs_muin=0.;
int stubs[2]={0};

std::array<double,2> chi_min;

    for(int j=0; j<tracks.size();j++)
    {
        if(tracks.at(j).sector()==0) sec0++;
        if(tracks.at(j).sector()==1) sec1++;
	}


int index=0;
double chi2_muin;
    for(int j=0; j<tracks.size();j++)
    {
        std::vector<MUonERecoOutputHit> hits_=tracks.at(j).hits();
	if(tracks.at(j).sector()==0 and sec0==1){
        stubs_muin=hits_.size();
	double th_inx=tracks.at(j).xSlope();
	double th_iny=tracks.at(j).ySlope();
	thin1.SetXYZ(th_inx,th_iny,1.0);
	thin1=thin1.Unit();
	chi2_muin=tracks.at(j).chi2perDegreeOfFreedom();
	th_muin=thin1.Theta(); 
	}

	if(tracks.at(j).sector()==1 and sec1==2){
	double th_inx=tracks.at(j).xSlope();
	double th_iny=tracks.at(j).ySlope();
	TVector3 p;
	p.SetXYZ(th_inx,th_iny,1.0);
	p=p.Unit();
	thin2_v.at(index)=p;
	chi_min.at(index)=tracks.at(j).chi2perDegreeOfFreedom();
	theta.at(index)=thin1.Angle(p);
	mx[index]=th_inx;
        my[index]=th_iny;
        x0[index]=tracks.at(j).x0();
        y0[index]=tracks.at(j).y0();
        stubs[index]=hits_.size();
	index++;
		 }

	}

double chi2_e,chi2_mu,id_e,id_mu;

if(sec0==1 and stubs_muin==6) mu_st0++;

if(sec0==1 and sec1==2 and stubs[0]==6 and stubs[1]==6 and stubs_muin==6){
                h_2d->Fill(th_el,th_mu);
                h_2d_vrtx->Fill(vrtx.electronTheta(),vrtx.muonTheta());

	if(theta.at(0)>theta.at(1)){id_e=1; id_mu=2; electron=thin2_v.at(0); muon=thin2_v.at(1); th_el=theta.at(0);
				    th_mu=theta.at(1); chi2_e=chi_min.at(0); chi2_mu=chi_min.at(1);}
	else{id_e=2; id_mu=1; electron=thin2_v.at(1); muon=thin2_v.at(0); th_el=theta.at(1); th_mu=theta.at(0); chi2_e=chi_min.at(1); chi2_mu=chi_min.at(0);}

						double dotProduct = muon.Dot(electron);
                                                TVector3 crossProduct = muon.Cross(electron);
                                                double T = thin1.Dot(crossProduct);
                                                TVector3 im= thin1.Cross(muon);
                                                TVector3 ie= thin1.Cross(electron);
                                                T = T>0? 1:-1;
                                                acoplanarity= T*(TMath::Pi()- acos( ((im).Dot(ie))/(im.Mag()*ie.Mag()) ));

	TVector3 p_muin(muin.xSlope(),muin.ySlope(),1.0);
        TVector3 p_mu(mu.xSlope(),mu.ySlope(),1.0);
        TVector3 p_e(e.xSlope(),e.ySlope(),1.0);


                                                double dotProduct_v = p_mu.Dot(p_e);
                                                TVector3 crossProduct_v = p_mu.Cross(p_e);
                                                double T_v = p_muin.Dot(crossProduct);
                                                TVector3 im_v= p_muin.Cross(p_mu);
                                                TVector3 ie_v= p_muin.Cross(p_e);
                                                T_v = T_v>0? 1:-1;
                                                acoplanarity_v= T_v*(TMath::Pi()- acos( ((im_v).Dot(ie_v))/(im_v.Mag()*ie_v.Mag()) ));



if(th_el<0.032){
 h_aco->Fill(acoplanarity);
 h_aco_v->Fill(acoplanarity_v);

if(chi>20 and th_mu<0.02){
elastic++;
 h_2d_chi->Fill(th_el,th_mu);
 h_2d_vrtx_chi->Fill(vrtx.electronTheta(),vrtx.muonTheta());
 h_aco_chi->Fill(acoplanarity);
 h_aco_chi_v->Fill(acoplanarity_v);
cout << "event " << i << endl;
//cout << " ReconstructedHitsMultiplicity " << ReconstructionOutput->reconstructedHitsMultiplicity() << endl;
vector<MUonERecoOutputHit> stubs=ReconstructionOutput->reconstructedHits();
double nstubs0=0.;
double nstubs1=0.;
h_stubs->Fill(ReconstructionOutput->reconstructedHitsMultiplicity());
for(int s=0; s<stubs.size(); s++){
if(stubs.at(s).stationID()==0)nstubs0++;
if(stubs.at(s).stationID()==1)nstubs1++;
}
h_stubs0->Fill(nstubs0);
h_stubs1->Fill(nstubs1);

			}//chi2
		}//th_mu<32mrad
	}//if(sec0==1 and sec1==2)

}

cout << mu_st0 << " muons reconstructed in station 0 over " << cbmsim->GetEntries() << ": " << 100*mu_st0/cbmsim->GetEntries() << "%" << endl;
cout << elastic << " elastic events over " << mu_st0 << ": " << 100*elastic/mu_st0 << "%" << endl;


    TCanvas n("n","n",700,700);
    n.Divide(2,2);
    n.cd(1);
    h_2d->Draw("COLZ");
    n.cd(2);
    h_2d_vrtx->Draw("COLZ");
    n.cd(3);
    h_2d_chi->Draw();
    n.cd(4);
    h_2d_vrtx_chi->Draw();
    n.SaveAs("2d_vtx_pileup.pdf");

    TCanvas n1("n1","n1",700,700);
    n1.Divide(1,2);
    n1.cd(1);
h_stubs0->Draw();
h_stubs1->SetLineColor(kRed);
h_stubs1->Draw("same");
    n1.cd(2);
h_stubs->Draw();
    n1.SaveAs("nstubs_pileup.pdf");


}


