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

double zcrosslines(double x0_0, double y0_0, double tanx_0, double tany_0,double x0_1, double y0_1, double tanx_1, double tany_1){

  double invn0norm=1./sqrt(tanx_0*tanx_0+tany_0*tany_0+1);
  double invn1norm=1./sqrt(tanx_1*tanx_1+tany_1*tany_1+1);
  double cs=(tanx_0*tanx_1+tany_0*tany_0+1)*invn0norm*invn1norm;
  double rr[]={x0_0-x0_1,y0_0-y0_1};

  double a0 = -1/(1-cs*cs)* ( (tanx_0*invn0norm-cs*tanx_1*invn1norm)*rr[0]+(tany_0*invn0norm-cs*tany_1*invn1norm)*rr[1] );
  double a1 = -1/(1-cs*cs)* ( (cs*tanx_0*invn0norm-tanx_1*invn1norm)*rr[0]+(cs*tany_0*invn0norm-tany_1*invn1norm)*rr[1] );

  //double xav = (t0.x0+t1.x0)/2. + (t0.tanx*invn0norm * a0 + + t1.tanx*invn1norm * a1)/2.
  //double yav = (t0.y0+t1.y0)/2. + (t0.tany*invn0norm * a0 + + t1.tany*invn1norm * a1)/2.
  double zav = (1.*invn0norm * a0 + 1.*invn1norm * a1)/2.;

  return zav;
}


TH1D* h_aco=new TH1D("aco","Acoplanarity of muone+electron from sig minbias",600,-3.14,3.14);
TH1D* h_chi=new TH1D("h_chi","chi2 of the kinematic vertex", 1000,0,1000);
TH1D* h_chi_tr1=new TH1D("h_chi1","chi2 of track1 e out", 100,0,100);
TH1D* h_chi_tr2=new TH1D("h_chi2","chi2 of track2 mu out", 100,0,100);
TH1D* h_chi_tr0=new TH1D("h_chi0","chi2 of track2 mu in", 100,0,100);

void RealDataAnalyzer(){

        TFile *inputfile = new TFile("/mnt/raid10/DATA/espedica/fairmu/dataReconstruction_3234-3235_try_shift.root");
//TRMesmer_3cm.root");
//dataReconstruction_3234-3235_filtered.root");//dataReconstruction_try_sigma.root");
//TRMesmer_3cm.root");
//dataReconstruction_try.root");
        TTree* cbmsim = (TTree*) inputfile->Get("cbmsim");

//TChain * cbmsim = new TChain("cbmsim");
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/dataReconstruction_3238_3239_new_alin.root");
//Add("TRMesmer_3cm.root");
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/dataReconstruction_3238_3239_new_alin.root");
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/dataReconstruction_merged_3238_3239_0hit_2.root");
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/dataReconstruction_merged_3238_3239_0hit_3.root");
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/dataReconstruction_merged_3238_3239_second_outlier_1.root");
//3232_3233_1_4M_outlier.root");
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/dataReconstruction_merged_3232_3233_small4M_outlier.root");

        MUonERecoOutput *ReconstructionOutput = 0;


   cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);

    TH2D *h_2d=new TH2D("h_2d","Electron VS muone angle" ,700,0.,0.07,90,0.,0.003);
    TH1D *h_Z=new TH1D("h_Z","Z of the vrtx found",800,800.,1000.);
    TH1D *h_mu_in=new TH1D("h_mu_in","Muon direction first station",100,-0.03,0.03);
    TH1D *h_prob_e=new TH1D("h_prob_e","Electron chi2 prob", 200,-0.5,1.5);
    TH1D *h_prob_mu=new TH1D("h_prob_mu","Muone chi2 prob", 200,-0.5,1.5);
    TH1D *h_prob_muin=new TH1D("h_prob_muin","Muone incoming chi2 prob", 200,-0.5,1.5);
    TH1D *h_chi_muin_vrtx=new TH1D("h_chi_muin_vrtx","Incoming muon chi2 best vrtx",1000,0.,1000.);
    TH1D *h_chi_e_vrtx=new TH1D("h_chi_e_vrtx","Outgoing electron chi2 best vrtx",1000,0.,1000.);
    TH1D *h_chi_mu_vrtx=new TH1D("h_chi_mu_vrtx","Outgoing muon chi2 best vrtx",1000,0.,1000.);
    TH2D *h_2d_vrtx=new TH2D("h_2d_vrtx","Electron VS muone angle from best vrtx",700,0.,0.07,90,0.,0.003);

// X or Y position on the track at a given Z
                auto pos_on_track = [](double q, double m, double z){return (q + m*z);};

double all=0.;
double e_c=0.;
double mu_c=0.;
double mu_in_c=0.;
double v=0.;

for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
		cbmsim->GetEntry(i);
		if(i%1000 == 0) cout<<"Entry "<<i<<endl;

vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
//vector<MUonERecoOutputVertex> vrtx = ReconstructionOutput->reconstructedVertices();
//vector<MUonERecoOutputAdaptiveFitterVertex> vrtx = ReconstructionOutput->adaptiveFitterVertices();

double chi=99.;
double Z_sig=0.;

/*cout <<"vrtx.size() " << vrtx.size() << endl;

  for(int j=0; j<vrtx.size();j++)
        {
	 if(j==0) {chi=vrtx.at(j).chi2(); h_chi->Fill(chi);
         cout << "vertex chi " << chi  << " at Z_sig " << vrtx.at(j).z() << endl;
	 if(vrtx.at(j).chi2()<50){h_Z->Fill(vrtx.at(j).z());}
	Z_sig=vrtx.at(j).z();
		}
        }*/

MUonERecoOutputVertex vrtx = ReconstructionOutput->bestVertex();
MUonERecoOutputTrack muin=vrtx.incomingMuon();
MUonERecoOutputTrack mu=vrtx.outgoingMuon();
MUonERecoOutputTrack e=vrtx.outgoingElectron();

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
double acoplanarity;
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
	th_muin=thin1.Theta(); h_mu_in->Fill(th_muin);
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

if(sec0==1 and sec1==2 and stubs[0]==6 and stubs[1]==6 and stubs_muin==6){


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

	double Z=zcrosslines(x0[0],y0[0],mx[0],my[0],x0[1],y0[1],mx[1],my[1]);

 double prob_e=TMath::Prob(chi2_e, 2);
 double prob_mu=TMath::Prob(chi2_mu, 2);
 double prob_muin=TMath::Prob(chi2_muin, 2);

if(th_el<0.032){
 h_aco->Fill(acoplanarity);
 h_prob_e->Fill(TMath::Prob(chi2_e, 2));
 h_prob_mu->Fill(TMath::Prob(chi2_mu, 2));
 h_prob_muin->Fill(prob_muin);
 h_chi_tr1->Fill(chi2_e);
 h_chi_tr2->Fill(chi2_mu);
 h_chi_tr0->Fill(chi2_muin);
 h_chi_muin_vrtx->Fill(chi2_muin_vrtx);
 h_chi_e_vrtx->Fill(chi2_e_vrtx);
 h_chi_mu_vrtx->Fill(chi2_mu_vrtx);


all++;
if(chi>20)v++;
if(chi2_e_vrtx>20)e_c++;
if(chi2_muin_vrtx>20)mu_in_c++;
if(chi2_mu_vrtx>20)mu_c++;

cout << "--------New event--------" << endl;
cout << "Incoming Muon chi2 before vrtx and after vrtxing: " << endl;
cout << "	" << chi2_muin << "  ->   " << chi2_muin_vrtx << endl;
cout << "Outgoing Muon before vrtx and after vrtxing : "  << endl;
cout << "       " << chi2_mu <<  "  ->   " << chi2_mu_vrtx << endl;
cout << "Outgoing Electron before vrtx and after vrtxing : "  << endl;
cout << "       " << chi2_e <<  "  ->   " << chi2_e_vrtx << endl;


double deltaX[3];
double deltaY[3];
    for(int j=0; j<tracks.size();j++){

          if(j==0) {    double x0=pos_on_track(tracks.at(j).x0(),tracks.at(j).xSlope(),z_v);
                        double y0=pos_on_track(tracks.at(j).y0(),tracks.at(j).ySlope(),z_v);
                        double x1=pos_on_track(muin.x0(),muin.xSlope(),0.);
                        double y1=pos_on_track(muin.y0(),muin.ySlope(),0.);

                        deltaX[0]=x0-x_v;
                        deltaY[0]=y0-y_v;
                        }
          if(j==id_mu) {
                        double x0=pos_on_track(tracks.at(j).x0(),tracks.at(j).xSlope(),z_v);
                        double y0=pos_on_track(tracks.at(j).y0(),tracks.at(j).ySlope(),z_v);
                        double x1=pos_on_track(mu.x0(),mu.xSlope(),0.);
                        double y1=pos_on_track(mu.y0(),mu.ySlope(),0.);

                        deltaX[1]=x0-x_v;
                        deltaY[1]=y0-y_v;
                        }
          if(j==id_e) {
                        double x0=pos_on_track(tracks.at(j).x0(),tracks.at(j).xSlope(),z_v);
                        double y0=pos_on_track(tracks.at(j).y0(),tracks.at(j).ySlope(),z_v);
                        double x1=pos_on_track(e.x0(),e.xSlope(),0.);
                        double y1=pos_on_track(e.y0(),e.ySlope(),0.);

			deltaX[2]=x0-x_v;
			deltaY[2]=y0-y_v;
                        }


 }

//th_mu>0.0002
	if(abs(acoplanarity)<1){/*and deltaX[0]<0.129 and deltaX[0]>-0.071
						and deltaX[1]<0.129 and deltaX[1]>-0.071
						and deltaX[2]<0.2 and deltaX[2]>-0.17
						and deltaY[0]<0.3 and deltaY[0]>0.
						and deltaY[1]<0.3 and deltaY[1]>0.
						and deltaY[2]<0.4 and deltaY[2]>-0.1){*/

	//if(th_el>0.01 and th_mu>0.0004) cout << "event " << i << endl;
                        h_2d->Fill(th_el,th_mu);
                        h_2d_vrtx->Fill(vrtx.electronTheta(),vrtx.muonTheta());
			}
		}//th_mu<32mrad
	}//if(sec0==1 and sec1==2)

}
cout << "Fraction of KF vrtx wth chi2>20 " << v/all << endl;
cout << "Fraction of muon_in tracks post-vrtx wth chi2>20 " << mu_in_c/all << endl;
cout << "Fraction of muon_out tracks post-vrtx wth chi2>20 " << mu_c/all << endl;
cout << "Fraction of electron_out tracks post-vrtx wth chi2>20 " << e_c/all << endl;

    TCanvas n("n","n",700,700);
    n.Divide(1,3);
    n.cd(1);
    h_2d->Draw();
    n.cd(2);
    h_2d_vrtx->Draw();
    n.cd(3);
    h_aco->Draw();
    n.SaveAs("info.pdf");


    TCanvas n1("n1","n1",1000,1000);
    n1.Divide(2,2);
    n1.cd(1);
h_prob_e->SetLineColor(kRed);
h_prob_e->Draw();
h_prob_muin->SetLineColor(kOrange);
h_prob_muin->Draw("same");
h_prob_mu->Draw("same");
gPad->SetLogy();
    n1.cd(2);
h_prob_e->SetLineColor(kRed);
h_prob_e->Draw();
h_prob_muin->SetLineColor(kOrange);
h_prob_muin->Draw("same");
h_prob_mu->Draw("same");
    n1.cd(3);
h_chi_tr1->SetLineColor(kRed);
h_chi_tr2->SetLineColor(kBlue);
h_chi_tr0->SetLineColor(kOrange);
h_chi_tr0->Draw();
h_chi_tr2->Draw("same");
h_chi_tr1->Draw("same");
//gPad->SetLogy();
    n1.cd(4);
h_chi_muin_vrtx->SetLineColor(kOrange);
h_chi_mu_vrtx->SetLineColor(kBlue);
h_chi_e_vrtx->SetLineColor(kRed);
h_chi_e_vrtx->Draw();
h_chi_mu_vrtx->Draw("same");
h_chi_muin_vrtx->Draw("same");
//gPad->SetLogy();
    n1.SaveAs("chi2.pdf");



}


