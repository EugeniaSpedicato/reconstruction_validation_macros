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

        TFile *inputfile = new TFile("/mnt/raid10/DATA/espedica/fairmu/TRMesmer_3cm_600k.root");//dataReconstruction_3234-3235_filtered.root");
//TRMesmer_3cm_600k.root");
//dataReconstruction_3234-3235_filtered.root");
//dataReconstruction_run3234_3235_1.root");//TRMesmer_3cm.root");//dataReconstruction_run3234_3235_1.root");
//run3234_3235_1.root");
//dataReconstruction_3234-3235_new.root");
//dataReconstruction_try_sigma.root");
//TRMesmer_3cm.root");
        TTree* cbmsim = (TTree*) inputfile->Get("cbmsim");

        MUonERecoOutput *ReconstructionOutput = 0;


   cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);

string title;


        std::vector<TH1D*> residual_mu_trk(4);
	residual_mu_trk.at(0)=new TH1D("residual_mu_trk0","Residual mu_out X track Vs mu_in X track",100,-0.05,0.05);
        residual_mu_trk.at(1)=new TH1D("residual_mu_trk1","Residual mu_out Y track Vs mu_in Y track",100,-0.05,0.05);
        residual_mu_trk.at(2)=new TH1D("residual_mu_trk2","Pull residual mu_out X track Vs mu_in X track",100,-5.,5.);//150,-0.3,0.3);
        residual_mu_trk.at(3)=new TH1D("residual_mu_trk3","Pull residual mu_out Y track Vs mu_in Y track",100,-5.,5.);//150,-0.3,0.3);


        std::vector<TH1D*> residual_e_trk(4);
        residual_e_trk.at(0)=new TH1D("residual_e_trk0","Residual e_out X track Vs mu_in X track",200,-0.8,0.8);
        residual_e_trk.at(1)=new TH1D("residual_e_trk1","Residual e_out Y track Vs mu_in Y track",200,-0.8,0.8);
        residual_e_trk.at(2)=new TH1D("residual_e_trk2","Pull residual e_out X track Vs mu_in X track",300,-30.,30.);//250,-1.,1.);
        residual_e_trk.at(3)=new TH1D("residual_e_trk3","Pull residual e_out Y track Vs mu_in Y track",300,-30.,30.);//250,-1.,1.);

        std::vector<TH1D*> error(6);
	error.at(0)=new TH1D("errorXmuin","Error x0 for incoming muon", 100,-0.05, 0.05);
        error.at(1)=new TH1D("errorYmuin","Error y0 for incoming muon", 100,-0.05, 0.05);
        error.at(2)=new TH1D("errorXmuout","Error x0 for outgoing muon", 100,-0.05, 0.05);
        error.at(3)=new TH1D("errorYmuout","Error y0 for outgoing muon", 100,-0.05, 0.05);
        error.at(4)=new TH1D("errorXel","Error x0 for outgoing electron", 100,-0.05, 0.05);
        error.at(5)=new TH1D("errorYel","Error y0 for outgoing electron", 100,-0.05, 0.05);


TH2D *h_2d_mu=new TH2D("h_2dmu","Electron VS muone angle with aco, chi muon coordinate cut" ,350,0.,0.035,90,0.,0.003);
TH2D *h_2d_all=new TH2D("h_2dall","Electron VS muone angle with aco, chi and mu, e coordinate cut" ,350,0.,0.035,90,0.,0.003);
TH2D *h_2d_chi=new TH2D("h_2dchi","Electron VS muone angle with aco and chi2 cut" ,350,0.,0.035,90,0.,0.003);
TH2D *h_2d_aco=new TH2D("h_2daco","Electron VS muone angle with aco cut" ,350,0.,0.035,90,0.,0.003);
TH2D *h_2d=new TH2D("h_2d","Electron VS muone angle (1 reco in- 2 reco out)" ,350,0.,0.035,90,0.,0.003);

TH2D *h_imp=new TH2D("h_imp", "Target impact point mu_in when there is an int", 500,-5.,5., 500,-5.,5.);

        TH1F* histo_nstubs_BX=new TH1F("histo_nstubs_BX", "nstubs per BX",50,0,50);

// X or Y position on the track at a given Z
                auto pos_on_track = [](double q, double m, double z){return (q + m*z);};

for(Long64_t i = 0; i <cbmsim->GetEntries(); i++) {
//1000000; i++) {//5000000; i++) {//cbmsim->GetEntries(); i++) {
		cbmsim->GetEntry(i);
		if(i%1000 == 0) cout<<"Entry "<<i<<endl;

vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();

double chi=99.;

MUonERecoOutputVertex vrtx = ReconstructionOutput->bestVertex();

chi=vrtx.chi2perDegreeOfFreedom();

MUonERecoOutputTrack muin=vrtx.incomingMuon();
MUonERecoOutputTrack mu=vrtx.outgoingMuon();
MUonERecoOutputTrack e=vrtx.outgoingElectron();

double chi2_muin_vrtx=vrtx.incomingMuon().chi2perDegreeOfFreedom();
double chi2_e_vrtx=vrtx.outgoingElectron().chi2perDegreeOfFreedom();
double chi2_mu_vrtx=vrtx.outgoingMuon().chi2perDegreeOfFreedom();

double x_v=0.;
double y_v=0.;
double z_v=0.;
if(chi!=99 and chi!=0) {x_v=vrtx.x(); y_v=vrtx.y(); z_v=vrtx.z();}

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
double x0_err[2];
double y0_err[2];
double mx[2];
double my[2];

UShort_t bx=0.;
UInt_t sid=0.;

std::array<double,2> chi_min;


    for(int j=0; j<tracks.size();j++)
    {
        if(tracks.at(j).sector()==0) sec0++;
        if(tracks.at(j).sector()==1) sec1++;
	}


int index=0;
double chi2_muin;
int stubs_muin=0.;
int stubs[2]={0};
int stubs_mu=0.;
int stubs_e=0.;
int stub_0=0.;
int stub_1=0.;
double x0_err_muin, y0_err_muin;
ROOT::Math::SMatrix<Double_t, 4> cov_1;
std::vector<ROOT::Math::SMatrix<Double_t, 4>> cov_2;
    for(int j=0; j<tracks.size();j++)
    {
	std::vector<MUonERecoOutputHit> hits_=tracks.at(j).hits();

        if(tracks.at(j).sector()==0) stub_0+=hits_.size();
        if(tracks.at(j).sector()==1) stub_1+=hits_.size();

	if(tracks.at(j).sector()==0 and sec0==1){
	stubs_muin=hits_.size();
	double th_inx=tracks.at(j).xSlope();
	double th_iny=tracks.at(j).ySlope();
	thin1.SetXYZ(th_inx,th_iny,1.0);
	thin1=thin1.Unit();
	chi2_muin=tracks.at(j).chi2perDegreeOfFreedom();
	th_muin=thin1.Theta();
	x0_err_muin=tracks.at(j).x0Error();
	y0_err_muin=tracks.at(j).y0Error();
	bx=hits_.at(0).bx();
	sid=hits_.at(0).superID();
	cov_1=tracks.at(j).linearFitCovarianceMatrix();
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
        x0_err[index]=tracks.at(j).x0Error();
        y0_err[index]=tracks.at(j).y0Error();
	stubs[index]=hits_.size();
	index++;
	cov_2.push_back(tracks.at(j).linearFitCovarianceMatrix());

		 }

	}

if(bx!=0 and sid!=0)
{ if(i<100)cout << "primi Bx " << bx << ", sid " << sid << endl;
  if(i>cbmsim->GetEntries()-100)cout << "utlimi Bx " << bx << ", sid " << sid << endl;
}

double chi2_e,chi2_mu;
int id_e,id_mu;
double x0_err_mu, x0_err_e, y0_err_mu, y0_err_e;
ROOT::Math::SMatrix<Double_t, 4> cov_mu,cov_e;

if(sec0==1 and sec1==2 and stubs_muin==6 and chi2_muin<5){//stubs[0]==6 and stubs[1]==6 and stubs_muin==6){


	if(theta.at(0)>theta.at(1)){id_e=1; id_mu=2; electron=thin2_v.at(0); muon=thin2_v.at(1); stubs_e=stubs[0]; stubs_mu=stubs[1];
				    th_el=theta.at(0); th_mu=theta.at(1); chi2_e=chi_min.at(0); chi2_mu=chi_min.at(1);
					x0_err_mu=x0_err[1]; x0_err_e=x0_err[0] ; y0_err_mu=y0_err[1] ; y0_err_e=y0_err[0];
					cov_mu=cov_2.at(1); cov_e=cov_2.at(0);}
	else{id_e=2; id_mu=1; electron=thin2_v.at(1); muon=thin2_v.at(0);  stubs_e=stubs[1]; stubs_mu=stubs[0];
	     th_el=theta.at(1); th_mu=theta.at(0); chi2_e=chi_min.at(1); chi2_mu=chi_min.at(0);
                                        x0_err_mu=x0_err[0]; x0_err_e=x0_err[1] ; y0_err_mu=y0_err[0] ; y0_err_e=y0_err[1];
                                        cov_mu=cov_2.at(0); cov_e=cov_2.at(1);}

						double dotProduct = muon.Dot(electron);
                                                TVector3 crossProduct = muon.Cross(electron);
                                                double T = thin1.Dot(crossProduct);
                                                TVector3 im= thin1.Cross(muon);
                                                TVector3 ie= thin1.Cross(electron);
                                                T = T>0? 1:-1;
                                                acoplanarity= T*(TMath::Pi()- acos( ((im).Dot(ie))/(im.Mag()*ie.Mag()) ));

double x_t[3]={0.};
double y_t[3]={0.};

    for(int j=0; j<tracks.size();j++){
std::vector<MUonERecoOutputHit> hits_=tracks.at(j).hits();

	  if(j==0) {x_t[0]=pos_on_track(tracks.at(j).x0(),tracks.at(j).xSlope(),z_v);
		    y_t[0]=pos_on_track(tracks.at(j).y0(),tracks.at(j).ySlope(),z_v);

			}
          if(j==id_mu) {x_t[1]=pos_on_track(tracks.at(j).x0(),tracks.at(j).xSlope(),z_v);
                        y_t[1]=pos_on_track(tracks.at(j).y0(),tracks.at(j).ySlope(),z_v);
			}

          if(j==id_e) {x_t[2]=pos_on_track(tracks.at(j).x0(),tracks.at(j).xSlope(),z_v);
                       y_t[2]=pos_on_track(tracks.at(j).y0(),tracks.at(j).ySlope(),z_v);
                       }

 }

h_imp->Fill(x_t[0],y_t[0]);

      //if(x_t[0]>-1.3 and x_t[0]<0.7 and y_t[0]>-0.5 and y_t[0]<0.5){
	if(x_t[0]>-1. and x_t[0]<1. and y_t[0]>-1. and y_t[0]<1.){
	     h_2d->Fill(th_el,th_mu);

      if(abs(acoplanarity)<1) h_2d_aco->Fill(th_el,th_mu);

      if(abs(acoplanarity)<1 and chi2_e<5 and chi2_mu<5 and th_el<0.032){// and th_mu>0.0002){//chi2_e<5

        h_2d_chi->Fill(th_el,th_mu);

        cout << "------ event " << i << " -------" << endl;

	double err1_x,err1_y,err2_x,err2_y,err3_x,err3_y;

	err1_x=sqrt((z_v*z_v*cov_1[2][2]) + cov_1[0][0] + 2*z_v*cov_1[0][2]);
        err1_y=sqrt((z_v*z_v*cov_1[3][3]) + cov_1[1][1] + 2*z_v*cov_1[1][3]);

        err2_x=sqrt((z_v*z_v*cov_mu[2][2]) + cov_mu[0][0] + 2*z_v*cov_mu[0][2]);
        err2_y=sqrt((z_v*z_v*cov_mu[3][3]) + cov_mu[1][1] + 2*z_v*cov_mu[1][3]);

        err3_x=sqrt((z_v*z_v*cov_e[2][2]) + cov_e[0][0] + 2*z_v*cov_e[0][2]);
        err3_y=sqrt((z_v*z_v*cov_e[3][3]) + cov_e[1][1] + 2*z_v*cov_e[1][3]);

			error.at(0)->Fill(err1_x);
                        error.at(1)->Fill(err1_y);

                        error.at(2)->Fill(err2_x);
                        error.at(3)->Fill(err2_y);

                        error.at(4)->Fill(err3_x);
                        error.at(5)->Fill(err3_y);

                        residual_mu_trk.at(0)->Fill(x_t[1]-x_t[0]);
                        residual_mu_trk.at(1)->Fill(y_t[1]-y_t[0]);
                        residual_mu_trk.at(2)->Fill( (x_t[1]-x_t[0])/( sqrt(err2_x*err2_x + err1_x*err1_x) ) );
                        residual_mu_trk.at(3)->Fill( (y_t[1]-y_t[0])/( sqrt(err2_y*err2_y + err1_y*err1_y) ) );
                        residual_e_trk.at(0)->Fill(x_t[2]-x_t[0]);
                        residual_e_trk.at(1)->Fill(y_t[2]-y_t[0]);
                        residual_e_trk.at(2)->Fill( (x_t[2]-x_t[0])/( sqrt(err3_x*err3_x + err1_x*err1_x) ) );
                        residual_e_trk.at(3)->Fill( (y_t[2]-y_t[0])/( sqrt(err3_y*err3_y + err1_y*err1_y) ) );

if( abs(x_t[1]-x_t[0])<0.02 and abs(y_t[1]-y_t[0])<0.02 and abs(x_t[2]-x_t[0])<0.3 and abs(y_t[2]-y_t[0])<0.3) h_2d_all->Fill(th_el,th_mu);
if( abs(x_t[1]-x_t[0])<0.02 and abs(y_t[1]-y_t[0])<0.02 )h_2d_mu->Fill(th_el,th_mu);

			}//abs(acoplanarity)<1
		}//abs(x_t[0])<1 and abs(y_t[0])<1
	}//if(sec0==1 and sec1==2)
if(stub_0>=5 and stub_1>=5 and (stub_1-stub_0)>=5){histo_nstubs_BX->Fill(stub_0+stub_1);}

  }


    TCanvas n("n","n",700,700);
    n.Divide(2,3);
    n.cd(1);
    h_2d->Draw();
    n.cd(2);
    h_2d_aco->Draw();
    n.cd(3);
    h_2d_chi->Draw();
    n.cd(4);
    h_2d_mu->Draw();
    n.cd(5);
    h_2d_all->Draw();
    //histo_nstubs_BX->Draw();
    n.SaveAs("2D_studio.pdf");

    TCanvas n1("n1","n1",700,700);
    n1.Divide(2,3);
   for(int m=0; m<6; m++){
    n1.cd(m+1);
    error.at(m)->Draw();
   }
    n1.SaveAs("error.pdf");

TCanvas n3("n3","n3",1000,1000);
n3.Divide(2,4);
n3.cd(1);
residual_mu_trk[0]->Draw();
//gPad->SetLogy(); 
n3.cd(3);
residual_mu_trk[1]->Draw();
//gPad->SetLogy();
n3.cd(5);
residual_e_trk[0]->Draw();
//gPad->SetLogy();
n3.cd(7);
residual_e_trk[1]->Draw();
//gPad->SetLogy();
n3.cd(2);
residual_mu_trk[2]->Draw();
//gPad->SetLogy(); 
n3.cd(4);
residual_mu_trk[3]->Draw();
//gPad->SetLogy();
n3.cd(6);
residual_e_trk[2]->Draw();
//gPad->SetLogy();
n3.cd(8);
residual_e_trk[3]->Draw();
//gPad->SetLogy();

n3.SaveAs("pull_res.pdf");


TCanvas n4 ("n4","n4",700,700);
h_imp->Draw();
n4.SaveAs("h_imp.pdf");


}


