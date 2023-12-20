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

        TFile *inputfile = new TFile("/mnt/raid10/DATA/espedica/fairmu/dataReconstruction_3234-3235_new12_st0_6_1shared_5M.root");
//minbias_1shared_5M.root");
//dataReconstruction_3234-3235_new12_st0_6_1shared_5M.root");
//dataReconstruction_3234-3235_new12_st0_6.root");
//minbias.root");
//dataReconstruction_3234-3235_new12_st0_6_1shared.root");
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

TH2D *h_2d=new TH2D("h_2d","Electron VS muone angle" ,700,0.,0.07,100,0.,0.005);
TH2D *h_2d_vrtx=new TH2D("h_2d_vrtx","Electron VS muone angle from best vrtx",700,0.,0.07,100,0.,0.005);

TH2D *h_2d_ip=new TH2D("h_2d_ip","Electron VS muone angle IP cut" ,700,0.,0.07,100,0.,0.005);
TH2D *h_2d_vrtx_ip=new TH2D("h_2d_vrtx_ip","Electron VS muone angle from best vrtx IP cut",700,0.,0.07,100,0.,0.005);

TH2D *h_2d_cut=new TH2D("h_2d_cut","Electron VS muone angle IP, th_mu>0.2mrad, |aco|<1, chi2<20 cuts" ,700,0.,0.07,100,0.,0.005);
TH2D *h_2d_vrtx_cut=new TH2D("h_2d_vrtx_cut","Electron VS muone angle from best vrtx IP, |aco|<1, chi2<20 cuts",700,0.,0.07,100,0.,0.005);

TH2D *h_2d_3=new TH2D("h_2d_3","Electron VS muone angle IP, tracks.size()>2" ,700,0.,0.07,100,0.,0.005);
TH2D *h_2d_3_v=new TH2D("h_2d_3_v","Electron VS muone angle from best vrtx IP, tracks.size()>2 " ,700,0.,0.07,100,0.,0.005);
TH2D *h_2d_3_cut=new TH2D("h_2d_3_cut","Electron VS muone angle from best vrtx IP, |aco|<1, chi2<20 cuts, tracks.size()>2",700,0.,0.07,100,0.,0.005);
TH2D *h_2d_3_cut_v=new TH2D("h_2d_3_cut_v","Electron VS muone angle from best vrtx IP, |aco|<1, chi2<20 cuts, tracks.size()>2",700,0.,0.07,100,0.,0.005);

TH1D* h_stubs=new TH1D("h_stubs","Number of stubs in st1 for events with golden muon in st0", 40,0,40);
TH1D* h_stubs0=new TH1D("h_stubs0","Number of st1 stubs for events with 2 tracks reco in st1", 50,0,50);
TH1D* h_stubs1=new TH1D("h_stubs1","Number of st1 stubs for events with >2 tracks reco in st1", 50,0,50);
TH1D* h_stubs2=new TH1D("h_stubs2","Number of st1 stubs for events with 1 tracks reco in st1", 50,0,50);
TH1D* h_stubs3=new TH1D("h_stubs3","Number of st1 stubs for events with 0 tracks reco in st1", 50,0,50);


TH1D* h_distX_p=new TH1D("h_distX_p","Impact Parameter sqrt(d_x^2+d_y^2) track 1 and track 2 at Z=target, theta<5mrad",700,-1.,6.);
TH2D* h_distY_p=new TH2D("h_distY_p","2D plot distance X vs distance Y of track 1 and track 2 at Z=target, theta<5mrad",700,-6.,6.,700,-6.,6.);
TH1D* h_distX=new TH1D("h_distX","Impact Parameter sqrt(d_x^2+d_y^2) track 1 and track 2 at Z=target, theta>5mrad",700,-1.,6.);
TH2D* h_distY=new TH2D("h_distY","2D plot distance X vs distance Y of track 1 and track 2 at Z=target, theta>5mrad",700,-6.,6.,700,-6.,6.);

TH2D* h_pos2d_gld=new TH2D("h_pos2d_gld","impact position of golden muons at Z=target",2000,-5.,5.,2000,-5.,5.);

TH1D* zv=new TH1D("zv","zv",1010,-10.,1000.);

TH1D* h_aco=new TH1D("aco","Acoplanarity of muone+electron with IP cut",600,-3.14,3.14);

TH1D* h_chi=new TH1D("chi2_h","Chi2 vertex with IP cut",1000,0.,1000.);

// X or Y position on the track at a given Z
                auto pos_on_track = [](double q, double m, double z){return (q + m*z);};

double all=0.;
double e_c=0.;
double mu_c=0.;
double mu_in_c=0.;
double v=0.;
double mu_st0=0.;
double mu_st0_a=0.;
double mu_st1=0.;
double mu_st1_s=0.;
double mu_st1_m=0.;
double mu_st1_0=0.;
double mu_st1_2=0.;
double elastic=0.;
double elastic3=0.;
double z_fix=912.7;

for(Long64_t i = 0; i < cbmsim->GetEntries(); i++){
		cbmsim->GetEntry(i);
		if(i%1000 == 0) cout<<"Entry "<<i<<endl;

std::vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();

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
if(chi!=99 and chi!=0) {x_v=vrtx.x(); y_v=vrtx.y(); z_v=vrtx.z();}

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
double th_inx,th_iny,x0_in,y0_in;

    for(int j=0; j<tracks.size();j++)
    {
        std::vector<MUonERecoOutputHit> hits_=tracks.at(j).hits();
	if(tracks.at(j).sector()==0 and sec0==1){
        stubs_muin=hits_.size();
	th_inx=tracks.at(j).xSlope();
	th_iny=tracks.at(j).ySlope();
	thin1.SetXYZ(th_inx,th_iny,1.0);
	thin1=thin1.Unit();
	x0_in=tracks.at(j).x0();
	y0_in=tracks.at(j).y0();
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

double posxIN_fake=pos_on_track(x0_in,th_inx,z_v);//fix);
double posyIN_fake=pos_on_track(y0_in,th_iny,z_v);

double posxIN=pos_on_track(x0_in,th_inx,z_fix);
double posyIN=pos_on_track(y0_in,th_iny,z_fix);

vector<MUonERecoOutputHit> stubs_=ReconstructionOutput->reconstructedHits();
double nstubs_0=0.;
double nstubs_1=0.;
for(int s=0; s<stubs_.size(); s++){
if(stubs_.at(s).stationID()==0)nstubs_0++;
if(stubs_.at(s).stationID()==1)nstubs_1++;
}

if(sec0==1 and stubs_muin==6) {mu_st0_a++;}

if(sec0==1 and stubs_muin==6 and abs(posxIN)<=1.5 and abs(posyIN)<=1.5 and chi2_muin<2){
h_pos2d_gld->Fill(posxIN,posyIN);
mu_st0++;
h_stubs->Fill(nstubs_1);
if(sec1==0) mu_st1_0++;
if(sec1==1 and (nstubs_1==6 or nstubs_1==5) ) mu_st1++;
if(sec1==1 and (nstubs_1!=6 and nstubs_1!=5) ) mu_st1_s++;
	if(sec1==2){//and stubs[0]==6 and stubs[1]==6){


	mu_st1_2++;
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
                                                double T_v = p_muin.Dot(crossProduct_v);
                                                TVector3 im_v= p_muin.Cross(p_mu);
                                                TVector3 ie_v= p_muin.Cross(p_e);
                                                T_v = T_v>0? 1:-1;
                                                acoplanarity_v= T_v*(TMath::Pi()- acos( ((im_v).Dot(ie_v))/(im_v.Mag()*ie_v.Mag()) ));

double x0_e,y0_e,x0_mu,y0_mu;
    for(int j=0; j<tracks.size();j++){
          if(j==id_mu) {
                        x0_mu=pos_on_track(tracks.at(j).x0(),tracks.at(j).xSlope(),z_fix);
                        y0_mu=pos_on_track(tracks.at(j).y0(),tracks.at(j).ySlope(),z_fix);
                        }
          if(j==id_e) {
                       	x0_e=pos_on_track(tracks.at(j).x0(),tracks.at(j).xSlope(),z_fix);
                        y0_e=pos_on_track(tracks.at(j).y0(),tracks.at(j).ySlope(),z_fix);
                        }
 }
double ip=sqrt( (x0_mu-x0_e)*(x0_mu-x0_e) + (y0_mu-y0_e)*(y0_mu-y0_e) );


 if(th_el>0.015 and th_el<=0.020){

                if(th_mu<0.005 and th_el>0.015 and th_el<=0.020 ){h_distX_p->Fill(ip);
                h_distY_p->Fill(x0_mu-x0_e,y0_mu-y0_e);}
                else{h_distX->Fill(ip);
                h_distY->Fill(x0_mu-x0_e,y0_mu-y0_e);}

		h_stubs0->Fill(nstubs_1);
                h_2d->Fill(th_el,th_mu);
                h_2d_vrtx->Fill(vrtx.electronTheta(),vrtx.muonTheta());

 if(ip<0.8){
  h_2d_ip->Fill(th_el,th_mu);
  h_2d_vrtx_ip->Fill(vrtx.electronTheta(),vrtx.muonTheta());

  h_aco->Fill(acoplanarity);
  h_chi->Fill(chi);

if( abs(acoplanarity)<1 and chi<20) {// and th_mu>0.0002
//cout << "event " << i << endl;
 elastic++;
 h_2d_cut->Fill(th_el,th_mu);
 h_2d_vrtx_cut->Fill(vrtx.electronTheta(),vrtx.muonTheta());
				}//aco and chi and th_mu
			}//ip
		}//th_e<32mrad
	}//if(sec1==2)
else if(sec1>2) {
h_stubs1->Fill(nstubs_1); mu_st1_m++;

int index0=0;
int index1=0;
vector<double> x0_t; x0_t.reserve(5);
vector<double> y0_t; y0_t.reserve(5);
vector<double> x1_t; x1_t.reserve(10);
vector<double> y1_t; y1_t.reserve(10);
vector<int> id; id.reserve(10);

  for(int j=0; j<tracks.size();j++){

			if(tracks.at(j).sector()==0){
			 x0_t.push_back(pos_on_track(tracks.at(j).x0(),tracks.at(j).xSlope(),z_fix));
                         y0_t.push_back(pos_on_track(tracks.at(j).y0(),tracks.at(j).ySlope(),z_fix));
			}
                        if(tracks.at(j).sector()==1){
                         x1_t.push_back(pos_on_track(tracks.at(j).x0(),tracks.at(j).xSlope(),z_fix));
                         y1_t.push_back(pos_on_track(tracks.at(j).y0(),tracks.at(j).ySlope(),z_fix));
			 id.push_back(j);
                        }
	}

vector<int> id_ip1; id_ip1.reserve(10);
vector<int> id_ip2; id_ip2.reserve(10);
vector<double> ip; ip.reserve(3);
int uno=99; int due=99;
for(int x=0; x<x1_t.size(); x++){
	for(int s=0; s<x1_t.size(); s++){
	 double ip_tmp=1.;
	 if(s!=x) ip_tmp=sqrt( (x1_t.at(x)-x1_t.at(s))*(x1_t.at(x)-x1_t.at(s)) + (y1_t.at(x)-y1_t.at(s))*(y1_t.at(x)-y1_t.at(s)) );
	 //questi sono gli ID delle particelle che hanno un IP<250micron
	 //cout << "("<<x<<","<<s<<")" << "ip " << ip_tmp << ", ID: " << id.at(x) << ", " << id.at(s) << endl;
	 if(ip_tmp<0.8 and uno!=s and due!=x){id_ip1.push_back(id.at(x)); id_ip2.push_back(id.at(s)); ip.push_back(ip_tmp); uno=x; due=s;}
		}
	}

if(ip.size()==1){

TVector3 p1,p2;
double theta1,theta2;
double stubs1,stubs2;
for(int j=0; j<tracks.size();j++){

	if(j==id_ip1.at(0)){
        p1.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0);
        p1=p1.Unit();
        theta1=thin1.Angle(p1);
        stubs1=tracks.at(j).hits().size();
	}
        if(j==id_ip2.at(0)){
        p2.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0);
        p2=p2.Unit();
        theta2=thin1.Angle(p2);
        stubs2=tracks.at(j).hits().size();
	}
                 }
         if(theta1>theta2){electron=p1; muon=p2; th_el=theta1;th_mu=theta2;}
         else{electron=p2; muon=p1; th_el=theta2;th_mu=theta1;}

                                                double dotProduct = muon.Dot(electron);
                                                TVector3 crossProduct = muon.Cross(electron);
                                                double T = thin1.Dot(crossProduct);
                                                TVector3 im= thin1.Cross(muon);
                                                TVector3 ie= thin1.Cross(electron);
                                                T = T>0? 1:-1;
                                                acoplanarity= T*(TMath::Pi()- acos( ((im).Dot(ie))/(im.Mag()*ie.Mag()) ));
 if(th_el>0.015 and th_el<=0.020){
 h_2d_3->Fill(th_el,th_mu);
 h_2d_3_v->Fill(vrtx.electronTheta(),vrtx.muonTheta());

if( abs(acoplanarity)<1 and chi<20) {

cout << "event " << i << endl;

 h_2d_3_cut->Fill(th_el,th_mu);
 h_2d_3_cut_v->Fill(vrtx.electronTheta(),vrtx.muonTheta());
 elastic3++;
		}//acoplanarity chi
	}//th_el>0.015 and th_el<=0.020

  }//ip.size()==1
}//sec1>2

else if(sec1<2) {if(sec1==1)h_stubs2->Fill(nstubs_1);
		 if(sec1==0)h_stubs3->Fill(nstubs_1);}
    }//mu_in conditions

}


cout << "Fraction of KF vrtx wth chi2>20 " << v/all << endl;
cout << "Fraction of muon_in tracks post-vrtx wth chi2>20 " << mu_in_c/all << endl;
cout << "Fraction of muon_out tracks post-vrtx wth chi2>20 " << mu_c/all << endl;
cout << "Fraction of electron_out tracks post-vrtx wth chi2>20 " << e_c/all << endl;
cout << mu_st0_a << " events with muons reconstructed in station 0 over " << cbmsim->GetEntries() << ": " << 100*mu_st0_a/cbmsim->GetEntries() << "%" << endl;
cout << mu_st0 << " events with golden muons reconstructed in station 0 and impacting the target |x,y|<1.5cm over " << cbmsim->GetEntries() << ": " << 100*mu_st0/cbmsim->GetEntries() << "%" << endl;
cout << mu_st1_0 << " events with 0 track reconstructed in st1 over golden muons " << mu_st0 << ": " << 100*mu_st1_0/mu_st0 << "%" << endl;
cout << mu_st1 << " events with 1 track reconstructed in st1 (with #stubs ==5,6) over golden muons " << mu_st0 << ": " << 100*mu_st1/mu_st0 << "%" << endl;
cout << mu_st1_s << " events with 1 track reconstructed in st1 (with #stubs !=5,6) over golden muons " << mu_st0 << ": " << 100*mu_st1_s/mu_st0 << "%" << endl;
cout << mu_st1_2 << " events with 2 track reconstructed in st1 over golden muons " << mu_st0 << ": " << 100*mu_st1_2/mu_st0 << "%" << endl;
cout << mu_st1_m << " events with >2 tracks reconstructed in st1 over golden muons " << mu_st0 << ": " << 100*mu_st1_m/mu_st0 << "%" << endl;

//sezioni d'urto LO SENZA TAGLIO IN TH_MU>0.2mrad
/*double s32=250.79;
double s25=150.91;
double s20=94.85;
double s15=51.60;
double s10=21.30;
double s5=4.28;*/

//sezioni d'urto LO con TAGLIO IN TH_MU>0.2mrad
double s32=249.23;
double s25=150.83;
double s20=94.85;
double s15=51.69;
double s10=21.29;
double s5=4.28;

//sezioni d'urto NLO con TAGLIO IN TH_MU>0.2mrad
/*double s32=248.36;
double s25=153.36;
double s20=97.88;
double s15=54.60;
double s10=23.67;
double s5=5.31;*/

double exp=(mu_st0*s32*5.5*1E+23*10*1E-30);
double error =sqrt((elastic/exp)*(1-(elastic/exp))/exp);
double error_0 =sqrt(elastic)/exp;
double error3 =sqrt((elastic3/exp)*(1-(elastic3/exp))/exp);
double error3_0 =sqrt(elastic3)/exp;


cout << elastic << " elastic events over the theoretical expectation " << exp << ": (" << 100*elastic/exp << " +- " << 100*error <<"), er1"<< endl;
cout << elastic << " elastic events over the theoretical expectation " << exp << ": (" << 100*elastic/exp << " +- " << 100*error_0 <<"), er2"<< endl;
cout << "which is the " << 100*elastic/mu_st0 << "% of the golden muons" << endl;
cout << elastic3 << " elastic events when tracks.size()>2 over the theoretical expectation " << exp << ": (" << 100*elastic3/exp << " +- " << 100*error3 <<"), er1"<< endl;
cout << elastic3 << " elastic events when tracks.size()>2 over the theoretical expectation " << exp << ": (" << 100*elastic3/exp << " +- " << 100*error3_0 <<"), er2"<< endl;

cout << "which is the " << 100*elastic3/mu_st0 << "% of the golden muons" << endl;

    TCanvas n("n","n",700,700);
    n.Divide(2,3);
    n.cd(1);
    h_2d->Draw();
    n.cd(2);
    h_2d_vrtx->Draw();
    n.cd(3);
    h_2d_ip->Draw();
    n.cd(4);
    h_2d_vrtx_ip->Draw();
    n.cd(5);
    h_2d_cut->Draw();
    n.cd(6);
    h_2d_vrtx_cut->Draw();
    n.SaveAs("2d_vtx.pdf");
/*
    TCanvas n4("n4","n4",700,700);
    n4.Divide(1,3);
    n4.cd(1);
h_stubs0->SetLineColor(kRed);
h_stubs1->SetLineColor(kOrange);
h_stubs2->SetLineColor(kBlue);
h_stubs3->SetLineColor(kGreen);
h_stubs2->Draw();
h_stubs1->Draw("same");
h_stubs0->Draw("same");
h_stubs3->Draw("same");
gPad->BuildLegend(0.35,0.15,0.35,0.15);
gPad->SetLogy();
    n4.cd(2);
//h_stubs->Draw();
//gPad->SetLogy();
zv->Draw();
gPad->SetLogy();
    n4.cd(3);
h_pos2d_gld->Draw("COLZ");
    n4.SaveAs("nstubs.pdf");

    TCanvas n5("n5","n5",700,700);
    n5.Divide(2,2);
    n5.cd(1);
h_distX_p->Draw();
    n5.cd(2);
h_distY_p->Draw("COLZ");
    n5.cd(3);
h_distX->Draw();
    n5.cd(4);
h_distY->Draw("COLZ");
    n5.SaveAs("distance.pdf");

TCanvas n6("n6","n6",700,700);
n6.Divide(2,2);
n6.cd(1);
h_2d_3->Draw();
n6.cd(2);
h_2d_3_v->Draw();
n6.cd(3);
h_2d_3_cut->Draw();
n6.cd(4);
h_2d_3_cut_v->Draw();
n6.SaveAs("2d_more.pdf");

TCanvas n7("n7","n7",700,700);
n7.Divide(1,2);
n7.cd(1);
h_chi->Draw();
gPad->SetLogy();
n7.cd(2);
h_aco->Draw();
n7.SaveAs("acochi.pdf");
*/


}


