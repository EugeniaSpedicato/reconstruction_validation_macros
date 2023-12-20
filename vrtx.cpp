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

        TFile *inputfile = new TFile("/mnt/raid10/DATA/espedica/fairmu/dataReconstruction_3234-3235_new12_5M.root");
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

TH2D *h_2d_chi=new TH2D("h_2d_chi","Electron VS muone angle chi2_vrtx cut" ,700,0.,0.07,100,0.,0.005);
TH2D *h_2d_vrtx_chi=new TH2D("h_2d_vrtx_chi","Electron VS muone angle from best vrtx chi2_vrtx cut",700,0.,0.07,100,0.,0.005);

TH2D *h_2d_aco1=new TH2D("h_2d_aco1","Electron VS muone angle |aco|<1 cut" ,700,0.,0.07,100,0.,0.005);
TH2D *h_2d_vrtx_aco1=new TH2D("h_2d_vrtx_aco1","Electron VS muone angle from best vrtx |aco|<1 cut",700,0.,0.07,100,0.,0.005);

TH2D *h_2d_aco2=new TH2D("h_2d_aco2","Electron VS muone angle |aco|<1 cut and chi2<20" ,700,0.,0.07,100,0.,0.005);
TH2D *h_2d_vrtx_aco2=new TH2D("h_2d_vrtx_aco2","Electron VS muone angle from best vrtx |aco|<1 cut and chi2<20",700,0.,0.07,100,0.,0.005);

TH2D *h_2d_shift=new TH2D("h_2d_shift","Electron VS muone angle PRE vrtx all cuts when d_th>1.5 mrad",700,0.,0.07,100,0.,0.005);
TH2D *h_2d_shift_v=new TH2D("h_2d_shift_v","Electron VS muone angle from best vrtx all cuts when d_th>1.5 mrad",700,0.,0.07,100,0.,0.005);

TH1D* h_aco=new TH1D("aco","Acoplanarity of muone+electron from sig minbias pre vrtx",600,-3.14,3.14);
TH1D* h_aco_v=new TH1D("aco_v","Acoplanarity of muone+electron from sig minbias after vrtx",600,-3.14,3.14);

TH1D* h_aco_chi=new TH1D("acoc","Acoplanarity of muone+electron from sig minbias pre vrtx with chi2<20 cut",600,-3.14,3.14);
TH1D* h_aco_chi_v=new TH1D("aco_vc","Acoplanarity of muone+electron from sig minbias after vrtx with chi2<20 cut",600,-3.14,3.14);

TH1D* h_chi_v=new TH1D("chi","Chi2 vertex",1000,0.,1000.);
TH1D* h_chi_v_aco1=new TH1D("chi_a1","Chi2 vertex with |aco|<1 cut",1000,0.,1000.);
TH1D* h_chi_v_aco2=new TH1D("chi_a2","Chi2 vertex with |aco|<1 cut",1000,0.,1000.);

TH1D* th_muin_d=new TH1D("th_muin_d","Difference theta mu_in pre and post vrtx",200,-0.0002,0.0002);
TH1D* th_mu_d=new TH1D("th_mu_d","Difference theta mu_out pre and post vrtx",200,-0.0002,0.0002);
TH1D* th_e_d=new TH1D("th_e_d","Difference theta e_out pre and post vrtx",100,-0.01,0.01);

TH2D *h_diff_th_el=new TH2D("h_diff_th_el", "Difference theta mu_in pre and post vrtx VS theta_el_pre" ,400,0.,0.04,100,-0.01,0.01);

TH1D* h_stubs=new TH1D("h_stubs","Number of stubs in strange events, where 12 are used", 40,0,40);
TH1D* h_stubs0=new TH1D("h_stubs0","Number of stubs in STATION0 (blue) and STATION1 (red)", 50,0,50);
TH1D* h_stubs1=new TH1D("h_stubs1","Number of stubs in STATION0 (blue) and STATION1 (red)", 50,0,50);

TH1D* h_distX_p=new TH1D("h_distX_p","Impact Parameter sqrt(d_x^2+d_y^2) track 1 and track 2 at Z=target, theta<5mrad",4000,-10.,10.);
TH2D* h_distY_p=new TH2D("h_distY_p","2D plot distance X vs distance Y of track 1 and track 2 at Z=target, theta<5mrad",4000,-10.,10.,4000,-10.,10.);
TH1D* h_distX=new TH1D("h_distX","Impact Parameter sqrt(d_x^2+d_y^2) track 1 and track 2 at Z=target, theta>5mrad",4000,-10.,10.);
TH2D* h_distY=new TH2D("h_distY","2D plot distance X vs distance Y of track 1 and track 2 at Z=target, theta>5mrad",4000,-10.,10.,4000,-10.,10.);

// X or Y position on the track at a given Z
                auto pos_on_track = [](double q, double m, double z){return (q + m*z);};

double all=0.;
double e_c=0.;
double mu_c=0.;
double mu_in_c=0.;
double v=0.;
double mu_st0=0.;
double mu_st1=0.;
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

double posxIN=pos_on_track(x0_in,th_inx,810);
double posyIN=pos_on_track(y0_in,th_iny,810);

vector<MUonERecoOutputHit> stubs_=ReconstructionOutput->reconstructedHits(); 
double nstubs_0=0.;
double nstubs_1=0.;
for(int s=0; s<stubs_.size(); s++){
if(stubs_.at(s).stationID()==0)nstubs_0++;
if(stubs_.at(s).stationID()==1)nstubs_1++;
}
h_stubs->Fill(ReconstructionOutput->reconstructedHitsMultiplicity());
h_stubs0->Fill(nstubs_0);
h_stubs1->Fill(nstubs_1);

//if(sec1==1 and nstubs_1==6)mu_st1++;
//if(sec0==1 and sec1==0 ) {cout << "event " << i << " nstubs_0,nstubs_1 " << nstubs_0 << ", " << nstubs_1 <<endl;}

if(sec0==1 and stubs_muin==6){ // and abs(posxIN)<=1.5 and abs(posyIN)<=1.5){
mu_st0++;

	if(sec1==2){// and stubs[0]==6 and stubs[1]==6){

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

double deltaX[3];
double deltaY[3];
double x0_e,y0_e,x0_mu,y0_mu;
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
                        x0_mu=x0; y0_mu=y0;
                        }
          if(j==id_e) {
                       	double x0=pos_on_track(tracks.at(j).x0(),tracks.at(j).xSlope(),z_v);
                        double y0=pos_on_track(tracks.at(j).y0(),tracks.at(j).ySlope(),z_v);
                        double x1=pos_on_track(e.x0(),e.xSlope(),0.);
                        double y1=pos_on_track(e.y0(),e.ySlope(),0.);

                        deltaX[2]=x0-x_v;
                        deltaY[2]=y0-y_v;
                        x0_e=x0; y0_e=y0;
                        }
 }
double ip=sqrt( (x0_mu-x0_e)*(x0_mu-x0_e) + (y0_mu-y0_e)*(y0_mu-y0_e) );

                if(th_mu<0.005 and th_el<0.005 ){h_distX_p->Fill(ip);
                h_distY_p->Fill(x0_mu-x0_e,y0_mu-y0_e);}
                else{h_distX->Fill(ip);
                h_distY->Fill(x0_mu-x0_e,y0_mu-y0_e);}



if(th_el<0.032){//and ip<0.250){
 h_aco->Fill(acoplanarity);
 h_aco_v->Fill(acoplanarity_v);
 h_chi_v->Fill(chi);

                h_2d->Fill(th_el,th_mu);
                h_2d_vrtx->Fill(vrtx.electronTheta(),vrtx.muonTheta());

if(chi<20){
//elastic++;
 h_2d_chi->Fill(th_el,th_mu);
 h_2d_vrtx_chi->Fill(vrtx.electronTheta(),vrtx.muonTheta());
 h_aco_chi->Fill(acoplanarity);
 h_aco_chi_v->Fill(acoplanarity_v);
	}

/*cout << "--------New event--------" << endl;
cout << "Incoming Muon chi2 before vrtx and after vrtxing: " << endl;
cout << "	" << chi2_muin << "  ->   " << chi2_muin_vrtx << endl;
cout << "Outgoing Muon before vrtx and after vrtxing : "  << endl;
cout << "       " << chi2_mu <<  "  ->   " << chi2_mu_vrtx << endl;
cout << "Outgoing Electron before vrtx and after vrtxing : "  << endl;
cout << "       " << chi2_e <<  "  ->   " << chi2_e_vrtx << endl;
*/

if( abs(acoplanarity)<1) {h_chi_v_aco2->Fill(chi);

  if(chi<20){ //h_2d->Fill(th_el,th_mu);

//cout << "event " << i << endl;
		if(vrtx.electronTheta()-th_el<-0.0015){h_2d_shift_v->Fill(vrtx.electronTheta(),vrtx.muonTheta());h_2d_shift->Fill(th_el,th_mu);}
elastic++;
                                 h_2d_aco2->Fill(th_el,th_mu);
                                 h_2d_vrtx_aco2->Fill(vrtx.electronTheta(),vrtx.muonTheta());
                th_muin_d->Fill(th_muin_v.Theta()-th_muin);
		th_mu_d->Fill(vrtx.muonTheta()-th_mu);
                th_e_d->Fill(vrtx.electronTheta()-th_el);
		h_diff_th_el->Fill(th_el,vrtx.electronTheta()-th_el);}

				}
	if( abs(acoplanarity)<1){/* and deltaX[0]<0.129 and deltaX[0]>-0.071
						and deltaX[1]<0.129 and deltaX[1]>-0.071
						and deltaX[2]<0.2 and deltaX[2]>-0.17
						and deltaY[0]<0.3 and deltaY[0]>0.
						and deltaY[1]<0.3 and deltaY[1]>0.
						and deltaY[2]<0.4 and deltaY[2]>-0.1){*/
//cout << "------ event " << i << "--------" << endl;

h_chi_v_aco1->Fill(chi);
all++;
if(chi>20)v++;
if(chi2_e_vrtx>20)e_c++;
if(chi2_muin_vrtx>20)mu_in_c++;
if(chi2_mu_vrtx>20)mu_c++;


	//if(th_el>0.01 and th_mu>0.0004) cout << "event " << i << endl;
                        h_2d_aco1->Fill(th_el,th_mu);
                        h_2d_vrtx_aco1->Fill(vrtx.electronTheta(),vrtx.muonTheta());

			}//aco
		}//th_mu<32mrad
	}//if(sec0==1 and sec1==2)
    }//mu_in conditions
//else{cout <<"event " << i << endl;}

}


cout << "Fraction of KF vrtx wth chi2>20 " << v/all << endl;
cout << "Fraction of muon_in tracks post-vrtx wth chi2>20 " << mu_in_c/all << endl;
cout << "Fraction of muon_out tracks post-vrtx wth chi2>20 " << mu_c/all << endl;
cout << "Fraction of electron_out tracks post-vrtx wth chi2>20 " << e_c/all << endl;
cout << mu_st0 << " muons reconstructed in station 0 over " << cbmsim->GetEntries() << ": " << 100*mu_st0/cbmsim->GetEntries() << "%" << endl;
cout << mu_st1 << " muons reconstructed in station 1 over " << cbmsim->GetEntries() << ": " << 100*mu_st1/cbmsim->GetEntries() << "%" << endl;
cout << elastic << " elastic events over " << mu_st0 << ": " << 100*elastic/mu_st0 << "%" << endl;

/*TCanvas n("n","n",700,700);
    n.Divide(1,2);
    n.cd(1);
    h_2d_aco2->Draw("COLZ");
    n.cd(2);
    h_2d_vrtx_aco2->Draw("COLZ");
    n.SaveAs("2d_vtx.pdf");*/
    TCanvas n("n","n",700,700);
    n.Divide(2,4);
    n.cd(1);
    h_2d->Draw();
    n.cd(2);
    h_2d_vrtx->Draw();
    n.cd(3);
    h_2d_chi->Draw();
    n.cd(4);
    h_2d_vrtx_chi->Draw();
    n.cd(5);
    h_2d_aco1->Draw();
    n.cd(6);
    h_2d_vrtx_aco1->Draw();
    n.cd(7);
    h_2d_aco2->Draw();
    n.cd(8);
    h_2d_vrtx_aco2->Draw();
    n.SaveAs("2d_vtx.pdf");
/*
    TCanvas n1("n1","n1",700,700);
    n1.Divide(1,5);
    n1.cd(1);
    h_aco->Draw();
    h_aco_v->SetLineColor(kRed);
    h_aco_v->Draw();
    h_aco->Draw("same");
    n1.cd(2);
h_aco_chi->Draw();
h_aco_chi_v->SetLineColor(kRed);
h_aco_chi_v->Draw("same");
    n1.cd(3);
h_chi_v->SetLineColor(kRed);
h_chi_v->Draw();
gPad->SetLogy();
    n1.cd(4);
h_chi_v_aco1->SetLineColor(kRed);
h_chi_v_aco1->Draw();
gPad->SetLogy();
    n1.cd(5);
h_chi_v_aco2->SetLineColor(kRed);
h_chi_v_aco2->Draw();
gPad->SetLogy();
    n1.SaveAs("aco_chi.pdf");

TCanvas n2("n2","n2",700,700);
n2.Divide(2,3);
n2.cd(1);
th_muin_d->Draw();
n2.cd(2);
th_mu_d->Draw();
n2.cd(3);
th_e_d->Draw();
n2.cd(4);
h_diff_th_el->Draw("COLZ");
n2.cd(5);
h_2d_shift->Draw();
h_2d_shift_v->SetMarkerColor(kRed);
h_2d_shift_v->Draw("same");
n2.SaveAs("diffth.pdf");
*/
    TCanvas n4("n4","n4",700,700);
    n4.Divide(2,2);
    n4.cd(1);
h_stubs1->SetLineColor(kRed);
h_stubs1->Draw();
h_stubs0->Draw("same");
gPad->SetLogy();
    n4.cd(2);
h_stubs->Draw();
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



}


