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

        TFile *inputfile = new TFile("/mnt/raid10/DATA/espedica/fairmu/dataReconstruction_3234-3235_filtered.root");
//dataReconstruction_3234-3235_filtered.root");
//run3234_3235_1.root");
//dataReconstruction_3234-3235_new.root");
//dataReconstruction_try_sigma.root");
//TRMesmer_3cm.root");
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

string title;

TH1D* h_aco=new TH1D("aco","Acoplanarity of muone+electron from sig minbias",100,-3.14,3.14);

        std::vector<TH1D*> residual_mu_trk(4);
        residual_mu_trk.at(0)=new TH1D("residual_mu_trk0","Residual mu_out X track Vs mu_in X track",100,-0.05,0.05);
        residual_mu_trk.at(1)=new TH1D("residual_mu_trk1","Residual mu_out Y track Vs mu_in Y track",100,-0.05,0.05);
        residual_mu_trk.at(2)=new TH1D("residual_mu_trk2","Pull residual mu_out X track Vs mu_in X track",150,-0.3,0.3);
        residual_mu_trk.at(3)=new TH1D("residual_mu_trk3","Pull residual mu_out Y track Vs mu_in Y track",150,-0.3,0.3);


        std::vector<TH1D*> residual_e_trk(4);
        residual_e_trk.at(0)=new TH1D("residual_e_trk0","Residual e_out X track Vs mu_in X track",200,-0.8,0.8);
        residual_e_trk.at(1)=new TH1D("residual_e_trk1","Residual e_out Y track Vs mu_in Y track",200,-0.8,0.8);
        residual_e_trk.at(2)=new TH1D("residual_e_trk2","Pull residual e_out X track Vs mu_in X track",250,-1.,1.);
        residual_e_trk.at(3)=new TH1D("residual_e_trk3","Pull residual e_out Y track Vs mu_in Y track",250,-1.,1.);


        std::vector<TH1D*> residual_muin_vrtx(4);
for(int m=0; m<4; m++){
        string name="residual_muin_vrtx_"+to_string(m);
                if(m%2==0) title="Residual mu_in X track Vs best vrtx ";
                if(m%2!=0) title="Residual mu_in Y track Vs best vrtx ";
	residual_muin_vrtx.at(m)=new TH1D(name.c_str(),title.c_str(),500,-0.5,0.5);//500,-0.005,0.005
                        }

        std::vector<TH1D*> residual_mu_vrtx(4);
for(int m=0; m<4; m++){
        string name="residual_mu_vrtx_"+to_string(m);
                if(m%2==0) title="Residual mu_out X track Vs best vrtx ";
                if(m%2!=0) title="Residual mu_out Y track Vs best vrtx ";
	residual_mu_vrtx.at(m)=new TH1D(name.c_str(),title.c_str(),500,-0.5,0.5);
                        }


        std::vector<TH1D*> residual_e_vrtx(4);
for(int m=0; m<4; m++){
        string name="residual_e_vrtx_"+to_string(m);
                if(m%2==0) title="Residual e_out X track Vs best vrtx ";
                if(m%2!=0) title="Residual e_out Y track Vs best vrtx ";
	residual_e_vrtx.at(m)=new TH1D(name.c_str(),title.c_str(),500,-0.5,0.5);
                        }




        std::vector<TH1D*> residual_muin(6);
for(int m=0; m<6; m++){
        string name="residual_muin_"+to_string(m);
                string title="Residual mu_in of module "+to_string(m);
if(m==2 or m==3) residual_muin.at(m)=new TH1D(name.c_str(),title.c_str(),600,-0.06,0.06);//500,-0.5,0.5)
else residual_muin.at(m)=new TH1D(name.c_str(),title.c_str(),600,-0.06,0.06);//500,-0.005,0.005
                        }

        std::vector<TH1D*> residual_mu(6);
for(int m=0; m<6; m++){
        string name="residual_mu_"+to_string(m);
                string title="Residual mu_out of module "+to_string(m);
if(m==2 or m==3) residual_mu.at(m)=new TH1D(name.c_str(),title.c_str(),600,-0.06,0.06);
else residual_mu.at(m)=new TH1D(name.c_str(),title.c_str(),600,-0.06,0.06);
                        }


        std::vector<TH1D*> residual_e(6);
for(int m=0; m<6; m++){
        string name="residual_e_"+to_string(m);
                string title="Residual e_out of module "+to_string(m);
if(m==2 or m==3) residual_e.at(m)=new TH1D(name.c_str(),title.c_str(),600,-0.06,0.06);
else residual_e.at(m)=new TH1D(name.c_str(),title.c_str(),600,-0.06,0.06);
                        }

TH2D *h_2d_mu=new TH2D("h_2dmu","Electron VS muone angle with aco, chi muon coordinate cut" ,700,0.,0.07,90,0.,0.003);
TH2D *h_2d_all=new TH2D("h_2dall","Electron VS muone angle with all events" ,320,0.,0.032,90,0.,0.003);
TH2D *h_2d_chi=new TH2D("h_2dchi","Electron VS muone angle with chi2 cut" ,320,0.,0.032,90,0.,0.003);
TH2D *h_2d_aco=new TH2D("h_2daco","Electron VS muone angle with  chi2 and aco cut" ,320,0.,0.032,90,0.,0.003);
TH2D *h_2d=new TH2D("h_2d","Electron VS muone angle (1 reco in- 2 reco out)" ,320,0.,0.032,90,0.,0.003);

TH1D *vrtx_x=new TH1D("vrtx_x","Kinemtic vertex x position",250,-5.,5.);
TH1D *vrtx_y=new TH1D("vrtx_y","Kinemtic vertex y position",250,-5.,5.);
TH1D *mu_in_x=new TH1D("mu_in_x","Muon_in x position at Z_t when interacts",250,-5.,5.);
TH1D *mu_in_y=new TH1D("mu_in_y","Muon_in y position at Z_t when interacts",250,-5.,5.);


// X or Y position on the track at a given Z
                auto pos_on_track = [](double q, double m, double z){return (q + m*z);};

for(Long64_t i = 0; i <cbmsim->GetEntries(); i++) {//5000000; i++) {//cbmsim->GetEntries(); i++) {
		cbmsim->GetEntry(i);
		if(i%1000 == 0) cout<<"Entry "<<i<<endl;

vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
//vector<MUonERecoOutputVertex> vrtx = ReconstructionOutput->reconstructedVertices();
//vector<MUonERecoOutputAdaptiveFitterVertex> vrtx = ReconstructionOutput->adaptiveFitterVertices();

double chi=99.;
double Z_sig=0.;

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
double x0_err[2];
double y0_err[2];

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
double x0_err_muin, y0_err_muin;

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
        x0_err_muin=tracks.at(j).x0Error();
        y0_err_muin=tracks.at(j).y0Error();
/*if(chi!=0 and chi!=99){
	cout << "------ event " << i << " -------" << endl;

        cout << "post vrtx: trackY in Z " << pos_on_track(muin.y0(),muin.ySlope(),0.) << " VS " << y_v << endl;
        cout << "pre vrtx: trackY in Z " << pos_on_track(tracks.at(j).y0(),tracks.at(j).ySlope(),z_v) << " VS " << y_v << endl;

        cout << "post vrtx: trackX in Z " << pos_on_track(muin.x0(),muin.xSlope(),0.) << " VS " << x_v << endl;
        cout << "pre vrtx: trackX in Z " << pos_on_track(tracks.at(j).x0(),tracks.at(j).xSlope(),z_v) << " VS " << x_v << endl;
}*/

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
		 }

	}

double chi2_e,chi2_mu;
int id_e,id_mu;
double x0_err_mu, x0_err_e, y0_err_mu, y0_err_e;

if(sec0==1 and sec1==2 and stubs[0]==6 and stubs[1]==6 and stubs_muin==6){


	if(theta.at(0)>theta.at(1)){id_e=1; id_mu=2; electron=thin2_v.at(0); muon=thin2_v.at(1); stubs_e=stubs[0]; stubs_mu=stubs[1];
				    th_el=theta.at(0); th_mu=theta.at(1); chi2_e=chi_min.at(0); chi2_mu=chi_min.at(1);
					x0_err_mu=x0_err[1]; x0_err_e=x0_err[0] ; y0_err_mu=y0_err[1] ; y0_err_e=y0_err[0];}

	else{id_e=2; id_mu=1; electron=thin2_v.at(1); muon=thin2_v.at(0);  stubs_e=stubs[1]; stubs_mu=stubs[0];
	     th_el=theta.at(1); th_mu=theta.at(0); chi2_e=chi_min.at(1); chi2_mu=chi_min.at(0);
		 x0_err_mu=x0_err[0]; x0_err_e=x0_err[1] ; y0_err_mu=y0_err[0] ; y0_err_e=y0_err[1];}

						double dotProduct = muon.Dot(electron);
                                                TVector3 crossProduct = muon.Cross(electron);
                                                double T = thin1.Dot(crossProduct);
                                                TVector3 im= thin1.Cross(muon);
                                                TVector3 ie= thin1.Cross(electron);
                                                T = T>0? 1:-1;
                                                acoplanarity= T*(TMath::Pi()- acos( ((im).Dot(ie))/(im.Mag()*ie.Mag()) ));

if(th_el<0.032 and th_mu>0.0002){

        h_2d_all->Fill(th_el,th_mu);
       if(chi2_e<5. and chi2_e<5. and chi2_muin<5.) {h_2d_chi->Fill(th_el,th_mu);  h_aco->Fill(acoplanarity);}
       if(chi2_e<5. and chi2_e<5. and chi2_muin<5. and abs(acoplanarity)<1){//and th_mu>0.0002){//chi2_e<5 and chi2_e<5 and chi2_muin<5){.]// abs(acoplanarity)<1)

        cout << "------ event " << i << " -------" << endl;

	h_2d_aco->Fill(th_el,th_mu);

double x_t[3]={0.};
double y_t[3]={0.};

    for(int j=0; j<tracks.size();j++){
std::vector<MUonERecoOutputHit> hits_=tracks.at(j).hits();

	  if(j==0) {for(int h=0;h<hits_.size();h++){residual_muin.at(hits_.at(h).moduleID())->Fill(hits_.at(h).perpendicularResiduum());}
			x_t[0]=pos_on_track(tracks.at(j).x0(),tracks.at(j).xSlope(),z_v);
			y_t[0]=pos_on_track(tracks.at(j).y0(),tracks.at(j).ySlope(),z_v);
                        double x1=pos_on_track(muin.x0(),muin.xSlope(),0.);
                        double y1=pos_on_track(muin.y0(),muin.ySlope(),0.);

       //cout << "Difference muon_in y_t[0]-y_v " << y_t[0]-y_v << " and x0-x_v " << x0-x_v << endl;


			residual_muin_vrtx.at(0)->Fill(x_t[0]-x_v);
                        residual_muin_vrtx.at(1)->Fill(y_t[0]-y_v);
                        residual_muin_vrtx.at(2)->Fill(x1-x_v);
                        residual_muin_vrtx.at(3)->Fill(y1-y_v);
			}
          if(j==id_mu) {for(int h=0;h<hits_.size();h++){residual_mu.at(hits_.at(h).moduleID())->Fill(hits_.at(h).perpendicularResiduum());}
			x_t[1]=pos_on_track(tracks.at(j).x0(),tracks.at(j).xSlope(),z_v);
                        y_t[1]=pos_on_track(tracks.at(j).y0(),tracks.at(j).ySlope(),z_v);
                        double x1=pos_on_track(mu.x0(),mu.xSlope(),0.);
                        double y1=pos_on_track(mu.y0(),mu.ySlope(),0.);

       //cout << "Difference muon y_t[0]-y_v " << y_t[0]-y_v << " and x0-x_v " << x0-x_v << endl;

                        residual_mu_vrtx.at(0)->Fill(x_t[1]-x_v);
                        residual_mu_vrtx.at(1)->Fill(y_t[1]-y_v);
                        residual_mu_vrtx.at(2)->Fill(x1-x_v);
                        residual_mu_vrtx.at(3)->Fill(y1-y_v);
			}

          if(j==id_e) {for(int h=0;h<hits_.size();h++){residual_e.at(hits_.at(h).moduleID())->Fill(hits_.at(h).perpendicularResiduum());}
                        x_t[2]=pos_on_track(tracks.at(j).x0(),tracks.at(j).xSlope(),z_v);
                        y_t[2]=pos_on_track(tracks.at(j).y0(),tracks.at(j).ySlope(),z_v);
                        double x1=pos_on_track(e.x0(),e.xSlope(),0.);
                        double y1=pos_on_track(e.y0(),e.ySlope(),0.);

       //cout << "Difference electron y_t[0]-y_v " << y_t[0]-y_v << " and x0-x_v " << x0-x_v << endl;

                        residual_e_vrtx.at(0)->Fill(x_t[2]-x_v);
                        residual_e_vrtx.at(1)->Fill(y_t[2]-y_v);
                        residual_e_vrtx.at(2)->Fill(x1-x_v);
                        residual_e_vrtx.at(3)->Fill(y1-y_v);
                        }

		 }//for()tracks



	                residual_mu_trk.at(0)->Fill(x_t[1]-x_t[0]);
                        residual_mu_trk.at(1)->Fill(y_t[1]-y_t[0]);
                        residual_mu_trk.at(2)->Fill( (x_t[1]-x_t[0])/( sqrt(x0_err_mu*x0_err_mu + x0_err_muin*x0_err_muin) ) );
                        residual_mu_trk.at(3)->Fill( (y_t[1]-y_t[0])/( sqrt(y0_err_mu*y0_err_mu + y0_err_muin*y0_err_muin) ) );
                        residual_e_trk.at(0)->Fill(x_t[2]-x_t[0]);
                        residual_e_trk.at(1)->Fill(y_t[2]-y_t[0]);
                        residual_e_trk.at(2)->Fill( (x_t[2]-x_t[0])/( sqrt(x0_err_e*x0_err_e + x0_err_muin*x0_err_muin) ) );
                        residual_e_trk.at(3)->Fill( (y_t[2]-y_t[0])/( sqrt(y0_err_e*y0_err_e + y0_err_muin*y0_err_muin) ) );

			vrtx_x->Fill(x_v);
			mu_in_x->Fill(x_t[0]);
                        vrtx_y->Fill(y_v);
                        mu_in_y->Fill(y_t[0]);

		}//th_el<32
	}//acop
    }//if(sec0==1 and sec1==2)
  }
/*
TCanvas n("n","n",700,700);
n.Divide(2,2);
n.cd(1);
h_2d_all->Draw();
n.cd(2);
h_2d_chi->Draw();
n.cd(3);
h_2d_aco->Draw();
n.SaveAs("h_2d.pdf");
*/
/*TCanvas n0("n0","n0",1000,1000);
n0.Divide(2,3);
for(int m=0; m<6; m++){
n0.cd(m+1);
residual_muin[m]->Draw();
gPad->SetLogy();}
n0.SaveAs("res_muin.pdf");

TCanvas n1("n1","n1",1000,1000);
n1.Divide(2,3);
for(int m=0; m<6; m++){
n1.cd(m+1);
residual_mu[m]->Draw();
gPad->SetLogy();}
n1.SaveAs("res_mu.pdf");

TCanvas n2("n2","n2",1000,1000);
n2.Divide(2,3);
for(int m=0; m<6; m++){
n2.cd(m+1);
residual_e[m]->Draw();
gPad->SetLogy();}
n2.SaveAs("res_e.pdf");
*/
/*
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


TCanvas n4("n4","n4",1000,1000);
n4.Divide(2,3);
n4.cd(1);
residual_muin_vrtx[0]->Draw();
residual_muin_vrtx[2]->SetLineColor(kRed);
//residual_muin_vrtx[2]->Draw("same");
//gPad->SetLogy();
n4.cd(2);
residual_muin_vrtx[1]->Draw();
residual_muin_vrtx[3]->SetLineColor(kRed);
//residual_muin_vrtx[3]->Draw("same");
//gPad->SetLogy();

n4.cd(3);
residual_mu_vrtx[0]->Draw();
residual_mu_vrtx[2]->SetLineColor(kRed);
//residual_mu_vrtx[2]->Draw("same");
//gPad->SetLogy(); 
n4.cd(4);
residual_mu_vrtx[1]->Draw();
residual_mu_vrtx[3]->SetLineColor(kRed);
//residual_mu_vrtx[3]->Draw("same");
//gPad->SetLogy();

n4.cd(5);
residual_e_vrtx[0]->Draw();
residual_e_vrtx[2]->SetLineColor(kRed);
//residual_e_vrtx[2]->Draw("same");
//gPad->SetLogy();
n4.cd(6);
residual_e_vrtx[1]->Draw();
residual_e_vrtx[3]->SetLineColor(kRed);
//residual_e_vrtx[3]->Draw("same");
//gPad->SetLogy();
n4.SaveAs("res_vrtx.pdf");


TCanvas n5("n5","n5",700,700);
h_aco->Draw();
n5.SaveAs("aco.pdf");*/

TCanvas n6("n6","n6",700,700);
n6.Divide(1,2);
n6.cd(1);
vrtx_x->SetLineColor(kRed);
vrtx_x->Draw();
mu_in_x->Draw("same");
n6.cd(2);
vrtx_y->SetLineColor(kRed);
vrtx_y->Draw();
mu_in_y->Draw("same");
n6.SaveAs("vrtx.pdf");

}


