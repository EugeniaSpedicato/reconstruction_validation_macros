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

void RealDataAnalyzer(){

// 	TFile *inputfile = new TFile("/mnt/raid10/DATA/espedica/fairmu/dataReconstruction_merged_3232_3233_1_4M_outlier.root");
//        TTree* cbmsim = (TTree*) inputfile->Get("cbmsim");

TChain * cbmsim = new TChain("cbmsim");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/dataReconstruction_merged_3238_3239_2hit.root");
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/dataReconstruction_merged_3238_3239_second_outlier_1.root");
//3232_3233_1_4M_outlier.root");
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/dataReconstruction_merged_3232_3233_small4M_outlier.root");

        MUonERecoOutput *ReconstructionOutput = 0;


   cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);

    TH2D *h_2d=new TH2D("h_2d","Electron VS muone angle" ,700,0.,0.07,90,0.,0.003);
    TH1D *h_Z=new TH1D("h_Z","Z of the vrtx found",800,800.,1000.);


// X or Y position on the track at a given Z
                auto pos_on_track = [](double q, double m, double z){return (q + m*z);};

for(Long64_t i = 0; i < 100000; i++) {//cbmsim->GetEntries(); i++) {
		cbmsim->GetEntry(i);
		if(i%1000 == 0) cout<<"Entry "<<i<<endl;

vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
vector<MUonERecoOutputVertex> vrtx = ReconstructionOutput->reconstructedVertices();
//vector<MUonERecoOutputAdaptiveFitterVertex> vrtx = ReconstructionOutput->adaptiveFitterVertices();

double chi=99.;
double Z_sig=0.;

  for(int j=0; j<vrtx.size();j++)
        {
	 if(j==0) {chi=vrtx.at(j).chi2perDegreeOfFreedom(); h_chi->Fill(chi);
        //if(vrtx.at(j).chi2()<50) h_Z->Fill(vrtx.at(j).z());
	 //if(vrtx.at(j).chi2perDegreeOfFreedom()<50){h_Z->Fill(vrtx.at(j).z());
	Z_sig=vrtx.at(j).z();
	 cout << "vertex chi " << chi  << " at Z_sig " << vrtx.at(j).z() << endl;
		}
        }
// }


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

std::array<double,2> chi_min;

    for(int j=0; j<tracks.size();j++)
    {
        if(tracks.at(j).sector()==0) sec0++;
        if(tracks.at(j).sector()==1) sec1++;
	}

int index=0;

    for(int j=0; j<tracks.size();j++)
    {
	if(tracks.at(j).sector()==0 and sec0==1){
	double th_inx=tracks.at(j).xSlope();
	double th_iny=tracks.at(j).ySlope();
	thin1.SetXYZ(th_inx,th_iny,1.0);
	thin1=thin1.Unit();
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
	index++;
		 }

	}


if(sec0==1 and sec1==2){

	if(theta.at(0)>theta.at(1)){electron=thin2_v.at(0); muon=thin2_v.at(1); th_el=theta.at(0); th_mu=theta.at(1);}
	else{electron=thin2_v.at(1); muon=thin2_v.at(0); th_el=theta.at(1); th_mu=theta.at(0);}

						double dotProduct = muon.Dot(electron);
                                                TVector3 crossProduct = muon.Cross(electron);
                                                double T = thin1.Dot(crossProduct);
                                                TVector3 im= thin1.Cross(muon);
                                                TVector3 ie= thin1.Cross(electron);
                                                T = T>0? 1:-1;
                                                acoplanarity= T*(TMath::Pi()- acos( ((im).Dot(ie))/(im.Mag()*ie.Mag()) ));

	double Z=zcrosslines(x0[0],y0[0],mx[0],my[0],x0[1],y0[1],mx[1],my[1]);
//tracks.at(1).x0(),tracks.at(1).y0(),tracks.at(1).xSlope(),tracks.at(1).ySlope(), tracks.at(2).x0(),tracks.at(2).y0(),tracks.at(2).xSlope(),tracks.at(2).ySlope());

  if(abs(Z-911.)<20. and chi_min.at(0)<50 and chi_min.at(1)<50) h_aco->Fill(acoplanarity);

//  if(abs(acoplanarity)<1) h_Z->Fill(Z);

//if(abs(Z-911.)<20. and
	if(abs(acoplanarity)<1 and chi<50){
                        h_2d->Fill(th_el,th_mu);
			//h_Z->Fill(Z);
		}

	}//if(sec0==1 and sec1==2)

}



    TCanvas n("n","n",700,700);
    n.Divide(1,3);
    n.cd(1);
    h_Z->Draw();
    n.cd(2);
    h_2d->Draw();
    n.cd(3);
    h_aco->Draw();
    n.cd(4);
    h_chi->Draw();
    n.SaveAs("info_2hit.pdf");

}


