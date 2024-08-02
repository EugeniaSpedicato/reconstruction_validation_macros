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

void pre_notReco_separate(int ix){

/*TChain * cbmsim = new TChain("cbmsim");
if(ix==0)cbmsim->Add("/home/espedica/fair_install/instFairRoot/share/MUonE/macros/prova_0_5.root");
*/

/*
TChain * cbmsim = new TChain("cbmsim");
if(ix==0)cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_0-5mrad_1M_2hitFirstModules_NOoutchi2_reassign.root");
else if(ix==1)cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_5-10mrad_1M_2hitFirstModules_NOoutchi2_reassign.root");
else if(ix==2)cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_10-15mrad_1M_2hitFirstModules_NOoutchi2_reassign.root");
else if(ix==3)cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_15-20mrad_1M_2hitFirstModules_NOoutchi2_reassign.root");
else if(ix==4)cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_20-25mrad_1M_2hitFirstModules_NOoutchi2_reassign.root");
else if(ix==5)cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_25-32mrad_1M_2hitFirstModules_NOoutchi2_reassign.root");
else if(ix==6)cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/theta_32-inf_mrad_1M_2hitFirstModules_NOoutchi2_reassign.root");
else cout << "ERROR"<< endl;
*/

/*TChain * cbmsim = new TChain("cbmsim");
if(ix==0)cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_0-5mrad_1M_1hit_NOoutchi2.root");
else if(ix==1)cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_5-10mrad_1M_1hit_NOoutchi2.root");
else if(ix==2)cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_10-15mrad_1M_1hit_NOoutchi2.root");
else if(ix==3)cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_15-20mrad_1M_1hit_NOoutchi2.root");
else if(ix==4)cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_20-25mrad_1M_1hit_NOoutchi2.root");
else if(ix==5)cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_25-32mrad_1M_1hit_NOoutchi2.root");
else cout << "ERROR"<< endl;*/



TString gen_filename;

TChain * cbmsim = new TChain("cbmsim");
if(ix==0){cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_commit_b892679b_0_5mrad_2hitFirstModules.root");
	  gen_filename="/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_b892679b_0_5mrad_v2.root";}

else if(ix==1){cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_commit_b892679b_5_10mrad_2hitFirstModules.root");
          gen_filename="/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_b892679b_5_10mrad.root";}

else if(ix==2){cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_commit_b892679b_10_15mrad_2hitFirstModules.root");
          gen_filename="/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_b892679b_10_15mrad.root";}

else if(ix==3){cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_commit_b892679b_15_20mrad_2hitFirstModules.root");
          gen_filename="/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_b892679b_15_20mrad.root";}

else if(ix==4){cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_commit_b892679b_20_25mrad_2hitFirstModules.root");
          gen_filename="/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_b892679b_20_25mrad.root";}

else if(ix==5){cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_commit_b892679b_25_32mrad_2hitFirstModules.root");
          gen_filename="/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_b892679b_25_32mrad.root";}

else cout << "ERROR"<< endl;


   TFile *f2 = new TFile(gen_filename);//"/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_b892679b_25_32mrad.root");
   TTree *t2 = (TTree*)f2->Get("cbmsim");

cout << "entries cbmsim " << cbmsim->GetEntries() << endl;
cout << "entries t2 " << t2->GetEntries() << endl;


   cbmsim->AddFriend(t2);
   cout << t2->GetEntries() << endl;
   cout << cbmsim->GetEntries() << endl;
   cbmsim->Print();




int index=ix;


        TClonesArray *MCTrack = 0;
        MuE::Event *MesmerEvent = 0;
        MUonERecoOutput *ReconstructionOutput = 0;

        cbmsim->SetBranchAddress("MCTrack", &MCTrack);
        cbmsim->SetBranchAddress("MesmerEvent", &MesmerEvent);
        cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);


double signal=0.;
double reco=0.;
double bin[7]={0.};

int yes2=0; int yes_v=0;


TH1::SetDefaultSumw2(kTRUE);

   const Int_t NBINS2 = 12;
   Double_t edges2[NBINS2 + 1] = {0.,0.0002,0.0004,0.0006,0.0008,0.001,0.00125,0.0015,0.00175,0.002,0.003,0.004,0.005};


                TH2D* h_2d=new TH2D(Form("h_2d%d",static_cast<char>(index)),"h_2d", 16, 0.,0.032,NBINS2,edges2);
                TH2D* h_2d_pre=new TH2D(Form("h_2d_pre%d",static_cast<char>(index)),"h_2d_pre", 16, 0.,0.032,NBINS2,edges2);
                TH1D* theta_e=new TH1D(Form("theta_e%d",static_cast<char>(index)),"theta_e",10,0.,0.035);
                TH1D* theta_mu=new TH1D(Form("theta_mu%d",static_cast<char>(index)),"theta_mu",20,0.,0.005);
		TH1D* h_opening=new TH1D(Form("h_opening%d",static_cast<char>(index)),"h_opening",35,0.,0.035);
		TH1D* d_aco=new TH1D(Form("d_aco%d",static_cast<char>(index)),"d_aco",600,-3.2,3.2);
                TH1D* theta_e_pre=new TH1D(Form("theta_e_pre%d",static_cast<char>(index)),"theta_e_pre",10,0.,0.035);
                TH1D* theta_mu_pre=new TH1D(Form("theta_mu_pre%d",static_cast<char>(index)),"theta_mu_pre",20,0.,0.005);
                TH1D* h_opening_pre=new TH1D(Form("h_opening_pre%d",static_cast<char>(index)),"h_opening_pre",35,0.,0.035);
                TH1D* d_aco_pre=new TH1D(Form("d_aco_pre%d",static_cast<char>(index)),"d_aco_pre",600,-3.2,3.2);


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
double stubs_muin;

double th_muin=0.;

int sec0=0;
int sec1=0;
    for(int j=0; j<tracks.size();j++)
    {
        if(tracks.at(j).sector()==0) sec0++;
        if(tracks.at(j).sector()==1) sec1++;
        }

TVector3 p_e,p_mu;
TVector3 p_muin;

double the_rec=-99.;
double thmu_rec=-99.;
        std::vector<MUonERecoOutputHit> hits_mu=vrtx.outgoingMuon().hits();
        std::vector<double> pos_mu; pos_mu.resize(6);
        for(int p=0;p<hits_mu.size();p++)pos_mu.push_back(hits_mu.at(p).position());
        std::vector<MUonERecoOutputHit> hits_e=vrtx.outgoingElectron().hits();
        std::vector<double> pos_e; pos_e.resize(6);
        for(int p=0;p<hits_e.size();p++)pos_e.push_back(hits_e.at(p).position());

for(int j=0; j<tracks.size();j++)
{

        std::vector<MUonERecoOutputHit> hits_=tracks.at(j).hits();


if( sec0==1 and tracks.at(j).sector()==0){
        th_inx=tracks.at(j).xSlope();
        th_iny=tracks.at(j).ySlope();
        x0_in=tracks.at(j).x0();
        y0_in=tracks.at(j).y0();
        chi2_muin=tracks.at(j).chi2perDegreeOfFreedom();
        stubs_muin=hits_.size();
        p_muin.SetXYZ(th_inx,th_iny,1.0);
        p_muin=p_muin.Unit();
	th_muin=p_muin.Theta();
                        }
if(tracks.at(j).sector()==1)
{
	std::vector<double> pos; pos.resize(6);
        for(int p=0;p<hits_.size();p++)pos.push_back(hits_.at(p).position());

         if(std::equal(pos.begin(),pos.end(),pos_mu.begin())){p_mu.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0); p_mu.Unit(); thmu_rec=p_mu.Angle(p_muin);}
        else if(std::equal(pos.begin(),pos.end(),pos_e.begin())){p_e.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0); p_e.Unit(); the_rec=p_e.Angle(p_muin);}

yes2++;

	}
}


double posxIN=99.;//pos_on_track(x0_in,th_inx,z_fix);
double posyIN=99.;//pos_on_track(y0_in,th_iny,z_fix);

std::vector<MUonERecoOutputHit> stubs=ReconstructionOutput->reconstructedHits();
int stub0 = 0;
int stub1 = 0;

for(int s=0; s<stubs.size(); s++){
if(stubs.at(s).stationID()==0){stub0++; if(stubs.at(s).moduleID()==4){posxIN=stubs.at(s).position();} else if(stubs.at(s).moduleID()==5){posyIN=stubs.at(s).position();}   }
if(stubs.at(s).stationID()==1)stub1++;
}


if(stubs_muin==6 and abs(posxIN)<=1.5 and abs(posyIN)<=1.5 and chi2_muin<=2 and stub0==6 and th_muin<0.004){//and stub1<=15){


bin[index]+=MesmerEvent->wgt_full;
signal+=MesmerEvent->wgt_full;


if(chi!=0 and the_rec!=-99 and thmu_rec!=-99){


                                                double dotProduct_v = p_mu.Dot(p_e);
                                                TVector3 crossProduct_v = p_mu.Cross(p_e);
                                                double T_v = p_muin.Dot(crossProduct_v);
                                                TVector3 im_v= p_muin.Cross(p_mu);
                                                TVector3 ie_v= p_muin.Cross(p_e);
                                                T_v = T_v>0? 1:-1;
                                                double acoplanarity_v= T_v*(TMath::Pi()- acos( ((im_v).Dot(ie_v))/(im_v.Mag()*ie_v.Mag()) ));

 d_aco_pre->Fill(acoplanarity_v,MesmerEvent->wgt_full);
 h_2d_pre->Fill(the_rec,thmu_rec,MesmerEvent->wgt_full);
 theta_mu_pre->Fill(thmu_rec,MesmerEvent->wgt_full);
 theta_e_pre->Fill(the_rec,MesmerEvent->wgt_full);
 h_opening_pre->Fill(p_mu.Angle(p_e),MesmerEvent->wgt_full);


 if(abs(acoplanarity_v)<=1 and chi<20 and thmu_rec>0.0002 and stub1<=15){ // and vrtx.zKinematicFit()<915. and vrtx.zKinematicFit()>907.){//and the_rec<=0.032 and stub1<=15){// and vrtx.electronThea()>=0.0005 and the_rec<=0.02){//the_rec<=0.032){

 if(yes2>=2){
  reco+=MesmerEvent->wgt_full;
  d_aco->Fill(acoplanarity_v,MesmerEvent->wgt_full);
  h_2d->Fill(the_rec,thmu_rec,MesmerEvent->wgt_full);
  theta_mu->Fill(thmu_rec,MesmerEvent->wgt_full);
  theta_e->Fill(the_rec,MesmerEvent->wgt_full);
  h_opening->Fill(p_mu.Angle(p_e),MesmerEvent->wgt_full);
 }


			}//aco e chi
		}//chi!=0
	}//chiusura mu_in

yes2=0;yes_v=0;

} //end of general for


cout << "bin0 " << bin[0] << endl;
cout << "bin1 " << bin[1] << endl;
cout << "bin2 " << bin[2] << endl;
cout << "bin3 " << bin[3] << endl;
cout << "bin4 " << bin[4] << endl;
cout << "bin5 " << bin[5] << endl;
cout << "bin6 " << bin[6] << endl;

cout << "ALL" << endl;
cout << "Su " << signal << " eventi di segnale, " << reco << " sono ricostruiti, con un rapporto del " << reco/signal*100 << "%"<< endl;




Int_t nx = h_2d->GetNbinsX();
Int_t ny = h_2d->GetNbinsY();
for (Int_t i=1; i<nx+1; i++) {
for (Int_t j=1; j<ny+1; j++) {
if (h_2d->GetBinContent(i,j)<=20) h_2d->SetBinContent(i,j,0);}}
h_2d->SaveAs(Form("comparison_RDMC/2D_MC_%d_new_pre2_separate_b892679b.root",static_cast<char>(index)));
Int_t nx_pre = h_2d_pre->GetNbinsX();
Int_t ny_pre = h_2d_pre->GetNbinsY();
for (Int_t i=1; i<nx_pre+1; i++) {
for (Int_t j=1; j<ny_pre+1; j++) {
if (h_2d_pre->GetBinContent(i,j)<=20) h_2d_pre->SetBinContent(i,j,0);}}
h_2d_pre->SaveAs(Form("comparison_RDMC/preCuts_2D_MC_%d_new_pre2_separate_b892679b.root",static_cast<char>(index)));


theta_e->SaveAs(Form("comparison_RDMC/theta_e_MC_%d_new_pre2_separate_b892679b.root",static_cast<char>(index)));
theta_mu->SaveAs(Form("comparison_RDMC/theta_mu_MC_%d_new_pre2_separate_b892679b.root",static_cast<char>(index)));
h_opening->SaveAs(Form("comparison_RDMC/opening_MC_%d_new_pre2_separate_b892679b.root",static_cast<char>(index)));
d_aco->SaveAs(Form("comparison_RDMC/d_aco_MC_%d_new_pre2_separate_b892679b.root",static_cast<char>(index))); 


theta_e_pre->SaveAs(Form("comparison_RDMC/preCuts_theta_e_MC_%d_new_pre2_separate_b892679b.root",static_cast<char>(index)));
theta_mu_pre->SaveAs(Form("comparison_RDMC/preCuts_theta_mu_MC_%d_new_pre2_separate_b892679b.root",static_cast<char>(index)));
h_opening_pre->SaveAs(Form("comparison_RDMC/preCuts_opening_MC_%d_new_pre2_separate_b892679b.root",static_cast<char>(index)));
d_aco_pre->SaveAs(Form("comparison_RDMC/preCuts_d_aco_MC_%d_new_pre2_separate_b892679b.root",static_cast<char>(index))); 



}


