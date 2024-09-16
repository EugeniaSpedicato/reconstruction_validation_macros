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

void single_pre_notReco(){


TChain * cbmsim = new TChain("cbmsim");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_commit_b892679b_0_infmrad_2hitFirstModules.root");

//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/single_sample_2hitFirstModules_NOoutchi2_MCcorrections_v2_0.root");
/*cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/single_sample_2hitFirstModules_NOoutchi2_MCcorrections_v2_1.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/single_sample_2hitFirstModules_NOoutchi2_MCcorrections_v2_2.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/single_sample_2hitFirstModules_NOoutchi2_MCcorrections_v2_3.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/single_sample_2hitFirstModules_NOoutchi2_MCcorrections_v2_4.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/single_sample_2hitFirstModules_NOoutchi2_MCcorrections_v2_5.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/single_sample_2hitFirstModules_NOoutchi2_MCcorrections_v2_6.root");
*/
/*cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/single_sample_2hitFirstModules_NOoutchi2_MCcorrections.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/single_sample_2hitFirstModules_NOoutchi2_MCcorrections1.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/single_sample_2hitFirstModules_NOoutchi2_MCcorrections2.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/single_sample_2hitFirstModules_NOoutchi2_MCcorrections3.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/single_sample_2hitFirstModules_NOoutchi2_MCcorrections4.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/single_sample_2hitFirstModules_NOoutchi2_MCcorrections5.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/single_sample_2hitFirstModules_NOoutchi2_MCcorrections6.root");
*/


        TClonesArray *MCTrack = 0;
        MuE::Event *MesmerEvent = 0;
        MUonERecoOutput *ReconstructionOutput = 0;

        cbmsim->SetBranchAddress("MCTrack", &MCTrack);
        cbmsim->SetBranchAddress("MesmerEvent", &MesmerEvent);
        cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);

                auto pos_on_track = [](double q, double m, double z){return (q + m*z);};

double signal=0.;

double reco=0.;


double error=0; double error1=0.;double error2=0.;double error3=0.;double error4=0.;double error5=0.;double error6=0.;

double reco_v=0.; double more_reco_v=0.; double reco0_v=0.;
double gen=0.;

int yes2=0; int yes_v=0;
int code_mu=-99; int code_e=-99; int code_mu_in=-99;
int TrackIdreco=-99;
double z_fix=912.7;

double sumW2=0.;
double sum_wgt=0.;

TH1::SetDefaultSumw2(kTRUE);
   const Int_t NBINS = 7;
   Double_t edges[NBINS + 1] = {0.0, 0.005, 0.010, 0.015, 0.020, 0.025, 0.032, 0.100};


   const Int_t NBINS2 = 12;
   Double_t edges2[NBINS2 + 1] = {0.,0.0002,0.0004,0.0006,0.0008,0.001,0.00125,0.0015,0.00175,0.002,0.003,0.004,0.005};


TH2D* h_2d=new TH2D("h2D","theta mu vs theta E with all cuts", 16, 0.,0.032,NBINS2,edges2);
TH2D* h_2d_pre=new TH2D("h2D_pre","theta mu vs theta E with all cuts", 16, 0.,0.032,NBINS2,edges2);

TH1::SetDefaultSumw2(kTRUE);

   TH1D* d_aco = new TH1D("d_aco_real", "Acoplanarity",600,-3.2,3.2);
   TH1D* theta_e = new TH1D("theta_e", "Electron scattering angles from MESMER",10,0.,0.035);
   TH1D* theta_mu = new TH1D("theta_mu", "Muon scattering angles from MESMER",20,0.,0.005);
   TH1D *h_opening=new TH1D("h_opening", "Opening angle reco events from MESMER",35,0.,0.035);


   TH1D* d_aco_pre = new TH1D("d_aco_real_pre", "Acoplanarity",600,-3.2,3.2);
   TH1D* theta_e_pre = new TH1D("theta_e_pre", "Electron scattering angles from MESMER",10,0.,0.035);
   TH1D* theta_mu_pre = new TH1D("theta_mu_pre", "Muon scattering angles from MESMER",20,0.,0.005);
   TH1D *h_opening_pre=new TH1D("h_opening_pre", "Opening angle reco events from MESMER",35,0.,0.035);



for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
		cbmsim->GetEntry(i);
		if(i%1000 == 0) cout<<"Entry "<<i<<endl;

        TVector3 p_muin_MC;
        TVector3 p_mu_MC;
        TVector3 p_e_MC;
	double the_gen, thmu_gen;
double emu=0.;



	for(int n = 0; n < MCTrack->GetEntries(); n++) {
	 const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
	 if(MCTr->interactionID()==0 and MCTr->pdgCode()==-13) {code_mu_in=n; p_muin_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); p_muin_MC.Unit();}
         if(MCTr->interactionID()==45 and MCTr->pdgCode()==11) {code_e=n; p_e_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); p_e_MC.Unit(); the_gen=p_muin_MC.Angle(p_e_MC);}
         if(MCTr->interactionID()==45 and MCTr->pdgCode()==-13) {code_mu=n; p_mu_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); p_mu_MC.Unit(); thmu_gen=p_muin_MC.Angle(p_mu_MC); emu=MCTr->energy();}
	}

 if(code_mu_in!=-99 and code_mu!=-99 and code_e!=-99){

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

//std::array<double,2> th_rec;

//std::array<TVector3,2> p;
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

if(tracks.at(j).processIDofLinkedTrack()==45 and tracks.at(j).sector()==1) TrackIdreco=tracks.at(j).linkedTrackID();
        std::vector<MUonERecoOutputHit> hits_=tracks.at(j).hits();


if(code_mu_in==tracks.at(j).linkedTrackID() and tracks.at(j).sector()==0){
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
if(tracks.at(j).processIDofLinkedTrack()==45 and tracks.at(j).sector()==1)
{
/*         TVector3 p1(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0); p1=p1.Unit();th_rec.at(yes2)=p1.Angle(p_muin); p.at(yes2)=p1;
         yes2++;
                 if(code_e==tracks.at(j).linkedTrackID()) {yes_e++;}
                 if(code_mu==tracks.at(j).linkedTrackID()) {yes_mu++;}*/
	std::vector<double> pos; pos.resize(6);
        for(int p=0;p<hits_.size();p++)pos.push_back(hits_.at(p).position());

         if(std::equal(pos.begin(),pos.end(),pos_mu.begin())){p_mu.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0); p_mu.Unit(); thmu_rec=p_mu.Angle(p_muin);}
        else if(std::equal(pos.begin(),pos.end(),pos_e.begin())){p_e.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0); p_e.Unit(); the_rec=p_e.Angle(p_muin);}

         yes2++;

	}
}


//estrapolo posizione negli ultimi due moduli della seconda stazione
/*double posxIN=pos_on_track(x0_in,th_inx,899.92180);
double posyIN=pos_on_track(y0_in,th_iny,903.76930);*/

// posizione locale negli ultimi due moduli della seconda stazione

double posxIN=99.;//pos_on_track(x0_in,th_inx,z_fix);
double posyIN=99.;//pos_on_track(y0_in,th_iny,z_fix);

//h_xy->Fill(posxIN,posyIN);
std::vector<MUonERecoOutputHit> stubs=ReconstructionOutput->reconstructedHits();

int stub0 = 0;
int stub1 = 0;

for(int s=0; s<stubs.size(); s++){
if(stubs.at(s).stationID()==0){stub0++; if(stubs.at(s).moduleID()==4){posxIN=stubs.at(s).position();} else if(stubs.at(s).moduleID()==5){posyIN=stubs.at(s).position();}   }
if(stubs.at(s).stationID()==1)stub1++;
}


if(stubs_muin==6 and abs(posxIN)<=1.5 and abs(posyIN)<=1.5 and chi2_muin<=2 and stub0==6 and th_muin<0.004){//and stub1<=15){


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
//first

sumW2+=(MesmerEvent->wgt_full*MesmerEvent->wgt_full);
sum_wgt+=MesmerEvent->wgt_full;

 reco+=MesmerEvent->wgt_full; error+=MesmerEvent->wgt_full*MesmerEvent->wgt_full;

 d_aco->Fill(acoplanarity_v,MesmerEvent->wgt_full);
 h_2d->Fill(the_rec,thmu_rec,MesmerEvent->wgt_full);

theta_mu->Fill(thmu_rec,MesmerEvent->wgt_full);
theta_e->Fill(the_rec,MesmerEvent->wgt_full);
h_opening->Fill(p_mu.Angle(p_e),MesmerEvent->wgt_full);



			}//aco e chi
		}//chi!=0
	}//chiusura mu_in
    }// chiusura if code_x!0-99
code_e=-99;code_mu=-99;code_mu_in=-99;
TrackIdreco=-99;
yes2=0;yes_v=0;
} //end of general for

cout << "gen VS reco " << gen << " - " << reco << endl;

cout << "ALL" << endl;
cout << "Su " << signal << " eventi di segnale, " << reco << " +- " << sqrt(error) << " sono ricostruiti, con un rapporto del " << reco/signal*100 << "%"<< endl;



Int_t nx = h_2d->GetNbinsX();
Int_t ny = h_2d->GetNbinsY();
for (Int_t i=1; i<nx+1; i++) {
for (Int_t j=1; j<ny+1; j++) {
if (h_2d->GetBinContent(i,j)<=20) h_2d->SetBinContent(i,j,0);}}
h_2d->SaveAs("comparison_RDMC/2D_MC_new_pre2_single_firstmod.root");
h_2d_pre->SaveAs("comparison_RDMC/preCuts_2D_MC_new_pre2_single_firstmod.root");

theta_e->SaveAs("comparison_RDMC/theta_e_MC_new_pre2_single_firstmod.root");
theta_mu->SaveAs("comparison_RDMC/theta_mu_MC_new_pre2_single_firstmod.root");
h_opening->SaveAs("comparison_RDMC/opening_MC_new_pre2_single_firstmod.root");
d_aco->SaveAs("comparison_RDMC/d_aco_MC_new_pre2_single_firstmod.root"); 


theta_e_pre->SaveAs("comparison_RDMC/preCuts_theta_e_MC_new_pre2_single_firstmod.root");
theta_mu_pre->SaveAs("comparison_RDMC/preCuts_theta_mu_MC_new_pre2_single_firstmod.root");
h_opening_pre->SaveAs("comparison_RDMC/preCuts_opening_MC_new_pre2_single_firstmod.root");
d_aco_pre->SaveAs("comparison_RDMC/preCuts_d_aco_MC_new_pre2_single_firstmod.root"); 


}

