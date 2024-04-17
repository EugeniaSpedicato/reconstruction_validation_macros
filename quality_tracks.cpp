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

void quality_tracks(){

TChain * cbmsim = new TChain("cbmsim");

//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/theta_0-32mrad_1M_1hit_NOoutchi2_1M.root");
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/theta_0-32mrad_1M_1hit_NOoutchi2.root");

/*
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_0-5mrad_1M_1hit_NOoutchi2.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_5-10mrad_1M_1hit_NOoutchi2.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_10-15mrad_1M_1hit_NOoutchi2.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_15-20mrad_1M_1hit_NOoutchi2.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_20-25mrad_1M_1hit_NOoutchi2.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_25-32mrad_1M_1hit_NOoutchi2.root");
*/

/*
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_0-5mrad_1M_2hit_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_5-10mrad_1M_2hit_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_10-15mrad_1M_2hit_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_15-20mrad_1M_2hit_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_20-25mrad_1M_2hit_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_25-32mrad_1M_2hit_NOoutchi2_1M.root");
*/


cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_0-5mrad_0hit_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_5-10mrad_0hit_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_10-15mrad_0hit_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_15-20mrad_0hit_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_20-25mrad_0hit_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_25-32mrad_0hit_NOoutchi2_1M.root");

/*
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_0-5mrad_1M_1hitFirstModules_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_5-10mrad_1M_1hitFirstModules_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_10-15mrad_1M_1hitFirstModules_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_15-20mrad_1M_1hitFirstModules_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_20-25mrad_1M_1hitFirstModules_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_25-32mrad_1M_1hitFirstModules_NOoutchi2_1M.root");
*/

/*cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_0-5mrad_1M_2hitFirstModules_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_5-10mrad_1M_2hitFirstModules_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_10-15mrad_1M_2hitFirstModules_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_15-20mrad_1M_2hitFirstModules_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_20-25mrad_1M_2hitFirstModules_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_25-32mrad_1M_2hitFirstModules_NOoutchi2_1M.root");
*/

        TClonesArray *MCTrack = 0;
        TClonesArray *TrackerStripDigis = 0;
        TClonesArray *TrackerPoints = 0;
        TClonesArray *TrackerStubs = 0;
        MuE::Event *MesmerEvent = 0;
        MUonERecoOutput *ReconstructionOutput = 0;

        cbmsim->SetBranchAddress("MCTrack", &MCTrack);
        cbmsim->SetBranchAddress("TrackerPoints", &TrackerPoints);
        cbmsim->SetBranchAddress("TrackerStripDigis", &TrackerStripDigis);
        cbmsim->SetBranchAddress("TrackerStubs", &TrackerStubs);
        cbmsim->SetBranchAddress("MesmerEvent", &MesmerEvent);
        cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);


                auto pos_on_track = [](double q, double m, double z){return (q + m*z);};
   const Int_t NBINS = 14;
   const Int_t NBINS_mu = 14;
   Double_t edges_el[NBINS + 1] = {0.,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.010, 0.015, 0.020, 0.025, 0.032};
   Double_t edges_mu[NBINS_mu + 1] = {0.,0.0001,0.0002,0.0003,0.0004,0.0005,0.0006,0.0007,0.0008,0.0009,0.001,0.002,0.003,0.004,0.005};

TH1D::SetDefaultSumw2(kTRUE);

TH1D *theta_e_clone=new TH1D("theta_e_clone", "Reco elastic event: Electron scattering reco angles from MESMER with clones",NBINS,edges_el);
TH1D *theta_mu_clone=new TH1D("theta_mu_clone", "Reco elastic event: Muon scattering reco angles from MESMER with clones",NBINS_mu,edges_mu);

TH1D *theta_e_single_clone=new TH1D("theta_single_e_clone", "Reco single electron: Electron scattering reco angles from MESMER with clones",NBINS,edges_el);
TH1D *theta_mu_single_clone=new TH1D("theta_single_mu_clone", "Reco single muon: Muon scattering reco angles from MESMER with clones",NBINS_mu,edges_mu);

TH1D *theta_e_clone_quality=new TH1D("theta_e_clone_quality", "Reco elastic event: Electron scattering reco angles from MESMER with clones quality check",NBINS,edges_el);
TH1D *theta_mu_clone_quality=new TH1D("theta_mu_clone_quality", "Reco elastic event: Muon scattering reco angles from MESMER with clones quality check",NBINS_mu,edges_mu);

TH1D *theta_e_single_clone_quality=new TH1D("theta_single_e_clone_quality", "Reco single electron: Electron scattering reco angles from MESMER with clones quality check",NBINS,edges_el);
TH1D *theta_mu_single_clone_quality=new TH1D("theta_single_mu_clone_quality", "Reco single muon: Muon scattering reco angles from MESMER with clones quality check",NBINS_mu,edges_mu);

TH1D *theta_e_gen=new TH1D("theta_e_gen", "Electron scattering generated angles from MESMER",NBINS,edges_el);
TH1D *theta_mu_gen=new TH1D("theta_mu_gen", "Muon scattering generated angles from MESMER",NBINS_mu,edges_mu);

TH1D *theta_fake=new TH1D("theta_fake", "Fakes tracks vs reconstructeds cattering angles",NBINS,edges_el);

TH1D *theta_reco=new TH1D("theta_reco", "Fakes tracks vs reconstructed scattering angles",NBINS,edges_el);
TH1D *theta_e_reco=new TH1D("theta_e_reco", "Reco tracks vs electron scattering reco angles from MESMER",NBINS,edges_el);
TH1D *theta_mu_reco=new TH1D("theta_mu_reco", "Reco tracks: Muon scattering reco angles from MESMER",NBINS_mu,edges_mu);

TH1D *theta_e_good=new TH1D("theta_e_good", "Good tracks vs electron scattering reco angles from MESMER",NBINS,edges_el);
TH1D *theta_mu_good=new TH1D("theta_mu_good", "Good tracks vs euon scattering reco angles from MESMER",NBINS_mu,edges_mu);



//wnorm 1M/1M NLO 300outlier/nooutlier NLO 2hitshared
	double r_wnorm[6]={20.786765103274643,33.313091221576336,42.396733790329876,50.584815206143652,61.828110400735824,106.88513370134392};

//1000chi2 1M single range
//      double r_wnorm[6]={169.10075957722117,169.10075957722117,169.10075957722117,169.10075957722117,169.10075957722117,169.10075957722117};


double sum_wgt=0.;
double sumW2[6]={0.};

double no_clones=0.;
double clones=0.;
double no_clones_q=0.;
double clones_q=0.;


double n_correct=0.;
double n_swapped=0.;
double n_wrong=0.;
double tot=0.;
double all=0.;
double n_novrtx=0.;
double n_yesvrtx=0.;
double n_yes_bestvrtx=0.;

for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
                cbmsim->GetEntry(i);
                if(i%100 == 0) cout<<"Entry "<<i<<endl;

int index=7;

Double_t code_mu_in=-99;
Double_t code_e=-99;
Double_t code_mu=-99;
        TVector3 p_muin_MC;
        TVector3 p_mu_MC;
        TVector3 p_e_MC;
        Double_t the_gen, thmu_gen,theX_gen,theY_gen,thmuX_gen,thmuY_gen;

           int hit_modXmu=0; int hit_modXe=0;
           int hit_modYmu=0; int hit_modYe=0;
           int stereo_mu=0; int stereo_e=0;
     int hit_modXmuin=0;
     int hit_modYmuin=0;
     int stereo_muin=0;

        for(int n = 0; n < MCTrack->GetEntries(); n++) {
         const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
         if(MCTr->interactionID()==0 and MCTr->pdgCode()==-13) {code_mu_in=n; p_muin_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); p_muin_MC.Unit();
        for(int s=0; s<TrackerPoints->GetEntries(); s++)
                         {const MUonETrackerPoint *TrackerPt = static_cast<const MUonETrackerPoint*>(TrackerPoints->At(s));
         if(TrackerPt->trackPDGCode()==-13 and TrackerPt->trackID()==n and TrackerPt->stationID()==0){
                                                                                                 if(TrackerPt->moduleID()==0 or TrackerPt->moduleID()==4) hit_modXmuin++;
                                                                                                 else if(TrackerPt->moduleID()==1 or TrackerPt->moduleID()==5) hit_modYmuin++;
                                                                                                 else if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_muin++;}
                }
	}
         if(MCTr->interactionID()==45 and MCTr->pdgCode()==11) {code_e=n; p_e_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); p_e_MC.Unit(); the_gen=p_muin_MC.Angle(p_e_MC);
                                                                                 theX_gen=MCTr->ax();
                                                                                 theY_gen=MCTr->ay();

        for(int s=0; s<TrackerPoints->GetEntries(); s++)
                         {const MUonETrackerPoint *TrackerPt = static_cast<const MUonETrackerPoint*>(TrackerPoints->At(s));
         if(TrackerPt->trackPDGCode()==11 and TrackerPt->trackID()==code_e and TrackerPt->stationID()==1){
                                                                                                 if(TrackerPt->moduleID()==0 or TrackerPt->moduleID()==4) hit_modXe++;
                                                                                                 else if(TrackerPt->moduleID()==1 or TrackerPt->moduleID()==5) hit_modYe++;
                                                                                                 else if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_e++;}
        	}
	 }
         if(MCTr->interactionID()==45 and MCTr->pdgCode()==-13) {code_mu=n; p_mu_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); p_mu_MC.Unit(); thmu_gen=p_muin_MC.Angle(p_mu_MC);
                                                                                 thmuX_gen=MCTr->ax();
                                                                                 thmuY_gen=MCTr->ay();
        for(int s=0; s<TrackerPoints->GetEntries(); s++)
                         {const MUonETrackerPoint *TrackerPt = static_cast<const MUonETrackerPoint*>(TrackerPoints->At(s));
	 if(TrackerPt->trackPDGCode()==-13 and TrackerPt->trackID()==code_mu and TrackerPt->stationID()==1){
                                                                                                 if(TrackerPt->moduleID()==0 or TrackerPt->moduleID()==4) hit_modXmu++;
                                                                                                 else if(TrackerPt->moduleID()==1 or TrackerPt->moduleID()==5) hit_modYmu++;
                                                                                                 else if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_mu++;}
		}
	 }
         i++;
	}

	if(code_mu_in!=-99 and code_e!=-99 and code_mu!=-99 and hit_modXmuin==4 and hit_modYmuin==4 and stereo_muin>1 and hit_modXe==4 and hit_modYe==4 and stereo_e>1 and hit_modXmu==4 and hit_modYmu==4 and stereo_mu>1){


double wnorm=0.;

if(the_gen>=0 and the_gen<0.005){index=0; wnorm=r_wnorm[0];}
if(the_gen>=0.005 and the_gen<0.010){index=1; wnorm=r_wnorm[1];}
if(the_gen>=0.010 and the_gen<0.015){index=2; wnorm=r_wnorm[2];}
if(the_gen>=0.015 and the_gen<0.020){index=3; wnorm=r_wnorm[3];}
if(the_gen>=0.020 and the_gen<0.025){index=4; wnorm=r_wnorm[4];}
if(the_gen>=0.025 and the_gen<0.032){index=5; wnorm=r_wnorm[5];}

all+=MesmerEvent->wgt_LO*wnorm;


double opening_angle=p_mu_MC.Angle(p_e_MC);


Int_t yes_mu=0;
Int_t yes_e=0;
Int_t yes2=0;
Double_t th_inx,th_iny,x0_in,y0_in;
Double_t chi2_muin=999.;
Int_t stubs_muin=0;
Int_t e=0;
Int_t mu=0;
TVector3 p_muin,p_e,p_mu;
Double_t the_rec,theX_rec,theY_rec,thmu_rec,thmuX_rec,thmuY_rec,the_rec_vrtx,thmu_rec_vrtx;
MUonERecoOutputVertex vrtx = ReconstructionOutput->bestVertex();

vector<MUonERecoOutputVertex> vrtx_all = ReconstructionOutput->reconstructedVertices();
Double_t chi=vrtx.chi2perDegreeOfFreedom();

if(vrtx_all.size()==0){n_novrtx+=MesmerEvent->wgt_LO*wnorm;}
else{n_yesvrtx+=MesmerEvent->wgt_LO*wnorm; if(chi!=0)n_yes_bestvrtx+=MesmerEvent->wgt_LO*wnorm;}

vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
int other=0;
int mu_in=0;

vector<double> quality_e; quality_e.reserve(5);
vector<double> quality_mu; quality_mu.reserve(5);

         for (auto&& track : tracks) {

        if(code_mu_in==track.linkedTrackID() and track.sector()==0){
	mu_in++;
         th_inx=track.xSlope();
         th_iny=track.ySlope();
         x0_in=track.x0();
         y0_in=track.y0();
         chi2_muin=track.chi2perDegreeOfFreedom();
         stubs_muin=track.hits().size();
        p_muin.SetXYZ(th_inx,th_iny,1.0);
        p_muin=p_muin.Unit();
                        }
      if(track.sector()==1){

                        theta_mu_reco->Fill(thmu_gen,MesmerEvent->wgt_LO*wnorm);
                        theta_e_reco->Fill(the_gen,MesmerEvent->wgt_LO*wnorm);
                        TVector3 p1(track.xSlope(),track.ySlope(),1.0);p1.Unit();
                        theta_reco->Fill(p1.Angle(p_muin),MesmerEvent->wgt_LO*wnorm);


        if(track.processIDofLinkedTrack()==45)
                {

                 if(code_e==track.linkedTrackID()) {yes2++; e++;
						    theX_rec=track.xSlope(); theY_rec=track.ySlope();
						    p_e.SetXYZ(theX_rec,theY_rec,1.0); p_e=p_e.Unit();the_rec=p_e.Angle(p_muin);
						    quality_e.push_back(track.fractionOfHitsSharedWithLinkedTrack());}
                 if(code_mu==track.linkedTrackID()) {yes2++; mu++;
						    thmuX_rec=track.xSlope(); thmuY_rec=track.ySlope();
                                                    p_mu.SetXYZ(thmuX_rec,thmuY_rec,1.0); p_mu=p_mu.Unit();thmu_rec=p_mu.Angle(p_muin);
                                                    quality_mu.push_back(track.fractionOfHitsSharedWithLinkedTrack());}
			if(track.fractionOfHitsSharedWithLinkedTrack()>=0.65){
                        theta_mu_good->Fill(thmu_gen,MesmerEvent->wgt_LO*wnorm);
                        theta_e_good->Fill(the_gen,MesmerEvent->wgt_LO*wnorm);
			}else{
			TVector3 p(track.xSlope(),track.ySlope(),1.0);p.Unit();
                        theta_fake->Fill(p.Angle(p_muin),MesmerEvent->wgt_LO*wnorm);}
                 }//ID 45
             else{
                        TVector3 p(track.xSlope(),track.ySlope(),1.0);p.Unit();
                        theta_fake->Fill(p.Angle(p_muin),MesmerEvent->wgt_LO*wnorm);}
		}//sector1

         }//for


//if(stubs_muin==6 and abs(posxIN)<=1.5 and abs(posyIN)<=1.5 and chi2_muin<2 and thmu_gen>0.0002){
if(stubs_muin>=5 and chi2_muin<2 and thmu_gen>0.0002){

tot+=MesmerEvent->wgt_LO*wnorm;
theta_mu_gen->Fill(thmu_gen,MesmerEvent->wgt_LO*wnorm);
theta_e_gen->Fill(the_gen,MesmerEvent->wgt_LO*wnorm);


if(e>=1){theta_e_single_clone->Fill(the_gen,MesmerEvent->wgt_LO*wnorm);}
if(mu>=1){theta_mu_single_clone->Fill(thmu_gen,MesmerEvent->wgt_LO*wnorm);}

if(e==1 and mu==1 and chi!=0) no_clones+=MesmerEvent->wgt_LO*wnorm;

if(e>=1 and mu>=1 and chi!=0){
clones+=MesmerEvent->wgt_LO*wnorm;
                  	theta_mu_clone->Fill(thmu_gen,MesmerEvent->wgt_LO*wnorm);
                        theta_e_clone->Fill(the_gen,MesmerEvent->wgt_LO*wnorm);}

//quality checks

if(find_if(quality_e.begin(),quality_e.end(),[](double i){return i>=0.65;})!=end(quality_e) and find_if(quality_mu.begin(),quality_mu.end(),[](double i){return i>=0.65;})!=end(quality_mu)){

if(e==1 and mu==1 and chi!=0){no_clones_q+=MesmerEvent->wgt_LO*wnorm;}
if(e>=1 and mu>=1 and chi!=0){
clones_q+=MesmerEvent->wgt_LO*wnorm;
			theta_mu_clone_quality->Fill(thmu_gen,MesmerEvent->wgt_LO*wnorm);
			theta_e_clone_quality->Fill(the_gen,MesmerEvent->wgt_LO*wnorm);}

			}//if quality
                }//if mu_in
        }//if generated

}//for

cout << "fractions elastic events without quality without clones " << no_clones/tot*100 << "%" << endl;
cout << "fractions elastic events without quality with clones " << clones/tot*100 << "%" << endl;

cout << "fractions elastic events with quality without clones " << no_clones_q/tot*100 << "%" << endl;
cout << "fractions elastic events with quality with clones " << clones_q/tot*100 << "%" << endl;


TH1D * h0_original = (TH1D*) theta_e_clone->Clone();
TH1D * h0gen = (TH1D*) theta_e_gen->Clone();
h0_original->Divide(h0_original,h0gen,1,1,"B");
TH1D * h1_original = (TH1D*) theta_mu_clone->Clone();
TH1D * h1gen = (TH1D*) theta_mu_gen->Clone();
h1_original->Divide(h1_original,h1gen,1,1,"B");

TH1D * h0 = (TH1D*) theta_e_clone_quality->Clone();
h0->Divide(h0,h0gen,1,1,"B");
TH1D * h1 = (TH1D*) theta_mu_clone_quality->Clone();
h1->Divide(h1,h1gen,1,1,"B");

TH1D * h0_quality = (TH1D*) theta_e_clone_quality->Clone();
TH1D * h0_no = (TH1D*) theta_e_clone->Clone();
h0_quality->Divide(h0_quality,h0_no,1,1,"B");
TH1D * h1_quality = (TH1D*) theta_mu_clone_quality->Clone();
TH1D * h1_no = (TH1D*) theta_mu_clone->Clone();
h1_quality->Divide(h1_quality,h1_no,1,1,"B");


TCanvas n7("n7","n7",700,700);
n7.Divide(2,2);
n7.cd(1);
TGaxis::SetMaxDigits(3);
h0_original->SetTitle("Reco eff. elastic event as a funtion of el angle");
h0_original->SetLineColor(kOrange+2);
h0_original->SetMinimum(0.6);
h0_original->Draw("E");
h0->Draw("E same");
h0->SaveAs("proposal/quality_prevertex_el_gen_0hit_LO.root");
gStyle->SetOptStat(0);
n7.cd(2);
TGaxis::SetMaxDigits(3);
h1_original->SetTitle("Reco eff. elastic event as a funtion of mu angle");
h1_original->SetLineColor(kOrange+2);
h1_original->SetMinimum(0.6);
h1_original->Draw("E");
h1->Draw("E same");
h1->SaveAs("proposal/quality_prevertex_mu_gen_0hit_LO.root");
gStyle->SetOptStat(0);
n7.cd(3);
h0_quality->SetTitle("Ratio between sample with quality cut/without quality track VS el angle");
TGaxis::SetMaxDigits(3);
h0_quality->SetMinimum(0.6);
h0_quality->Draw("E");
h0_quality->SaveAs("proposal/quality_prevertex_el_clone_0hit_LO.root");
gStyle->SetOptStat(0);
n7.cd(4);
h1_quality->SetTitle("Ratio between sample with quality cut/without quality track VS el angle");
TGaxis::SetMaxDigits(3);
h1_quality->SetMinimum(0.6);
h1_quality->Draw("E");
h1_quality->SaveAs("proposal/quality_prevertex_mu_clone_0hit_LO.root");
gStyle->SetOptStat(0);
n7.SaveAs("prova_0hit_LO.pdf");


TH1D * h2 = (TH1D*) theta_e_reco->Clone();
TH1D * h3 = (TH1D*) theta_mu_reco->Clone();
TH1D * h4 = (TH1D*) theta_reco->Clone();

TH1D * h2_good = (TH1D*) theta_e_good->Clone();
h2_good->Divide(h2_good,h2,1,1,"B");

TH1D * h3_good = (TH1D*) theta_mu_good->Clone();
h3_good->Divide(h3_good,h3,1,1,"B");

TH1D * h3_fake = (TH1D*) theta_fake->Clone();
h3_fake->Divide(h3_fake,h4,1,1,"B");


TCanvas n8("n8","n8",2100,700);
n8.Divide(3,1);
TGaxis::SetMaxDigits(3);
n8.cd(1);
h2_good->Draw("E");
h2_good->SetMinimum(0.6);
h2_good->SaveAs("proposal/good_el_0hit_LO.root");
n8.cd(2);
h3_good->Draw("E");
h3_good->SetMinimum(0.6);
h3_good->SaveAs("proposal/good_mu_0hit_LO.root");
n8.cd(3);
h3_fake->Draw("E");
h3_fake->SaveAs("proposal/fake_0hit_LO.root");
n8.SaveAs("proposal/good_fake_0hit_LO.pdf");

}

