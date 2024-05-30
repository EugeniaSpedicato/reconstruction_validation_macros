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

void fake_rates_and_quality_prevrtx(int nhits){

TChain * cbmsim = new TChain("cbmsim");


if(nhits==0){
cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_0-5mrad_%dhit_NOoutchi2_1M.root",static_cast<char>(nhits)));
cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_5-10mrad_%dhit_NOoutchi2_1M.root",static_cast<char>(nhits)));
cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_10-15mrad_%dhit_NOoutchi2_1M.root",static_cast<char>(nhits)));
cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_15-20mrad_%dhit_NOoutchi2_1M.root",static_cast<char>(nhits)));
cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_20-25mrad_%dhit_NOoutchi2_1M.root",static_cast<char>(nhits)));
cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_25-32mrad_%dhit_NOoutchi2_1M.root",static_cast<char>(nhits)));
}
else{
cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_0-5mrad_1M_%dhitFirstModules_NOoutchi2_1M.root",static_cast<char>(nhits)));
cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_5-10mrad_1M_%dhitFirstModules_NOoutchi2_1M.root",static_cast<char>(nhits)));
cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_10-15mrad_1M_%dhitFirstModules_NOoutchi2_1M.root",static_cast<char>(nhits)));
cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_15-20mrad_1M_%dhitFirstModules_NOoutchi2_1M.root",static_cast<char>(nhits)));
cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_20-25mrad_1M_%dhitFirstModules_NOoutchi2_1M.root",static_cast<char>(nhits)));
cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_25-32mrad_1M_%dhitFirstModules_NOoutchi2_1M.root",static_cast<char>(nhits)));
}

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


   const Int_t NBINS = 14;
   const Int_t NBINS_mu = 14;
   Double_t edges_el[NBINS + 1] = {0.,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.010, 0.015, 0.020, 0.025, 0.032};
   Double_t edges_mu[NBINS_mu + 1] = {0.,0.0001,0.0002,0.0003,0.0004,0.0005,0.0006,0.0007,0.0008,0.0009,0.001,0.002,0.003,0.004,0.005};

TH1D::SetDefaultSumw2(kTRUE);

TH1D *theta_e_clone=new TH1D("theta_e_clone", "Reco elastic event: Electron scattering reco angles from MESMER with clones",NBINS,edges_el);
TH1D *theta_mu_clone=new TH1D("theta_mu_clone", "Reco elastic event: Muon scattering reco angles from MESMER with clones",NBINS_mu,edges_mu);

TH1D *theta_e_clone_quality=new TH1D("theta_e_clone_quality", "Reco elastic event: Electron scattering reco angles from MESMER with clones quality check",NBINS,edges_el);
TH1D *theta_mu_clone_quality=new TH1D("theta_mu_clone_quality", "Reco elastic event: Muon scattering reco angles from MESMER with clones quality check",NBINS_mu,edges_mu);

TH1D *theta_e_gen=new TH1D("theta_e_gen", "Electron scattering generated angles from MESMER",NBINS,edges_el);
TH1D *theta_mu_gen=new TH1D("theta_mu_gen", "Muon scattering generated angles from MESMER",NBINS_mu,edges_mu);

TH1D *theta_fake=new TH1D("theta_fake", "Fakes tracks vs reconstructeds cattering angles",NBINS,edges_el);

TH1D *theta_reco=new TH1D("theta_reco", "Fakes tracks vs reconstructed scattering angles",NBINS,edges_el);
TH1D *theta_e_reco=new TH1D("theta_e_reco", "Reco tracks vs electron scattering reco angles from MESMER",NBINS,edges_el);
TH1D *theta_mu_reco=new TH1D("theta_mu_reco", "Reco tracks: Muon scattering reco angles from MESMER",NBINS_mu,edges_mu);

TH1D *theta_e_good=new TH1D("theta_e_good", "Good tracks vs electron scattering reco angles from MESMER",NBINS,edges_el);
TH1D *theta_mu_good=new TH1D("theta_mu_good", "Good tracks vs euon scattering reco angles from MESMER",NBINS_mu,edges_mu);


	double r_wnorm[6]={20.786765103274643,33.313091221576336,42.396733790329876,50.584815206143652,61.828110400735824,106.88513370134392};



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

// Checking if in the MCTracks container there are elastic particles and if they are in acceptance (4 hits in X modules (== 2 stubs, as TrackerPoints give the hit per sensor) + 4 hits in Y modules + at least 2 hits in UV)

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
	}

// Look at reconstruction if events are reconstructible (all three particles with necessary hits to be potentially reconstructed)

	if(code_mu_in!=-99 and code_e!=-99 and code_mu!=-99 and hit_modXmuin==4 and hit_modYmuin==4 and stereo_muin>1 and hit_modXe==4 and hit_modYe==4 and stereo_e>1 and hit_modXmu==4 and hit_modYmu==4 and stereo_mu>1){

// wnorm is a number needed for normalization when we use different mesmer sample together (like in this case, 6 sample in differen kinematic region of the electron)

double wnorm=99.;

if(the_gen>=0 and the_gen<0.005){index=0; wnorm=r_wnorm[0];}
if(the_gen>=0.005 and the_gen<0.010){index=1; wnorm=r_wnorm[1];}
if(the_gen>=0.010 and the_gen<0.015){index=2; wnorm=r_wnorm[2];}
if(the_gen>=0.015 and the_gen<0.020){index=3; wnorm=r_wnorm[3];}
if(the_gen>=0.020 and the_gen<0.025){index=4; wnorm=r_wnorm[4];}
if(the_gen>=0.025 and the_gen<=0.032){index=5; wnorm=r_wnorm[5];}



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

// Define best vertex and its chi2

MUonERecoOutputVertex vrtx = ReconstructionOutput->bestVertex();
Double_t chi=vrtx.chi2perDegreeOfFreedom();


// Look at reconstructed track container to see if outgoing mu and e are reconstructed
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
		//if the elastic track has a good quality, fill the histos, otherwise count it as fake even if it is an elastic particle
			if(track.fractionOfHitsSharedWithLinkedTrack()>=0.65){
                        theta_mu_good->Fill(thmu_gen,MesmerEvent->wgt_LO*wnorm);
                        theta_e_good->Fill(the_gen,MesmerEvent->wgt_LO*wnorm);
			}else{
			TVector3 p(track.xSlope(),track.ySlope(),1.0);p.Unit();
                        theta_fake->Fill(p.Angle(p_muin),MesmerEvent->wgt_LO*wnorm);}

                 }// end of ID==45 if
             else{
	    // if the particle is not elastic count it has fake
                        TVector3 p(track.xSlope(),track.ySlope(),1.0);p.Unit();
                        theta_fake->Fill(p.Angle(p_muin),MesmerEvent->wgt_LO*wnorm);}
		}// ef of sector1 if

         }//end of for cycle


// Check that the incoming muon is well reconstructed
if(stubs_muin>=5 and chi2_muin<2){//and thmu_gen>0.0002){

theta_mu_gen->Fill(thmu_gen,MesmerEvent->wgt_LO*wnorm);
theta_e_gen->Fill(the_gen,MesmerEvent->wgt_LO*wnorm);


// Plot to build elastic event reco efficiency without quality requirement
if(e>=1 and mu>=1 and chi!=0){
                  	theta_mu_clone->Fill(thmu_gen,MesmerEvent->wgt_LO*wnorm);
                        theta_e_clone->Fill(the_gen,MesmerEvent->wgt_LO*wnorm);}


// Plot to build elastic event reco efficiency with quality requirement

if(find_if(quality_e.begin(),quality_e.end(),[](double i){return i>=0.65;})!=end(quality_e) and find_if(quality_mu.begin(),quality_mu.end(),[](double i){return i>=0.65;})!=end(quality_mu)){

if(e>=1 and mu>=1 and chi!=0){
			theta_mu_clone_quality->Fill(thmu_gen,MesmerEvent->wgt_LO*wnorm);
			theta_e_clone_quality->Fill(the_gen,MesmerEvent->wgt_LO*wnorm);}

			}//if quality
                }//if mu_in
        }//if generated

}//for


// Elastic events reconstruction efficiency with and without quality cuts
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
h0_original->SetTitle("Elastic event reco #epsilon as a funtion of e- angle");
h0_original->SetLineColor(kOrange+2);
h0_original->SetMinimum(0.6);
h0_original->Draw("E");
h0->Draw("E same");
h0->SaveAs(Form("/home/espedica/macros_fairmu/clean_codes/validation/quality_prevertex_el_gen_%dhit_LO.root",static_cast<char>(nhits)));
gStyle->SetOptStat(0);
n7.cd(2);
TGaxis::SetMaxDigits(3);
h1_original->SetTitle("Elastic event reco #epsilon as a funtion of #mu angle");
h1_original->SetLineColor(kOrange+2);
h1_original->SetMinimum(0.6);
h1_original->Draw("E");
h1->Draw("E same");
h1->SaveAs(Form("/home/espedica/macros_fairmu/clean_codes/validation/quality_prevertex_mu_gen_%dhit_LO.root",static_cast<char>(nhits)));
gStyle->SetOptStat(0);
n7.cd(3);
h0_quality->SetTitle("Ratio between (quality cut/without quality track) VS e- angle");
TGaxis::SetMaxDigits(3);
h0_quality->SetMinimum(0.6);
h0_quality->Draw("E");
h0_quality->SaveAs(Form("/home/espedica/macros_fairmu/clean_codes/validation/quality_prevertex_el_clone_%dhit_LO.root",static_cast<char>(nhits)));
gStyle->SetOptStat(0);
n7.cd(4);
h1_quality->SetTitle("Ratio between (quality cut/without quality track) VS #mu angle");
TGaxis::SetMaxDigits(3);
h1_quality->SetMinimum(0.6);
h1_quality->Draw("E");
h1_quality->SaveAs(Form("/home/espedica/macros_fairmu/clean_codes/validation/quality_prevertex_mu_clone_%dhit_LO.root",static_cast<char>(nhits)));
gStyle->SetOptStat(0);


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
h2_good->SaveAs(Form("/home/espedica/macros_fairmu/clean_codes/validation/good_el_%dhit_LO.root",static_cast<char>(nhits)));
n8.cd(2);
h3_good->Draw("E");
h3_good->SetMinimum(0.6);
h3_good->SaveAs(Form("/home/espedica/macros_fairmu/clean_codes/validation/good_mu_%dhit_LO.root",static_cast<char>(nhits)));
n8.cd(3);
h3_fake->Draw("E");
h3_fake->SaveAs(Form("/home/espedica/macros_fairmu/clean_codes/validation/fake_%dhit_LO.root",static_cast<char>(nhits)));
n8.SaveAs(Form("/home/espedica/macros_fairmu/clean_codes/validation/good_fake_%dhit_LO.pdf",static_cast<char>(nhits)));

}

