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

void vrtx_wrong_correct(int nhits){

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


                auto pos_on_track = [](double q, double m, double z){return (q + m*z);};
   const Int_t NBINS = 14;
   const Int_t NBINS_mu = 14;
   Double_t edges_el[NBINS + 1] = {0.,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.010, 0.015, 0.020, 0.025, 0.032};
   Double_t edges_mu[NBINS_mu + 1] = {0.,0.0001,0.0002,0.0003,0.0004,0.0005,0.0006,0.0007,0.0008,0.0009,0.001,0.002,0.003,0.004,0.005};

TH1D::SetDefaultSumw2(kTRUE);


TH1D *h_opening_clone=new TH1D("h_opening_clone", "Opening angle reco events from MESMER with clones",33,0.,0.033);

TH1D *h_opening_clone_vrtx=new TH1D("h_opening_clone_vrtx", "True opening angle reco events from MESMER with clones when vrtx is created with good quality",33,0.,0.033);

TH1D *h_opening_wrong_vrtx=new TH1D("h_opening_wrong_vrtx", "True opening angle reco events from MESMER with clones when vrtx is created with wrong association",33,0.,0.033);



	double r_wnorm[6]={20.786765103274643,33.313091221576336,42.396733790329876,50.584815206143652,61.828110400735824,106.88513370134392};

double n_wrong=0.;

double sum_wgt=0.;
double sumW2[6]={0.};

double n_clones_el=0.;
double n_clones_mu=0.;
double total=0.;
double no_muin=0.;
double all=0.;
double n_el2=0.;
double n_el3=0.;
double n_cl=0.;
double n_one_el=0.;
double n_other=0.;
double n_zero=0.;
double n_both=0.;
double tot=0.;

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
     double E_e=0.;


// Checking if in the MCTracks container there are elastic particles and if they are in acceptance (4 hits in X modules (== 2 stubs, as TrackerPoints give the hit per sensor) + 4 hits in Y modules + at least 2 hits in UV)

        for(int n = 0; n < MCTrack->GetEntries(); n++) {
         const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
         if(MCTr->interactionID()==0 and MCTr->pdgCode()==-13) {code_mu_in=n; p_muin_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); p_muin_MC.Unit();E_e=MCTr->energy();
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

double wnorm=99.;

// wnorm is a number needed for normalization when we use different mesmer sample together (like in this case, 6 sample in differen kinematic region of the electron)
if(the_gen>=0 and the_gen<0.005){index=0; wnorm=r_wnorm[0];}
if(the_gen>=0.005 and the_gen<0.010){index=1; wnorm=r_wnorm[1];}
if(the_gen>=0.010 and the_gen<0.015){index=2; wnorm=r_wnorm[2];}
if(the_gen>=0.015 and the_gen<0.020){index=3; wnorm=r_wnorm[3];}
if(the_gen>=0.020 and the_gen<0.025){index=4; wnorm=r_wnorm[4];}
if(the_gen>=0.025 and the_gen<=0.032){index=5; wnorm=r_wnorm[5];}

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
Double_t the_rec,theX_rec,theY_rec,thmu_rec,thmuX_rec,thmuY_rec;

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
        if(track.processIDofLinkedTrack()==45 and track.sector()==1)
                {
                 if(code_e==track.linkedTrackID()) {yes2++; e++; theX_rec=track.xSlope(); theY_rec=track.ySlope();p_e.SetXYZ(theX_rec,theY_rec,1.0); p_e=p_e.Unit();the_rec=p_e.Angle(p_muin);
                                                quality_e.push_back(track.fractionOfHitsSharedWithLinkedTrack());}
                 if(code_mu==track.linkedTrackID()) {yes2++; mu++;thmuX_rec=track.xSlope(); thmuY_rec=track.ySlope();p_mu.SetXYZ(thmuX_rec,thmuY_rec,1.0); p_mu=p_mu.Unit();thmu_rec=p_mu.Angle(p_muin);
                                                quality_mu.push_back(track.fractionOfHitsSharedWithLinkedTrack());}
                }
        if(track.processIDofLinkedTrack()!=45 and track.sector()==1){other++;}

         }// end for cycle


// Check that the incoming muon is well reconstructed

if(stubs_muin>=5 and chi2_muin<2){//and thmu_gen>0.0002){//and E_e>1){



// Check that at least one muon and one electron are reconstructed with a good quality -0.65 means that at least 4 hits are shared with the linked MC track-

 if(e>=1 and mu>=1 and chi!=0 and find_if(quality_e.begin(),quality_e.end(),[](double i){return i>=0.65;})!=end(quality_e) and find_if(quality_mu.begin(),quality_mu.end(),[](double i){return i>=0.65;})!=end(quality_mu)){


 h_opening_clone->Fill(opening_angle,MesmerEvent->wgt_LO*wnorm);

// If the particles associeted to the vrtx tracks are 1 mu + 1 e (also with swapped identity), then the vrtx did it well, otherwise count it as wrong
  if(((vrtx.outgoingMuon().linkedTrackID()==code_e and vrtx.outgoingElectron().linkedTrackID()==code_mu) or (vrtx.outgoingMuon().linkedTrackID()==code_mu and vrtx.outgoingElectron().linkedTrackID()==code_e))){// and vrtx.outgoingMuon().fractionOfHitsSharedWithLinkedTrack()>=0.65 and vrtx.outgoingElectron().fractionOfHitsSharedWithLinkedTrack()>=0.65){
	 h_opening_clone_vrtx->Fill(opening_angle,MesmerEvent->wgt_LO*wnorm);
	}
  else{   n_wrong+=MesmerEvent->wgt_LO*wnorm;
          h_opening_wrong_vrtx->Fill(opening_angle,MesmerEvent->wgt_LO*wnorm);}

 }

                }//if mu_in

        }//if generated

}//for



TH1D * h2vrtx = (TH1D*) h_opening_clone_vrtx->Clone();
TH1D * h2_clone = (TH1D*) h_opening_clone->Clone();
h2vrtx->Divide(h2vrtx,h2_clone,1,1,"B");

TH1D * h2vrtx_w = (TH1D*) h_opening_wrong_vrtx->Clone();
h2vrtx_w->Divide(h2vrtx_w,h2_clone,1,1,"B");



TCanvas a("a","a",1400,700);
a.Divide(1,2);
a.cd(1);
h2vrtx->Draw();
a.cd(2);
h2vrtx_w->Draw();
a.SaveAs(Form("/home/espedica/macros_fairmu/clean_codes/validation/vrtx_LO_%dhitFirstModules_quality.pdf",static_cast<char>(nhits)));

h2vrtx->SaveAs(Form("/home/espedica/macros_fairmu/clean_codes/validation/vrtx_correct_LO_%dhitFirstModules_true_clone_quality.root",static_cast<char>(nhits)));
h2vrtx_w->SaveAs(Form("/home/espedica/macros_fairmu/clean_codes/validation/vrtx_wrong_LO_%dhitFirstModules_true_clone_quality.root",static_cast<char>(nhits)));


}

