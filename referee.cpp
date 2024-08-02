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

void referee(){


TChain * cbmsim = new TChain("cbmsim");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/gen_digi_reco/WiP_v0140_commit_258faf6b_0_infmrad_1M.root");

   TH1D* theta_e_bv=new TH1D("theta_e_bv", "Electron scattering reco before vertex angles",175,0.,0.035);
   TH1D* theta_mu_bv=new TH1D("theta_mu_bv", "Muon scattering reco before vertex angles",250,0.,0.005);
   TH1D* theta_e_gen=new TH1D("theta_e", "Electron scattering gen angles",175,0.,0.035);
   TH1D* theta_mu_gen=new TH1D("theta_mu", "Muon scattering gen angles",250,0.,0.005);
   TH1D* theta_e_v=new TH1D("theta_e_v", "Electron scattering reco vrtx angles",175,0.,0.035);
   TH1D* theta_mu_v=new TH1D("theta_mu_v", "Muon scattering reco vrtx angles",250,0.,0.005);

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


TH1D::SetDefaultSumw2(kTRUE);



for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
                cbmsim->GetEntry(i);
                if(i%10000 == 0) cout<<"Entry "<<i<<endl;

int index=99;
int index_mu=99;

int code_mu_in=-99;
int code_e=-99;
int code_mu=-99;
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
	}

	if(code_mu_in!=-99 and code_e!=-99 and code_mu!=-99 and hit_modXmuin==4 and hit_modYmuin==4 and stereo_muin>1 and hit_modXe==4 and hit_modYe==4 and stereo_e>1 and hit_modXmu==4 and hit_modYmu==4 and stereo_mu>1){


double wnorm=99.;


Int_t yes2=0;
Double_t th_inx,th_iny,x0_in,y0_in;
Double_t chi2_muin=999.;
Int_t stubs_muin=0;
Int_t e=0;
Int_t mu=0;
TVector3 p_muin,p_e,p_mu;
Double_t theX_rec,theY_rec,thmuX_rec,thmuY_rec,the_rec_vrtx,thmu_rec_vrtx;
double thmu_rec=99.;
double the_rec=99;
MUonERecoOutputVertex vrtx = ReconstructionOutput->bestVertex();

Double_t chi=vrtx.chi2perDegreeOfFreedom();


vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
int other=0;
int mu_in=0;


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
        else if(track.processIDofLinkedTrack()==45 and track.sector()==1 and chi!=0)
                {
                 if(code_e==track.linkedTrackID() and code_e==vrtx.outgoingElectron().linkedTrackID() and track.index()==vrtx.outgoingElectron().index()) {
                 //if(track.index()==vrtx.outgoingElectron().index()) {
			yes2++; e++; theX_rec=track.xSlope(); theY_rec=track.ySlope();p_e.SetXYZ(theX_rec,theY_rec,1.0); p_e=p_e.Unit();
							the_rec=p_e.Angle(p_muin);
			theta_e_bv->Fill(the_rec,MesmerEvent->wgt_full);theta_e_gen->Fill(the_gen,MesmerEvent->wgt_full);theta_e_v->Fill(vrtx.electronTheta(),MesmerEvent->wgt_full);}
		else if(code_e==track.linkedTrackID() and code_e==vrtx.outgoingMuon().linkedTrackID() and track.index()==vrtx.outgoingMuon().index()) {
                        yes2++; e++; theX_rec=track.xSlope(); theY_rec=track.ySlope();p_e.SetXYZ(theX_rec,theY_rec,1.0); p_e=p_e.Unit();
                                                        the_rec=p_e.Angle(p_muin);
                        theta_e_bv->Fill(the_rec,MesmerEvent->wgt_full);theta_e_gen->Fill(the_gen,MesmerEvent->wgt_full);theta_e_v->Fill(vrtx.muonTheta(),MesmerEvent->wgt_full);}


                 if(code_mu==track.linkedTrackID() and code_mu==vrtx.outgoingMuon().linkedTrackID() and track.index()==vrtx.outgoingMuon().index()) {
                 //if(track.index()==vrtx.outgoingMuon().index()) {
			yes2++; mu++;thmuX_rec=track.xSlope(); thmuY_rec=track.ySlope();p_mu.SetXYZ(thmuX_rec,thmuY_rec,1.0); p_mu=p_mu.Unit();
							thmu_rec=p_mu.Angle(p_muin);
			theta_mu_bv->Fill(thmu_rec,MesmerEvent->wgt_full);theta_mu_gen->Fill(thmu_gen,MesmerEvent->wgt_full);theta_mu_v->Fill(vrtx.muonTheta(),MesmerEvent->wgt_full);}
		else if(code_mu==track.linkedTrackID() and code_mu==vrtx.outgoingElectron().linkedTrackID() and track.index()==vrtx.outgoingElectron().index()){
                        yes2++; mu++;thmuX_rec=track.xSlope(); thmuY_rec=track.ySlope();p_mu.SetXYZ(thmuX_rec,thmuY_rec,1.0); p_mu=p_mu.Unit();
                                                        thmu_rec=p_mu.Angle(p_muin);
			theta_mu_bv->Fill(thmu_rec,MesmerEvent->wgt_full);theta_mu_gen->Fill(thmu_gen,MesmerEvent->wgt_full);theta_mu_v->Fill(vrtx.electronTheta(),MesmerEvent->wgt_full);}

		}
         }//for

        }//if generated
}//for


TCanvas n0("n0","n0",2100,3500);
n0.Divide(1,2);
n0.cd(1);
theta_e_bv->SetLineColor(kRed);
theta_e_gen->SetLineColor(kBlack);
theta_e_v->SetLineColor(kOrange);
theta_e_bv->Draw("hist");
theta_e_gen->Draw("hist same");
theta_e_v->Draw("hist same");
n0.cd(2);
theta_mu_bv->SetLineColor(kRed);
theta_mu_gen->SetLineColor(kBlack);
theta_mu_v->SetLineColor(kOrange);
theta_mu_bv->Draw("hist");
theta_mu_gen->Draw("hist same");
theta_mu_v->Draw("hist same");
gPad->SetLogy();
n0.SaveAs("theta_full.pdf");

}

