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

void fake_rates_and_quality_prevrtx(const char* reco_filename, const char* gen_filename, int nhits, const char* index, const char* info, string path){

   TFile *f1 = new TFile(reco_filename);
   TFile *f2 = new TFile(gen_filename);
   TTree *cbmsim = (TTree*)f1->Get("cbmsim");
   TTree *t2 = (TTree*)f2->Get("cbmsim");

        t2->SetEntries(cbmsim->GetEntries());

cout << "entries cbmsim " << cbmsim->GetEntries() << endl;
cout << "entries t2 " << t2->GetEntries() << endl;


   cbmsim->AddFriend(t2);
   cout << t2->GetEntries() << endl;
   cout << cbmsim->GetEntries() << endl;




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


        double r_wnorm[6]={18.558632798881270,31.156709183434657,41.306609022934722,50.678758808390178,61.782157006784871,106.98964295859898};



for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
                cbmsim->GetEntry(i);
                if(i%100 == 0) cout<<"Entry "<<i<<endl;

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

size_t sz=sizeof(index);
if( strncmp(index,"0_5",sz) ){wnorm=r_wnorm[0];}
else if( strncmp(index,"5_10",sz) ){wnorm=r_wnorm[1];}
else if( strncmp(index,"10_15",sz) ){wnorm=r_wnorm[2];}
else if( strncmp(index,"15_20",sz) ){wnorm=r_wnorm[3];}
else if( strncmp(index,"20_25",sz) ){wnorm=r_wnorm[4];}
else if( strncmp(index,"25_32",sz) ){wnorm=r_wnorm[5];}

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


theta_e_reco->SaveAs(Form("%s/reco_el_%dhit_LO_%s_%s.root",path.c_str(),static_cast<char>(nhits),index,info));
theta_mu_reco->SaveAs(Form("%s/reco_mu_%dhit_LO_%s_%s.root",path.c_str(),static_cast<char>(nhits),index,info));
theta_reco->SaveAs(Form("%s/reco_%dhit_LO_%s_%s.root",path.c_str(),static_cast<char>(nhits),index,info));

theta_e_good->SaveAs(Form("%s/good_el_%dhit_LO_%s_%s.root",path.c_str(),static_cast<char>(nhits),index,info));
theta_mu_good->SaveAs(Form("%s/good_mu_%dhit_LO_%s_%s.root",path.c_str(),static_cast<char>(nhits),index,info));
theta_fake->SaveAs(Form("%s/fake_%dhit_LO_%s_%s.root",path.c_str(),static_cast<char>(nhits),index,info));


}

