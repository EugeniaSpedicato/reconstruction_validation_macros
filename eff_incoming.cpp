
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

	TFile *inputfile = new TFile("test.root");//TRMesmer_box_100k.root");//TRMesmer_beamProfile_mu-.root");
        TTree* cbmsim = (TTree*) inputfile->Get("cbmsim");

        TClonesArray *MCTrack = 0;
        TClonesArray *SignalTracks = 0;
        TClonesArray *TrackerStripDigis = 0;
        TClonesArray *TrackerPoints = 0;
        TClonesArray *TrackerStubs = 0;
        MuE::Event *MesmerEvent = 0;
        MUonERecoOutput *ReconstructionOutput = 0;

        cbmsim->SetBranchAddress("MCTrack", &MCTrack);
        cbmsim->SetBranchAddress("SignalTracks", &SignalTracks);
        cbmsim->SetBranchAddress("TrackerPoints", &TrackerPoints);
        cbmsim->SetBranchAddress("TrackerStripDigis", &TrackerStripDigis);
        cbmsim->SetBranchAddress("TrackerStubs", &TrackerStubs);
        cbmsim->SetBranchAddress("MesmerEvent", &MesmerEvent);
        cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);


double thmu_in_gen=0;
double thmu_in_rec=0;
double thmu_in_gen_x=0; double thmu_in_gen_y=0;
double thmu_in_rec_x=0; double thmu_in_rec_y=0;

double signal=0.;
double signal_el=0.;
double reco=0.; double reco_el=0.;
int yes_mu_in_g=0; int yes_mu_g=0; int yes_e_g=0;
int point_mu_in=0;int point_mu=0;int point_el=0;
int code_mu=-99;int code_e=-99;

double X,Y,mx,my,p_in;

TH2D* pos_in=new TH2D("pos_in","Position incoming muon at Z=target in cm",200,-10,10,200,-10,10);
TH2D* pos_in_not=new TH2D("pos_in_not","Position incoming muon at Z=target in cm when not reco",200,-10,10,200,-10,10);

TH2D* posx_mx=new TH2D("posx_mx","Position X versus slope X incoming muon at Z=target in cm",200,-10,10,80,-0.004,0.004);
TH2D* posx_mx_not=new TH2D("posx_mx_not","Position X versus slope X incoming muon at Z=target in cm when not reco",200,-10,10,80,-0.004,0.004);

TH1D* energy=new TH1D("energy","Energy of incoming muon in GeV",300,0,200);
TH1D* energy_not=new TH1D("energy_not","Energy of incoming muon in GeV when not reco",360,0,180);


for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
		cbmsim->GetEntry(i);
		if(i%1000 == 0) cout<<"Entry "<<i<<endl;

vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();

	TVector3 pmuin,pe,pmu;

	double Emu_in=0;
	TVector3 muin;
	TVector3 pmuin_dir=pmuin.Unit();

	for(int n = 0; n < MCTrack->GetEntries(); n++) {
	 const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
	 if(MCTr->interactionID()==0 and MCTr->pdgCode()==13){
	 pmuin.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz());
         energy->Fill(pmuin.Mag());//,MesmerEvent->wgt_full);

	 pmuin=pmuin.Unit();
	thmu_in_gen = pmuin.Theta();
	 Emu_in = MCTr->energy();
	 p_in=MCTr->p();
         thmu_in_gen_x=MCTr->ax();
         thmu_in_gen_y=MCTr->ay();
         muin.SetXYZ(thmu_in_gen_x,thmu_in_gen_y,1.0);


           mx=MCTr->ax();
           my=MCTr->ay();
           double qx=MCTr->bx();
           double qy=MCTr->by();

		X=qx+(mx*1031.7);
                Y=qy+(my*1031.7);

	 }
         if(MCTr->interactionID()==45 and MCTr->pdgCode()==11) {code_e=n;}
         if(MCTr->interactionID()==45 and MCTr->pdgCode()==13) {code_mu=n;}

	}

           int last_modXmu_in=0;
           int last_modYmu_in=0;
	   int stereo_mu_in=0;

                         for(int s=0; s<TrackerPoints->GetEntries(); s++)
                         {const MUonETrackerPoint *TrackerPt = static_cast<const MUonETrackerPoint*>(TrackerPoints->At(s));
                          if(TrackerPt->trackPDGCode()==13 and TrackerPt->trackID()==0 and TrackerPt->stationID()==0){ point_mu_in++;
												 if(TrackerPt->moduleID()==4) last_modXmu_in++;
                                                                                                 if(TrackerPt->moduleID()==5) last_modYmu_in++;
												 if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_mu_in++;}
                         }


if(last_modXmu_in==2 and last_modYmu_in==2 and stereo_mu_in>1){
		yes_mu_in_g=1;
			}

	 if(yes_mu_in_g!=1) {pos_in_not->Fill(X,Y,MesmerEvent->wgt_full);
			     posx_mx_not->Fill(X,mx,MesmerEvent->wgt_full);
		             //energy_not->Fill(Emu_in,MesmerEvent->wgt_full);
}

//energy->Fill(Emu_in,MesmerEvent->wgt_full);

int yes_mu_in=0;
	if( yes_mu_in_g==1){
	   signal+=1;
	   pos_in->Fill(X,Y,MesmerEvent->wgt_full);
	   posx_mx->Fill(X,mx,MesmerEvent->wgt_full);
	   //energy->Fill(Emu_in,MesmerEvent->wgt_full);

TVector3 in;
for(int j=0; j<tracks.size();j++)
{
if(tracks.at(j).processIDofLinkedTrack()==0 and tracks.at(j).sector()==0){
thmu_in_rec_x=tracks.at(j).xSlope();
thmu_in_rec_y=tracks.at(j).ySlope();
in.SetXYZ(thmu_in_rec_x,thmu_in_rec_y,1.0);
thmu_in_rec=in.Theta();
	}
}

for(int j=0; j<tracks.size();j++)
{
 if(tracks.at(j).processIDofLinkedTrack()==0 and tracks.at(j).sector()==0)
	{

		 if(0==tracks.at(j).linkedTrackID()) { yes_mu_in=1;
		 TVector3 mu_outv(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0);
		TVector3 z(0.,0.,1.);
                mu_outv=mu_outv.Unit();
		thmu_in_rec_x=tracks.at(j).xSlope();
                thmu_in_rec_y=tracks.at(j).ySlope();
                 thmu_in_rec=sqrt(thmu_in_rec_x*thmu_in_rec_x+thmu_in_rec_y*thmu_in_rec_y);
					}
	}
}


if(yes_mu_in==1){reco+=1;}

}

	for(int n = 0; n < SignalTracks->GetEntries(); n++) {
        const MUonETrack *SigTracks = static_cast<const MUonETrack*>(SignalTracks->At(n));
         if(SignalTracks->GetEntries()>1 and SigTracks->interactionID()==45 )
         {
           int last_modXmu=0; int last_modXe=0;
           int last_modYmu=0; int last_modYe=0;
           int stereo_mu=0; int stereo_e=0;

                         for(int s=0; s<TrackerPoints->GetEntries(); s++)
                         {const MUonETrackerPoint *TrackerPt = static_cast<const MUonETrackerPoint*>(TrackerPoints->At(s));
                          if(TrackerPt->trackPDGCode()==11 and SigTracks->pdgCode()==11 and TrackerPt->trackID()==code_e and TrackerPt->stationID()==1){ point_el++;
                                                                                                 if(TrackerPt->moduleID()==4) last_modXe++;
                                                                                                 if(TrackerPt->moduleID()==5)last_modYe++;
                                                                                                 if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_e++;
                                                                                                }
                          if(TrackerPt->trackPDGCode()==13 and SigTracks->pdgCode()==13 and TrackerPt->trackID()==code_mu and TrackerPt->stationID()==1){ point_mu++;
                                                                                                 if(TrackerPt->moduleID()==4) last_modXmu++;
                                                                                                 if(TrackerPt->moduleID()==5) last_modYmu++;
                                                                                                 if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_mu++;
                                                                                                }
                         }

    if(SigTracks->pdgCode()==11 and last_modXe==2 and last_modYe==2 and stereo_e>1){yes_e_g=1;}
    if(SigTracks->pdgCode()==13 and last_modXmu==2 and last_modYmu==2 and stereo_mu>1){yes_mu_g=1;}
         }
	}

        if(yes_e_g==1 and yes_mu_g==1 and yes_mu_in==1 ){
           signal_el+=MesmerEvent->wgt_full;

int yes_e=0;
int yes_mu=0;

for(int j=0; j<tracks.size();j++)
{
 if(tracks.at(j).processIDofLinkedTrack()==45 and tracks.size()>=3 and tracks.at(j).sector()==1 and tracks.at(j).percentageOfHitsSharedWithLinkedTrack()>=0)
        {

                 if(code_e==tracks.at(j).linkedTrackID()) {yes_e++;}

                 if(code_mu==tracks.at(j).linkedTrackID()) {yes_mu++;}
        }
}


if(yes_mu>=1 and yes_e>=1 ){reco_el+=MesmerEvent->wgt_full;}


		}


yes_mu_in_g=0;
code_mu=-99;
code_e=-99;
yes_mu_g=0;
yes_e_g=0;
point_mu_in=0;
point_el=0;
point_mu=0;


} //end of general for

double ratio =reco/signal;

cout << "Su " << cbmsim->GetEntries() << " muoni entranti, " << signal << " sono ricostruibili: " << (signal/cbmsim->GetEntries())*100 << "%"<< endl;
cout << "Su " << signal << " muoni entranti ricostruibili, " << reco << " sono ricostruiti: " << (reco/signal)*100 << "%"<< endl;
cout << "Su " << reco << " muoni entranti ricostruiti, " << signal_el << " sono gli eventi elastici di segnale ricostruibili: " << (signal_el/reco)*100 << "%"<< endl;
cout << "Su " << reco << " muoni entranti ricostruiti, " << reco_el << " sono gli eventi elastici di segnale ricostruiti: " << (reco_el/reco)*100 << "%"<< endl;
cout << "Su " << signal_el << " eventi di segnale ricostruibili, " << reco_el << " sono ricostruiti: " << (reco_el/signal_el)*100 << "%"<< endl;

Int_t nxTh = pos_in->GetNbinsX();
Int_t nyTh = pos_in->GetNbinsY();
for (Int_t i=1; i<nxTh+1; i++) {
for (Int_t j=1; j<nyTh+1; j++) {
if (pos_in->GetBinContent(i,j)<1) pos_in->SetBinContent(i,j,0);}}

Int_t nxTh_not = pos_in_not->GetNbinsX();
Int_t nyTh_not = pos_in_not->GetNbinsY();
for (Int_t i=1; i<nxTh_not+1; i++) {
for (Int_t j=1; j<nyTh_not+1; j++) {
if (pos_in_not->GetBinContent(i,j)<1) pos_in_not->SetBinContent(i,j,0);}}


Int_t nxTh1 = posx_mx->GetNbinsX();
Int_t nyTh1 = posx_mx->GetNbinsY();
for (Int_t i=1; i<nxTh1+1; i++) {
for (Int_t j=1; j<nyTh1+1; j++) {
if (posx_mx->GetBinContent(i,j)<1) posx_mx->SetBinContent(i,j,0);}}

Int_t nxTh_not1 = posx_mx_not->GetNbinsX();
Int_t nyTh_not1 = posx_mx_not->GetNbinsY();
for (Int_t i=1; i<nxTh_not1+1; i++) {
for (Int_t j=1; j<nyTh_not1+1; j++) {
if (posx_mx_not->GetBinContent(i,j)<1) posx_mx_not->SetBinContent(i,j,0);}}

TCanvas a("pos","pos",700,700);
a.Divide(1,3);
a.cd(1);
pos_in->Draw();
pos_in_not->SetMarkerColor(kRed);
pos_in_not->Draw("same");
a.cd(2);
posx_mx->Draw();
posx_mx_not->SetMarkerColor(kRed);
posx_mx_not->Draw("same");
a.cd(3);
energy->Draw("hist");
energy_not->SetLineColor(kRed);
energy_not->Draw("hist same");
//gPad->SetLogy();
a.SaveAs("ingo_muin.pdf");
}
