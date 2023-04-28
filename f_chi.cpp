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

 	TFile *inputfile = new TFile("TRMesmer_box_100k_2GeV.root");
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


TH2D *h_res=new TH2D("res", "(the_rec-the_true) VS energy Ee<5 GeV POST-VRTX",25,0,5,400,-0.04,0.04);
TH2D *h_res1=new TH2D("res1", "(the_rec-the_true) VS energy 5<E<15 GeV POST-VRTX",20,5,15,80,-0.004,0.004);
TH2D *h_res2=new TH2D("res2", "(the_rec-the_true) VS energy 15<E<25 GeV POST-VRTX",10,15,20,20,-0.001,0.001);
TH2D *h_res3=new TH2D("res3", "(the_rec-the_true) VS energy Ee>25 GeV POST-VRTX",6,25,65,20,-0.001,0.001);

TH2D *h_res_mu=new TH2D("h_res_mu", "(thmu_rec-thmu_true) VS energy Emu(GeV) POST-VRTX",60,157,160,22,-0.0006,0.0006);//30,140,170,22,-0.0006,0.0006);
TH2D *h_res_mu_in=new TH2D("h_res_mu_in", "(thmu_in_rec-thmu_in_true) VS energy Emu(GeV) POST-VRTX",4,160,160.2,22,-0.0006,0.0006);



TH2D *h_resTR=new TH2D("resTR", "(the_rec-the_true) VS energy Ee<5 GeV PRE-VRTX",25,0,5,400,-0.04,0.04);
TH2D *h_res1TR=new TH2D("res1TR", "(the_rec-the_true) VS energy 5<E<15 GeV PRE-VRTX",20,5,15,80,-0.004,0.004);
TH2D *h_res2TR=new TH2D("res2TR", "(the_rec-the_true) VS energy 15<E<25 GeV PRE-VRTX",10,15,20,20,-0.001,0.001);
TH2D *h_res3TR=new TH2D("res3TR", "(the_rec-the_true) VS energy Ee>25 GeV PRE-VRTX",6,25,65,20,-0.001,0.001);

TH2D *h_res_muTR=new TH2D("h_res_muTR", "(thmu_rec-thmu_true) VS energy Emu(GeV) PRE-VRTX",40,150,165,22,-0.0006,0.0006);
TH2D *h_res_mu_inTR=new TH2D("h_res_mu_inTR", "(thmu_in_rec-thmu_in_true) VS energy Emu(GeV) PRE-VRTX",4,160,160.2,22,-0.0006,0.0006);


TH2D *h_echi_e=new TH2D("h_chidoe", "chi2 per dof vs energy electron Ee<5 GeV POST-VRTX",50,0,5,300,0,100);
TH2D *h_echi_mu=new TH2D("h_chidomu", "chi2 per dof vs energy muon Ee<5 GeV POST-VRTX",50,0,5,300,0,100);

TH2D *h_echi_e1=new TH2D("h_chidoe1", "chi2 per dof vs energy electron 5<E<15 GeV POST-VRTX",20,5,15,300,0,100);
TH2D *h_echi_mu1=new TH2D("h_chidomu1", "chi2 per dof vs energy muon 5<E<15 GeV POST-VRTX",20,5,15,300,0,100);

TH2D *h_echi_e2=new TH2D("h_chidoe2", "chi2 per dof vs energy electron 15<E<25 GeV POST-VRTX",20,15,25,300,0,100);
TH2D *h_echi_mu2=new TH2D("h_chidomu2", "chi2 per dof vs energy muon 15<E<25 GeV POST-VRTX",20,15,25,300,0,100);

TH2D *h_echi_e3=new TH2D("h_chidoe3", "chi2 per dof vs energy electron Ee>25 GeV POST-VRTX",6,25,65,300,0,100);
TH2D *h_echi_mu3=new TH2D("h_chidomu3", "chi2 per dof vs energy muon Ee>25 GeV POST-VRTX",6,25,65,300,0,100);


TH2D *h_echi_eTR=new TH2D("h_chidoeTR", "chi2 per dof vs energy electron Ee<5 GeV PRE-VRTX",50,0,5,300,0,100);
TH2D *h_echi_muTR=new TH2D("h_chidomuTR", "chi2 per dof vs energy muon Ee<5 GeV PRE-VRTX",50,0,5,300,0,100);

TH2D *h_echi_e1TR=new TH2D("h_chidoe1TR", "chi2 per dof vs energy electron 5<E<15 GeV PRE-VRTX",20,5,15,300,0,100);
TH2D *h_echi_mu1TR=new TH2D("h_chidomu1TR", "chi2 per dof vs energy muon 5<E<15 GeV PRE-VRTX",20,5,15,300,0,100);

TH2D *h_echi_e2TR=new TH2D("h_chidoe2TR", "chi2 per dof vs energy electron 15<E<25 GeV PRE-VRTX",20,15,25,300,0,100);
TH2D *h_echi_mu2TR=new TH2D("h_chidomu2TR", "chi2 per dof vs energy muon 15<E<25 GeV PRE-VRTX",20,15,25,300,0,100);

TH2D *h_echi_e3TR=new TH2D("h_chidoe3TR", "chi2 per dof vs energy electron Ee>25 GeV PRE-VRTX",6,25,65,300,0,100);
TH2D *h_echi_mu3TR=new TH2D("h_chidomu3TR", "chi2 per dof vs energy muon Ee>25 GeV PRE-VRTX",6,25,65,300,0,100);



TH1D *h_chi_e=new TH1D("chidoe", "tracks chi2 per dof electron Ee<5 GeV",300,0,100);
TH1D *h_chi_mu=new TH1D("chidomu", "tracks chi2 per dof muon Ee<5 GeV",300,0,100);

TH1D *h_chi_e1=new TH1D("chidoe1", "tracks chi2 per dof electron 5<E<15 GeV",300,0,100);
TH1D *h_chi_mu1=new TH1D("chidomu1", "tracks chi2 per dof  muon 5<E<15 GeV",300,0,100);

TH1D *h_chi_e2=new TH1D("chidoe2", "tracks chi2 per dof electron 15<E<25 GeV",300,0,100);
TH1D *h_chi_mu2=new TH1D("chidomu2", "tracks chi2 per dof muon 15<E<25 GeV",300,0,100);

TH1D *h_chi_e3=new TH1D("chidoe3", "tracks chi2 per dof electron Ee>25 GeV",300,0,100);
TH1D *h_chi_mu3=new TH1D("chidomu3", "tracks chi2 per dof muon Ee>25 GeV",300,0,100);

TH1D *h_chi_eTOT=new TH1D("chidoeTOT", "tracks chi2 per dof electron all Ee",300,0,100);
TH1D *h_chi_muTOT=new TH1D("chidomuTOT", "tracks chi2 per dof muon all Ee",300,0,100);

TH1D *h_chi_eTR=new TH1D("chidoeTR", "tracks chi2 per dof electron Ee<5 GeV",300,0,100);
TH1D *h_chi_muTR=new TH1D("chidomuTR", "tracks chi2 per dof muon Ee<5 GeV",300,0,100);

TH1D *h_chi_e1TR=new TH1D("chidoe1TR", "tracks chi2 per dof electron 5<E<15 GeV",300,0,100);
TH1D *h_chi_mu1TR=new TH1D("chidomu1TR", "tracks chi2 per dof  muon 5<E<15 GeV",300,0,100);

TH1D *h_chi_e2TR=new TH1D("chidoe2TR", "tracks chi2 per dof electron 15<E<25 GeV",300,0,100);
TH1D *h_chi_mu2TR=new TH1D("chidomu2TR", "tracks chi2 per dof muon 15<E<25 GeV",300,0,100);

TH1D *h_chi_e3TR=new TH1D("chidoe3TR", "tracks chi2 per dof electron Ee>25 GeV",300,0,100);
TH1D *h_chi_mu3TR=new TH1D("chidomu3TR", "tracks chi2 per dof muon Ee>25 GeV",300,0,100);


TH1D *angle_e=new TH1D("angle_e","Electron angle when not correctly reconstructed",1000,0.,0.03);
TH1D *angle_mu=new TH1D("angle_mu","Muon angle when not correctly reconstructed",200,0.,0.005);
TH1D *energy_e=new TH1D("energy_e","Electron enrgy when not correctly reconstructed",180,0.,60.);
TH2D *h_2d_e=new TH2D("h_2d_e","Electron angle Vs Ee when not correctly reconstructed",1000,0.,0.03,500,0.,100.);
TH2D *h_2d_emu=new TH2D("h_2d_emu","Electron angle Vs muon angle when not correctly reconstructed",1000,0.,0.03,200,0.,0.005);



double the_gen=0; double thmu_gen=0; double th_in_gen=0;
double the_rec_vrtx=0; double thmu_rec_vrtx=0; double th_in_rec_vrtx=0;
double the_rec_tracks=0; double thmu_rec_tracks=0; double th_in_rec_tracks=0;

double signal=0.; double reco=0.; double reco1=0.; double more_reco=0.; double reco0=0.; double reco3=0.;
double reco_v=0.; double more_reco_v=0.; double reco0_v=0.; double one=0;
double e=0.; double mu=0.;

int yes_e=0;int yes_mu=0; int yes2=0; int yes_v=0;
int point_mu=0; int point_el=0;
int code_mu=-99; int code_e=-99;
int TrackIdreco=-99;


for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
		cbmsim->GetEntry(i);
		if(i%1000 == 0) cout<<"Entry "<<i<<endl;
	TVector3 pmuin,pe,pmu;

	double Emu_in=0;
        double thmu_sdr,the_sdr;
        double Ee=0;
        double Emu=0;
	TVector3 muin;
	TVector3 pmuin_dir=pmuin.Unit();
        TVector3 pmu_dir,pe_dir;


	for(int n = 0; n < MCTrack->GetEntries(); n++) {
	 const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
	 if(MCTr->interactionID()==0 and MCTr->pdgCode()==-13){
	 pmuin.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz());
	 th_in_gen = pmuin.Theta();
	 Emu_in = MCTr->energy();
         double mx=MCTr->ax();
         double my=MCTr->ay();
         muin.SetXYZ(mx,my,1.0);
	 }

         if(MCTr->interactionID()==45 and MCTr->pdgCode()==11) {code_e=n;}
         if(MCTr->interactionID()==45 and MCTr->pdgCode()==-13) {code_mu=n;}
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
                                                                                                 if(TrackerPt->moduleID()==5) last_modYe++;
												 if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_e++;}
                          if(TrackerPt->trackPDGCode()==-13 and SigTracks->pdgCode()==-13 and TrackerPt->trackID()==code_mu and TrackerPt->stationID()==1){ point_mu++;
												 if(TrackerPt->moduleID()==4) last_modXmu++;
                                                                                                 if(TrackerPt->moduleID()==5) last_modYmu++;
												 if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_mu++;}
                         }

           if(SigTracks->pdgCode()==11 and last_modXe==2 and last_modYe==2 and stereo_e>1){
		pe.SetXYZ(SigTracks->px(),SigTracks->py(),SigTracks->pz()); the_gen=pe.Theta();
		Ee=SigTracks->energy();
		double mx1=SigTracks->ax();
		double my1=SigTracks->ay();
		pe_dir=pe.Unit();
		TVector3 e(mx1,my1,1.0);
                the_sdr=muin.Angle(e);
		//the_sdr=acos(pmuin_dir.Dot(pe_dir));
		yes_e=1;}
 	   if(SigTracks->pdgCode()==-13 and last_modXmu==2 and last_modYmu==2 and stereo_mu>1){
		Emu=SigTracks->energy();
		pmu.SetXYZ(SigTracks->px(),SigTracks->py(),SigTracks->pz());
		//thmu_gen=pmu.Theta();
		double mx1=SigTracks->ax();
		double my1=SigTracks->ay();
                pmu_dir=pmu.Unit();
		TVector3 mu(mx1,my1,1.0);
                thmu_gen=mu.Theta();
		thmu_sdr=muin.Angle(mu);
		//thmu_sdr=acos(pmuin_dir.Dot(pmu_dir));
		yes_mu=1;}
	 }
	}

 if(yes_e!=1 or yes_mu!=1) cout << "NOT RECONSTRUCTIBLE" << endl;

 if(yes_e==1 and yes_mu==1){
cout << "RECONSTRUCTIBLE" << endl;

	   signal+=MesmerEvent->wgt_full;

vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
vector<MUonERecoOutputVertex> vrtx = ReconstructionOutput->reconstructedVertices();



for(int j=0; j<vrtx.size();j++)
{
if(vrtx.at(j).stationIndex()==1)
{
 MUonERecoOutputTrack mu_in = vrtx.at(j).incomingMuon();
 MUonERecoOutputTrack mu_out = vrtx.at(j).outgoingMuon();
 MUonERecoOutputTrack e_out = vrtx.at(j).outgoingElectron();

TVector3 in1(mu_in.xSlope(),mu_in.ySlope(),1.0);
th_in_rec_vrtx=in1.Theta();


// scelgo il vertice con miglior chi^2
 if(j==0){
// quante volte sbaglio associazione ID traccia vertice con ID traccia MC?
if(code_e!=e_out.linkedTrackID()) {e+=MesmerEvent->wgt_full;
					cout << "el " << code_e << ", " << e_out.linkedTrackID() << " on process id " << e_out.processIDofLinkedTrack() << endl;
				  }
if(code_mu!=mu_out.linkedTrackID()) {mu+=MesmerEvent->wgt_full;
					cout << "mu " << code_mu << ", " << mu_out.linkedTrackID() << " on process id " << mu_out.processIDofLinkedTrack() << endl;
				  }

//associazione corretta e un certo numero di tracce riscostruite
if(code_mu==mu_out.linkedTrackID() and code_e==e_out.linkedTrackID() and tracks.size()==3){
	yes_v++;

TVector3 mu1(mu_out.xSlope(),mu_out.ySlope(),1.0);
TVector3 e1(e_out.xSlope(),e_out.ySlope(),1.0);


	the_rec_vrtx = e1.Theta();
        thmu_rec_vrtx = mu1.Theta();


	if( Ee<5 ){
	h_echi_e->Fill(Ee,e_out.chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);
        h_echi_mu->Fill(Ee,mu_out.chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);

        h_chi_e->Fill(e_out.chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);
        h_chi_mu->Fill(mu_out.chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);
	}
        if( Ee<15 and Ee>5){
	h_echi_e1->Fill(Ee,e_out.chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);
        h_echi_mu1->Fill(Ee,mu_out.chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);

	h_chi_e1->Fill(e_out.chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);
        h_chi_mu1->Fill(mu_out.chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);
	}
        if( Ee<25 and Ee>15){
	h_echi_e2->Fill(Ee,e_out.chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);
        h_echi_mu2->Fill(Ee,mu_out.chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);

        h_chi_e2->Fill(e_out.chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);
        h_chi_mu2->Fill(mu_out.chi2perDegreeOfFreedom(),MesmerEvent->wgt_full); 
	}
        if( Ee>25){
	h_echi_e3->Fill(Ee,e_out.chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);
        h_echi_mu3->Fill(Ee,mu_out.chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);

	h_chi_e3->Fill(e_out.chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);
        h_chi_mu3->Fill(mu_out.chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);
	}

		}

// quante volte sbaglia contemporaneamente l'associazione, quindi scambia e per mu?
  else if(code_mu!=mu_out.linkedTrackID() and code_e!=e_out.linkedTrackID()){one+=MesmerEvent->wgt_full;}

		}
	}
}

TVector3 in;
for(int j=0; j<tracks.size();j++)
{
if(tracks.at(j).processIDofLinkedTrack()==0 and tracks.at(j).sector()==0){
double th_inx=tracks.at(j).xSlope();
double th_iny=tracks.at(j).ySlope();
in.SetXYZ(th_inx,th_iny,1.0);
th_in_rec_tracks=in.Theta();
}
}

int yes_er=0;
int yes_mur=0;

for(int j=0; j<tracks.size();j++)
{

if(tracks.at(j).sector()==1 and tracks.at(j).processIDofLinkedTrack()==0) cout << "c'e' un muone 0!! " << endl;
if(tracks.at(j).processIDofLinkedTrack()==45 and tracks.at(j).sector()==1) TrackIdreco=tracks.at(j).linkedTrackID();
if(tracks.at(j).processIDofLinkedTrack()==45 and tracks.size()>=3 and tracks.at(j).sector()==1 and tracks.at(j).percentageOfHitsSharedWithLinkedTrack()>=0) 
{ yes2++; cout << "tracks.size " << tracks.size() << endl;

 int sum = tracks.at(j).numberOfXProjectionHits() + tracks.at(j).numberOfYProjectionHits() + tracks.at(j).numberOfStereoHits();
 cout << "TrackerStubs " << TrackerStubs->GetEntries() << endl;
 cout << "TrackID " << tracks.at(j).linkedTrackID() << endl;
 cout << "#%hitsshared " <<tracks.at(j).percentageOfHitsSharedWithLinkedTrack() << " and sum of hits " << sum << endl;
 cout << "chi2perDegreeOfFreedom " << tracks.at(j).chi2perDegreeOfFreedom() << endl;

//if(tracks.size()==3)
//{
        if( Ee<5 ){
        if(code_e==tracks.at(j).linkedTrackID()) {        h_echi_eTR->Fill(Ee,tracks.at(j).chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);
							  h_chi_eTR->Fill(tracks.at(j).chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);}
        if(code_mu==tracks.at(j).linkedTrackID()) {       h_echi_muTR->Fill(Ee,tracks.at(j).chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);
                                                          h_chi_muTR->Fill(tracks.at(j).chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);}
        }
	if( Ee<15 and Ee>5){
        if(code_e==tracks.at(j).linkedTrackID()) {        h_echi_e1TR->Fill(Ee,tracks.at(j).chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);
                                                          h_chi_e1TR->Fill(tracks.at(j).chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);}
        if(code_mu==tracks.at(j).linkedTrackID()) {       h_echi_mu1TR->Fill(Ee,tracks.at(j).chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);
                                                          h_chi_mu1TR->Fill(tracks.at(j).chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);}
        }
	if( Ee<25 and Ee>15){
        if(code_e==tracks.at(j).linkedTrackID()) {        h_echi_e2TR->Fill(Ee,tracks.at(j).chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);
                                                          h_chi_e2TR->Fill(tracks.at(j).chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);}
        if(code_mu==tracks.at(j).linkedTrackID()) {       h_echi_mu2TR->Fill(Ee,tracks.at(j).chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);
                                                          h_chi_mu2TR->Fill(tracks.at(j).chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);}
        }
	if( Ee>25){
        if(code_e==tracks.at(j).linkedTrackID()) {        h_echi_e3TR->Fill(Ee,tracks.at(j).chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);
                                                          h_chi_e3TR->Fill(tracks.at(j).chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);}
        if(code_mu==tracks.at(j).linkedTrackID()) {       h_echi_mu3TR->Fill(Ee,tracks.at(j).chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);
                                                          h_chi_mu3TR->Fill(tracks.at(j).chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);}
        }

if(code_e==tracks.at(j).linkedTrackID()) {
	yes_er++;
	h_chi_eTOT->Fill(tracks.at(j).chi2(),MesmerEvent->wgt_full);
	TVector3 e_outv(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0);
        the_rec_tracks=e_outv.Theta();
					}

if(code_mu==tracks.at(j).linkedTrackID()) {
        yes_mur++;
	h_chi_muTOT->Fill(tracks.at(j).chi2(),MesmerEvent->wgt_full);
	TVector3 mu_outv(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0);
	thmu_rec_tracks=mu_outv.Theta();
					}
//}


	}
}

cout << "Numero elettroni " << yes_er << ", numero muoni " << yes_mur << endl;

if(yes_er==1 and yes_mur==1 and tracks.size()==3) {reco3+=MesmerEvent->wgt_full; cout << "reco3" << endl;}
//and tracks.size()==3){reco3+=MesmerEvent->wgt_full; cout << "reco3" << endl;}

if(yes2>=2){reco+=MesmerEvent->wgt_full;
 cout << "reco" << endl;
if(yes_v>=1){
double res = the_rec_vrtx-the_gen;
double res_mu = thmu_rec_vrtx-thmu_gen;
double res_muin = th_in_rec_vrtx-th_in_gen;

double resTR = the_rec_tracks-the_gen;
double res_muTR = thmu_rec_tracks-thmu_gen;
double res_muinTR = th_in_rec_tracks-th_in_gen;

if( Ee<5 )h_res->Fill(Ee,res,MesmerEvent->wgt_full);
if( Ee<15 and Ee>5)h_res1->Fill(Ee,res,MesmerEvent->wgt_full);
if( Ee<25 and Ee>15)h_res2->Fill(Ee,res,MesmerEvent->wgt_full);
if( Ee>25)h_res3->Fill(Ee,res,MesmerEvent->wgt_full);
h_res_mu->Fill(Emu,res_mu,MesmerEvent->wgt_full);
h_res_mu_in->Fill(Emu_in,res_muin,MesmerEvent->wgt_full);


if( Ee<5 )h_resTR->Fill(Ee,resTR,MesmerEvent->wgt_full);
if( Ee<15 and Ee>5)h_res1TR->Fill(Ee,resTR,MesmerEvent->wgt_full);
if( Ee<25 and Ee>15)h_res2TR->Fill(Ee,resTR,MesmerEvent->wgt_full);
if( Ee>25)h_res3TR->Fill(Ee,resTR,MesmerEvent->wgt_full);

h_res_muTR->Fill(Emu,res_muTR,MesmerEvent->wgt_full);
h_res_mu_inTR->Fill(Emu_in,res_muinTR,MesmerEvent->wgt_full);


 }
}
else{reco0+=MesmerEvent->wgt_full;
     angle_e->Fill(the_sdr,MesmerEvent->wgt_full);
     angle_mu->Fill(thmu_sdr,MesmerEvent->wgt_full);
     energy_e->Fill(Ee,MesmerEvent->wgt_full);
     h_2d_emu->Fill(the_sdr,thmu_sdr,MesmerEvent->wgt_full);
     h_2d_e->Fill(the_sdr,Ee,MesmerEvent->wgt_full);}

if(yes2>2){more_reco+=MesmerEvent->wgt_full;}

//if(yes2==0){reco0+=MesmerEvent->wgt_full;
//cout <<"NOT RECONSTRUCTED"<<endl;
//	}

if(yes2<2 and TrackIdreco!=-99){reco1+=MesmerEvent->wgt_full;
cout <<"NOT RECONSTRUCTED"<<endl;
	}

if(yes_v>=1){reco_v+=MesmerEvent->wgt_full;
        }

if(yes_v>1){more_reco_v+=MesmerEvent->wgt_full;
        }

if(yes_v<1){reco0_v+=MesmerEvent->wgt_full;
cout <<"NOT RECONSTRUCTED vertex"<<endl;
        }



}
cout << "---------------------"<<endl;
yes_e=0;yes_mu=0;
code_e=-99;code_mu=-99;
point_el=0;point_mu=0;
TrackIdreco=-99;
yes2=0;yes_v=0;
} //end of general for

double ratio =reco/signal;
double ratio1 =reco1/signal;
double ratio0 =reco0/signal;
double ratioM =more_reco/signal;

cout << "Su " << signal << " eventi di segnale, " << reco << " sono ricostruiti, con un rapporto del " << ratio*100 << "%"<< endl;
cout << "Su " << signal << " eventi di segnale, " << reco3 << " sono ricostruiti con 3 tracce, con un rapporto del " << (reco3/signal)*100 << "%"<< endl;
cout << "Su " << signal << " eventi di segnale con piu tracce (stesso id), " << more_reco << " sono ricostruiti, con un rapporto del " << ratioM*100 << "%"<< endl;
cout << "Su " << signal << " eventi di segnale con 0 tracce di segnale reco, " << reco0 << ", con un rapporto del " << ratio0*100 << "%"<< endl;
cout << "Su " << signal << " eventi di con 1 sola traccia di segnale reco, " << reco1 << ", con un rapporto del " << ratio1*100 << "%"<< endl;

double ratio_v =reco_v/signal;
double ratio0_v =reco0_v/signal;
double ratioM_v =more_reco_v/signal;

cout << "Su " << signal << " eventi di segnale, " << reco_v << " sono ricostruiti con almeno 1 vertice, con un rapporto del " << ratio_v*100 << "%"<< endl;
cout << "Su " << signal << " eventi di segnale, " << more_reco_v << " sono ricostruiti con piu' di unun vertice, con un rapporto del " << ratioM_v*100 << "%"<< endl;
cout << "Su " << signal << " eventi di segnale, " << reco0_v << " hanno 0 vertici, con un rapporto del " << ratio0_v*100 << "%"<< endl;

cout << "Su " << reco_v << " eventi di segnale ricostruiti, " << e << " hanno vrtxID vs linkID sbagliate per elettrone: " << (e*100)/reco_v << "%" << endl;
cout << "Su " << reco_v << " eventi di segnale ricostruiti, " << mu << " hanno vrtxID vs linkID sbagliate per muone: " << (mu*100)/reco_v << "%" << endl; 
cout << "Su " << reco_v << " eventi di segnale ricostruiti, " << one << " hanno entrambi i link sbagliati: " << (one*100)/reco_v << "%" << endl;

Int_t nxTh = h_2d_e->GetNbinsX();
Int_t nyTh = h_2d_e->GetNbinsY();
for (Int_t i=1; i<nxTh+1; i++) {
for (Int_t j=1; j<nyTh+1; j++) {
if (h_2d_e->GetBinContent(i,j)<1) h_2d_e->SetBinContent(i,j,0);}}

Int_t nxTh1 = h_2d_emu->GetNbinsX();
Int_t nyTh1 = h_2d_emu->GetNbinsY();
for (Int_t i=1; i<nxTh1+1; i++) {
for (Int_t j=1; j<nyTh1+1; j++) {
if (h_2d_emu->GetBinContent(i,j)<1) h_2d_emu->SetBinContent(i,j,0);}}

TCanvas s("s","s",700,700);
s.Divide(2,3);
s.cd(1);
angle_e->Draw("hist");
s.cd(2);
angle_mu->Draw("hist");
s.cd(3);
energy_e->Draw("hist");
gPad->SetLogy();
s.cd(4);
h_2d_emu->Draw("COLZ");
s.cd(5);
h_2d_e->Draw("COLZ");
s.SaveAs("not_reco.pdf");

/*
Int_t nxTh = h_echi_e->GetNbinsX();
Int_t nyTh = h_echi_e->GetNbinsY();
for (Int_t i=1; i<nxTh+1; i++) {
for (Int_t j=1; j<nyTh+1; j++) {
if (h_echi_e->GetBinContent(i,j)<1) h_echi_e->SetBinContent(i,j,0);}}


Int_t nxTh0 = h_echi_mu->GetNbinsX();
Int_t nyTh0 = h_echi_mu->GetNbinsY();
for (Int_t i=1; i<nxTh0+1; i++) {
for (Int_t j=1; j<nyTh0+1; j++) {
if (h_echi_mu->GetBinContent(i,j)<1) h_echi_mu->SetBinContent(i,j,0);}}

Int_t nxTh1 = h_echi_e1->GetNbinsX();
Int_t nyTh1 = h_echi_e1->GetNbinsY();
for (Int_t i=1; i<nxTh1+1; i++) {
for (Int_t j=1; j<nyTh1+1; j++) {
if (h_echi_e1->GetBinContent(i,j)<1) h_echi_e1->SetBinContent(i,j,0);}}


Int_t nxTh2 = h_echi_mu1->GetNbinsX();
Int_t nyTh2 = h_echi_mu1->GetNbinsY();
for (Int_t i=1; i<nxTh2+1; i++) {
for (Int_t j=1; j<nyTh2+1; j++) {
if (h_echi_mu->GetBinContent(i,j)<1) h_echi_mu->SetBinContent(i,j,0);}}

Int_t nxTh3 = h_echi_e2->GetNbinsX();
Int_t nyTh3 = h_echi_e2->GetNbinsY();
for (Int_t i=1; i<nxTh3+1; i++) {
for (Int_t j=1; j<nyTh3+1; j++) {
if (h_echi_e2->GetBinContent(i,j)<1) h_echi_e2->SetBinContent(i,j,0);}}


Int_t nxTh4 = h_echi_mu2->GetNbinsX();
Int_t nyTh4 = h_echi_mu2->GetNbinsY();
for (Int_t i=1; i<nxTh4+1; i++) {
for (Int_t j=1; j<nyTh4+1; j++) {
if (h_echi_mu2->GetBinContent(i,j)<1) h_echi_mu2->SetBinContent(i,j,0);}}

Int_t nxTh5 = h_echi_e3->GetNbinsX();
Int_t nyTh5 = h_echi_e3->GetNbinsY();
for (Int_t i=1; i<nxTh5+1; i++) {
for (Int_t j=1; j<nyTh5+1; j++) {
if (h_echi_e3->GetBinContent(i,j)<1) h_echi_e3->SetBinContent(i,j,0);}}

Int_t nxTh6 = h_echi_mu3->GetNbinsX();
Int_t nyTh6 = h_echi_mu3->GetNbinsY();
for (Int_t i=1; i<nxTh6+1; i++) {
for (Int_t j=1; j<nyTh6+1; j++) {
if (h_echi_mu3->GetBinContent(i,j)<1) h_echi_mu3->SetBinContent(i,j,0);}}



Int_t nxThTR = h_echi_eTR->GetNbinsX();
Int_t nyThTR = h_echi_eTR->GetNbinsY();
for (Int_t i=1; i<nxThTR+1; i++) {
for (Int_t j=1; j<nyThTR+1; j++) {
if (h_echi_eTR->GetBinContent(i,j)<1) h_echi_eTR->SetBinContent(i,j,0);}}

Int_t nxTh0TR = h_echi_muTR->GetNbinsX();
Int_t nyTh0TR = h_echi_muTR->GetNbinsY();
for (Int_t i=1; i<nxTh0TR+1; i++) {
for (Int_t j=1; j<nyTh0TR+1; j++) {
if (h_echi_muTR->GetBinContent(i,j)<1) h_echi_muTR->SetBinContent(i,j,0);}}

Int_t nxTh1TR = h_echi_e1TR->GetNbinsX();
Int_t nyTh1TR = h_echi_e1TR->GetNbinsY();
for (Int_t i=1; i<nxTh1TR+1; i++) {
for (Int_t j=1; j<nyTh1TR+1; j++) {
if (h_echi_e1TR->GetBinContent(i,j)<1) h_echi_e1TR->SetBinContent(i,j,0);}}


Int_t nxTh2TR = h_echi_mu1TR->GetNbinsX();
Int_t nyTh2TR = h_echi_mu1TR->GetNbinsY();
for (Int_t i=1; i<nxTh2TR+1; i++) {
for (Int_t j=1; j<nyTh2TR+1; j++) {
if (h_echi_muTR->GetBinContent(i,j)<1) h_echi_muTR->SetBinContent(i,j,0);}}

Int_t nxTh3TR = h_echi_e2TR->GetNbinsX();
Int_t nyTh3TR = h_echi_e2TR->GetNbinsY();
for (Int_t i=1; i<nxTh3TR+1; i++) {
for (Int_t j=1; j<nyTh3TR+1; j++) {
if (h_echi_e2TR->GetBinContent(i,j)<1) h_echi_e2TR->SetBinContent(i,j,0);}}


Int_t nxTh4TR = h_echi_mu2TR->GetNbinsX();
Int_t nyTh4TR = h_echi_mu2TR->GetNbinsY();
for (Int_t i=1; i<nxTh4TR+1; i++) {
for (Int_t j=1; j<nyTh4TR+1; j++) {
if (h_echi_mu2TR->GetBinContent(i,j)<1) h_echi_mu2TR->SetBinContent(i,j,0);}}

Int_t nxTh5TR = h_echi_e3TR->GetNbinsX();
Int_t nyTh5TR = h_echi_e3TR->GetNbinsY();
for (Int_t i=1; i<nxTh5TR+1; i++) {
for (Int_t j=1; j<nyTh5TR+1; j++) {
if (h_echi_e3TR->GetBinContent(i,j)<1) h_echi_e3TR->SetBinContent(i,j,0);}}

Int_t nxTh6TR = h_echi_mu3TR->GetNbinsX();
Int_t nyTh6TR = h_echi_mu3TR->GetNbinsY();
for (Int_t i=1; i<nxTh6TR+1; i++) {
for (Int_t j=1; j<nyTh6TR+1; j++) {
if (h_echi_mu3TR->GetBinContent(i,j)<1) h_echi_mu3TR->SetBinContent(i,j,0);}}




Int_t nxTh7 = h_res->GetNbinsX();
Int_t nyTh7 = h_res->GetNbinsY();
for (Int_t i=1; i<nxTh7+1; i++) {
for (Int_t j=1; j<nyTh7+1; j++) {
if (h_res->GetBinContent(i,j)<1) h_res->SetBinContent(i,j,0);}}

Int_t nxTh8 = h_res1->GetNbinsX();
Int_t nyTh8 = h_res1->GetNbinsY();
for (Int_t i=1; i<nxTh8+1; i++) {
for (Int_t j=1; j<nyTh8+1; j++) {
if (h_res1->GetBinContent(i,j)<1) h_res1->SetBinContent(i,j,0);}}

Int_t nxTh9 = h_res2->GetNbinsX();
Int_t nyTh9 = h_res2->GetNbinsY();
for (Int_t i=1; i<nxTh9+1; i++) {
for (Int_t j=1; j<nyTh9+1; j++) {
if (h_res2->GetBinContent(i,j)<1) h_res2->SetBinContent(i,j,0);}}

Int_t nxTh10 = h_res3->GetNbinsX();
Int_t nyTh10 = h_res3->GetNbinsY();
for (Int_t i=1; i<nxTh10+1; i++) {
for (Int_t j=1; j<nyTh10+1; j++) {
if (h_res3->GetBinContent(i,j)<1) h_res3->SetBinContent(i,j,0);}}

Int_t nxTh11 = h_res_mu->GetNbinsX();
Int_t nyTh11 = h_res_mu->GetNbinsY();
for (Int_t i=1; i<nxTh11+1; i++) {
for (Int_t j=1; j<nyTh11+1; j++) {
if (h_res_mu->GetBinContent(i,j)<1) h_res_mu->SetBinContent(i,j,0);}}

Int_t nxTh12 = h_res_mu_in->GetNbinsX();
Int_t nyTh12 = h_res_mu_in->GetNbinsY();
for (Int_t i=1; i<nxTh12+1; i++) {
for (Int_t j=1; j<nyTh12+1; j++) {
if (h_res_mu_in->GetBinContent(i,j)<1) h_res_mu_in->SetBinContent(i,j,0);}}



Int_t nxTh7TR = h_resTR->GetNbinsX();
Int_t nyTh7TR = h_resTR->GetNbinsY();
for (Int_t i=1; i<nxTh7TR+1; i++) {
for (Int_t j=1; j<nyTh7TR+1; j++) {
if (h_resTR->GetBinContent(i,j)<1) h_resTR->SetBinContent(i,j,0);}}

Int_t nxTh8TR = h_res1TR->GetNbinsX();
Int_t nyTh8TR = h_res1TR->GetNbinsY();
for (Int_t i=1; i<nxTh8TR+1; i++) {
for (Int_t j=1; j<nyTh8TR+1; j++) {
if (h_res1TR->GetBinContent(i,j)<1) h_res1TR->SetBinContent(i,j,0);}}

Int_t nxTh9TR = h_res2TR->GetNbinsX();
Int_t nyTh9TR = h_res2TR->GetNbinsY();
for (Int_t i=1; i<nxTh9TR+1; i++) {
for (Int_t j=1; j<nyTh9TR+1; j++) {
if (h_res2TR->GetBinContent(i,j)<1) h_res2TR->SetBinContent(i,j,0);}}

Int_t nxTh10TR = h_res3TR->GetNbinsX();
Int_t nyTh10TR = h_res3TR->GetNbinsY();
for (Int_t i=1; i<nxTh10TR+1; i++) {
for (Int_t j=1; j<nyTh10TR+1; j++) {
if (h_res3TR->GetBinContent(i,j)<1) h_res3TR->SetBinContent(i,j,0);}}

Int_t nxTh11TR = h_res_muTR->GetNbinsX();
Int_t nyTh11TR = h_res_muTR->GetNbinsY();
for (Int_t i=1; i<nxTh11TR+1; i++) {
for (Int_t j=1; j<nyTh11TR+1; j++) {
if (h_res_muTR->GetBinContent(i,j)<1) h_res_muTR->SetBinContent(i,j,0);}}

Int_t nxTh12TR = h_res_mu_inTR->GetNbinsX();
Int_t nyTh12TR = h_res_mu_inTR->GetNbinsY();
for (Int_t i=1; i<nxTh12TR+1; i++) {
for (Int_t j=1; j<nyTh12TR+1; j++) {
if (h_res_mu_inTR->GetBinContent(i,j)<1) h_res_mu_inTR->SetBinContent(i,j,0);}}

TCanvas a("pa","p",700,700);
a.Divide(2,4);
a.cd(1);
auto p1 = h_echi_e->ProfileX();
p1->SetLineColor(40);

auto p1TR = h_echi_eTR->ProfileX();
p1->Draw();
p1TR->Draw("same");

a.cd(2);
auto p2 = h_echi_mu->ProfileX();
p2->SetLineColor(40);

auto p2TR = h_echi_muTR->ProfileX();
p2->Draw();
p2TR->Draw("same");

a.cd(3);
auto p3 = h_echi_e1->ProfileX();
p3->SetLineColor(46);

auto p3TR = h_echi_e1TR->ProfileX();
p3->Draw();
p3TR->Draw("same");

a.cd(4);
auto p4 = h_echi_mu1->ProfileX();
p4->SetLineColor(46);

auto p4TR = h_echi_mu1TR->ProfileX();
p4->Draw();
p4TR->Draw("same");

a.cd(5);
auto p5 = h_echi_e2->ProfileX();
p5->SetLineColor(42);

auto p5TR = h_echi_e2TR->ProfileX();
p5->Draw();
p5TR->Draw("same");

a.cd(6);
auto p6 = h_echi_mu2->ProfileX();
p6->SetLineColor(42);

auto p6TR = h_echi_mu2TR->ProfileX();
p6->Draw();
p6TR->Draw("same");

a.cd(7);
auto p7 = h_echi_e3->ProfileX();
p7->SetLineColor(30);

auto p7TR = h_echi_e3TR->ProfileX();
p7->Draw();
p7TR->Draw("same");

a.cd(8);
auto p8 = h_echi_mu3->ProfileX();
p8->SetLineColor(30);

auto p8TR = h_echi_mu3TR->ProfileX();
p8->Draw();
p8TR->Draw("same");

a.SaveAs("chiVSen.pdf");


TCanvas b("b","b",700,700);
b.Divide(3,4);
b.cd(1);
h_res->Draw("COLZ");
b.cd(2);
TObjArray aSlices;
h_res->SetLineColor(40);
auto q2 = h_res->ProfileX();
q2->SetLineColor(40);
q2->Draw();
h_res->FitSlicesY(0,0,-1,0,"QLW", &aSlices);
b.cd(3);
//aSlices[2]->SetLineColor(40);
aSlices[2]->Draw();

b.cd(4);
h_res1->Draw("COLZ");
b.cd(5);
TObjArray aSlices1;
auto q3 = h_res1->ProfileX();
h_res1->FitSlicesY(0,0,-1,0,"QLW", &aSlices1);
q3->SetLineColor(46);
q3->Draw();
b.cd(6);
//aSlices1[2]->SetLineColor(46);
aSlices1[2]->Draw();

b.cd(7);
h_res2->Draw("COLZ");
b.cd(8);
TObjArray aSlices2;
auto q4 = h_res2->ProfileX();
 h_res2->FitSlicesY(0,0,-1,0,"QLW", &aSlices2);
q4->SetLineColor(42);
q4->Draw();
b.cd(9);
//aSlices2[2]->SetLineColor(42);
aSlices2[2]->Draw();

b.cd(10);
h_res3->Draw("COLZ");
b.cd(11);
TObjArray aSlices3;
auto q1 = h_res3->ProfileX();
h_res3->FitSlicesY(0,0,-1,0,"QLW", &aSlices3);
q1->SetLineColor(30);
q1->Draw();
b.cd(12);
//aSlices3[2]->SetLineColor(30);
aSlices3[2]->Draw();

b.SaveAs("res_postVrtx.pdf");


TCanvas b1("b1","b1",700,700);
b1.Divide(3,4);
b1.cd(1);
h_resTR->Draw("COLZ");
b1.cd(2);
TObjArray aSlicesTR;
h_resTR->SetLineColor(40);
auto q2TR = h_resTR->ProfileX();
q2TR->SetLineColor(40);
q2TR->Draw();
h_resTR->FitSlicesY(0,0,-1,0,"QLW", &aSlicesTR);
b1.cd(3);
//aSlices[2]->SetLineColor(40);
aSlicesTR[2]->Draw();

b1.cd(4);
h_res1TR->Draw("COLZ");
b1.cd(5);
TObjArray aSlices1TR;
auto q3TR = h_res1TR->ProfileX();
h_res1TR->FitSlicesY(0,0,-1,0,"QLW", &aSlices1TR);
q3TR->SetLineColor(46);
q3TR->Draw();
b1.cd(6);
//aSlices1[2]->SetLineColor(46);
aSlices1TR[2]->Draw();

b1.cd(7);
h_res2TR->Draw("COLZ");
b1.cd(8);
TObjArray aSlices2TR;
auto q4TR = h_res2TR->ProfileX();
 h_res2TR->FitSlicesY(0,0,-1,0,"QLW", &aSlices2TR);
q4TR->SetLineColor(42);
q4TR->Draw();
b1.cd(9);
//aSlices2[2]->SetLineColor(42);
aSlices2TR[2]->Draw();

b1.cd(10);
h_res3TR->Draw("COLZ");
b1.cd(11);
TObjArray aSlices3TR;
auto q1TR = h_res3TR->ProfileX();
h_res3TR->FitSlicesY(0,0,-1,0,"QLW", &aSlices3TR);
q1TR->SetLineColor(30);
q1TR->Draw();
b1.cd(12);
//aSlices3[2]->SetLineColor(30);
aSlices3TR[2]->Draw();

b1.SaveAs("res_preVrtx.pdf");

TCanvas a1("pa1","p",700,700);
a1.Divide(2,5);
a1.cd(1);
h_chi_e->SetLineColor(40);
h_chi_eTR->Draw("hist");
h_chi_e->Draw("hist same");
a1.cd(2);
h_chi_mu->SetLineColor(40);
h_chi_muTR->Draw("hist");
h_chi_mu->Draw("hist same");
a1.cd(3);
h_chi_e1->SetLineColor(46);
h_chi_e1TR->Draw("hist");
h_chi_e1->Draw("hist same");
a1.cd(4);
h_chi_mu1->SetLineColor(46);
h_chi_mu1TR->Draw("hist");
h_chi_mu1->Draw("hist same");
a1.cd(5);
h_chi_e2->SetLineColor(42);
h_chi_e2TR->Draw("hist");
h_chi_e2->Draw("hist same");
a1.cd(6);
h_chi_mu2->SetLineColor(42);
h_chi_mu2TR->Draw("hist");
h_chi_mu2->Draw("hist same");
a1.cd(7);
h_chi_e3->SetLineColor(30);
h_chi_e3TR->Draw("hist");
h_chi_e3->Draw("hist same ");
a1.cd(8);
h_chi_mu3->SetLineColor(30);
h_chi_mu3TR->Draw("hist");
h_chi_mu3->Draw("hist same");
a1.cd(9);
h_chi_eTOT->Draw("hist");
a1.cd(10);
h_chi_muTOT->Draw("hist");
a1.SaveAs("chi.pdf");


TCanvas d("d","d",700,700);
d.Divide(2,3);
d.cd(1);
h_res_mu->Draw("COLZ");
d.cd(3);
TObjArray aSlices_mu;
auto m1 = h_res_mu->ProfileX();
h_res_mu->FitSlicesY(0,0,-1,0,"QLW", &aSlices_mu);
m1->SetLineColor(30);
m1->Draw();
d.cd(5);
//aSlices3[2]->SetLineColor(30);
aSlices_mu[2]->Draw();

d.cd(2);
h_res_mu_in->Draw("COLZ");
d.cd(4);
TObjArray aSlices_mu_in;
auto m1_in = h_res_mu_in->ProfileX();
h_res_mu_in->FitSlicesY(0,0,-1,0,"QLW", &aSlices_mu_in);
m1_in->SetLineColor(35);
m1_in->Draw();
d.cd(6);
//aSlices3[2]->SetLineColor(30);
aSlices_mu_in[2]->Draw();

d.SaveAs("res_mu_POSTvrtx.pdf");


TCanvas d1("d1","d1",700,700);
d1.Divide(2,3);
d1.cd(1);
h_res_muTR->Draw("COLZ");
d1.cd(3);
TObjArray aSlices_muTR;
auto m1TR = h_res_muTR->ProfileX();
h_res_muTR->FitSlicesY(0,0,-1,0,"QLW", &aSlices_muTR);
m1TR->SetLineColor(30);
m1TR->Draw();
d1.cd(5);
//aSlices3[2]->SetLineColor(30);
aSlices_muTR[2]->Draw();
d1.cd(2);
h_res_mu_inTR->Draw("COLZ");
d1.cd(4);
TObjArray aSlices_mu_inTR;
auto m1_inTR = h_res_mu_inTR->ProfileX();
h_res_mu_inTR->FitSlicesY(0,0,-1,0,"QLW", &aSlices_mu_inTR);
m1_inTR->SetLineColor(35);
m1_inTR->Draw();
d1.cd(6);
//aSlices3[2]->SetLineColor(30);
aSlices_mu_inTR[2]->Draw();

d1.SaveAs("res_mu_PREvrtx.pdf");
*/

}


