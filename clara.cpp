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

 	TFile *inputfile = new TFile("TRMesmer_100k_box.root");//beamprofile/TRMesmer_1M.root");
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

TH1D *h_E1 = new TH1D("ene1","generated electron energy when event is reco with the_el_rec<10" ,175,0,45);
TH1D *h_E2 = new TH1D("ene2","generated electron energy when event is reco with 10<=the_el_rec<30" ,75,0,15);
TH1D *h_E3 = new TH1D("ene3","generated electron energy when event is reco with the_el_rec>=30" ,75,0,15);

TH1D *h_E1_LO = new TH1D("ene1lo","generated electron energy when event is reco with the_el_rec<10 LEADING ORDER" ,75,0,15);
TH1D *h_E2_LO = new TH1D("ene2lo","generated electron energy when event is reco with 10<=the_el_rec<30 LEADING ORDER" ,75,0,15);
TH1D *h_E3_LO = new TH1D("ene3lo","generated electron energy when event is reco with the_el_rec>=30 LEADING ORDER" ,75,0,15);

TH1D *h_E1_ph = new TH1D("ene1_ph","generated electron energy when event is reco with the_el_rec<10 and there is a photon in MCtracks" ,175,0,45);
TH1D *h_E2_ph = new TH1D("ene2_ph","generated electron energy when event is reco with 10<=the_el_rec<30 and there is a photon in MCtracks" ,75,0,15);
TH1D *h_E3_ph = new TH1D("ene3_ph","generated electron energy when event is reco with the_el_rec>=30 and there is a photon in MCtracks" ,75,0,15);


TH1D *h_E1_gen = new TH1D("ene1_gen","generated electron energy when event is reco with the_el_gen<10" ,175,0,45);
TH1D *h_E2_gen = new TH1D("ene2_gen","generated electron energy when event is reco with 10<=the_el_gen<30" ,75,0,15);
TH1D *h_E3_gen = new TH1D("ene3_gen","generated electron energy when event is reco with the_el_gen>=30" ,75,0,15);

TH1D *h_E1_LO_gen = new TH1D("ene1lo_gen","generated electron energy when event is reco with the_el_gen<10 LEADING ORDER" ,75,0,15);
TH1D *h_E2_LO_gen = new TH1D("ene2lo_gen","generated electron energy when event is reco with 10<=the_el_gen<30 LEADING ORDER" ,75,0,15);
TH1D *h_E3_LO_gen = new TH1D("ene3lo_gen","generated electron energy when event is reco with the_el_gen>=30 LEADING ORDER" ,75,0,15);

TH1D *h_E1_ph_gen = new TH1D("ene1_ph_gen","generated electron energy when event is reco with the_el_gen<10 and there is a photon in MCtracks" ,175,0,45);
TH1D *h_E2_ph_gen = new TH1D("ene2_ph_gen","generated electron energy when event is reco with 10<=the_el_gen<30 and there is a photon in MCtracks" ,75,0,15);
TH1D *h_E3_ph_gen = new TH1D("ene3_ph_gen","generated electron energy when event is reco with the_el_gen>=30 and there is a photon in MCtracks" ,75,0,15);





TH1D *pid1 = new TH1D("pid1","ID of the MCTracks in the case where 10<=the_el_rec<30 mrad",50,-25,25);
TH1D *pid2 = new TH1D("pid2","ID of the MCTracks in the case where 10<=the_el_rec<30 mrad",50,-25,25);
TH1D *pid3 = new TH1D("pid3","ID of the MCTracks in the case where 10<=the_el_rec<30 mrad",50,-25,25);

TH1D *pr1 = new TH1D("pr1","process ID of the gamma in the case where 10<=the_el_rec<30 mrad",50,0,50);
TH1D *pr2 = new TH1D("pr2","process ID of the gamma in the case where 10<=the_el_rec<30 mrad",50,0,50);
TH1D *pr3 = new TH1D("pr3","process ID of the gamma in the case where 10<=the_el_rec<30 mrad",50,0,50);

TH1D *mot1 = new TH1D("mot1","mother ID of the MCTracks in the case where 10<=the_el_rec<30 mrad",50,-25,25);
TH1D *mot2 = new TH1D("mot2","mother ID of the MCTracks in the case where 10<=the_el_rec<30 mrad",50,-25,25);
TH1D *mot3 = new TH1D("mot3","motherID of the MCTracks in the case where 10<=the_el_rec<30 mrad",50,-25,25);


double the_gen=0; double thmu_gen=0; double th_in_gen=0;
double the_rec_vrtx=0; double thmu_rec_vrtx=0; double th_in_rec_vrtx=0;
double the_rec_tracks=0; double thmu_rec_tracks=0; double th_in_rec_tracks=0;

double signal=0.; double reco=0.; double reco1=0.; double more_reco=0.; double reco0=0.;
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

int no=0;int noM=0;
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
         //if(MCTr->pdgCode()==22) {no=1;}
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
		thmu_gen=pmu.Theta();
		double mx1=SigTracks->ax();
		double my1=SigTracks->ay();
                pmu_dir=pmu.Unit();
		TVector3 mu(mx1,my1,1.0);
		thmu_sdr=muin.Angle(mu);
		//thmu_sdr=acos(pmuin_dir.Dot(pmu_dir));
		yes_mu=1;}
	 }
	}

 if(yes_e!=1 or yes_mu!=1) cout << "NOT RECONSTRUCTIBLE" << endl;

 if(yes_e==1 and yes_mu==1){// and no==0){
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

			}
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

for(int j=0; j<tracks.size();j++)
{

if(tracks.at(j).processIDofLinkedTrack()==45 and tracks.at(j).sector()==1) TrackIdreco=tracks.at(j).linkedTrackID();
if(tracks.at(j).processIDofLinkedTrack()==45 and tracks.size()>=3 and tracks.at(j).sector()==1)
{ yes2++; cout << "tracks.size " << tracks.size() << endl;


if(tracks.size()==3)
{

if(code_e==tracks.at(j).linkedTrackID()) {

        TVector3 e_outv(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0);
        the_rec_tracks=in.Angle(e_outv);


// PRIMO BLOCCO DI ELETTRONI
if(the_rec_tracks<0.010) {
	for(int n = 0; n < MCTrack->GetEntries(); n++) {
         const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
cout << n <<  ") IntID : " <<MCTr->interactionID() << " and pdf " << MCTr->pdgCode() << " mother " << MCTr->motherID() << endl;
pid1->Fill(MCTr->pdgCode(),MesmerEvent->wgt_full);
if(MCTr->pdgCode()==22){pr1->Fill(MCTr->interactionID(),MesmerEvent->wgt_full);

        if(MCTr->motherID()!=-1){no=1;
         const MUonETrack *mother = static_cast<const MUonETrack*>(MCTrack->At(MCTr->motherID()));
                                mot1->Fill(mother->pdgCode(),MesmerEvent->wgt_full);}
        if(MCTr->motherID()==-1){noM=1;mot1->Fill(-1,MesmerEvent->wgt_full);}
                        }

						}

if(no==1)h_E1_ph->Fill(Ee,MesmerEvent->wgt_full);
if(noM==1)h_E1_LO->Fill(Ee,MesmerEvent->wgt_full);
if(no==0 and noM==0)h_E1->Fill(Ee,MesmerEvent->wgt_full);
no=0;noM=0;
}

if(the_sdr<0.010) {
        for(int n = 0; n < MCTrack->GetEntries(); n++) {
         const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
if(MCTr->pdgCode()==22){if(MCTr->interactionID()!=45)no=1;if(MCTr->interactionID()==45)noM=1;}
                                                }

if(no==1)h_E1_ph_gen->Fill(Ee,MesmerEvent->wgt_full);
if(noM==1)h_E1_LO_gen->Fill(Ee,MesmerEvent->wgt_full);
if(no==0 and noM==0)h_E1_gen->Fill(Ee,MesmerEvent->wgt_full);
no=0;noM=0;
}

// SECONDO BLOCCO DI ELETTRONI

if(the_rec_tracks>=0.010 and the_rec_tracks<0.030) {
        for(int n = 0; n < MCTrack->GetEntries(); n++) {
         const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
cout << n <<  ") IntID : " <<MCTr->interactionID() << " and pdf " << MCTr->pdgCode() << " mother " << MCTr->motherID() << endl;
pid2->Fill(MCTr->pdgCode(),MesmerEvent->wgt_full);
if(MCTr->pdgCode()==22){pr2->Fill(MCTr->interactionID(),MesmerEvent->wgt_full);
	if(MCTr->motherID()!=-1){no=1;
         const MUonETrack *mother = static_cast<const MUonETrack*>(MCTrack->At(MCTr->motherID()));
				mot2->Fill(mother->pdgCode(),MesmerEvent->wgt_full);}
	if(MCTr->motherID()==-1){noM=1;mot2->Fill(-1,MesmerEvent->wgt_full);}
			}
                                                }
if(no==1)h_E2_ph->Fill(Ee,MesmerEvent->wgt_full);
if(noM==1)h_E2_LO->Fill(Ee,MesmerEvent->wgt_full);
if(no==0 and noM==0)h_E2->Fill(Ee,MesmerEvent->wgt_full);
no=0;noM=0;
}

if(the_sdr>=0.010 and the_sdr<0.030) {
        for(int n = 0; n < MCTrack->GetEntries(); n++) {
         const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
if(MCTr->pdgCode()==22){if(MCTr->motherID()!=-1)no=1;if(MCTr->motherID()==-1)noM=1;}
                                                }
if(no==1)h_E2_ph_gen->Fill(Ee,MesmerEvent->wgt_full);
if(noM==1)h_E2_LO_gen->Fill(Ee,MesmerEvent->wgt_full);
if(no==0 and noM==0)h_E2_gen->Fill(Ee,MesmerEvent->wgt_full);
no=0;noM=0;
}

// TERZO BLOCCO DI ELETTRONI

if(the_rec_tracks>=0.030) {h_E3->Fill(Ee,MesmerEvent->wgt_full);
        for(int n = 0; n < MCTrack->GetEntries(); n++) {
         const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
cout << n <<  ") IntID : " <<MCTr->interactionID() << " and pdf " << MCTr->pdgCode() << " mother " << MCTr->motherID() << endl;
pid3->Fill(MCTr->pdgCode(),MesmerEvent->wgt_full);
if(MCTr->pdgCode()==22){pr3->Fill(MCTr->interactionID(),MesmerEvent->wgt_full);	 
        if(MCTr->motherID()!=-1){no=1;
         const MUonETrack *mother = static_cast<const MUonETrack*>(MCTrack->At(MCTr->motherID()));
                                mot3->Fill(mother->pdgCode(),MesmerEvent->wgt_full);}
        if(MCTr->motherID()==-1){noM=1;mot3->Fill(-1,MesmerEvent->wgt_full);}
                        }
                                                }
if(no==1)h_E3_ph->Fill(Ee,MesmerEvent->wgt_full);
if(noM==1)h_E3_LO->Fill(Ee,MesmerEvent->wgt_full);
if(no==0 and noM==0)h_E3->Fill(Ee,MesmerEvent->wgt_full);
no=0;noM=0;
					}
if(the_sdr>=0.030) {h_E3->Fill(Ee,MesmerEvent->wgt_full);
        for(int n = 0; n < MCTrack->GetEntries(); n++) {
         const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
if(MCTr->pdgCode()==22){if(MCTr->motherID()!=-1)no=1;if(MCTr->motherID()==-1)noM=1;}
                                                }
if(no==1)h_E3_ph_gen->Fill(Ee,MesmerEvent->wgt_full);
if(noM==1)h_E3_LO_gen->Fill(Ee,MesmerEvent->wgt_full);
if(no==0 and noM==0)h_E3_gen->Fill(Ee,MesmerEvent->wgt_full);
no=0;noM=0;
                                        }

	}

if(code_mu==tracks.at(j).linkedTrackID()) {

	TVector3 mu_outv(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0);
	thmu_rec_tracks=mu_outv.Theta();
					}

		}
	}
}


if(yes2==2 ){reco+=MesmerEvent->wgt_full;

}


if(yes2>2){more_reco+=MesmerEvent->wgt_full;}

if(yes2<2 and TrackIdreco==-99){reco0+=MesmerEvent->wgt_full;
cout <<"NOT RECONSTRUCTED"<<endl;
	}

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


TCanvas a("a","a",700,700);
a.Divide(1,3);
a.cd(1);
h_E1_LO->SetLineColor(kOrange);
h_E1_LO->Draw("hist");
h_E1->Draw("hist same");
h_E1_ph->SetLineColor(kPink);
h_E1_ph->Draw("hist same");
//gPad->SetLogy();
a.cd(2);
h_E2->Draw("hist");
h_E2_LO->SetLineColor(kOrange);
h_E2_LO->Draw("hist same");
h_E2_ph->SetLineColor(kPink);
h_E2_ph->Draw("hist same");
gPad->SetLogy();
a.cd(3);
h_E3->Draw("hist");
h_E3_LO->SetLineColor(kOrange);
h_E3_LO->Draw("hist same");
h_E3_ph->SetLineColor(kPink);
h_E3_ph->Draw("hist same");
gPad->SetLogy();
a.SaveAs("electron_e.pdf");

TCanvas a1("a1","a1",700,700);
a1.Divide(1,3);
a1.cd(1);
h_E1_LO_gen->SetLineColor(kOrange);
h_E1_LO_gen->Draw("hist");
h_E1_gen->Draw("hist same");
h_E1_ph_gen->SetLineColor(kPink);
h_E1_ph_gen->Draw("hist same");
//gPad->SetLogy();
a1.cd(2);
h_E2_gen->Draw("hist");
h_E2_LO_gen->SetLineColor(kOrange);
h_E2_LO_gen->Draw("hist same");
h_E2_ph_gen->SetLineColor(kPink);
h_E2_ph_gen->Draw("hist same");
gPad->SetLogy();
a1.cd(3);
h_E3->Draw("hist");
h_E3_LO_gen->SetLineColor(kOrange);
h_E3_LO_gen->Draw("hist same");
h_E3_ph_gen->SetLineColor(kPink);
h_E3_ph_gen->Draw("hist same");
gPad->SetLogy();
a1.SaveAs("gen_electron_e.pdf");

TCanvas b("b","b",700,700);
b.Divide(3,3);
b.cd(1);
pid1->Draw("hist");
b.cd(2);
pr1->Draw("hist");
b.cd(3);
mot1->Draw("hist");
b.cd(4);
pid2->Draw("hist");
b.cd(5);
pr2->Draw("hist");
b.cd(6);
mot2->Draw("hist");
b.cd(7);
pid3->Draw("hist");
b.cd(8);
pr3->Draw("hist");
b.cd(9);
mot3->Draw("hist");
b.SaveAs("pid.pdf");

}


