
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

	TFile *inputfile = new TFile("TRMesmer_box_100k_2GeV_div.root");
	// 	TFile *inputfile = new TFile("TRMesmer_ohit.root");
	//TRMesmer_thin_5k.root");//TRMesmer_100k_box.root");//beamprofile/TRMesmer_1M.root");
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



/*TH2D *h_resTR=new TH2D("resTR", "(the_rec-the_true) VS energy Ee<5 GeV PRE-VRTX",25,0,5,400,-0.04,0.04);
TH2D *h_res1TR=new TH2D("res1TR", "(the_rec-the_true) VS energy 5<E<15 GeV PRE-VRTX",20,5,15,80,-0.004,0.004);
TH2D *h_res2TR=new TH2D("res2TR", "(the_rec-the_true) VS energy 15<E<25 GeV PRE-VRTX",10,15,20,20,-0.001,0.001);
TH2D *h_res3TR=new TH2D("res3TR", "(the_rec-the_true) VS energy Ee>25 GeV PRE-VRTX",6,25,65,20,-0.001,0.001);*/

TH1D *h_resTR=new TH1D("resTR", "(the_rec-the_true) Ee<5 GeV PRE-VRTX",100,-0.01,0.01);
TH1D *h_res1TR=new TH1D("res1TR", "(the_rec-the_true) 5<E<15 GeV PRE-VRTX",50,-0.0025,0.0025);
TH1D *h_res2TR=new TH1D("res2TR", "(the_rec-the_true) 15<E<25 GeV PRE-VRTX",40,-0.001,0.001);
TH1D *h_res3TR=new TH1D("res3TR", "(the_rec-the_true) Ee>25 GeV PRE-VRTX",40,-0.001,0.001);

TH2D *h_res_muTR=new TH2D("h_res_muTR", "(thmu_rec-thmu_true) VS energy Emu(GeV) PRE-VRTX",30,155,160,120,-0.0006,0.0006);//30,140,170,120,-0.0006,0.0006);
TH2D *h_res_mu_inTR=new TH2D("h_res_mu_inTR", "(thmu_in_rec-thmu_in_true) VS energy Emu(GeV) PRE-VRTX",10,150,160,22,-0.0006,0.0006);

TH1D *h_res_muTR1x=new TH1D("h_res_muTR1x", "(thmu_rec-thmu_true) proj. X  155<Emu<157(GeV) PRE-VRTX",120,-0.0006,0.0006);
TH1D *h_res_muTR2x=new TH1D("h_res_muTR2x", "(thmu_rec-thmu_true) proj. X 157<Emu<158(GeV) PRE-VRTX",120,-0.0006,0.0006);
TH1D *h_res_muTR3x=new TH1D("h_res_muTR3x", "(thmu_rec-thmu_true) proj. X 158<Emu<160(GeV) PRE-VRTX",120,-0.0006,0.0006);

TH1D *h_res_muTR1y=new TH1D("h_res_muTR1y", "(thmu_rec-thmu_true) proj. Y 155<Emu<157(GeV) PRE-VRTX",120,-0.0006,0.0006);
TH1D *h_res_muTR2y=new TH1D("h_res_muTR2y", "(thmu_rec-thmu_true) proj. Y 157<Emu<158(GeV) PRE-VRTX",120,-0.0006,0.0006);
TH1D *h_res_muTR3y=new TH1D("h_res_muTR3y", "(thmu_rec-thmu_true) proj. Y 158<Emu<160(GeV) PRE-VRTX",120,-0.0006,0.0006);

TH1D *h_res_muTR1=new TH1D("h_res_muTR1", "(thmu_rec-thmu_true) VS energy 155<Emu<157(GeV) PRE-VRTX",180,-0.0006,0.0006);
TH1D *h_res_muTR2=new TH1D("h_res_muTR2", "(thmu_rec-thmu_true) VS energy 157<Emu<158(GeV) PRE-VRTX",180,-0.0006,0.0006);
TH1D *h_res_muTR3=new TH1D("h_res_muTR3", "(thmu_rec-thmu_true) VS energy 158<Emu<160(GeV) PRE-VRTX",180,-0.0006,0.0006);
TH1D *h_res_mu_inTR1=new TH1D("h_res_mu_inTR1", "(thmu_in_rec-thmu_in_true)  PRE-VRTX",80,-0.0002,0.0002);

TH1D *theta1 =new TH1D("theta1" , "theta muon tail residuum" , 250,0.,0.005);
TH1D *theta2 =new TH1D("theta2" , "theta muon peak residuum" , 250,0.,0.005);
TH1D *theta1g =new TH1D("theta1g" , "theta muon tail residuum" , 250,0.,0.005);
TH1D *theta2g =new TH1D("theta2g" , "theta muon peak residuum" , 250,0.,0.005);

double the_gen=0; double thmu_gen=0; double th_in_gen=0;
double the_rec_vrtx=0; double thmu_rec_vrtx=0; double th_in_rec_vrtx=0;
double the_rec_tracks=0; double thmu_rec_tracks=0; double th_in_rec_tracks=0;

double thmu_gen_x=0; double thmu_gen_y=0;
double thmu_rec_x=0; double thmu_rec_y=0;

double signal=0.; double reco=0.; double reco1=0.; double more_reco=0.; double reco0=0.;
double reco_v=0.; double more_reco_v=0.; double reco0_v=0.; double one=0;
double e=0.; double mu=0.;

int yes_e=0;int yes_mu=0; int yes2=0; int yes_v=0;
int point_mu=0; int point_el=0;
int code_mu=-99; int code_e=-99;
int TrackIdreco=-99;


TH1D *th_mu44g = new TH1D("th_mu44g"," angle mu 44 wrt incoming muon (generated)",200,0,0.005);
TH1D *th_mu44 = new TH1D("th_mu44"," angle mu 44 wrt incoming muon (reconstructed)",200,0,0.005);
TH1D *th_e44g = new TH1D("th_e44e"," angle e 44 wrt incoming muon (generated)",200,0,0.05);
TH1D *th_e44 = new TH1D("th_mu44_r"," angle e 44 wrt incoming muon (reconstructed)",200,0,0.05);




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
	 pmuin=pmuin.Unit();
//	TVector3 z (0.,0.,1.);
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
		pe=pe.Unit();
		TVector3 z(0.,0.,1.);
                //the_gen=pe.Theta();
		//the_gen=acos(pe.Dot(z));
		yes_e=1;}
 	   if(SigTracks->pdgCode()==-13 and last_modXmu==2 and last_modYmu==2 and stereo_mu>1){
		Emu=SigTracks->energy();
		pmu.SetXYZ(SigTracks->px(),SigTracks->py(),SigTracks->pz());
	        TVector3 z(0.,0.,1.);
		pmu=pmu.Unit();
		thmu_gen=pmu.Theta();
         double mx=SigTracks->ax();
         double my=SigTracks->ay();
                thmu_gen_x=mx;
		thmu_gen_y=my;
		yes_mu=1;}
	 }
	}

	// if(yes_e!=1 or yes_mu!=1) cout << "NOT RECONSTRUCTIBLE" << endl;

	if(yes_e==1 and yes_mu==1){
	  //cout << "RECONSTRUCTIBLE" << endl;


	   signal+=MesmerEvent->wgt_full;

vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();

TVector3 in;
for(int j=0; j<tracks.size();j++)
{
if(tracks.at(j).processIDofLinkedTrack()==0 and tracks.at(j).sector()==0){
double th_inx=tracks.at(j).xSlope();
double th_iny=tracks.at(j).ySlope();
in.SetXYZ(th_inx,th_iny,1.0);
//TVector3 z(0.,0.,1.0);
th_in_rec_tracks=in.Theta();
	}
}

int yes_mu=0;
int yes_e=0;

for(int j=0; j<tracks.size();j++)
{
 if(tracks.at(j).processIDofLinkedTrack()==45 and tracks.size()==3 and tracks.at(j).sector()==1 and tracks.at(j).percentageOfHitsSharedWithLinkedTrack()>=100)
	{yes2++;

	 	 if(code_e==tracks.at(j).linkedTrackID()) { yes_e=1;
	 	 TVector3 e_outv(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0);
		TVector3 z(0.,0.,1.);
		e_outv=e_outv.Unit();
        	 the_rec_tracks=e_outv.Theta();
					}

		 if(code_mu==tracks.at(j).linkedTrackID()) { yes_mu=1;
		 TVector3 mu_outv(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0);
		TVector3 z(0.,0.,1.);
                mu_outv=mu_outv.Unit();
                //thmu_rec_tracks=acos(mu_outv.Dot(z));
		thmu_rec_x=tracks.at(j).xSlope();
                thmu_rec_y=tracks.at(j).ySlope();

		 thmu_rec_tracks=mu_outv.Theta();
					}
	}
}


if(yes_mu==1 and yes_e==1 ){reco+=MesmerEvent->wgt_full;

 th_mu44->Fill(thmu_rec_tracks);
 th_e44->Fill(the_rec_tracks);
 th_mu44g->Fill(thmu_gen);
 th_e44g->Fill(the_gen);

double resTR = the_rec_tracks-the_gen;
double res_muTR = thmu_rec_tracks-thmu_gen;
double res_muinTR = th_in_rec_tracks-th_in_gen;
double resx=thmu_rec_x-thmu_gen_x;
double resy=thmu_rec_y-thmu_gen_y;

cout << "th_in_rec_tracks " << th_in_rec_tracks << endl;
cout << "th_in_gen " << th_in_gen << endl;
cout << "res_muinTR " << res_muinTR << endl;

if( Ee<5) h_resTR->Fill(resTR,MesmerEvent->wgt_full);
if( Ee<15 and Ee>5) h_res1TR->Fill(resTR,MesmerEvent->wgt_full);
if( Ee<25 and Ee>15) h_res2TR->Fill(resTR,MesmerEvent->wgt_full);
if( Ee>25) h_res3TR->Fill(resTR,MesmerEvent->wgt_full);

if(Emu>150 and Emu<=157) {h_res_muTR1->Fill(res_muTR,MesmerEvent->wgt_full);h_res_muTR1x->Fill(resx,MesmerEvent->wgt_full); h_res_muTR1y->Fill(resy,MesmerEvent->wgt_full);}
if(Emu>157 and Emu<=160) {h_res_muTR2->Fill(res_muTR,MesmerEvent->wgt_full);h_res_muTR2x->Fill(resx,MesmerEvent->wgt_full); h_res_muTR2y->Fill(resy,MesmerEvent->wgt_full);}
//if(Emu>158 and Emu<=160) {h_res_muTR3->Fill(res_muTR,MesmerEvent->wgt_full);h_res_muTR3x->Fill(resx,MesmerEvent->wgt_full); h_res_muTR3y->Fill(resy,MesmerEvent->wgt_full);}
h_res_muTR->Fill(Emu,res_muTR,MesmerEvent->wgt_full);
h_res_mu_inTR->Fill(Emu_in,res_muinTR,MesmerEvent->wgt_full);
h_res_mu_inTR1->Fill(res_muinTR,MesmerEvent->wgt_full);
if(res_muTR>0.1e-03){theta1->Fill(thmu_rec_tracks,MesmerEvent->wgt_full);theta1g->Fill(thmu_gen,MesmerEvent->wgt_full);}
if(res_muTR<0.1e-03 and res_muTR>-0.1e-03){theta2->Fill(thmu_rec_tracks,MesmerEvent->wgt_full);theta2g->Fill(thmu_gen,MesmerEvent->wgt_full);}


}

}
//cout << "---------------------"<<endl;
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


/*
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


TCanvas d1("d1","d1",700,700);
d1.Divide(2,3);
d1.cd(1);
h_res_muTR->Draw("COLZ");
d1.cd(2);
TObjArray aSlices_muTR;
auto m1TR = h_res_muTR->ProfileX();
h_res_muTR->FitSlicesY(0,0,-1,0,"QLW", &aSlices_muTR);
m1TR->SetLineColor(30);
m1TR->Draw();
d1.cd(3);
aSlices_muTR[0]->Draw();
d1.cd(4);
aSlices_muTR[1]->Draw();
d1.cd(5);
aSlices_muTR[2]->Draw();
d1.cd(6);
aSlices_muTR[3]->Draw();
gPad->SetLogy();
d1.SaveAs("res_mu_PREvrtx.pdf");

*/


TF1 *a1 = new TF1("a1", "[2]*TMath::Gaus(x,[0],[1])");
a1->SetParameters(5.3e-05,110e-05,1);

TF1 *a2 = new TF1("a2", "[2]*TMath::Gaus(x,[0],[1])");
a2->SetParameters(2.6e-05,50e-05,1);


TCanvas d1("d1","d1",700,700);
d1.Divide(2,2);
d1.cd(1);
h_resTR->Fit("a1","R","",-0.006,+0.006);
h_resTR->Draw("hist same");
gStyle->SetOptFit(1111);
d1.cd(2);
h_res1TR->Fit("a2","R","",-0.002,+0.002);
h_res1TR->Draw("hist same");
gStyle->SetOptFit(1111);
d1.cd(3);
h_res2TR->Fit("a2","R","",-0.001,+0.001);
h_res2TR->Draw("hist same");
gStyle->SetOptFit(1111);
d1.cd(4);
h_res3TR->Fit("a2","R","",-0.001,+0.001);
h_res3TR->Draw("hist same");
gStyle->SetOptFit(1111);
d1.SaveAs("electron_res.pdf");

/*
TF1 *f1 = new TF1("f1", "[2]*TMath::Gaus(x,[0],[1])");
f1->SetParameters(7.7e-06,30e-06,1);

 TF1 *f2 = new TF1("f2", "[2]*TMath::Gaus(x,[0],[1])");
f2->SetParameters(3.8e-05,25e-06,1);

TF1 *f3 = new TF1("f3", "[2]*TMath::Gaus(x,[0],[1])");
f3->SetParameters(2e-05,25e-06,1);

TCanvas d2("d2","d2",700,700);
d2.Divide(2,3);
d2.cd(1);
h_res_muTR1->Draw();
h_res_muTR1->Draw("hist same");
h_res_muTR1->Fit("f1","R","",-0.0001,+0.0001);
//h_res_muTR1->Fit("gaus","R","",-0.0001,+0.0001);
 gStyle->SetOptFit(1111);
d2.cd(2);
h_res_muTR2->Draw();
h_res_muTR2->Draw("hist same");
h_res_muTR2->Fit("f2","R","",-0.0001,+0.0001);
//h_res_muTR2->Fit("gaus","R","",-0.0001,+0.0001);
gStyle->SetOptFit(1111);
d2.cd(3);
theta1g->SetLineColor(kOrange);
theta1g->Draw("hist");
theta1->Draw("hist same");
d2.cd(4);
theta2g->SetLineColor(kOrange);
theta2g->Draw("hist");
theta2->Draw("hist same");
d2.cd(5);
h_res_mu_inTR1->Draw();
h_res_mu_inTR1->Fit("f1","R","",-0.0001,+0.0001);
h_res_mu_inTR1->Draw("hist same");
gStyle->SetOptFit(1111);
d2.SaveAs("mu_out_peren_GA.pdf");


TF1 *f4 = new TF1("f4", "[2]*TMath::Gaus(x,[0],[1])");
f4->SetParameters(-1.9e-05,30e-06,1);

TF1 *f5 = new TF1("f5", "[2]*TMath::Gaus(x,[0],[1])");
f5->SetParameters(-1.2e-05,28e-06,1);

TF1 *f6 = new TF1("f6", "[2]*TMath::Gaus(x,[0],[1])");
f6->SetParameters(-7e-05,28e-06,1);

TCanvas d3("d3","d3",700,700);
d3.Divide(2,2);
d3.cd(1);
h_res_muTR1x->Draw();
h_res_muTR1x->Draw("hist same");
h_res_muTR1x->Fit("f4","R","",-0.0001,+0.0001);
gStyle->SetOptFit(1111);
d3.cd(3);
h_res_muTR2x->Draw();
h_res_muTR2x->Draw("hist same");
h_res_muTR2x->Fit("f5","R","",-0.0001,+0.0001);
gStyle->SetOptFit(1111);
d3.cd(2);
h_res_muTR1y->Draw();
h_res_muTR1y->Draw("hist same");
h_res_muTR1y->Fit("f5","R","",-0.0001,+0.0001);
gStyle->SetOptFit(1111);
d3.cd(4);
h_res_muTR2y->Draw();
h_res_muTR2y->Draw("hist same");
h_res_muTR2y->Fit("f6","R","",-0.0001,+0.0001);
gStyle->SetOptFit(1111);
d3.SaveAs("proj_mu_GA.pdf");

TCanvas b1("b1","b1",700,700);
b1.Divide(1,2);
b1.cd(1);
th_mu44g->Draw("hist");
th_mu44->SetLineColor(kViolet);
th_mu44->Draw("hist same");
gPad->SetLogy();
b1.cd(2);
th_e44g->Draw("hist");
th_e44->SetLineColor(kViolet);
th_e44->Draw("hist same");
b1.SaveAs("angles.pdf");
*/

}


