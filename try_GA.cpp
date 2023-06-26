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

	TFile *inputfile = new TFile("TRMesmer_box_offset_100k_2GeV_0.root");
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
TH1D *h_res1TR=new TH1D("res1TR", "(the_rec-the_true) 5<E<15 GeV PRE-VRTX",150,-0.0025,0.0025);
TH1D *h_res2TR=new TH1D("res2TR", "(the_rec-the_true) 15<E<25 GeV PRE-VRTX",40,-0.001,0.001);
TH1D *h_res3TR=new TH1D("res3TR", "(the_rec-the_true) Ee>25 GeV PRE-VRTX",40,-0.001,0.001);

TH2D *h_res_muTR=new TH2D("h_res_muTR", "(thmu_rec-thmu_true) VS energy Emu(GeV) PRE-VRTX",30,155,160,120,-0.0006,0.0006);//30,140,170,120,-0.0006,0.0006);
TH2D *h_res_mu_inTR=new TH2D("h_res_mu_inTR", "(thmu_in_rec-thmu_in_true) VS energy Emu(GeV) PRE-VRTX",10,150,160,22,-0.0006,0.0006);

TH1D *h_res_muTR1x=new TH1D("h_res_muTR1x", "(thmu_rec-thmu_true) proj. X  155<Emu<157(GeV) PRE-VRTX",120,-0.0006,0.0006);
TH1D *h_res_muTR2x=new TH1D("h_res_muTR2x", "(thmu_rec-thmu_true) proj. X 157<Emu<160(GeV) PRE-VRTX",120,-0.0006,0.0006);
TH1D *h_res_muTR3x=new TH1D("h_res_muTR3x", "(thmu_rec-thmu_true) proj. X 158<Emu<160(GeV) PRE-VRTX",120,-0.0006,0.0006);

TH1D *h_res_muTR1y=new TH1D("h_res_muTR1y", "(thmu_rec-thmu_true) proj. Y 155<Emu<157(GeV) PRE-VRTX",120,-0.0006,0.0006);
TH1D *h_res_muTR2y=new TH1D("h_res_muTR2y", "(thmu_rec-thmu_true) proj. Y 157<Emu<160(GeV) PRE-VRTX",120,-0.0006,0.0006);
TH1D *h_res_muTR3y=new TH1D("h_res_muTR3y", "(thmu_rec-thmu_true) proj. Y 158<Emu<160(GeV) PRE-VRTX",120,-0.0006,0.0006);

TH1D *h_res_muTR1=new TH1D("h_res_muTR1", "(thmu_rec-thmu_true) VS energy 155<Emu<157(GeV) PRE-VRTX",180,-0.0006,0.0006);
TH1D *h_res_muTR2=new TH1D("h_res_muTR2", "(thmu_rec-thmu_true) VS energy 157<Emu<160(GeV) PRE-VRTX",180,-0.0006,0.0006);
TH1D *h_res_muTR3=new TH1D("h_res_muTR3", "(thmu_rec-thmu_true) VS energy 158<Emu<160(GeV) PRE-VRTX",180,-0.0006,0.0006);
TH1D *h_res_mu_inTR1=new TH1D("h_res_mu_inTR1", "(thmu_in_rec-thmu_in_true)  PRE-VRTX",80,-0.0002,0.0002);

TH1D *theta1 =new TH1D("theta1" , "theta of the muon for events in the peak of angular residuum" , 150,0.,0.005);
TH1D *theta2 =new TH1D("theta2" , "theta of the muon for events in the tail of angular residuum" , 150,0.,0.005);
TH1D *theta1g =new TH1D("theta1g" , "theta of the muon for events in the peak of angular residuum" , 150,0.,0.005);
TH1D *theta2g =new TH1D("theta2g" , "theta of the muon for events in the tail of angular residuum" , 150,0.,0.005);

TH1D *theta1x =new TH1D("theta1x" , "theta x of the muon for events in the peak of angular residuum" , 150,-0.0025,0.0025);//1150,-0.0025,0.0025
TH1D *theta2x =new TH1D("theta2x" , "theta x of the muon for events in the tail of angular residuum" , 150,-0.0025,0.0025);
TH1D *theta1gx =new TH1D("theta1gx" , "theta x of the muon for events in the peak of angular residuum" , 150,-0.0025,0.0025);
TH1D *theta2gx =new TH1D("theta2gx" , "theta x of the muon for events in the tail of angular residuum" , 150,-0.0025,0.0025);

TH1D *theta1y =new TH1D("theta1y" , "theta y of the muon for events in the peak of angular residuum" , 150,-0.0025,0.0025);
TH1D *theta2y =new TH1D("theta2y" , "theta y of the muon for events in the tail of angular residuum" , 150,-0.0025,0.0025);
TH1D *theta1gy =new TH1D("theta1gy" , "theta y of the muon for events in the peak of angular residuum" , 150,-0.0025,0.0025);
TH1D *theta2gy =new TH1D("theta2gy" , "theta y of the muon for events in the tail of angular residuum" , 150,-0.0025,0.0025);


TH1D *h_phiin =new TH1D("h_phiin"," entering muon phi angle",180,-180,180);

TH1D *h_phi =new TH1D("h_phi"," outgoing muon phi angle for events in the tails of the angular residuum",180,-180,180);
TH1D *h_phiP =new TH1D("h_phip"," outgoing muon phi angle for events in the peaks of the angular residuum",180,-180,180);

TH1D *h_phie =new TH1D("h_phie"," outgoing electron phi angle for events in the tails of the angular residuum",180,-180,180);
TH1D *h_phiPe =new TH1D("h_phipe"," outgoing electron phi angle for events in the peaks of the angular residuum",180,-180,180);

TH1D *h_phitot =new TH1D("h_phitot"," outgoing muon phi angle",180,-180,180);
TH1D *h_phietot =new TH1D("h_phietot"," outgoing muon phi angle",180,-180,180);

TH1D* nstubs=new TH1D("nstubs","nstubs for well reconstructed events",25,5,30);
TH1D* nstubs90=new TH1D("nstubs90","nstubs when 88<phi_mu<93",25,5,30);
TH1D* nstubs180=new TH1D("nstubs180","nnstubs when 178<phi_mu<183",25,5,30);

TH2D *h_thphi=new TH2D("resTR", "theta-phi for events in the tail of angular residuum",180,-180,180,1150,-0.0025,0.0025);

TH1D *h_op0 = new TH1D("op0","opening angle for events in the tail of angular residuum",400,0,0.04);


double the_gen=0; double thmu_gen=0; double th_in_gen=0;
double the_rec_vrtx=0; double thmu_rec_vrtx=0; double th_in_rec_vrtx=0;
double the_rec_tracks=0; double thmu_rec_tracks=0; double th_in_rec_tracks=0;

double thmu_gen_x=0; double thmu_gen_y=0;
double thmu_rec_x=0; double thmu_rec_y=0;

double the_gen_x=0; double the_gen_y=0;

double signal=0.; double reco=0.; double reco1=0.; double more_reco=0.; double reco0=0.;
double reco_v=0.; double more_reco_v=0.; double reco0_v=0.; double one=0;
double e=0.; double mu=0.;
double tailx=0; double peakx=0;
double taily=0; double peaky=0;

int yes_e_g=0;int yes_mu_g=0; int yes2=0; int yes_v=0;
int point_mu=0; int point_el=0;
int code_mu=-99; int code_e=-99;
int TrackIdreco=-99;

double phi=0;
double phie=0;
double phi_in=0.;


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
         phi_in=pmuin.Phi()*(180/M_PI);
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
                                                                                                 if(TrackerPt->moduleID()==5)last_modYe++;
												 if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_e++;
												}
                          if(TrackerPt->trackPDGCode()==-13 and SigTracks->pdgCode()==-13 and TrackerPt->trackID()==code_mu and TrackerPt->stationID()==1){ point_mu++;
												 if(TrackerPt->moduleID()==4) last_modXmu++;
                                                                                                 if(TrackerPt->moduleID()==5) last_modYmu++;
												 if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_mu++;
												}
                         }

           if(SigTracks->pdgCode()==11 and last_modXe==2 and last_modYe==2 and stereo_e>1){
		pe.SetXYZ(SigTracks->px(),SigTracks->py(),SigTracks->pz()); the_gen=pe.Theta();
		Ee=SigTracks->energy();
		pe=pe.Unit();
		TVector3 z(0.,0.,1.);
                //the_gen=pe.Theta();
		//the_gen=acos(pe.Dot(z));
         double mx=SigTracks->ax();
         double my=SigTracks->ay();
                the_gen_x=mx;
                the_gen_y=my;
                phie=pe.Phi()*(180/M_PI);
		yes_e_g=1;}
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
		yes_mu_g=1;
                phi=pmu.Phi()*(180/M_PI);;
			}
	 }
	}

	// if(yes_e!=1 or yes_mu!=1) cout << "NOT RECONSTRUCTIBLE" << endl;

	if(yes_e_g==1 and yes_mu_g==1){
	  //cout << "RECONSTRUCTIBLE" << endl;
	   signal+=MesmerEvent->wgt_full;

double opening = acos(pe.Dot(pmu));

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

if( (phi<91 and phi>89) or (phi>-91 and phi<-89) ) nstubs90->Fill(TrackerStubs->GetEntries()-6,MesmerEvent->wgt_full);
if( (phi<181 and phi>179) or (phi>-181 and phi<-179) ) nstubs180->Fill(TrackerStubs->GetEntries()-6,MesmerEvent->wgt_full);

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
                 //thmu_rec_tracks=sqrt(thmu_rec_x*thmu_rec_x+thmu_rec_y*thmu_rec_y);
		 thmu_rec_tracks=mu_outv.Theta();
					}
	}
}


if(yes_mu==1 and yes_e==1 ){reco+=MesmerEvent->wgt_full;

h_phitot->Fill(phi,MesmerEvent->wgt_full); h_phietot->Fill(phie,MesmerEvent->wgt_full);

h_phiin->Fill(phi_in,MesmerEvent->wgt_full);

nstubs->Fill(TrackerStubs->GetEntries()-6,MesmerEvent->wgt_full);

double resTR = the_rec_tracks-the_gen;
double res_muTR = thmu_rec_tracks-thmu_gen;
double res_muinTR = th_in_rec_tracks-th_in_gen;
double resx=thmu_rec_x-thmu_gen_x;
double resy=thmu_rec_y-thmu_gen_y;


if( Ee<5) h_resTR->Fill(resTR,MesmerEvent->wgt_full);
if( Ee<15 and Ee>5) h_res1TR->Fill(resTR,MesmerEvent->wgt_full);
if( Ee<25 and Ee>15) h_res2TR->Fill(resTR,MesmerEvent->wgt_full);
if( Ee>25) h_res3TR->Fill(resTR,MesmerEvent->wgt_full);



//if(Emu>155 and Emu<=157) {
h_res_muTR1->Fill(res_muTR,MesmerEvent->wgt_full);h_res_muTR1x->Fill(resx,MesmerEvent->wgt_full); h_res_muTR1y->Fill(resy,MesmerEvent->wgt_full);
//if(Emu>157 and Emu<=160) {h_res_muTR2->Fill(res_muTR,MesmerEvent->wgt_full);h_res_muTR2x->Fill(resx,MesmerEvent->wgt_full); h_res_muTR2y->Fill(resy,MesmerEvent->wgt_full);}

h_res_muTR->Fill(Emu,res_muTR,MesmerEvent->wgt_full);
h_res_mu_inTR->Fill(Emu_in,res_muinTR,MesmerEvent->wgt_full);
h_res_mu_inTR1->Fill(res_muinTR,MesmerEvent->wgt_full);

if(abs(res_muTR)<0.1e-03){theta1->Fill(thmu_rec_tracks,MesmerEvent->wgt_full);theta1g->Fill(thmu_gen,MesmerEvent->wgt_full);}
if(abs(res_muTR)>0.1e-03){theta2->Fill(thmu_rec_tracks,MesmerEvent->wgt_full);theta2g->Fill(thmu_gen,MesmerEvent->wgt_full);}

if(abs(resx)<0.1e-03) peakx+=MesmerEvent->wgt_full;
if(abs(resy)<0.1e-03) peaky+=MesmerEvent->wgt_full;

if(abs(resx)<0.1e-03 and abs(resy)<0.1e-03) {theta1x->Fill(thmu_rec_x,MesmerEvent->wgt_full);theta1gx->Fill(thmu_gen_x,MesmerEvent->wgt_full);
						theta1y->Fill(thmu_rec_y,MesmerEvent->wgt_full);theta1gy->Fill(thmu_gen_y,MesmerEvent->wgt_full);}

if(abs(resx)>0.1e-03){tailx+=MesmerEvent->wgt_full; 

h_thphi->Fill(phi,thmu_gen_x,MesmerEvent->wgt_full);h_op0->Fill(opening,MesmerEvent->wgt_full);

                                   theta2x->Fill(thmu_rec_x,MesmerEvent->wgt_full);theta2gx->Fill(thmu_gen_x,MesmerEvent->wgt_full);}
if(abs(resy)>0.1e-03){taily+=MesmerEvent->wgt_full;
                                   theta2y->Fill(thmu_rec_y,MesmerEvent->wgt_full);theta2gy->Fill(thmu_gen_y,MesmerEvent->wgt_full);}

     if(abs(resx)>0.1e-03 or abs(resy)>0.1e-03) {h_phi->Fill(phi,MesmerEvent->wgt_full);h_phie->Fill(phie,MesmerEvent->wgt_full);

}
     else{h_phiP->Fill(phi,MesmerEvent->wgt_full); h_phiPe->Fill(phie,MesmerEvent->wgt_full);}

}

}
//cout << "---------------------"<<endl;
yes_e_g=0;yes_mu_g=0;
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
cout << "Su " << signal << " eventi di segnale, eventi di con proiezione x nella coda " << (tailx/reco)*100 << "%"<< endl;
cout << "Su " << signal << " eventi di segnale, eventi di con proiezione y nella coda " << (taily/reco)*100 << "%"<< endl;

cout << "Su " << signal << " eventi di segnale, eventi di con proiezione x nella picco " << (peakx/reco)*100 << "%"<< endl;
cout << "Su " << signal << " eventi di segnale, eventi di con proiezione y nella picco " << (peaky/reco)*100 << "%"<< endl;

Int_t nxThTR = h_thphi->GetNbinsX();
Int_t nyThTR = h_thphi->GetNbinsY();
for (Int_t i=1; i<nxThTR+1; i++) {
for (Int_t j=1; j<nyThTR+1; j++) {
if (h_thphi->GetBinContent(i,j)<1) h_thphi->SetBinContent(i,j,0);}}



/*TCanvas b1("b1","b1",700,700);
h_phiin->Draw("hist");
b1.SaveAs("pdf_try/phi_in.pdf");*/



TCanvas b1("b1","b1",700,700);
b1.Divide(2,3);

b1.cd(1);
h_phitot->Draw("hist");
gStyle->SetOptStat(222001111);
b1.cd(2);
h_phietot->SetLineColor(kRed);
h_phietot->Draw("hist");
gStyle->SetOptStat(222001111);
b1.cd(3);
h_phiP->Draw("hist");
gStyle->SetOptStat(222001111);
b1.cd(4);
h_phi->Draw("hist");
gStyle->SetOptStat(222001111);
b1.cd(5);
h_phiPe->SetLineColor(kRed);
h_phiPe->Draw("hist");
gStyle->SetOptStat(222001111);
b1.cd(6);
h_phie->SetLineColor(kRed);
gStyle->SetOptStat(222001111);
h_phie->Draw("hist");
b1.SaveAs("pdf_try/phi.pdf");



TCanvas b2("b1","b1",700,700);
b2.Divide(2,2);
b2.cd(1);
nstubs->Draw("hist");
nstubs90->SetLineColor(kRed);
nstubs90->Draw("hist same");
//gPad->SetLogy();
b2.cd(2);
nstubs->Draw("hist");
nstubs180->SetLineColor(kRed);
nstubs180->Draw("hist same");
//gPad->SetLogy();
b2.cd(3);
h_thphi->Draw("COLZ");
gStyle->SetOptStat(222001111);
b2.cd(4);
h_op0->Draw("hist");
gStyle->SetOptStat(222001111);
b2.SaveAs("pdf_try/2Dth_phi.pdf");


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
gStyle->SetOptStat(222001111);
//gPad->SetLogy();
d1.cd(2);
h_res1TR->Fit("a2","R","",-0.002,+0.002);
h_res1TR->Draw("hist same");
gStyle->SetOptFit(1111);
gStyle->SetOptStat(222001111);
//gPad->SetLogy();
d1.cd(3);
h_res2TR->Fit("a2","R","",-0.001,+0.001);
h_res2TR->Draw("hist same");
gStyle->SetOptFit(1111);
gStyle->SetOptStat(222001111);
//gPad->SetLogy();
d1.cd(4);
h_res3TR->Fit("a2","R","",-0.0001,+0.0001);
h_res3TR->Draw("hist same");
gStyle->SetOptFit(1111);
gStyle->SetOptStat(222001111);
//gPad->SetLogy();
d1.SaveAs("pdf_try/electron_res.pdf");

TF1 *f1 = new TF1("f1", "[2]*TMath::Gaus(x,[0],[1])");
f1->SetParameters(7.7e-06,30e-06,1);

 TF1 *f2 = new TF1("f2", "[2]*TMath::Gaus(x,[0],[1])");
f2->SetParameters(3.8e-05,25e-06,1);

TF1 *f3 = new TF1("f3", "[2]*TMath::Gaus(x,[0],[1])");
f3->SetParameters(2e-05,25e-06,1);

TCanvas d2("d2","d2",700,700);
d2.Divide(2,2);
d2.cd(1);
h_res_muTR1->Draw();
h_res_muTR1->Draw("hist same");
h_res_muTR1->Fit("f1","R","",-0.0001,+0.0001);
//h_res_muTR1->Fit("gaus","R","",-0.0001,+0.0001);
 gStyle->SetOptFit(1111);
gStyle->SetOptStat(222001111);
//gPad->SetLogy();
d2.cd(2);
theta1g->SetLineColor(kOrange);
theta1g->Draw("hist");
gStyle->SetOptStat(222001111);
theta1->Draw("hist same");
//gPad->SetLogy();
d2.cd(3);
theta2g->SetLineColor(kOrange);
theta2g->Draw("hist");
theta2->Draw("hist same");
gStyle->SetOptStat(222001111);
//gPad->SetLogy();
d2.cd(4);
h_res_mu_inTR1->Draw();
h_res_mu_inTR1->Fit("f1","R","",-0.0001,+0.0001);
h_res_mu_inTR1->Draw("hist same");
gStyle->SetOptFit(1111);
gStyle->SetOptStat(222001111);
//gPad->SetLogy();
d2.SaveAs("pdf_try/mu_out_peren_GA.pdf");


TF1 *f4 = new TF1("f4", "[2]*TMath::Gaus(x,[0],[1])");
f4->SetParameters(-1.9e-05,30e-06,1);

TF1 *f5 = new TF1("f5", "[2]*TMath::Gaus(x,[0],[1])");
f5->SetParameters(-1.2e-05,28e-06,1);

TF1 *f6 = new TF1("f6", "[2]*TMath::Gaus(x,[0],[1])");
f6->SetParameters(-7e-05,28e-06,1);

TCanvas d3("d3","d3",700,700);
d3.Divide(2,2);
d3.cd(3);
h_res_muTR1x->SetMinimum(1.);
h_res_muTR1x->Draw();
h_res_muTR1x->Draw("hist same");
gPad->SetLogy();
gStyle->SetOptStat(222001111);
d3.cd(4);
h_res_muTR1y->SetMinimum(1.);
h_res_muTR1y->Draw();
h_res_muTR1y->Draw("hist same");
gPad->SetLogy();
gStyle->SetOptStat(222001111);
d3.cd(1);
h_res_muTR1x->Draw();
h_res_muTR1x->Draw("hist same");
h_res_muTR1x->Fit("f4","R","",-0.0001,+0.0001);
gStyle->SetOptFit(1111);
gStyle->SetOptStat(222001111);
//gPad->SetLogy();
d3.cd(2);
h_res_muTR1y->Draw();
h_res_muTR1y->Draw("hist same");
h_res_muTR1y->Fit("f5","R","",-0.0001,+0.0001);
gStyle->SetOptFit(1111);
gStyle->SetOptStat(222001111);
//gPad->SetLogy();
d3.SaveAs("pdf_try/proj_mu_GA.pdf");


TCanvas d3a("d3a","d3a",700,700);
d3a.Divide(2,3);
d3a.cd(1);
theta1gx->SetLineColor(kOrange);
theta1gx->Draw("hist");
theta1x->Draw("hist same");
gStyle->SetOptStat(222001111);
d3a.cd(2);
theta1gy->SetLineColor(kOrange);
theta1gy->Draw("hist");
theta1y->Draw("hist same");
gStyle->SetOptStat(222001111);
d3a.cd(3);
theta2gx->SetLineColor(kOrange);
theta2gx->Draw("hist");
theta2x->Draw("hist same");
gStyle->SetOptStat(222001111);
d3a.cd(4);
theta2gy->SetLineColor(kOrange);
theta2gy->Draw("hist");
theta2y->Draw("hist same");
gStyle->SetOptStat(222001111);
d3a.SaveAs("pdf_try/th_proj_mu_GA.pdf");

}


