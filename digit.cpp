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

	TFile *inputfile = new TFile("TRMesmer_box_100k_2GeV_0_noTilt.root");
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


TH1D *h_res_muTR2r=new TH1D("h_res_muTR2r", "(thmu_rec-thmu_true) VS energy 157<Emu<160(GeV) PRE-VRTX",180,-0.0006,0.0006);


TH1D *h_x1=new TH1D("h_x1","residuum (x_g4-x_digit) module 0 hit peak angular res. (blue) and tail angular res. (red)",100,-0.05,0.05);//100,-0.05,0.05);
TH1D *h_x1_g=new TH1D("h_x1g","residuum (x_g4-x_digit) module 0 hit peak angular res. (blue) and tail angular res. (red)",100,-0.05,0.05);//100,-0.05,0.05);
TH1D *h_x2=new TH1D("h_x2","residuum (x_g4-x_digit) module 4 hit peak angular res. (blue) and tail angular res. (red)",100,-0.05,0.05);//100,-0.05,0.05);
TH1D *h_x2_g=new TH1D("h_x2g","residuum (x_g4-y_digit) module 4 hit peak angular res. (blue) and tail angular res. (red)",100,-0.05,0.05);//100,-0.05,0.05);
TH1D *h_y1=new TH1D("h_y1","residuum (y_g4-y_digit) module 1 hit peak angular res. (blue) and tail angular res. (red)",100,-0.05,0.05);//100,-0.05,0.05);
TH1D *h_y1_g=new TH1D("h_y1g","residuum (y_g4-y_digit) module 1 hit peak angular res. (blue) and tail angular res. (red)",100,-0.05,0.05);//100,-0.05,0.05);
TH1D *h_y2=new TH1D("h_y2","residuum (y_g4-y_digit) module 5 hit peak angular res. (blue) and tail angular res. (red)",100,-0.05,0.05);//100,-0.05,0.05);
TH1D *h_y2_g=new TH1D("h_y2g","residuum (y_g4-y_digit) module 5 hit peak angular res. (blue) and tail angular res. (red)",100,-0.05,0.05);//100,-0.05,0.05);

for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) { //for(Long64_t i = numb; i < numb+1; i++) {
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
		yes_mu=1;

		}
	 }
	}

	// if(yes_e!=1 or yes_mu!=1) cout << "NOT RECONSTRUCTIBLE" << endl;

	if(yes_e==1 and yes_mu==1){

vector<double> Xmu; Xmu.reserve(4);
vector<double> Ymu; Ymu.reserve(4);
vector<double> Zmu; Zmu.reserve(8);

int ok=0;
double X1,X2, Y1,Y2, Z1,Z2,Z3,Z4;
double dX1=-99; double dX2=-99; double dY1=-99; double dY2=-99;

	for(int s=0; s<TrackerPoints->GetEntries(); s++)
                         {const MUonETrackerPoint *TrackerPt = static_cast<const MUonETrackerPoint*>(TrackerPoints->At(s));

if(TrackerPt->trackPDGCode()==-13 and TrackerPt->trackID()==code_mu and TrackerPt->stationID()==1)
//if(TrackerPt->trackPDGCode()==-13 and TrackerPt->trackID()==0 and TrackerPt->stationID()==0)
//if(TrackerPt->trackPDGCode()==-13 and TrackerPt->trackID()==code_mu and TrackerPt->stationID()==1)
				{ ok=1;
				  TVector3 entering=TrackerPt->enteringPositionLocalCoordinates();
				  TVector3 exiting=TrackerPt->exitingPositionLocalCoordinates();

				  if(TrackerPt->moduleID()==0 or TrackerPt->moduleID()==4) Xmu.push_back(exiting.X());
				  if(TrackerPt->moduleID()==1 or TrackerPt->moduleID()==5) Ymu.push_back(exiting.Y());
				  if(TrackerPt->moduleID()!=2 and TrackerPt->moduleID()!=3) Zmu.push_back(exiting.Z());

				}
			 }

if(ok==1){
      X1=(Xmu.at(0)+Xmu.at(1))/2;
        X2=(Xmu.at(2)+Xmu.at(3))/2;
        Y1=(Ymu.at(0)+Ymu.at(1))/2;
        Y2=(Ymu.at(2)+Ymu.at(3))/2;

        Z1=(Zmu.at(0)+Zmu.at(1))/2;
        Z2=(Zmu.at(2)+Zmu.at(3))/2;
        Z3=(Zmu.at(4)+Zmu.at(5))/2;
        Z4=(Zmu.at(6)+Zmu.at(7))/2;
        }

vector<double> digiXmu1; digiXmu1.reserve(4);
vector<double> digiYmu1; digiYmu1.reserve(4);
vector<double> digiXmu2; digiXmu2.reserve(4);
vector<double> digiYmu2; digiYmu2.reserve(4);
vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();

/*
for(int j=0; j<tracks.size();j++)
{
 if(tracks.at(j).processIDofLinkedTrack()==45 and tracks.size()==3 and tracks.at(j).sector()==1 and tracks.at(j).percentageOfHitsSharedWithLinkedTrack()>=100)
        {
vector<MUonERecoOutputTrackHit> hits_=tracks.at(j).hits();
if(code_mu==tracks.at(j).linkedTrackID()){yes2++;
for(int h=0;h<hits_.size();h++){
if(hits_.at(h).moduleID()==0 or hits_.at(h).moduleID()==4) {digiXmu.push_back(hits_.at(h).positionPerpendicular());}
if(hits_.at(h).moduleID()==1 or hits_.at(h).moduleID()==5) {digiYmu.push_back(hits_.at(h).positionPerpendicular());}
                }
        }
   }
}
*/


if(tracks.size()==3){
        for(int t=0; t<TrackerStubs->GetEntries(); t++)
                         {const MUonETrackerStub *stubs = static_cast<const MUonETrackerStub*>(TrackerStubs->At(t));
	if(stubs->stationID()==1){ double stub=(stubs->seedClusterCenterStrip() + 0.5 + 0.5 * stubs->bend()) * 9.144 / 1016 - 0.5 * 9.144;
			  if(stubs->moduleID()==0) {digiXmu1.push_back(stub);}
                          if( stubs->moduleID()==4) {digiXmu2.push_back(stub);}
                          if(stubs->moduleID()==1) {digiYmu1.push_back(stub*sin(90));}
			  if(stubs->moduleID()==5) {digiYmu2.push_back(stub*sin(90));}
			}
		}
}

vector<double> diff1,diff2;
for(int d=0; d<digiXmu1.size(); d++){
diff1.push_back(abs(digiXmu1.at(d)-X1));
}

for(int d=0; d<digiXmu2.size(); d++){
diff2.push_back(abs(digiXmu2.at(d)-X2));
}

if(digiXmu1.size()!=0 and digiXmu2.size()!=0){
auto it = min_element(begin(diff1),end(diff1));
auto it2 = min_element(begin(diff2),end(diff2));

dX1=digiXmu1.at( distance(begin(diff1), it) );
dX2=digiXmu2.at( distance(begin(diff2), it2) );

 }

vector<double> diff1Y,diff2Y;
for(int d=0; d<digiYmu1.size(); d++){
diff1Y.push_back(abs(digiYmu1.at(d)-Y1));
}

for(int d=0; d<digiYmu2.size(); d++){
diff2Y.push_back(abs(digiYmu2.at(d)-Y2));
}

if(digiYmu1.size()!=0 and digiYmu2.size()!=0){
auto itY = min_element(begin(diff1Y),end(diff1Y));
auto it2Y = min_element(begin(diff2Y),end(diff2Y));

dY1=digiYmu1.at( distance(begin(diff1Y), itY) );
dY2=digiYmu2.at( distance(begin(diff2Y), it2Y) );

 }

int ok_mu=0;
int ok_e=0;
double thmu_rec_tracks2=0.;
for(int j=0; j<tracks.size();j++)
{
 if(tracks.at(j).processIDofLinkedTrack()==45 and tracks.size()==3 and tracks.at(j).sector()==1 and tracks.at(j).percentageOfHitsSharedWithLinkedTrack()>=100)
        {

               if(code_e==tracks.at(j).linkedTrackID()) { ok_e=1;
                 TVector3 e_outv(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0);
                TVector3 z(0.,0.,1.);
                e_outv=e_outv.Unit();
                 the_rec_tracks=e_outv.Theta();
                                        }

                 if(code_mu==tracks.at(j).linkedTrackID()) { ok_mu=1;
                 TVector3 mu_outv(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0);
                TVector3 z(0.,0.,1.);
                mu_outv=mu_outv.Unit();
                //thmu_rec_tracks=acos(mu_outv.Dot(z));
                thmu_rec_x=tracks.at(j).xSlope();
                thmu_rec_y=tracks.at(j).ySlope();

                 //thmu_rec_tracks=mu_outv.Theta();
                 thmu_rec_tracks2=sqrt(thmu_rec_x*thmu_rec_x+thmu_rec_y*thmu_rec_y);
                                        }


        }
}

double res_muTR2=-999;
if(ok_mu==1 and ok_e==1 ){

 res_muTR2 = thmu_rec_tracks2-thmu_gen;
}

if(ok==1 and ok_mu==1 and dX1!=-99 and dX2!=-99 and dY1!=-99 and dY2!=-99 and abs(res_muTR2)<0.1e-03){

h_x1->Fill(X1-dX1,MesmerEvent->wgt_full);// h_x1_g->Fill(digiXmu.at(0),MesmerEvent->wgt_full);
h_x2->Fill(X2-dX2,MesmerEvent->wgt_full); //h_x2_g->Fill(digiXmu.at(1),MesmerEvent->wgt_full);
h_y1->Fill(Y1-dY1,MesmerEvent->wgt_full); //h_y1_g->Fill(digiYmu.at(0),MesmerEvent->wgt_full);
h_y2->Fill(Y2-dY2,MesmerEvent->wgt_full); //h_y2_g->Fill(digiYmu.at(1),MesmerEvent->wgt_full);
}

if(ok==1 and ok_mu==1 and dX1!=-99 and dX2!=-99 and dY1!=-99 and dY2!=-99 and abs(res_muTR2)>0.1e-03){

h_x1_g->Fill(X1-dX1,MesmerEvent->wgt_full);// h_x1_g->Fill(digiXmu.at(0),MesmerEvent->wgt_full);
h_x2_g->Fill(X2-dX2,MesmerEvent->wgt_full); //h_x2_g->Fill(digiXmu.at(1),MesmerEvent->wgt_full);
h_y1_g->Fill(Y1-dY1,MesmerEvent->wgt_full); //h_y1_g->Fill(digiYmu.at(0),MesmerEvent->wgt_full);
h_y2_g->Fill(Y2-dY2,MesmerEvent->wgt_full); //h_y2_g->Fill(digiYmu.at(1),MesmerEvent->wgt_full);
}

ok=0;

	}//end reconstuibility
yes_e=0;yes_mu=0;
code_e=-99;code_mu=-99;
point_el=0;point_mu=0;
TrackIdreco=-99;
yes2=0;yes_v=0;
} //end of general for


TCanvas d3a("d3a","d3a",700,700);
d3a.Divide(2,2);
d3a.cd(1);
h_x1->Draw("hist");
h_x1_g->SetLineColor(kRed);
h_x1_g->Draw("hist same");

d3a.cd(2);
h_x2->Draw("hist");
h_x2_g->SetLineColor(kRed);
h_x2_g->Draw("hist same");

d3a.cd(3);
h_y1->Draw("hist");
h_y1_g->SetLineColor(kRed);
h_y1_g->Draw("hist same");


d3a.cd(4);
h_y2->Draw("hist");
h_y2_g->SetLineColor(kRed);
h_y2_g->Draw("hist same");


d3a.SaveAs("h_digi.pdf");

}



