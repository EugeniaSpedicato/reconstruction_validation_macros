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
double resx=0.; double resy=0.;
double thmu_gen_x=0; double thmu_gen_y=0;
double thmu_rec_x=0; double thmu_rec_y=0;
double the_gen_x=0; double the_gen_y=0;
double the_rec_x=0; double the_rec_y=0;
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

TH1D *h_phi =new TH1D("h_phi"," outgoing muon phi angle for events in the tails of the angular residuum",180,-180,180);
TH1D *h_phiP =new TH1D("h_phip"," outgoing muon phi angle for events in the peaks of the angular residuum",180,-180,180);

TH1D *h_phie =new TH1D("h_phie"," outgoing electron phi angle for events in the tails of the angular residuum",180,-180,180);
TH1D *h_phiPe =new TH1D("h_phipe"," outgoing electron phi angle for events in the peaks of the angular residuum",180,-180,180);

TH1D *h_phitot =new TH1D("h_phitot"," outgoing muon phi angle",180,-180,180);
TH1D *h_phietot =new TH1D("h_phietot"," outgoing muon phi angle",180,-180,180);

TH1D *theta1x =new TH1D("theta1x" , "theta x of the muon for events in the peak of angular residuum" , 150,-0.0025,0.0025);
TH1D *theta2x =new TH1D("theta2x" , "theta x of the muon for events in the tail of angular residuum" , 150,-0.0025,0.0025);
TH1D *theta1gx =new TH1D("theta1gx" , "theta x of the muon for events in the peak of angular residuum" , 150,-0.0025,0.0025);
TH1D *theta2gx =new TH1D("theta2gx" , "theta x of the muon for events in the tail of angular residuum" , 150,-0.0025,0.0025);

TH1D *theta1y =new TH1D("theta1y" , "theta y of the muon for events in the peak of angular residuum" , 150,-0.0025,0.0025);
TH1D *theta2y =new TH1D("theta2y" , "theta y of the muon for events in the tail of angular residuum" , 150,-0.0025,0.0025);

TH1D *theta1gy =new TH1D("theta1gy" , "theta y of the muon for events in the peak of angular residuum" , 150,-0.0025,0.0025);
TH1D *theta2gy =new TH1D("theta2gy" , "theta y of the muon for events in the tail of angular residuum" , 150,-0.0025,0.0025);

TH1D *h_cooX =new TH1D("X1","Reconstructed X1_mu-X1_e for events in tail of angular residuum",40,-0.2,0.2);
TH1D *h_cooX_g =new TH1D("X1g","Generated X1_mu-X1_e for events in tail of angular residuum",40,-0.2,0.2);
TH1D *h_cooXp =new TH1D("X1p","Reconstructed X1_mu-X1_e for events in peak of angular residuum",200,-1,1);
TH1D *h_cooX_gp =new TH1D("X1gp","Generated X1_mu-X1_e for events in peak of angular residuum",200,-1,1);

TH1D *h_stripX =new TH1D("sX1","strip_mu-strip_e for events in tail of angular residuum",40,-20,20);
TH1D *h_stripXp =new TH1D("sX1p","strip_mu-strip_e for events in peak of angular residuum",200,-100,100);


TH1D *h_seede=new TH1D("h_seede","Seed cluster size electron events in the tail",6,0.,6.);
TH1D *h_seedPe=new TH1D("h_seedPe","Seed cluster size electron events in the peak",6,0.,6.);
TH1D *h_seedmu=new TH1D("h_seedmu","Seed cluster size muon events in the tail",6,0.,6.);
TH1D *h_seedPmu=new TH1D("h_seedPmu","Seed cluster size muon events in the peak",6,0.,6.);

TH1D *h_corre=new TH1D("h_corre","corr cluster size electron events in the tail",6,0.,6.);
TH1D *h_corrPe=new TH1D("h_corrPe","corr cluster size electron events in the peak",6,0.,6.);

TH1D *h_corrmu=new TH1D("h_corrmu","corr cluster size muon events in the tail",6,0.,6.);
TH1D *h_corrPmu=new TH1D("h_corrPmu","corr cluster size muon events in the peak",6,0.,6.);

double phi=0;
double phie=0;

double tailx=0.;
double peakx=0.;
double taily=0.;
double peaky=0.;


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
         double mx=SigTracks->ax();
         double my=SigTracks->ay();
                the_gen_x=mx;
                the_gen_y=my;
                phie=pe.Phi()*(180/M_PI);;
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
                phi=pmu.Phi()*(180/M_PI);;

		}
	 }
	}

	// if(yes_e!=1 or yes_mu!=1) cout << "NOT RECONSTRUCTIBLE" << endl;

	if(yes_e==1 and yes_mu==1){

vector<double> Xmu; Xmu.reserve(4);
vector<double> Ymu; Ymu.reserve(4);
vector<double> Zmu; Zmu.reserve(8);

vector<double> Xe;

int ok=0;int oke=0;
double X1,X2, Y1,Y2, Z1,Z2,Z3,Z4, X1e;
double dX1=-99; double dX2=-99; double dY1=-99; double dY2=-99; double dX1e=-99;
double stripX=-999; double stripXe=-999;

	for(int s=0; s<TrackerPoints->GetEntries(); s++)
                         {const MUonETrackerPoint *TrackerPt = static_cast<const MUonETrackerPoint*>(TrackerPoints->At(s));

if(TrackerPt->trackPDGCode()==-13 and TrackerPt->trackID()==code_mu and TrackerPt->stationID()==1)
//if(TrackerPt->trackPDGCode()==-13 and TrackerPt->trackID()==0 and TrackerPt->stationID()==0)
//if(TrackerPt->trackPDGCode()==-13 and TrackerPt->trackID()==code_mu and TrackerPt->stationID()==1)
				{ ok=1;
				  TVector3 entering=TrackerPt->enteringPositionGlobalCoordinates();
				  TVector3 exiting=TrackerPt->exitingPositionGlobalCoordinates();

				  if(TrackerPt->moduleID()==0 or TrackerPt->moduleID()==4) Xmu.push_back(exiting.X());
				  if(TrackerPt->moduleID()==1 or TrackerPt->moduleID()==5) Ymu.push_back(exiting.Y());
				  if(TrackerPt->moduleID()!=2 and TrackerPt->moduleID()!=3) Zmu.push_back(exiting.Z());

				}
if(TrackerPt->trackPDGCode()==11 and TrackerPt->trackID()==code_e and TrackerPt->stationID()==1)
{ oke=1;
                                  TVector3 exiting=TrackerPt->exitingPositionGlobalCoordinates();

				  if(TrackerPt->moduleID()==0) Xe.push_back(exiting.X());

}


			 }


if(ok==1 and oke==1){
      X1=(Xmu.at(0)+Xmu.at(1))/2;
        X2=(Xmu.at(2)+Xmu.at(3))/2;
        Y1=(Ymu.at(0)+Ymu.at(1))/2;
        Y2=(Ymu.at(2)+Ymu.at(3))/2;

        Z1=(Zmu.at(0)+Zmu.at(1))/2;
        Z2=(Zmu.at(2)+Zmu.at(3))/2;
        Z3=(Zmu.at(4)+Zmu.at(5))/2;
        Z4=(Zmu.at(6)+Zmu.at(7))/2;

      X1e=(Xe.at(0)+Xe.at(1))/2;

        }



vector<double> digiXmu1; digiXmu1.reserve(4);
vector<double> digiYmu1; digiYmu1.reserve(4);
vector<double> digiXmu2; digiXmu2.reserve(4);
vector<double> digiYmu2; digiYmu2.reserve(4);

vector<double> digistripX1; digistripX1.reserve(4);
vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();

double seed_clustersize_e;
double seed_clustersize_mu;
double corr_clustersize_e;
double corr_clustersize_mu;

for(int j=0; j<tracks.size();j++)
{
 if(tracks.at(j).processIDofLinkedTrack()==45 and tracks.size()==3 and tracks.at(j).sector()==1 and tracks.at(j).percentageOfHitsSharedWithLinkedTrack()>=100)
        {
vector<MUonERecoOutputTrackHit> hits_=tracks.at(j).hits();
if(code_mu==tracks.at(j).linkedTrackID()){yes2++;
for(int h=0;h<hits_.size();h++){
if(hits_.at(h).moduleID()==0){dX1=hits_.at(h).positionPerpendicular(); stripX= hits_.at(h).seedClusterCenterStrip();
seed_clustersize_mu=hits_.at(h).seedClusterWidth();corr_clustersize_mu=hits_.at(h).correlationClusterWidth();}
if(hits_.at(h).moduleID()==4)dX2=hits_.at(h).positionPerpendicular();
if(hits_.at(h).moduleID()==1)dY1=hits_.at(h).positionPerpendicular();
if(hits_.at(h).moduleID()==5)dY2=hits_.at(h).positionPerpendicular();
                }
	}
if(code_e==tracks.at(j).linkedTrackID()){
for(int h=0;h<hits_.size();h++){if(hits_.at(h).moduleID()==0)
{dX1e=hits_.at(h).positionPerpendicular();seed_clustersize_e=hits_.at(h).seedClusterWidth();corr_clustersize_e=hits_.at(h).correlationClusterWidth();
						 stripXe= hits_.at(h).seedClusterCenterStrip();}	
		}

	}
   }
}

/*
if(tracks.size()==3){
        for(int t=0; t<TrackerStubs->GetEntries(); t++)
                         {const MUonETrackerStub *stubs = static_cast<const MUonETrackerStub*>(TrackerStubs->At(t));
	if(stubs->stationID()==1){ double stub=(stubs->seedClusterCenterStrip() + 0.5 + 0.5 * stubs->bend()) * 9.144 / 1016 - 0.5 * 9.144;
			  if(stubs->moduleID()==0) {digiXmu1.push_back(stub); digistripX1.push_back(stubs->seedClusterCenterStrip());}
                          if( stubs->moduleID()==4) {digiXmu2.push_back(stub);}
                          if(stubs->moduleID()==1) {digiYmu1.push_back(stub*sin(90));}
			  if(stubs->moduleID()==5) {digiYmu2.push_back(stub*sin(90));}
			}
		}
}

vector<double> diff1,diff2,diffe1;
for(int d=0; d<digiXmu1.size(); d++){
diff1.push_back(abs(digiXmu1.at(d)-X1));
diffe1.push_back(abs(digiXmu1.at(d)-X1e));
}

for(int d=0; d<digiXmu2.size(); d++){
diff2.push_back(abs(digiXmu2.at(d)-X2));
}


if(digiXmu1.size()!=0 and digiXmu2.size()!=0){
auto it = min_element(begin(diff1),end(diff1));
auto it2 = min_element(begin(diff2),end(diff2));
auto ite = min_element(begin(diffe1),end(diffe1));

dX1=digiXmu1.at( distance(begin(diff1), it) );
dX2=digiXmu2.at( distance(begin(diff2), it2) );
stripX = digistripX1.at( distance(begin(diff1), it) );

dX1e=digiXmu1.at( distance(begin(diffe1), ite) );
stripXe = digistripX1.at( distance(begin(diffe1), ite) );

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

*/

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
                the_rec_x=tracks.at(j).xSlope();
                the_rec_y=tracks.at(j).ySlope();
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



vector<double> digistrip;
vector<double> digibend;
vector<double> dig1;
double bend, bende;
if(tracks.size()==3){
        for(int t=0; t<TrackerStubs->GetEntries(); t++)
                         {const MUonETrackerStub *stubs = static_cast<const MUonETrackerStub*>(TrackerStubs->At(t));
        if(stubs->stationID()==1){ double stub=(stubs->seedClusterCenterStrip() + 0.5 ) * 9.144 / 1016 - 0.5 * 9.144;
                          if(stubs->moduleID()==0) {dig1.push_back(stub); digistrip.push_back(stubs->seedClusterCenterStrip());digibend.push_back(stubs->bend());}
		}
	}
}
vector<double> diff1,diffe1;
double stub;
double stube;
for(int d=0; d<digistrip.size(); d++){
diff1.push_back(abs(digistrip.at(d)-stripX));
diffe1.push_back(abs(digistrip.at(d)-stripXe));
}
if(digistrip.size()!=0){
auto it = min_element(begin(diff1),end(diff1));
auto ite = min_element(begin(diffe1),end(diffe1));

stub=dig1.at( distance(begin(diff1), it) );
stube=dig1.at( distance(begin(diffe1), ite) );
bend=digibend.at( distance(begin(diff1), it) );
bende=digibend.at( distance(begin(diffe1), ite) );
 }







double res_muTR2=-999;
if(ok_mu==1 and ok_e==1 ){reco+=MesmerEvent->wgt_full;

 res_muTR2 = thmu_rec_tracks2-thmu_gen;
double resx=thmu_rec_x-thmu_gen_x;
double resy=thmu_rec_y-thmu_gen_y;
/*double resx=the_rec_x-the_gen_x;
double resy=the_rec_y-the_gen_y;*/


h_phitot->Fill(phi,MesmerEvent->wgt_full);h_phietot->Fill(phie,MesmerEvent->wgt_full);

if(abs(dX1-X1)>0.006 or abs(dX2-X2)>0.006){
				   tailx+=MesmerEvent->wgt_full;h_phi->Fill(phi,MesmerEvent->wgt_full);h_phie->Fill(phie,MesmerEvent->wgt_full);
                                   theta2x->Fill(thmu_rec_x,MesmerEvent->wgt_full);theta2gx->Fill(thmu_gen_x,MesmerEvent->wgt_full);
					}
if(abs(dY1-Y1)>0.006 or abs(dY2-Y2)>0.006){
				   taily+=MesmerEvent->wgt_full;h_phi->Fill(phi,MesmerEvent->wgt_full);h_phie->Fill(phie,MesmerEvent->wgt_full);
                                   theta2y->Fill(thmu_rec_y,MesmerEvent->wgt_full);theta2gy->Fill(thmu_gen_y,MesmerEvent->wgt_full);
					}
if(abs(dX1-X1)<0.006 and abs(dX2-X2)<0.006 and abs(dY1-Y1)<0.006 and abs(dY2-Y2)<0.006) {
h_phiP->Fill(phi,MesmerEvent->wgt_full); h_phiPe->Fill(phie,MesmerEvent->wgt_full);
theta1x->Fill(thmu_rec_x,MesmerEvent->wgt_full);theta1gx->Fill(thmu_gen_x,MesmerEvent->wgt_full);
theta1y->Fill(thmu_rec_y,MesmerEvent->wgt_full);theta1gy->Fill(thmu_gen_y,MesmerEvent->wgt_full);
}

if(abs(dX1-X1)<0.006 and abs(dX2-X2)<0.006) peakx+=MesmerEvent->wgt_full;
if(abs(dY1-Y1)<0.006 and abs(dY2-Y2)<0.006)peaky+=MesmerEvent->wgt_full;

if( dX1!=-99 and dX2!=-99 and dY1!=-99 and dY2!=-99 and abs(resx)<0.2e-03 and abs(resy)<0.1e-03){

h_x1->Fill(stub-X1,MesmerEvent->wgt_full);// h_x1_g->Fill(digiXmu.at(0),MesmerEvent->wgt_full);
h_x2->Fill(dX2-X2,MesmerEvent->wgt_full); //h_x2_g->Fill(digiXmu.at(1),MesmerEvent->wgt_full);
h_y1->Fill(dY1-Y1,MesmerEvent->wgt_full); //h_y1_g->Fill(digiYmu.at(0),MesmerEvent->wgt_full);
h_y2->Fill(dY2-Y2,MesmerEvent->wgt_full); //h_y2_g->Fill(digiYmu.at(1),MesmerEvent->wgt_full);
}

if(abs(dX1e-X1e)>0.006 or abs(dX1-X1)>0.006){
cout <<"---------"<< endl;

 double corre= stripXe + (bende-3);
 double corr= stripX + (bend-3);
 cout << "T) corr electron " << corre << endl;
 cout << "T) corr muon " << corr << endl;

if(abs(dX1e-X1e)>0.006){h_seede->Fill(seed_clustersize_e,MesmerEvent->wgt_full);h_corre->Fill(corr_clustersize_e,MesmerEvent->wgt_full);
 }
if(abs(dX1-X1)>0.006){h_seedmu->Fill(seed_clustersize_mu,MesmerEvent->wgt_full);h_corrmu->Fill(corr_clustersize_mu,MesmerEvent->wgt_full);
 }
h_cooX->Fill(dX1-dX1e,MesmerEvent->wgt_full);
h_cooX_g->Fill(X1-X1e,MesmerEvent->wgt_full);
h_stripX->Fill(stripX-stripXe,MesmerEvent->wgt_full);
}

if(abs(dX1e-X1e)<0.006 and abs(dX1-X1)<0.006){

 double corre= stripXe + (bende-3);
 double corr= stripX + (bend-3);
 cout << "P) corr electron " << corre << endl;
 cout << "P) corr muon " << corr << endl;

if(abs(dX1e-X1e)<0.006){ h_seedPe->Fill(seed_clustersize_e,MesmerEvent->wgt_full);h_corrPe->Fill(corr_clustersize_e,MesmerEvent->wgt_full);}
if(abs(dX1-X1)<0.006){ h_seedPmu->Fill(seed_clustersize_mu,MesmerEvent->wgt_full);h_corrPmu->Fill(corr_clustersize_mu,MesmerEvent->wgt_full);}

h_cooXp->Fill(dX1-dX1e,MesmerEvent->wgt_full);
h_cooX_gp->Fill(X1-X1e,MesmerEvent->wgt_full);
h_stripXp->Fill(stripX-stripXe,MesmerEvent->wgt_full);
}


if(dX1!=-99 and dX2!=-99 and dY1!=-99 and dY2!=-99 and (abs(resx)>0.2e-03 or abs(resy)>0.1e-03)){

h_x1_g->Fill(dX1-X1,MesmerEvent->wgt_full);// h_x1_g->Fill(digiXmu.at(0),MesmerEvent->wgt_full);
h_x2_g->Fill(dX2-X2,MesmerEvent->wgt_full); //h_x2_g->Fill(digiXmu.at(1),MesmerEvent->wgt_full);
h_y1_g->Fill(dY1-Y1,MesmerEvent->wgt_full); //h_y1_g->Fill(digiYmu.at(0),MesmerEvent->wgt_full);
h_y2_g->Fill(dY2-Y2,MesmerEvent->wgt_full); //h_y2_g->Fill(digiYmu.at(1),MesmerEvent->wgt_full);
}


}

ok=0;
oke=0;
	}//end reconstuibility
yes_e=0;yes_mu=0;
code_e=-99;code_mu=-99;
point_el=0;point_mu=0;
TrackIdreco=-99;
yes2=0;yes_v=0;
} //end of general for

cout << "Su " << signal << " eventi di segnale, eventi con residui di posizioneX nella coda " << (tailx/reco)*100 << "%"<< endl;

cout << "Su " << signal << " eventi di segnale, eventi con residui di posizioneX nel picco " << (peakx/reco)*100 << "%"<< endl;
cout << "Su " << signal << " eventi di segnale, eventi con residui di posizioneY  nella coda " << (taily/reco)*100 << "%"<< endl;

cout << "Su " << signal << " eventi di segnale, eventi con residui di posizioneY nel picco " << (peaky/reco)*100 << "%"<< endl;

TCanvas f("f","f",700,700);
f.Divide(2,2);
f.cd(1);
h_cooXp->Draw("hist");
h_cooX_gp->SetLineColor(kOrange);
h_cooX_gp->Draw("hist same");
f.cd(2);
h_cooX->Draw("hist");
h_cooX_g->SetLineColor(kOrange);
h_cooX_g->Draw("hist same");
f.cd(3);
h_stripXp->Draw("hist");
f.cd(4);
h_stripX->Draw("hist");
f.SaveAs("X1coo.pdf");

TCanvas f1("f1","f1",700,700);
f1.Divide(2,2);
f1.cd(1);
h_seedPe->Draw("hist");
h_seede->SetLineColor(kViolet);
h_seede->Draw("hist same");
f1.cd(2);
h_seedPmu->Draw("hist");
h_seedmu->SetLineColor(kViolet);
h_seedmu->Draw("hist same");
f1.cd(3);
h_corrPe->Draw("hist");
h_corre->SetLineColor(kViolet);
h_corre->Draw("hist same");
f1.cd(4);
h_corrPmu->Draw("hist");
h_corrmu->SetLineColor(kViolet);
h_corrmu->Draw("hist same");
f1.SaveAs("cluster_width.pdf");


TCanvas d3a("d3a","d3a",700,700);
d3a.Divide(2,2);
d3a.cd(1);
h_x1->Draw("hist");
h_x1_g->SetLineColor(kRed);
gPad->SetLogy();
h_x1->SetMinimum(1.0);
h_x1_g->Draw("hist same");
gStyle->SetOptStat(222001111);
d3a.cd(2);
h_x2->Draw("hist");
h_x2_g->SetLineColor(kRed);
gPad->SetLogy();
h_x2->SetMinimum(1.0);
h_x2_g->Draw("hist same");
gStyle->SetOptStat(222001111);
d3a.cd(3);
h_y1->Draw("hist");
h_y1_g->SetLineColor(kRed);
gPad->SetLogy();
h_y1_g->Draw("hist same");
h_y1->SetMinimum(1.0);
gStyle->SetOptStat(222001111);
d3a.cd(4);
h_y2->Draw("hist");
h_y2_g->SetLineColor(kRed);
h_y2_g->Draw("hist same");
gPad->SetLogy();
h_y2->SetMinimum(1.0);
gStyle->SetOptStat(222001111);
d3a.SaveAs("h_digi.pdf");

/*TCanvas b1("b1","b1",700,700);
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
h_phie->Draw("hist");
gStyle->SetOptStat(222001111);
b1.SaveAs("phi_digit.pdf");

TCanvas d3aa("d3aa","d3aa",700,700);
d3aa.Divide(2,3);
d3aa.cd(1);
theta1gx->SetLineColor(kOrange);
theta1gx->Draw("hist");
theta1x->Draw("hist same");
gStyle->SetOptStat(222001111);
d3aa.cd(2);
theta1gy->SetLineColor(kOrange);
theta1gy->Draw("hist");
theta1y->Draw("hist same");
gStyle->SetOptStat(222001111);
d3aa.cd(3);
theta2gx->SetLineColor(kOrange);
theta2gx->Draw("hist");
theta2x->Draw("hist same");
gStyle->SetOptStat(222001111);
d3aa.cd(4);
theta2gy->SetLineColor(kRed);
theta2gy->Draw("hist");
theta2y->Draw("hist same");
gStyle->SetOptStat(222001111);
d3aa.SaveAs("th_proj_mu_GA_digit.pdf");
*/

}



