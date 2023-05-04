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

	TFile *inputfile = new TFile("TRMesmer_box_100k_2GeV_0.root");
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

TH1D *h_res_muTR1x=new TH1D("h_res_muTR1x", "(thmu_rec-thmu_true) proj. X  155<Emu<157(GeV) PRE-VRTX",120,-0.0006,0.0006);
TH1D *h_res_muTR2x=new TH1D("h_res_muTR2x", "(thmu_rec-thmu_true) proj. X 157<Emu<160(GeV) PRE-VRTX",120,-0.0006,0.0006);
TH1D *h_res_muTR3x=new TH1D("h_res_muTR3x", "(thmu_rec-thmu_true) proj. X 158<Emu<160(GeV) PRE-VRTX",120,-0.0006,0.0006);

TH1D *h_res_muTR1y=new TH1D("h_res_muTR1y", "(thmu_rec-thmu_true) proj. Y 155<Emu<157(GeV) PRE-VRTX",120,-0.0006,0.0006);
TH1D *h_res_muTR2y=new TH1D("h_res_muTR2y", "(thmu_rec-thmu_true) proj. Y 157<Emu<160(GeV) PRE-VRTX",120,-0.0006,0.0006);
TH1D *h_res_muTR3y=new TH1D("h_res_muTR3y", "(thmu_rec-thmu_true) proj. Y 158<Emu<160(GeV) PRE-VRTX",120,-0.0006,0.0006);

TH1D *h_res_muTR1=new TH1D("h_res_muTR1", "(thmu_g4-thmu_true) VS energy 155<Emu<157(GeV) PRE-VRTX",180,-0.0006,0.0006);
TH1D *h_res_muTR2=new TH1D("h_res_muTR2", "(thmu_g4-thmu_true) VS energy 157<Emu<160(GeV) PRE-VRTX",180,-0.0006,0.0006);


TH1D *h_res_muTR1r=new TH1D("h_res_muTR1r", "(thmu_rec-thmu_true) VS energy 155<Emu<157(GeV) PRE-VRTX",180,-0.0006,0.0006);
TH1D *h_res_muTR2r=new TH1D("h_res_muTR2r", "(thmu_rec-thmu_true) VS energy 157<Emu<160(GeV) PRE-VRTX",180,-0.0006,0.0006);


TH1D *theta1 =new TH1D("theta1" , "theta muon peak residuum" , 150,0.,0.005);
TH1D *theta2 =new TH1D("theta2" , "theta muon tail residuum" , 150,0.,0.005);
TH1D *theta1g =new TH1D("theta1g" , "theta muon peak residuum" , 150,0.,0.005);
TH1D *theta2g =new TH1D("theta2g" , "theta muon tail residuum" , 150,0.,0.005);

TH1D *theta1x =new TH1D("theta1x" , "theta x muon peak residuum" , 150,-0.0025,0.0025);
TH1D *theta2x =new TH1D("theta2x" , "theta x muon tail residuum" , 150,-0.0025,0.0025);
TH1D *theta1gx =new TH1D("theta1gx" , "theta x muon peak residuum" , 150,-0.0025,0.0025);
TH1D *theta2gx =new TH1D("theta2gx" , "theta x muon tail residuum" , 150,-0.0025,0.0025);

TH1D *theta1y =new TH1D("theta1y" , "theta y muon peak residuum" , 150,-0.0025,0.0025);
TH1D *theta2y =new TH1D("theta2y" , "theta y muon tail residuum" , 150,-0.0025,0.0025);
TH1D *theta1gy =new TH1D("theta1gy" , "theta y muon peak residuum" , 150,-0.0025,0.0025);
TH1D *theta2gy =new TH1D("theta2gy" , "theta y muon tail residuum" , 150,-0.0025,0.0025);

TH1D *theta2x_if =new TH1D("theta2x_if" , "theta x muon residuum if theta y in the tail" , 150,-0.0025,0.0025);
TH1D *theta2gx_if =new TH1D("theta2gx_if" , "theta x muon residuum if theta y in the tail" , 150,-0.0025,0.0025);
TH1D *theta2y_if =new TH1D("theta2y_if" , "theta y muon residuum if theta x in the tail" , 150,-0.0025,0.0025);
TH1D *theta2gy_if =new TH1D("theta2gy_if" , "theta y muon peak residuum if theta x in the tail" , 150,-0.0025,0.0025);


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

                cout << "thmu_gen " << thmu_gen << endl;
                cout << "mx " << thmu_gen_x << endl;
                cout << "my " << thmu_gen_y << endl;

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
double startZ,startX,startY;

	for(int n = 0; n < SignalTracks->GetEntries(); n++) {
        const MUonETrack *SigTracks = static_cast<const MUonETrack*>(SignalTracks->At(n));
	if(SigTracks->interactionID()==45 and SigTracks->pdgCode()==-13 ){startZ=SigTracks->startZ(); startX=SigTracks->startX(); startY=SigTracks->startY();}

	for(int s=0; s<TrackerPoints->GetEntries(); s++)
                         {const MUonETrackerPoint *TrackerPt = static_cast<const MUonETrackerPoint*>(TrackerPoints->At(s));
			  if(TrackerPt->trackPDGCode()==-13 and SigTracks->interactionID()==45 and TrackerPt->trackID()==code_mu and TrackerPt->stationID()==1)
				{ ok=1;
				  TVector3 entering=TrackerPt->enteringPositionGlobalCoordinates();
				  TVector3 exiting=TrackerPt->exitingPositionGlobalCoordinates();

				  if(TrackerPt->moduleID()==0 or TrackerPt->moduleID()==4) Xmu.push_back(exiting.X());
				  if(TrackerPt->moduleID()==1 or TrackerPt->moduleID()==5) Ymu.push_back(exiting.Y());
				  if(TrackerPt->moduleID()!=2 and TrackerPt->moduleID()!=3) Zmu.push_back(exiting.Z());

				}
			 }
		}

vector<double> noXmu; noXmu.reserve(4);
vector<double> noYmu; noYmu.reserve(4);

vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
for(int j=0; j<tracks.size();j++)
{
 if(tracks.at(j).processIDofLinkedTrack()==45 and tracks.size()==3 and tracks.at(j).sector()==1 and tracks.at(j).percentageOfHitsSharedWithLinkedTrack()>=100)
        {
std::vector<MUonERecoOutputTrackHit> hits_=tracks.at(j).hits();
if(code_mu==tracks.at(j).linkedTrackID()){yes2++;
for(int h=0;h<hits_.size();h++){
if(hits_.at(h).moduleID()==0 or hits_.at(h).moduleID()==4) {noXmu.push_back(hits_.at(h).positionPerpendicular());}
if(hits_.at(h).moduleID()==1 or hits_.at(h).moduleID()==5) {noYmu.push_back(hits_.at(h).positionPerpendicular());}
                }
        }
   }
}

if(ok==1 and yes2==1 and noXmu.size()==2 and noYmu.size()==2){
      X1=(Xmu.at(0)+Xmu.at(1))/2;
        X2=(Xmu.at(2)+Xmu.at(3))/2;
        Y1=(Ymu.at(0)+Ymu.at(1))/2;
        Y2=(Ymu.at(2)+Ymu.at(3))/2;

        Z1=(Zmu.at(0)+Zmu.at(1))/2;
        Z2=(Zmu.at(2)+Zmu.at(3))/2;
        Z3=(Zmu.at(4)+Zmu.at(5))/2;
        Z4=(Zmu.at(6)+Zmu.at(7))/2;

auto gx = new TGraphErrors();
auto gy = new TGraphErrors();

gx->SetPoint(0,startZ,startX);
gy->SetPoint(0,startZ,startY);
gx->SetPoint(1,Z1,X1);
gy->SetPoint(1,Z2,Y1);
gx->SetPoint(2,Z3,X2);
gy->SetPoint(2,Z4,Y2);
gx->SetPointError(0,0.,0.0012);
gx->SetPointError(1,0.,0.0012);
gx->SetPointError(2,0.,0.0012);
gy->SetPointError(0,0.,0.0012);
gy->SetPointError(1,0.,0.0012);
gy->SetPointError(2,0.,0.0012);


TF1 *linx = new TF1("linx", "[0]+[1]*x");
TF1 *liny = new TF1("liny", "[0]+[1]*x");
gx->Fit("linx");
gy->Fit("liny");

double mxf=linx->GetParameter(1);
double myf=liny->GetParameter(1);
thmu_rec_tracks=sqrt(mxf*mxf+myf*myf);

                cout << "thmu_rec " << thmu_rec_tracks << endl;
                cout << "mxr " << mxf << endl;
                cout << "myr " << myf << endl;

double res_muTR = thmu_rec_tracks-thmu_gen;
double resx=mxf-thmu_gen_x;
double resy=myf-thmu_gen_y;

if(Emu>150 and Emu<=157)
{h_res_muTR1->Fill(res_muTR,MesmerEvent->wgt_full);h_res_muTR1x->Fill(resx,MesmerEvent->wgt_full); h_res_muTR1y->Fill(resy,MesmerEvent->wgt_full);}
if(Emu>157 and Emu<=160)
{h_res_muTR2->Fill(res_muTR,MesmerEvent->wgt_full);h_res_muTR2x->Fill(resx,MesmerEvent->wgt_full); h_res_muTR2y->Fill(resy,MesmerEvent->wgt_full);}


if(res_muTR>0.1e-03){theta1->Fill(thmu_rec_tracks,MesmerEvent->wgt_full);theta1g->Fill(thmu_gen,MesmerEvent->wgt_full);}
if(res_muTR<0.1e-03 and res_muTR>-0.1e-03){theta2->Fill(thmu_rec_tracks,MesmerEvent->wgt_full);theta2g->Fill(thmu_gen,MesmerEvent->wgt_full);}


if(abs(resx)<0.1e-03 and abs(resy)<0.1e-03) {theta1x->Fill(mxf,MesmerEvent->wgt_full);theta1gx->Fill(thmu_gen_x,MesmerEvent->wgt_full);
                                                theta1y->Fill(myf,MesmerEvent->wgt_full);theta1gy->Fill(thmu_gen_y,MesmerEvent->wgt_full);}

if(abs(resx)>0.1e-03){theta2x->Fill(mxf,MesmerEvent->wgt_full);theta2gx->Fill(thmu_gen_x,MesmerEvent->wgt_full);
                                   theta2y_if->Fill(myf,MesmerEvent->wgt_full);theta2gy_if->Fill(thmu_gen_y,MesmerEvent->wgt_full);}
if(abs(resy)>0.1e-03){theta2y->Fill(myf,MesmerEvent->wgt_full);theta2gy->Fill(thmu_gen_y,MesmerEvent->wgt_full);
                                   theta2x_if->Fill(mxf,MesmerEvent->wgt_full);theta2gx_if->Fill(thmu_gen_x,MesmerEvent->wgt_full);}

delete gx,gy,linx,liny;

        }
ok=0;
yes2=0;

//vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();

int ok_mu=0;
int ok_e=0;
double thmu_rec_tracks2=0.;
for(int j=0; j<tracks.size();j++)
{
 if(tracks.at(j).processIDofLinkedTrack()==45 and tracks.size()==3 and tracks.at(j).sector()==1 and tracks.at(j).percentageOfHitsSharedWithLinkedTrack()>=100)
        {yes2++;

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

if(ok_mu==1 and ok_e==1 ){

double res_muTR2 = thmu_rec_tracks2-thmu_gen;
double resx2=thmu_rec_x-thmu_gen_x;
double resy2=thmu_rec_y-thmu_gen_y;

                cout << "thmu_rec " << thmu_rec_tracks2 << endl;
                cout << "thmu_rec_x " << thmu_rec_x << endl;
                cout << "thmu_rec_y " << thmu_rec_y << endl;

if(Emu>150 and Emu<=157) {h_res_muTR1r->Fill(res_muTR2,MesmerEvent->wgt_full);}
//h_res_muTR1xr->Fill(resx2,MesmerEvent->wgt_full); h_res_muTR1yr->Fill(resy2,MesmerEvent->wgt_full);}
if(Emu>157 and Emu<=160) {h_res_muTR2r->Fill(res_muTR2,MesmerEvent->wgt_full);}
//h_res_muTR2xr->Fill(resx2,MesmerEvent->wgt_full); h_res_muTR2yr->Fill(resy2,MesmerEvent->wgt_full);}
}



	}//end reconstuibility
yes_e=0;yes_mu=0;
code_e=-99;code_mu=-99;
point_el=0;point_mu=0;
TrackIdreco=-99;
yes2=0;yes_v=0;
} //end of general for


/*TF1 *lin = new TF1("lin", "[0]+[1]*x");

TCanvas a("a","a",700,700);
a.Divide(1,2);
a.cd(1);
gx->SetMinimum(-6);
gx->SetMaximum(6);
gx->SetMarkerColor(kRed);
gx->SetTitle("x projection");
gx->Fit("lin");
gx->Draw("A* same");
gStyle->SetOptFit(1111);
a.cd(2);
gy->SetMinimum(-6);
gy->SetMaximum(6);
gy->SetMarkerColor(kRed);
gy->SetTitle("y projection");
gy->Fit("lin");
gy->Draw("A* same");
gStyle->SetOptFit(1111);
a.SaveAs("g_manual.pdf");
*/

TF1 *f1 = new TF1("f1", "[2]*TMath::Gaus(x,[0],[1])");
f1->SetParameters(7.7e-06,30e-06,1);

 TF1 *f2 = new TF1("f2", "[2]*TMath::Gaus(x,[0],[1])");
f2->SetParameters(3.8e-05,25e-06,1);

TF1 *f1r = new TF1("f1r", "[2]*TMath::Gaus(x,[0],[1])");
f1r->SetParameters(7.7e-06,30e-06,1);

 TF1 *f2r = new TF1("f2r", "[2]*TMath::Gaus(x,[0],[1])");
f2r->SetParameters(3.8e-05,25e-06,1);

TCanvas d2("d2","d2",700,700);
d2.Divide(2,2);
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
h_res_muTR1r->Draw();
h_res_muTR1r->Draw("hist same");
h_res_muTR1r->Fit("f1r","R","",-0.0001,+0.0001);
//h_res_muTR1->Fit("gaus","R","",-0.0001,+0.0001);
 gStyle->SetOptFit(1111);
d2.cd(4);
h_res_muTR2r->Draw();
h_res_muTR2r->Draw("hist same");
h_res_muTR2r->Fit("f2r","R","",-0.0001,+0.0001);
//h_res_muTR2->Fit("gaus","R","",-0.0001,+0.0001);
gStyle->SetOptFit(1111);
d2.SaveAs("g4_mu_out_manual.pdf");


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
d3.SaveAs("g4_proj_mu_GA_manual.pdf");


TCanvas d3a("d3a","d3a",700,700);
d3a.Divide(2,3);
d3a.cd(1);
theta1gx->SetLineColor(kOrange);
theta1gx->Draw("hist");
theta1x->Draw("hist same");
d3a.cd(2);
theta1gy->SetLineColor(kOrange);
theta1gy->Draw("hist");
theta1y->Draw("hist same");
d3a.cd(3);
theta2gx->SetLineColor(kOrange);
theta2gx->Draw("hist");
theta2x->Draw("hist same");
d3a.cd(4);
theta2gy_if->SetLineColor(kOrange);
theta2gy_if->Draw("hist");
theta2y_if->Draw("hist same");
d3a.cd(5);
theta2gx_if->SetLineColor(kOrange);
theta2gx_if->Draw("hist");
theta2x_if->Draw("hist same");
d3a.cd(6);
theta2gy->SetLineColor(kRed);
theta2gy->Draw("hist");
theta2y->Draw("hist same");
d3a.SaveAs("g4_manual_th_proj_mu_GA.pdf");

}


