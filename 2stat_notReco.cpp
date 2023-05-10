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

	TFile *inputfile = new TFile("TRMesmer_box_100k.root");//narrowBP_20GeV_1shared.root");//trPoints,20GeV.root");//exampleProductionJob.root");narrowBP.root
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

TH1D* inter= new TH1D("inter","interactionID reco", 50,0,50);
TH1D *tracksize = new TH1D("tr","reco_tracks.size() event reco", 10,0,10);

TH1D* inter0= new TH1D("inter0","interactionID when event noreco and 0 signal particle reco", 50,0,50);
TH1D *tracksize0 = new TH1D("tr0","reco_tracks.size() event noreco and 0 signal particle reco", 10,0,10);
TH1D *h_anglee0 = new TH1D("e0","electron angle when not reco in event noreco and 0 signal particle reco",400,0,0.04);
TH1D *h_anglemu0 = new TH1D("mu0","muon angle when not reco in event noreco and 0 signal particle reco",50,0,0.005);
TH2D *h_angen_el = new TH2D("emu0","electron angle-energy when not reco",400,0,0.04,120,0,60);
TH1F* h_pte0=new TH1F("ynr","Tracker pointes electron event noreco and 0 signal particle reco", 20,0,20);
TH1F* h_ptmu0=new TH1F("xnr","Tracker pointes muon event noreco and 0 signal particle reco", 20,0,20);
TH1F *h_trackerStubs0=new TH1F("trs0","#TrackerStubs event noreco and 0 signal particle reco",35,0,35);
TH1D *h_op0 = new TH1D("op0","opening angle when not reco in event noreco and 0 signal particle reco",400,0,0.04);

TH1D* inter1= new TH1D("inter1","interactionID when event noreco and 1 signal particle reco", 50,0,50);
TH1D *tracksize1 = new TH1D("tr1","reco_tracks.size() event noreco and 1 signal particle reco", 10,0,10);
TH1D *h_anglee1 = new TH1D("e1","electron angle when not reco in event noreco and 1 signal particle reco",400,0,0.04);
TH1D *h_anglemu1 = new TH1D("mu1","muon angle when not reco in event noreco and 1 signal particle reco",50,0,0.005);
TH2D *h_angle1 = new TH2D("emu1","electron-muon angle when not reco in event noreco and 1 signal particle reco",400,0,0.04,50,0,0.005);
TH1F* h_pte1=new TH1F("ynr1","Tracker pointes electron event noreco and 1 signal particle reco", 20,0,20);
TH1F* h_ptmu1=new TH1F("xnr1","Tracker pointes muon event noreco and 1 signal particle reco", 20,0,20);
TH1F *h_trackerStubs1=new TH1F("trs1","#TrackerStubs event noreco and 1 signal particle reco",35,0,35);
TH1D *h_op1 = new TH1D("op1","opening angle when not reco in event noreco and 1 signal particle reco",400,0,0.04);

TH1D* interM= new TH1D("interM","interactionID when reco and yes>2", 50,0,50);
TH1D *tracksizeM = new TH1D("trM","reco_tracks.size() event reco and yes>2", 10,0,10);

TH1D *h_part=new TH1D("h_part","reconstructed particle when event is NOT reco", 5,10,15);

TH1D *h_chidofM=new TH1D("h_chidofM", "chi2 per dof when I have also fake signal tracks",200,0,200);
TH1D *h_chidof2=new TH1D("h_chidof2", "chi2 per dof when I have only 2 signal tracks",200,0,200);

TH1D *h_shared=new TH1D("h_shared","moduels where shared hits are",6,0,6);
TH1D *h_op = new TH1D("op","opening angle when reco",400,0,0.04);

           auto trackatZ = [](double q, double m,double z) {

                return q + (z ) * m;
            };

double the_gen=0; double thmu_gen=0; double thmu=0; double the=0;
double signal=0; double reco=0; double reco1=0; double more_reco=0; double reco0=0;
int link0; int link1=0;
vector<int> link;
int yes_e=0;int yes_mu=0; int yes2=0;
int point_mu=0; int point_el=0;
int code_mu=-99; int code_e=-99;
int TrackIdreco=-99;
for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
		cbmsim->GetEntry(i);
		if(i%1000 == 0) cout<<"Entry "<<i<<endl;
	TVector3 pmuin,pe,pmu;

	for(int n = 0; n < MCTrack->GetEntries(); n++) {
	 const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
	 if(MCTr->interactionID()==0 and MCTr->pdgCode()==-13){
	 pmuin.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz());
	 double energy=MCTr->energy();
	 double mx=MCTr->ax();
	 double my=MCTr->ay();
	 double qx=MCTr->bx();
         double qy=MCTr->by();
	 double x1= qx+mx*(118.0218);
         double y1= qy+my*(121.8693);

	 }
         if(MCTr->interactionID()==45 and MCTr->pdgCode()==11) code_e=n;
         if(MCTr->interactionID()==45 and MCTr->pdgCode()==-13) code_mu=n;
	}

        TVector3 pmuin_dir=pmuin.Unit(); 
	TVector3 pmu_dir,pe_dir;
double thmu_sdr,the_sdr;
double Ee=0;
	for(int n = 0; n < SignalTracks->GetEntries(); n++) {
	const MUonETrack *SigTracks = static_cast<const MUonETrack*>(SignalTracks->At(n));
	 if(SignalTracks->GetEntries()>1 and SigTracks->interactionID()==45 )
	 {
           double mx=SigTracks->ax();
           double my=SigTracks->ay();
           double qx=SigTracks->bx();
           double qy=SigTracks->by();

           int last_modXmu=0; int last_modXe=0; //= qx+mx*(189.9218);
           int last_modYmu=0; int last_modYe=0; // = qy+my*(193.7693);
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

           if(SigTracks->pdgCode()==11 and last_modXe==2 and last_modYe==2 and stereo_e>1){ //and abs(SigTracks->startX())<4 and abs(SigTracks->startY())<4 and abs(last_modX)<4 and abs(last_modY)<4){ 
		pe.SetXYZ(SigTracks->px(),SigTracks->py(),SigTracks->pz()); the_gen=pe.Theta();
		Ee=SigTracks->energy();
		double mx=SigTracks->ax();
		double my=SigTracks->ay();
		the=sqrt(mx*mx+my*my);
		pe_dir=pe.Unit();    the_sdr=acos(pmuin_dir.Dot(pe_dir));
		yes_e=1;}//if(the_sdr<0.035)yes_e=1;}
 	   if(SigTracks->pdgCode()==-13 and last_modXmu==2 and last_modYmu==2 and stereo_mu>1){//and abs(SigTracks->startX())<4 and abs(SigTracks->startY())<4 and abs(last_modX)<4 and abs(last_modY)<4){ 
		pmu.SetXYZ(SigTracks->px(),SigTracks->py(),SigTracks->pz());
		thmu_gen=pmu.Theta();double mx=SigTracks->ax();
		double my=SigTracks->ay();
		thmu=sqrt(mx*mx+my*my);
                pmu_dir=pmu.Unit();  thmu_sdr=acos(pmuin_dir.Dot(pmu_dir));
		yes_mu=1;}//if(thmu_sdr>0.0001)yes_mu=1;}
	 }
	}

 if(yes_e!=1 or yes_mu!=1) cout << "NOT RECONSTRUCTIBLE" << endl;
 if(yes_e==1 and yes_mu==1){

cout << "RECONSTRUCTIBLE" << endl;

std::array<std::vector<double>,6> position_mu;//={{{0.,0.},{0.,0.},{0.,0.},{0.,0.},{0.,0.},{0.,0.}}};
std::array<std::vector<double>,6> position_e;


double opening = acos(pe_dir.Dot(pmu_dir));
	   signal+=MesmerEvent->wgt_full;
vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
vector<MUonERecoOutputVertex> vrtx = ReconstructionOutput->reconstructedVertices();

int sig=0;
for(int j=0; j<vrtx.size();j++)
{
if(vrtx.at(j).stationIndex()==1) //and tracks.at(0).processIDofLinkedTrack()==45 and tracks.at(0).linkedTrackID()!=tracks.at(1).linkedTrackID()){
{
 MUonERecoOutputTrack mu_in = vrtx.at(j).incomingMuon();
 MUonERecoOutputTrack mu_out = vrtx.at(j).outgoingMuon();
 MUonERecoOutputTrack e_out = vrtx.at(j).outgoingElectron();

if(vrtx.size()>1){h_chidofM->Fill(vrtx.at(j).chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);

		  double distancex = trackatZ(mu_out.x0(),mu_out.xSlope(),0.) - trackatZ(e_out.x0(),e_out.xSlope(),0.) ;
		  double distancey = trackatZ(mu_out.y0(),mu_out.ySlope(),0.) - trackatZ(e_out.y0(),e_out.ySlope(),0.) ;
		  double d = sqrt(distancex*distancex + distancey*distancey);
			cout << "      vertex " << j << ", chi^2 " << vrtx.at(j).chi2perDegreeOfFreedom() << endl;
		  cout << "total distance " << d << endl;
		}
if(vrtx.size()==1)h_chidof2->Fill(vrtx.at(j).chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);

if(mu_out.processIDofLinkedTrack()==45 and e_out.processIDofLinkedTrack()==45) yes2++;

std::vector<MUonERecoOutputTrackHit> hits_mu=mu_out.hits();

for(int h=0;h<hits_mu.size();h++){
position_mu.at(hits_mu.at(h).moduleID()).push_back(hits_mu.at(h).position());
		}

std::vector<MUonERecoOutputTrackHit> hits_e=e_out.hits();

for(int h=0;h<hits_e.size();h++){
position_e.at(hits_e.at(h).moduleID()).push_back(hits_e.at(h).position());
                }

inter->Fill(mu_out.processIDofLinkedTrack(),MesmerEvent->wgt_full);
inter->Fill(e_out.processIDofLinkedTrack(),MesmerEvent->wgt_full);
	}
}

std::vector<double> x_e;
std::vector<double> x_mu;
std::vector<double> y_e;
std::vector<double> y_mu;
std::array<std::vector<double>,6> position;

for(int j=0; j<tracks.size();j++)
{ 

if(tracks.at(j).processIDofLinkedTrack()==45 and tracks.at(j).sector()==1) TrackIdreco=tracks.at(j).linkedTrackID();
if(tracks.at(j).processIDofLinkedTrack()==45 and tracks.size()>=3 and tracks.at(j).sector()==1 and tracks.at(j).percentageOfHitsSharedWithLinkedTrack()>=100) //and tracks.at(0).processIDofLinkedTrack()==45 and tracks.at(0).linkedTrackID()!=tracks.at(1).linkedTrackID()){
{cout << "tracks.size " << tracks.size() << endl;

if(code_e==tracks.at(j).linkedTrackID()) {x_e.push_back(trackatZ(tracks.at(j).x0(),tracks.at(j).xSlope(),1031.7)); y_e.push_back(trackatZ(tracks.at(j).y0(),tracks.at(j).ySlope(),1031.7));}
if(code_mu==tracks.at(j).linkedTrackID()) {x_mu.push_back(trackatZ(tracks.at(j).x0(),tracks.at(j).xSlope(),1031.7)); y_mu.push_back(trackatZ(tracks.at(j).y0(),tracks.at(j).ySlope(),1031.7));}


std::vector<MUonERecoOutputTrackHit> hits_=tracks.at(j).hits();

for(int h=0;h<hits_.size();h++){
position.at(hits_.at(h).moduleID()).push_back(hits_.at(h).position());
                }
        }
}

if(x_e.size()!=0 and y_e.size()!=0 and x_mu.size()!=0 and y_mu.size()!=0){
for(int e=0; e<x_e.size(); e++){
	for(int m=0; m<x_mu.size(); m++){

		double distancex = x_e.at(e) - x_mu.at(m);
		double distancey = y_e.at(e) - y_mu.at(m);
                double d = sqrt(distancex*distancex + distancey*distancey);
		cout << "-------------------" << endl;
		cout << "electron " << e << ", muon " << m << ", trks total distance at Z=trgt " << d << endl;

		}
	}
}
if(yes2>=1){reco+=MesmerEvent->wgt_full;
//        h_angen_el->Fill(the_sdr,Ee,MesmerEvent->wgt_full);
        tracksize->Fill(vrtx.size(),MesmerEvent->wgt_full);
	h_trackerStubs0->Fill(tracks.size(),MesmerEvent->wgt_full);
	h_op->Fill(opening,MesmerEvent->wgt_full);
	}

if(yes2>1){more_reco+=MesmerEvent->wgt_full;
cout << " event " << i << endl;
        tracksizeM->Fill(vrtx.size(),MesmerEvent->wgt_full);

        for(int j=0; j<tracks.size();j++)
        {inter1->Fill(tracks.at(j).processIDofLinkedTrack(),MesmerEvent->wgt_full);}

	}

if(yes2<1){reco0+=MesmerEvent->wgt_full;
cout <<"NOT RECONSTRUCTED"<<endl;
	tracksize0->Fill(vrtx.size(),MesmerEvent->wgt_full);
	h_anglee0->Fill(the_sdr,MesmerEvent->wgt_full); h_pte0->Fill(point_el,MesmerEvent->wgt_full);
	h_angle1->Fill(the_sdr,thmu_sdr,MesmerEvent->wgt_full);
        h_anglemu0->Fill(thmu_sdr,MesmerEvent->wgt_full); h_ptmu0->Fill(point_mu,MesmerEvent->wgt_full);

	h_op1->Fill(opening,MesmerEvent->wgt_full);
	//if(tracks.size()==0) h_trackerStubs0->Fill(TrackerStubs->GetEntries(),MesmerEvent->wgt_full);

	for(int j=0; j<tracks.size();j++)
	{inter0->Fill(tracks.at(j).processIDofLinkedTrack(),MesmerEvent->wgt_full);}
	}

}

cout << "---------------------"<<endl;
yes_e=0;yes_mu=0;
code_e=-99;code_mu=-99;
point_el=0;point_mu=0;
TrackIdreco=-99;
yes2=0;
link.clear();
} //end of general for

double ratio =reco/signal;
double ratio1 =reco1/signal;
double ratio0 =reco0/signal;
double ratioM =more_reco/signal;

cout << "Su " << signal << " eventi di segnale, " << reco << " sono ricostruiti, con un rapporto del " << ratio*100 << "%"<< endl;
cout << "Su " << signal << " eventi di segnale con piu' di un vertice, " << more_reco << " sono ricostruiti, con un rapporto del " << ratioM*100 << "%"<< endl;
cout << "Su " << signal << " eventi di segnale con 0 vertici di segnale reco, " << reco0 << ", con un rapporto del " << ratio0*100 << "%"<< endl;



/*
TCanvas o("o","o",700,700);
h_op->Draw("hist");
o.SaveAs("op.pdf");

TCanvas b("b","b",700,700);
h_shared->Draw("hist");
b.SaveAs("shared.pdf");

TCanvas a("pa","p",700,700);
gPad->SetLogy();
h_chidof2->Draw("hist");
h_chidofM->SetLineColor(kRed);
h_chidofM->Draw("hist same");
a.SaveAs("chi.pdf");

TCanvas a3("p3","p3",700,700);
a3.Divide(2,2);
a3.cd(1);
gPad->SetLogy();
h_pte0->SetLineColor(kRed);
h_pte0->Draw("hist");
h_ptmu0->Draw("hist same");
a3.cd(2);
gPad->SetLogy();
h_pte1->SetLineColor(kRed);
h_pte1->Draw("hist");
h_ptmu1->Draw("hist same");
a3.cd(3);
gPad->SetLogy();
h_trackerStubs0->Draw("hist");
a3.cd(4);
gPad->SetLogy();
h_trackerStubs1->Draw("hist");
a3.SaveAs("detector_pt.pdf");

TCanvas a4("size","size",700,700);
a4.Divide(2,4);
a4.cd(1);
tracksize->Draw("hist");
a4.cd(2);
inter->Draw("hist");
a4.cd(3);
tracksizeM->Draw("hist");
a4.cd(4);
interM->Draw("hist");
a4.cd(5);
tracksize0->Draw("hist");
a4.cd(6);
inter0->Draw("hist");
a4.cd(7);
tracksize1->Draw("hist");
a4.cd(8);
inter1->Draw("hist");
a4.SaveAs("trackint.pdf");

TCanvas a1("p","p",700,700);
h_part->Draw("hist");
a1.SaveAs("hpart.pdf");

Int_t nxTh = h_angen_el->GetNbinsX();
Int_t nyTh = h_angen_el->GetNbinsY();
for (Int_t i=1; i<nxTh+1; i++) {
for (Int_t j=1; j<nyTh+1; j++) {
if (h_angen_el->GetBinContent(i,j)<1) h_angen_el->SetBinContent(i,j,0);}}

Int_t nxTh1 = h_angle1->GetNbinsX();
Int_t nyTh1 = h_angle1->GetNbinsY();
for (Int_t i=1; i<nxTh1+1; i++) {
for (Int_t j=1; j<nyTh1+1; j++) {
if (h_angle1->GetBinContent(i,j)<1) h_angle1->SetBinContent(i,j,0);}}


TCanvas a2("a2","a2",1400,1400);
a2.Divide(2,4);
a2.cd(1);
h_anglee0->Draw("hist");
a2.cd(2);
h_anglee1->Draw("hist");
a2.cd(3);
h_anglemu0->Draw("hist");
a2.cd(4);
h_anglemu1->Draw("hist");
a2.cd(5);
h_angen_el->Draw("COLZ");
a2.cd(6);
h_angle1->Draw("COLZ");
a2.cd(7);
h_op0->Draw("hist");
a2.cd(8);
h_op1->Draw("hist");
a2.SaveAs("angle.pdf");

//-----------------
Int_t nxTh = h_angle->GetNbinsX();
Int_t nyTh = h_angle->GetNbinsY();
for (Int_t i=1; i<nxTh+1; i++) {
for (Int_t j=1; j<nyTh+1; j++) {
if (h_angle->GetBinContent(i,j)<1) h_angle->SetBinContent(i,j,0);}}

Int_t nxTh1 = h_angle1->GetNbinsX();
Int_t nyTh1 = h_angle1->GetNbinsY();
for (Int_t i=1; i<nxTh1+1; i++) {
for (Int_t j=1; j<nyTh1+1; j++) {
if (h_angle1->GetBinContent(i,j)<1) h_angle1->SetBinContent(i,j,0);}}

Int_t nxTh2 = h_angle3->GetNbinsX();
Int_t nyTh2 = h_angle3->GetNbinsY();
for (Int_t i=1; i<nxTh2+1; i++) {
for (Int_t j=1; j<nyTh2+1; j++) {
if (h_angle3->GetBinContent(i,j)<1) h_angle3->SetBinContent(i,j,0);}}

Int_t nxTh3 = h_angle_noreco->GetNbinsX();
Int_t nyTh3 = h_angle_noreco->GetNbinsY();
for (Int_t i=1; i<nxTh3+1; i++) {
for (Int_t j=1; j<nyTh3+1; j++) {
if (h_angle_noreco->GetBinContent(i,j)<1) h_angle_noreco->SetBinContent(i,j,0);}}

Int_t nxTh4 = h_angle_noreco0->GetNbinsX();
Int_t nyTh4 = h_angle_noreco0->GetNbinsY();
for (Int_t i=1; i<nxTh4+1; i++) {
for (Int_t j=1; j<nyTh4+1; j++) {
if (h_angle_noreco0->GetBinContent(i,j)<1) h_angle_noreco0->SetBinContent(i,j,0);}}

Int_t nxTh5 = h_angle_noreco1->GetNbinsX();
Int_t nyTh5 = h_angle_noreco1->GetNbinsY();
for (Int_t i=1; i<nxTh5+1; i++) {
for (Int_t j=1; j<nyTh5+1; j++) {
if (h_angle_noreco1->GetBinContent(i,j)<1) h_angle_noreco1->SetBinContent(i,j,0);}}

h_anglee->Sumw2();
h_anglee_noreco->Sumw2();
h_anglemu->Sumw2();
h_anglemu_noreco->Sumw2();

TH1D *h_ratio_e = (TH1D*)h_anglee_noreco->Clone("h_ratio_e");
h_ratio_e->Divide(h_anglee_noreco,h_anglee);

TH1D *h_ratio_mu = (TH1D*)h_anglemu_noreco->Clone("h_ratio_mu");
h_ratio_mu->Divide(h_anglemu_noreco,h_anglemu);

TCanvas a1("p","p",700,700);
a1.Divide(2,3);
a1.cd(1);
h_energy->Fit("gaus");
h_energy->Draw();
gStyle->SetOptFit(1011);
a1.cd(2);
h_mx->Fit("gaus");
h_mx->Draw();
gStyle->SetOptFit(1011);
a1.cd(3);
h_my->Fit("gaus");
h_my->Draw();
gStyle->SetOptFit(1011);
a1.cd(4);
h_x->Fit("gaus");
h_x->Draw();
gStyle->SetOptFit(1011);
a1.cd(5);
h_y->Fit("gaus");
h_y->Draw();
gStyle->SetOptFit(1011);
a1.cd(6);
h_xy->Draw("COLZ");

a1.SaveAs("propr_narrow.pdf");

TCanvas a2("a2","a2",1400,1400);
a2.Divide(2,3);
a2.cd(1);
//h_angle->Draw("COLZ");
h_anglemu->SetLineColor(kRed);
h_anglemu->Draw("hist");
h_anglemu_noreco->Draw("hist same");
a2.cd(2);
//h_angle1->Draw("COLZ");
h_anglee->SetLineColor(kRed);
h_anglee->Draw("hist");
h_anglee_noreco->Draw("hist same");
a2.cd(5);
h_angle3->Draw("COLZ");
a2.cd(6);
h_angle_noreco->Draw("COLZ");
a2.cd(3);
h_ratio_e->Draw();
a2.cd(4);
h_ratio_mu->Draw();
a2.SaveAs("angles.pdf");

TCanvas a3("p3","p3",700,700);
a3.Divide(1,2);
a3.cd(1);
gPad->SetLogy();
h_pte_noreco->SetLineColor(kRed);
h_pte_noreco->Draw("hist");
h_ptmu_noreco->Draw("hist same");
a3.cd(2);
gPad->SetLogy();
h_pte->SetLineColor(kRed);
h_pte->Draw("hist");
h_ptmu->Draw("hist same");
a3.SaveAs("detector_pt.pdf");

TCanvas a4("size","size",700,700);
a4.Divide(1,3);
a4.cd(1);
tracksize_noreco->Draw("hist");
a4.cd(2);
tracksize->Draw("hist");
a4.cd(3);
inter->Draw("hist");
a4.SaveAs("tracksize.pdf");

TCanvas a5("a5","a5",1400,1400);
a5.Divide(1,3);
a5.cd(1);
h_anglemu_noreco0->Draw("hist");
a5.cd(2);
h_anglee_noreco0->Draw("hist");
a5.cd(3);
h_angle_noreco0->Draw("hist");
a5.SaveAs("angles_0trk.pdf");

TCanvas a6("a6","a6",1400,1400);
gPad->SetLogy();
h_trackerStubs_noreco0->SetLineColor(kGreen);
h_trackerStubs_noreco1->SetLineColor(kRed);
h_trackerStubs->Draw("hist");
h_trackerStubs_noreco0->Draw("hist same");
h_trackerStubs_noreco1->Draw("hist same");
a6.BuildLegend();
a6.SaveAs("trackerStubs.pdf");

TCanvas a7("a7","a7",1400,1400);
a7.Divide(1,3);
a7.cd(1);
h_anglemu_noreco1->Draw("hist");
a7.cd(2);
h_anglee_noreco1->Draw("hist");
a7.cd(3);
h_angle_noreco1->Draw("hist");
a7.SaveAs("angles_1trk.pdf");
*/
}


