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

 	TFile *inputfile = new TFile("TRMesmer_box_offset/TRMesmer_boxdiv_offset_2.root");
//      TFile *inputfile = new TFile("TRMesmer_box_nobend_100k_2.root");
//TRMesmer_box_nobend_parallel_20k_2GeV_0.root");

int fn=2;
//TRMesmer_box_offset_100k_2GeV_2.root");//TRMesmer_boxdiv_offset_2.root");

//TChain * cbmsim = new TChain("cbmsim");
//cbmsim->Add("TRPP_1k.root");
//cbmsim->Add("TRMesmer_1k_unw.root");
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

TH2D *h_res=new TH2D("res", "(thmu_rec-the_true) VS energy Ee<5 GeV",50,0,5,400,-0.2,0.2);
TH2D *h_res1=new TH2D("res1", "(thmu_rec-the_true) VS energy 5<E<15 GeV",20,5,15,400,-0.2,0.2);
TH2D *h_res2=new TH2D("res2", "(thmu_rec-the_true) VS energy 15<E<25 GeV",20,15,25,400,-0.2,0.2);
TH2D *h_res3=new TH2D("res3", "(thmu_rec-the_true) VS energy Ee>25 GeV",6,25,65,400,-0.2,0.2);

TH1D *h_chidofM=new TH1D("h_chidofM", "chi2 per dof when I have MORE THAN 2 signal tracks",101,0,101);
TH1D *h_chidof2=new TH1D("h_chidof2", "chi2 per dof when I have only 2 signal tracks",101,0,101);

TH1D *h_shared=new TH1D("h_shared","multiplicity for well reconstructed elastic events",20,0,20);
TH1D *h_op = new TH1D("op","opening angle when reco",400,0,0.04);

TH1D *h_vrtx=new TH1D("vrt","number of vertices when ghosts tracks",12,0,12);

TH1D *h_quality=new TH1D("q","Quality of the well reconstructed events",101,0,101);
TH1D *h_quality_more=new TH1D("qmore","Quality of the well reconstructed events with ghosts",101,0,101);
TH1D *vrtx_chi=new TH1D("chie","Chi2 per DOF of the kinematic vrtx for signal",40,0,200);//300,0,6000);
TH1D *vrtx_chi_after=new TH1D("chie_after","Chi2 per DOF of the kinematic vrtx for signal after selection",40,0,200);

TH2D *rid2D=new TH2D("emu1","electron-muon angle cut by th_mu and chi2 cuts",400,0,0.04,143,0,0.005);
TH2D *pre2D=new TH2D("preemu1","electron-muon angle before  th_mu and chi2 cuts",400,0,0.04,143,0,0.005);
TH2D *post2D=new TH2D("postemu1","electron-muon angle postth_mu and chi2 cuts",400,0,0.04,143,0,0.005);

TH1D* h_aco=new TH1D("aco","Acoplanarity of muone+electron from sig mesmer",600,-3.14,3.14);
TH1D* h_aco_02=new TH1D("aco02","Acoplanarity of muone+electron from sig mesmer after th_mu>0.2mrad",600,-3.14,3.14);


//TH2D *h_eff=new TH2D("h_eff","efficiency VS electron energy",15,0,150,200,0,1);

double the_gen=0; double thmu_gen=0; double thmu=0; double the=0;
double thmu_rec=0;double the_rec=0;
double signal=0.; double reco=0.; double reco1=0.; double more_reco=0.; double reco0=0.;double reco3=0.;
double reco_v=0.; double more_reco_v=0.; double reco0_v=0.;
int link0; int link1=0;
vector<int> link;
int yes_e_g=0;int yes_mu_g=0; int yes2=0; int yes_v=0;
int point_mu=0; int point_el=0;
int code_mu=-99; int code_e=-99;
int TrackIdreco=-99;
std::array<double,14> eff;


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
		yes_e_g=1;}//if(the_sdr<0.035)yes_e=1;}
 	   if(SigTracks->pdgCode()==-13 and last_modXmu==2 and last_modYmu==2 and stereo_mu>1){//and abs(SigTracks->startX())<4 and abs(SigTracks->startY())<4 and abs(last_modX)<4 and abs(last_modY)<4){ 
		pmu.SetXYZ(SigTracks->px(),SigTracks->py(),SigTracks->pz());
		thmu_gen=pmu.Theta();double mx=SigTracks->ax();
		double my=SigTracks->ay();
		thmu=sqrt(mx*mx+my*my);
                pmu_dir=pmu.Unit();  thmu_sdr=acos(pmuin_dir.Dot(pmu_dir));
		yes_mu_g=1;}//if(thmu_sdr>0.0001)yes_mu=1;}
	 }
	}

 if(yes_e_g!=1 or yes_mu_g!=1) cout << "NOT RECONSTRUCTIBLE" << endl;
 if(yes_e_g==1 and yes_mu_g==1){

cout << "RECONSTRUCTIBLE" << endl;

std::array<std::vector<double>,6> position;//={{{0.,0.},{0.,0.},{0.,0.},{0.,0.},{0.,0.},{0.,0.}}};
double opening = acos(pe_dir.Dot(pmu_dir));
	   signal+=MesmerEvent->wgt_full;
vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
vector<MUonERecoOutputVertex> vrtx = ReconstructionOutput->reconstructedVertices();

double chi;

for(int j=0; j<vrtx.size();j++)
{
if(vrtx.at(j).stationIndex()==1) //and tracks.at(0).processIDofLinkedTrack()==45 and tracks.at(0).linkedTrackID()!=tracks.at(1).linkedTrackID()){
{
 MUonERecoOutputTrack mu_in = vrtx.at(j).incomingMuon();
 MUonERecoOutputTrack mu_out = vrtx.at(j).outgoingMuon();
 MUonERecoOutputTrack e_out = vrtx.at(j).outgoingElectron();

//if(mu_out.processIDofLinkedTrack()==45 and e_out.processIDofLinkedTrack()==45) yes_v++;
yes_v++;
if(j==0){chi=vrtx.at(j).chi2perDegreeOfFreedom();}

}
}


int sig=0;
int yes_mu=0;
int yes_e=0;

TVector3 in;
vector<TVector3> p_muout;
vector<TVector3> p_eout;


for(int j=0; j<tracks.size();j++)
{
if(tracks.at(j).processIDofLinkedTrack()==0 and tracks.at(j).sector()==0){
double th_inx=tracks.at(j).xSlope();
double th_iny=tracks.at(j).ySlope();
in.SetXYZ(th_inx,th_iny,1.0);
in=in.Unit();
 }
}
vector<double> thmu_rec_vec,chi_min_mu,the_rec_vec,chi_min_e;
chi_min_mu.reserve(5);thmu_rec_vec.reserve(5);chi_min_e.reserve(5);the_rec_vec.reserve(5);

for(int j=0; j<tracks.size();j++)
{

if(tracks.at(j).processIDofLinkedTrack()==45 and tracks.at(j).sector()==1) TrackIdreco=tracks.at(j).linkedTrackID();
if(tracks.at(j).processIDofLinkedTrack()==45 and tracks.size()>=3 and tracks.at(j).sector()==1 and  tracks.at(j).percentageOfHitsSharedWithLinkedTrack()>=0) //and tracks.at(0).processIDofLinkedTrack()==45 and tracks.at(0).linkedTrackID()!=tracks.at(1).linkedTrackID()){
{yes2++; cout << "tracks.size " << tracks.size() << endl;

 int sum = tracks.at(j).numberOfXProjectionHits() + tracks.at(j).numberOfYProjectionHits() + tracks.at(j).numberOfStereoHits();

                 if(code_e==tracks.at(j).linkedTrackID()) { yes_e++;
TVector3 pe_rec;
 pe_rec.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.);
 pe_rec= pe_rec.Unit();
 the_rec_vec.push_back(acos(in.Dot(pe_rec)));
 chi_min_e.push_back(tracks.at(j).chi2perDegreeOfFreedom());
 p_eout.push_back(pe_rec);
                                        }

                 if(code_mu==tracks.at(j).linkedTrackID()) { yes_mu++;

TVector3 pmu_rec;
 pmu_rec.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.);
 pmu_rec= pmu_rec.Unit();
 thmu_rec_vec.push_back(acos(in.Dot(pmu_rec)));
 chi_min_mu.push_back(tracks.at(j).chi2perDegreeOfFreedom());
 p_muout.push_back(pmu_rec);
                                        }

	}
}
TVector3 p_erec,p_murec;
if( thmu_rec_vec.size()!=0){
auto it = min_element(chi_min_mu.begin(),chi_min_mu.end()); thmu_rec = thmu_rec_vec.at(std::distance(chi_min_mu.begin(), it));
 p_murec= p_muout.at(std::distance(chi_min_mu.begin(), it));
        }
if( the_rec_vec.size()!=0){
auto it = min_element(chi_min_e.begin(),chi_min_e.end()); the_rec = the_rec_vec.at(std::distance(chi_min_e.begin(), it));
 p_erec= p_eout.at(std::distance(chi_min_e.begin(), it));

        }

double acoplanarity;

if(yes_e>=1 and yes_mu>=1){//and thmu_rec>0.0002 and chi>0){

						double dotProduct = p_murec.Dot(p_erec);
                                                TVector3 crossProduct = p_murec.Cross(p_erec);
                                                double T = in.Dot(crossProduct);
                                                TVector3 im= in.Cross(p_murec);
                                                TVector3 ie= in.Cross(p_erec);
                                                T = T>0? 1:(-1);
                                                acoplanarity= T*(TMath::Pi()-acos( ((im).Dot(ie))/(im.Mag()*ie.Mag()) ));

double res_mu = thmu_rec-thmu_sdr;
double res_e = the_rec-the_sdr;
reco+=MesmerEvent->wgt_full;
if(thmu_rec<0.0002 or chi>100) rid2D->Fill(the_rec,thmu_rec,MesmerEvent->wgt_full);
if(thmu_rec>0.0002 and chi<100) post2D->Fill(the_rec,thmu_rec,MesmerEvent->wgt_full);
pre2D->Fill(the_rec,thmu_rec,MesmerEvent->wgt_full);

if(the_rec<0.035){
vrtx_chi->Fill(chi,MesmerEvent->wgt_full);
h_aco->Fill(acoplanarity);
if(thmu_rec>0.0002){vrtx_chi_after->Fill(chi,MesmerEvent->wgt_full);h_aco_02->Fill(acoplanarity);}
}

//if(thmu_rec>0.0002){reco+=MesmerEvent->wgt_full;}
//else rid2D->Fill(the_sdr,thmu_sdr,MesmerEvent->wgt_full);


/*if(E_e>=0.2 and E_e<=2) eff.at(0)+=MesmerEvent->wgt_full;
else if(E_e>=2 and E_e<=20) eff.at(1)+=MesmerEvent->wgt_full;
else if(E_e>=20 and E_e<=30) eff.at(2)+=MesmerEvent->wgt_full;
else if(E_e>=30 and E_e<=40) eff.at(3)+=MesmerEvent->wgt_full;
else if(E_e>=40 and E_e<=50) eff.at(4)+=MesmerEvent->wgt_full;
else if(E_e>=50 and E_e<=60) eff.at(5)+=MesmerEvent->wgt_full;
else if(E_e>=60 and E_e<=70) eff.at(6)+=MesmerEvent->wgt_full;
else if(E_e>=70 and E_e<=80) eff.at(7)+=MesmerEvent->wgt_full;
else if(E_e>=80 and E_e<=90) eff.at(8)+=MesmerEvent->wgt_full;
else if(E_e>=90 and E_e<=100) eff.at(9)+=MesmerEvent->wgt_full;
else if(E_e>=100 and E_e<=110) eff.at(10)+=MesmerEvent->wgt_full;
else if(E_e>=110 and E_e<=120) eff.at(11)+=MesmerEvent->wgt_full;
else if(E_e>=120 and E_e<=130) eff.at(12)+=MesmerEvent->wgt_full;
else if(E_e>=130 and E_e<=140) eff.at(13)+=MesmerEvent->wgt_full;*/


        tracksize->Fill(tracks.size(),MesmerEvent->wgt_full);
	h_trackerStubs0->Fill(TrackerStubs->GetEntries(),MesmerEvent->wgt_full);
	h_op->Fill(opening,MesmerEvent->wgt_full);
        for(int j=0; j<tracks.size();j++)
        {inter->Fill(tracks.at(j).processIDofLinkedTrack(),MesmerEvent->wgt_full);}}

if(yes_e==1 and yes_mu==1 and tracks.size()==3 and thmu_rec>0.0002 and chi>0) {reco3+=MesmerEvent->wgt_full;}
else{h_angle1->Fill(the_sdr,thmu_sdr,MesmerEvent->wgt_full);}

if(yes_e==1 and yes_mu==1 and thmu_rec>0.0002 and chi>0){
h_shared->Fill(tracks.size(),MesmerEvent->wgt_full);

        for(int j=0; j<tracks.size();j++)
        {
	if(tracks.at(j).processIDofLinkedTrack()==45) h_chidof2->Fill(tracks.at(j).chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);
if(tracks.at(j).processIDofLinkedTrack()==45 and tracks.size()>=3 and tracks.at(j).sector()==1 and  tracks.at(j).percentageOfHitsSharedWithLinkedTrack()>=0) 
	{h_quality->Fill(tracks.at(j).percentageOfHitsSharedWithLinkedTrack(),MesmerEvent->wgt_full);}
	}
}

if((yes_e>1 or yes_mu>1) and yes_e!=0 and yes_mu!=0 and thmu_rec>0.0002 and chi>0){more_reco+=MesmerEvent->wgt_full;

        for(int j=0; j<tracks.size();j++)
        {
if(tracks.at(j).processIDofLinkedTrack()==45 and tracks.size()>=3 and tracks.at(j).sector()==1 and  tracks.at(j).percentageOfHitsSharedWithLinkedTrack()>=0) 
        {h_quality_more->Fill(tracks.at(j).percentageOfHitsSharedWithLinkedTrack(),MesmerEvent->wgt_full);}
        }

h_vrtx->Fill(vrtx.size(),MesmerEvent->wgt_full);
cout << "event " << i << endl;
int sig=0;
int sig_1=0;
int sig_2=0;
for(int s=0;s<6;s++)
{sort(position.at(s).begin(), position.at(s).end());
 for(int j=1; j<position.at(s).size();j++)
        {
         if(position.at(s).at(j)==position.at(s).at(j-1))sig=1;
        }
 if(s<4 and sig==1)sig_1=1;
 if(s>=4 and sig==1)sig_2=1; sig=0;
}

//if(sig_1==1 and sig_2==0)more_reco+=MesmerEvent->wgt_full;
//else if(sig_1==0 and sig_2==0)more_reco+=MesmerEvent->wgt_full;

        tracksizeM->Fill(tracks.size(),MesmerEvent->wgt_full);
        for(int j=0; j<tracks.size();j++)
	{//if(tracks.at(j).processIDofLinkedTrack()==45) h_chidofM->Fill(tracks.at(j).chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);
	 interM->Fill(tracks.at(j).processIDofLinkedTrack(),MesmerEvent->wgt_full);}}

if(yes2<2 and TrackIdreco==-99){reco0+=MesmerEvent->wgt_full;
cout <<"NOT RECONSTRUCTED"<<endl;
        if(TrackIdreco==code_e) h_part->Fill(11,MesmerEvent->wgt_full);
        if(TrackIdreco==code_mu) h_part->Fill(-13,MesmerEvent->wgt_full);
	tracksize0->Fill(tracks.size(),MesmerEvent->wgt_full);
	h_anglee0->Fill(the_sdr,MesmerEvent->wgt_full); h_pte0->Fill(point_el,MesmerEvent->wgt_full);
	h_angle1->Fill(the_sdr,thmu_sdr,MesmerEvent->wgt_full);
        h_anglemu0->Fill(thmu_sdr,MesmerEvent->wgt_full); h_ptmu0->Fill(point_mu,MesmerEvent->wgt_full);
        //h_angen_el->Fill(the_sdr,thmu_sdr,MesmerEvent->wgt_full);
h_op1->Fill(opening,MesmerEvent->wgt_full);//h_op0
	//if(tracks.size()==0) h_trackerStubs0->Fill(TrackerStubs->GetEntries(),MesmerEvent->wgt_full);

	for(int j=0; j<tracks.size();j++)
	{inter0->Fill(tracks.at(j).processIDofLinkedTrack(),MesmerEvent->wgt_full);}
	}

if(yes2<2 and TrackIdreco!=-99){reco1+=MesmerEvent->wgt_full;
cout <<"NOT RECONSTRUCTED"<<endl;
	if(TrackIdreco==code_e) h_part->Fill(11,MesmerEvent->wgt_full);
        if(TrackIdreco==code_mu) h_part->Fill(-13,MesmerEvent->wgt_full);
        tracksize1->Fill(tracks.size(),MesmerEvent->wgt_full);
        if(TrackIdreco!=code_e){h_anglee1->Fill(the_sdr,MesmerEvent->wgt_full); h_pte1->Fill(point_el,MesmerEvent->wgt_full);}
        if(TrackIdreco!=code_mu){h_anglemu1->Fill(thmu_sdr,MesmerEvent->wgt_full); h_ptmu1->Fill(point_mu,MesmerEvent->wgt_full);}
        h_angle1->Fill(the_sdr,thmu_sdr,MesmerEvent->wgt_full);
	h_angen_el->Fill(the_sdr,Ee,MesmerEvent->wgt_full);
        if(tracks.size()==1) h_trackerStubs1->Fill(TrackerStubs->GetEntries(),MesmerEvent->wgt_full);
        h_op1->Fill(opening,MesmerEvent->wgt_full);

        for(int j=0; j<tracks.size();j++)
        {inter1->Fill(tracks.at(j).processIDofLinkedTrack(),MesmerEvent->wgt_full);}
	}

if(yes_v>=1 and yes_e>=1 and yes_mu>=1 and thmu_rec>0.0002){reco_v+=MesmerEvent->wgt_full;
        }

if(yes_v>1 and yes_e>=1 and yes_mu>=1 and thmu_rec>0.0002){more_reco_v+=MesmerEvent->wgt_full;
        }

if(yes_v<1 and yes_e>=1 and yes_mu>=1 and thmu_rec>0.0002){reco0_v+=MesmerEvent->wgt_full;
cout <<"NOT RECONSTRUCTED vertex"<<endl;
        }



}
cout << "---------------------"<<endl;
yes_e_g=0;yes_mu_g=0;
code_e=-99;code_mu=-99;
point_el=0;point_mu=0;
TrackIdreco=-99;
yes2=0;yes_v=0;
link.clear();
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

double ratio_v =reco_v/reco;
double ratio0_v =reco0_v/reco;
double ratioM_v =more_reco_v/reco;

cout << "Su " << signal << " eventi di segnale, " << reco_v << " sono ricostruiti con almeno 1 vertice, con un rapporto del " << ratio_v*100 << "%"<< endl;
cout << "Su " << signal << " eventi di segnale, " << more_reco_v << " sono ricostruiti con piu' di unun vertice, con un rapporto del " << ratioM_v*100 << "%"<< endl;
cout << "Su " << signal << " eventi di segnale, " << reco0_v << " hanno 0 vertici, con un rapporto del " << ratio0_v*100 << "%"<< endl;

/*Int_t nxTh1 = rid2D->GetNbinsX();
Int_t nyTh1 = rid2D->GetNbinsY();
for (Int_t i=1; i<nxTh1+1; i++) {
for (Int_t j=1; j<nyTh1+1; j++) {
if (rid2D->GetBinContent(i,j)<1) rid2D->SetBinContent(i,j,0);}}

Int_t nxTh2 = pre2D->GetNbinsX();
Int_t nyTh2 = pre2D->GetNbinsY();
for (Int_t i=1; i<nxTh2+1; i++) {
for (Int_t j=1; j<nyTh2+1; j++) {
if (pre2D->GetBinContent(i,j)<1) pre2D->SetBinContent(i,j,0);}}

Int_t nxTh3 = post2D->GetNbinsX();
Int_t nyTh3 = post2D->GetNbinsY();
for (Int_t i=1; i<nxTh3+1; i++) {
for (Int_t j=1; j<nyTh3+1; j++) {
if (post2D->GetBinContent(i,j)<1) post2D->SetBinContent(i,j,0);}}

TCanvas t("t","theta",700,700);
t.Divide(1,3);
t.cd(1);
pre2D->Draw("COLZ");
t.cd(2);
post2D->Draw("COLZ");
t.cd(3);
rid2D->Draw("COLZ");
t.SaveAs("2d.pdf");

TCanvas b7("b7","b7",700,700);
b7.Divide(2,2);
b7.cd(1);
gPad->SetLogy();
int nBins1 = vrtx_chi->GetNbinsX()+1;
TH1F *h1s = (TH1F*)(vrtx_chi->Clone("h1s"));
h1s->Scale(1/h1s->Integral(0,nBins1));
h1s->Draw("hist");
gPad->SetLogy();
gPad->SetLogy();
b7.cd(2);
int nBins2 = vrtx_chi_after->GetNbinsX()+1;
TH1F *h2s = (TH1F*)(vrtx_chi_after->Clone("h2s"));
h2s->Scale(1/h2s->Integral(0,nBins2));
h2s->Draw("hist");
gPad->SetLogy();
b7.cd(3);
    int nx    = h1s->GetNbinsX()+1;
    double *xbins= new double[nx+1];
    for (int i=0;i<nx;i++)
        xbins[i]=h1s->GetBinLowEdge(i+1);
    xbins[nx]=xbins[nx-1]+h1s->GetBinWidth(nx);
    //book a temporary histogram having extra bins for overflows
    TH1F *htmp_s = new TH1F("htmp_s", h1s->GetTitle(), nx, xbins);
    htmp_s->Sumw2();
    //fill the new histogram including the overflows
    for (int i=1; i<=nx; i++) {
        htmp_s->SetBinContent(htmp_s->FindBin(htmp_s->GetBinCenter(i)),h1s->GetBinContent(i));
        htmp_s->SetBinError(htmp_s->FindBin(htmp_s->GetBinCenter(i)),h1s->GetBinError(i));
    }
    htmp_s->SetBinContent(htmp_s->FindBin(h1s->GetBinLowEdge(1)-1), h1s->GetBinContent(0));
    htmp_s->SetBinError(htmp_s->FindBin(h1s->GetBinLowEdge(1)-1), h1s->GetBinError(0));
    // Restore the number of entries
    htmp_s->SetEntries(h1s->GetEffectiveEntries());

htmp_s->Draw("hist");
gPad->SetLogy();

b7.cd(4);
    int nx1    = h2s->GetNbinsX()+1;
    double *xbins1= new double[nx1+1];
    for (int i=0;i<nx1;i++)
        xbins1[i]=h1s->GetBinLowEdge(i+1);
    xbins1[nx1]=xbins1[nx1-1]+h1s->GetBinWidth(nx1);
    //book a temporary histogram having extra bins for overflows
    TH1F *htmp_s2 = new TH1F(h2s->GetName(), h2s->GetTitle(), nx1, xbins1);
    htmp_s2->Sumw2();
    //fill the new histogram including the overflows
    for (int i=1; i<=nx1; i++) {
        htmp_s2->SetBinContent(htmp_s2->FindBin(htmp_s2->GetBinCenter(i)),h2s->GetBinContent(i));
        htmp_s2->SetBinError(htmp_s2->FindBin(htmp_s2->GetBinCenter(i)),h2s->GetBinError(i));
    }

htmp_s2->Draw("hist");
gPad->SetLogy();

b7.SaveAs("vrtx_chi_sig.pdf");
h1s->SaveAs("vrtx_chi_sig2gev_norm.root");
h2s->SaveAs("vrtx_chi_sig2gev_norm_after02.root");
htmp_s->SaveAs("vrtx_chi_sig2gev_norm_ovrf.root");
htmp_s2->SaveAs("vrtx_chi_sig2gev_norm_after02_ovrf.root");
*/

TCanvas ac("ac","ac",700,700);
ac.Divide(1,2);
ac.cd(1);
TH1F *h_acoN = (TH1F*)(h_aco->Clone("h_acoN"));
h_acoN->Scale(1/h_acoN->Integral());
h_acoN->Draw("hist");
ac.cd(2);
TH1F *h_acoN_02 = (TH1F*)(h_aco_02->Clone("h_acoN_02"));
h_acoN_02->Scale(1/h_acoN_02->Integral());
h_acoN_02->Draw("hist");
if(fn==0){
ac.SaveAs("h_aco_sig_02gev.pdf");
h_acoN->SaveAs("h_aco_sig_02gev.root");
h_acoN_02->SaveAs("h_aco_sig_02gev_after02.root");}
else{
ac.SaveAs("h_aco_sig_2gev.pdf");
h_acoN->SaveAs("h_aco_sig_2gev.root");
h_acoN_02->SaveAs("h_aco_sig_2gev_after02.root");}

/*
TCanvas o1("o1","o1",700,700);
o1.Divide(1,2);
o1.cd(1);
rid2D->Draw("COLZ");
o1.cd(2);
vrtx_chi->Draw("hist");
vrtx_chi_after->SetLineColor(kOrange);
vrtx_chi_after->Draw("hist same");
//vrtx_chi->SetMinimum(0);
gPad->SetLogy();
//vrtx_chi->SetMinimum(0);
o1.SaveAs("vrtx_chi_sig.pdf");
TCanvas o("o","o",700,700);
h_op->Draw("hist");
o.SaveAs("pdf_notReco/op.pdf");

TCanvas b("b","b",700,700);
h_shared->Draw("hist");
gPad->SetLogy();
b.SaveAs("pdf_notReco/shared.pdf");


TCanvas a("pa","p",700,700);
gPad->SetLogy();
h_chidof2->Draw("hist");
h_chidofM->SetLineColor(kRed);
h_chidofM->Draw("hist same");
a.SaveAs("pdf_notReco/chi.pdf");

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
a3.SaveAs("pdf_notReco/detector_pt.pdf");

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
a4.SaveAs("pdf_notReco/trackint.pdf");

TCanvas a1("p","p",700,700);
h_part->Draw("hist");
a1.SaveAs("pdf_notReco/hpart.pdf");

Int_t nxTh = h_angen_el->GetNbinsX();
Int_t nyTh = h_angen_el->GetNbinsY();
for (Int_t i=1; i<nxTh+1; i++) {
for (Int_t j=1; j<nyTh+1; j++) {
if (h_angen_el->GetBinContent(i,j)<1) h_angen_el->SetBinContent(i,j,0);}}
*/
/*Int_t nxTh1 = h_angle1->GetNbinsX();
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
a2.SaveAs("pdf_notReco/angle.pdf");
*/
/*
TCanvas m("m","m",1400,1400);
h_vrtx->Draw("hist");
m.SaveAs("pdf_notReco/vrtx.pdf");


------------------------cut giu
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

a1.SaveAs("pdf_notReco/propr_narrow.pdf");

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
a2.SaveAs("pdf_notReco/angles.pdf");

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
a3.SaveAs("pdf_notReco/detector_pt.pdf");

TCanvas a4("size","size",700,700);
a4.Divide(1,3);
a4.cd(1);
tracksize_noreco->Draw("hist");
a4.cd(2);
tracksize->Draw("hist");
a4.cd(3);
inter->Draw("hist");
a4.SaveAs("pdf_notReco/tracksize.pdf");

TCanvas a5("a5","a5",1400,1400);
a5.Divide(1,3);
a5.cd(1);
h_anglemu_noreco0->Draw("hist");
a5.cd(2);
h_anglee_noreco0->Draw("hist");
a5.cd(3);
h_angle_noreco0->Draw("hist");
a5.SaveAs("pdf_notReco/angles_0trk.pdf");

TCanvas a6("a6","a6",1400,1400);
gPad->SetLogy();
h_trackerStubs_noreco0->SetLineColor(kGreen);
h_trackerStubs_noreco1->SetLineColor(kRed);
h_trackerStubs->Draw("hist");
h_trackerStubs_noreco0->Draw("hist same");
h_trackerStubs_noreco1->Draw("hist same");
a6.BuildLegend();
a6.SaveAs("pdf_notReco/trackerStubs.pdf");

TCanvas a7("a7","a7",1400,1400);
a7.Divide(1,3);
a7.cd(1);
h_anglemu_noreco1->Draw("hist");
a7.cd(2);
h_anglee_noreco1->Draw("hist");
a7.cd(3);
h_angle_noreco1->Draw("hist");
a7.SaveAs("pdf_notReco/angles_1trk.pdf");
*/
}


