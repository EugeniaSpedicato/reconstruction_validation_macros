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

//TFile *inputfile = new TFile("TRMesmer_box_offset_100k_2GeV_2.root");//TRMesmer_box_offset_100k_2GeV_2.root");//TRMesmer_boxdiv_offset_2.root");

TChain * cbmsim = new TChain("cbmsim");
cbmsim->Add("energy_files/TRMesmer_boxdiv_offset_30k_02-2GeV_2.root");
cbmsim->Add("energy_files/TRMesmer_boxdiv_offset_30k_100-110GeV_2.root");
cbmsim->Add("energy_files/TRMesmer_boxdiv_offset_30k_120-130GeV_2.root");
cbmsim->Add("energy_files/TRMesmer_boxdiv_offset_30k_130-140GeV_2.root");
cbmsim->Add("energy_files/TRMesmer_boxdiv_offset_30k_20-30GeV_2.root");
cbmsim->Add("energy_files/TRMesmer_boxdiv_offset_30k_2-20GeV_2.root");
cbmsim->Add("energy_files/TRMesmer_boxdiv_offset_30k_30-40GeV_2.root");
cbmsim->Add("energy_files/TRMesmer_boxdiv_offset_30k_40-50GeV_2.root");
cbmsim->Add("energy_files/TRMesmer_boxdiv_offset_30k_50-60GeV_2.root");
cbmsim->Add("energy_files/TRMesmer_boxdiv_offset_30k_60-70GeV_2.root");
cbmsim->Add("energy_files/TRMesmer_boxdiv_offset_30k_70-80GeV_2.root");
cbmsim->Add("energy_files/TRMesmer_boxdiv_offset_30k_80-90GeV_2.root");
cbmsim->Add("energy_files/TRMesmer_boxdiv_offset_30k_90-100GeV_2.root");



        //TTree* cbmsim = (TTree*) inputfile->Get("cbmsim");

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
std::array<double,14> sig;

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

signal+=MesmerEvent->wgt_full;

if(Ee>=0.2 and Ee<=2) reco1+=MesmerEvent->wgt_full;
else if(Ee>=2 and Ee<=20) sig.at(1)+=MesmerEvent->wgt_full;
else if(Ee>=20 and Ee<=30) sig.at(2)+=MesmerEvent->wgt_full;
else if(Ee>=30 and Ee<=40) sig.at(3)+=MesmerEvent->wgt_full;
else if(Ee>=40 and Ee<=50) sig.at(4)+=MesmerEvent->wgt_full;
else if(Ee>=50 and Ee<=60) sig.at(5)+=MesmerEvent->wgt_full;
else if(Ee>=60 and Ee<=70) sig.at(6)+=MesmerEvent->wgt_full;
else if(Ee>=70 and Ee<=80) sig.at(7)+=MesmerEvent->wgt_full;
else if(Ee>=80 and Ee<=90) sig.at(8)+=MesmerEvent->wgt_full;
else if(Ee>=90 and Ee<=100) sig.at(9)+=MesmerEvent->wgt_full;
else if(Ee>=100 and Ee<=110) sig.at(10)+=MesmerEvent->wgt_full;
else if(Ee>=110 and Ee<=120) sig.at(11)+=MesmerEvent->wgt_full;
else if(Ee>=120 and Ee<=130) sig.at(12)+=MesmerEvent->wgt_full;
else if(Ee>=130 and Ee<=140) sig.at(13)+=MesmerEvent->wgt_full;

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

int yes_mu=0;
int yes_e=0;

TVector3 in;
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
{yes2++; 

 int sum = tracks.at(j).numberOfXProjectionHits() + tracks.at(j).numberOfYProjectionHits() + tracks.at(j).numberOfStereoHits();

                 if(code_e==tracks.at(j).linkedTrackID()) { yes_e++;
TVector3 pe_rec;
 pe_rec.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.);
 pe_rec= pe_rec.Unit();
 the_rec_vec.push_back(acos(in.Dot(pe_rec)));
 chi_min_e.push_back(tracks.at(j).chi2perDegreeOfFreedom());
                                        }

                 if(code_mu==tracks.at(j).linkedTrackID()) { yes_mu++;

TVector3 pmu_rec;
 pmu_rec.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.);
 pmu_rec= pmu_rec.Unit();
 thmu_rec_vec.push_back(acos(in.Dot(pmu_rec)));
 chi_min_mu.push_back(tracks.at(j).chi2perDegreeOfFreedom());
                                        }


std::vector<MUonERecoOutputTrackHit> hits_=tracks.at(j).hits();

	}
}

if( thmu_rec_vec.size()!=0){
auto it = min_element(chi_min_mu.begin(),chi_min_mu.end()); thmu_rec = thmu_rec_vec.at(std::distance(chi_min_mu.begin(), it));
        }
if( the_rec_vec.size()!=0){
auto it = min_element(chi_min_e.begin(),chi_min_e.end()); the_rec = the_rec_vec.at(std::distance(chi_min_e.begin(), it));
        }

if(yes_e>=1 and yes_mu>=1 and thmu_rec>0.0002 and chi<100){ 
reco+=MesmerEvent->wgt_full; //chi<100 and thmu_rec>0.0002
double res_mu = thmu_rec-thmu_sdr;
double res_e = the_rec-the_sdr;


if(Ee>=0.2 and Ee<=2) more_reco+=MesmerEvent->wgt_full;
else if(Ee>=2 and Ee<=20) eff.at(1)+=MesmerEvent->wgt_full;
else if(Ee>=20 and Ee<=30) eff.at(2)+=MesmerEvent->wgt_full;
else if(Ee>=30 and Ee<=40) eff.at(3)+=MesmerEvent->wgt_full;
else if(Ee>=40 and Ee<=50) eff.at(4)+=MesmerEvent->wgt_full;
else if(Ee>=50 and Ee<=60) eff.at(5)+=MesmerEvent->wgt_full;
else if(Ee>=60 and Ee<=70) eff.at(6)+=MesmerEvent->wgt_full;
else if(Ee>=70 and Ee<=80) eff.at(7)+=MesmerEvent->wgt_full;
else if(Ee>=80 and Ee<=90) eff.at(8)+=MesmerEvent->wgt_full;
else if(Ee>=90 and Ee<=100) eff.at(9)+=MesmerEvent->wgt_full;
else if(Ee>=100 and Ee<=110) eff.at(10)+=MesmerEvent->wgt_full;
else if(Ee>=110 and Ee<=120) eff.at(11)+=MesmerEvent->wgt_full;
else if(Ee>=120 and Ee<=130) eff.at(12)+=MesmerEvent->wgt_full;
else if(Ee>=130 and Ee<=140) eff.at(13)+=MesmerEvent->wgt_full;
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

cout << sig.at(0) << endl;
cout << eff.at(0) << endl;

cout << "Su " << signal << " eventi di segnale con, energia tra 0.2 e 140  GeV, "
<< reco << " sono ricostruiti, con un rapporto del " << (reco/signal)*100 << "%"<< endl;


cout << "Su " << reco1 << " eventi di segnale con, energia tra" << 0.2 << " e " << 2 << " GeV, "
<< more_reco << " sono ricostruiti, con un rapporto del " << (more_reco/reco1)*100 << "%"<< endl;

for(int i=1; i<14; i++){
cout << "Su " << sig.at(i) << " eventi di segnale con, energia tra" << i*10 << " e " << (i+1)*10 << " GeV, "
<< eff.at(i) << " sono ricostruiti, con un rapporto del " << (eff.at(i)/sig.at(i))*100 << "%"<< endl;
 }

}


