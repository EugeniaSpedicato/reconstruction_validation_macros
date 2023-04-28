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

	TFile *inputfile = new TFile("TRMesmer_100k_box.root");//narrowBP_20GeV_1shared.root");//trPoints,20GeV.root");//exampleProductionJob.root");narrowBP.root
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


TH1D *h_anglee0 = new TH1D("e0","polar reconstructed electron angle when event is reco",700,0,0.07);
TH2D *h_anE_el = new TH2D("eth0","reconstructed electron angle-energy when event is reco",700,0,0.07,300,0,60);
TH2D *h_angle_MC_rec = new TH2D("thrmc0","polar electron angle reco(X) VS generated(Y)",700,0,0.07,700,0,0.07);
TH1D *h_En_e = new TH1D("ene0"," generated electron energy when event is reco" ,300,0,60);
TH1D *h_anglee1 = new TH1D("e1","electron angle when event isn't reco",700,0,0.07);

TH1D *h_E1 = new TH1D("ene1","generated electron energy when event is reco with the_el_rec<10" ,450,0,70);
TH1D *h_E2 = new TH1D("ene2","generated electron energy when event is reco with 10<=the_el_rec<30" ,450,0,70);
TH1D *h_E3 = new TH1D("ene3","generated electron energy when event is reco with the_el_rec>=30" ,450,0,70);
TH1D *h_E4 = new TH1D("ene4","generated electron energy when event is NOT reco" ,450,0,70);

TH1D *h_E1MC = new TH1D("ene1MC","generated electron energy when event is reco with the_el_gen<10" ,450,0,70);
TH1D *h_E2MC = new TH1D("ene2MC","generated electron energy when event is reco with 10<=the_el_gen<30" ,450,0,70);
TH1D *h_E3MC = new TH1D("ene3MC","generated electron energy when event is reco with the_el_gen>=30" ,450,0,70);


TH2D *h_angle_MC1 = new TH2D("thmc1","polar angle electron vs muon generated the_el_gen<10",200,0,0.01,100,0,0.005);
TH2D *h_angle_rec1 = new TH2D("thr1","polar angle electron vs muon reconstructed the_el_rec<10",200,0.,0.01,100,0,0.005);

TH2D *h_angle_MC2 = new TH2D("thmc2","polar angle electron vs muon generated 10<=the_el_gen<30",400,0.01,0.03,100,0,0.005);
TH2D *h_angle_rec2 = new TH2D("thr2","polar angle electron vs muon reconstructed 10<=the_el_rec<30",400,0.01,0.03,100,0,0.005);

TH2D *h_angle_MC3 = new TH2D("thmc3","polar angle electron vs muon generated the_el_gen>=30",800,0.03,0.04,100,0,0.005);
TH2D *h_angle_rec3 = new TH2D("thr3","polar angle electron vs muon reconstructed the_el_rec>=30",800,0.03,0.04,100,0,0.005);


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
int no=0;
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
         if(MCTr->interactionID()==5) no=1;
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
                          if(TrackerPt->trackPDGCode()==11 and SigTracks->pdgCode()==11 and TrackerPt->trackID()==code_e){ point_el++; 
                                                                                                 if(TrackerPt->moduleID()==4) last_modXe++; 
                                                                                                 if(TrackerPt->moduleID()==5) last_modYe++;
												 if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_e++;}
                          if(TrackerPt->trackPDGCode()==-13 and SigTracks->pdgCode()==-13 and TrackerPt->trackID()==code_mu){ point_mu++;
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

 if(yes_e!=1 or yes_mu!=1 ) cout << "NOT RECONSTRUCTIBLE" << endl;

double the_rec;
double thmu_rec;

 if(yes_e==1 and yes_mu==1 and no==0){

cout << "RECONSTRUCTIBLE" << endl;

std::array<std::vector<double>,6> position;//={{{0.,0.},{0.,0.},{0.,0.},{0.,0.},{0.,0.},{0.,0.}}};
double opening = acos(pe_dir.Dot(pmu_dir));
	   signal+=MesmerEvent->wgt_full;
vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
int vtx = ReconstructionOutput->adaptiveFitterVerticesMultiplicity();

int sig=0;
vector<double> thmu_rec_v;
thmu_rec_v.reserve(5);
vector<double> the_rec_v;
the_rec_v.reserve(5);
vector<double> chi_min_e;
chi_min_e.reserve(5);
vector<double> chi_min_mu;
chi_min_mu.reserve(5);


TVector3 in; double th_in_rec_tracks;
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

if(tracks.at(j).processIDofLinkedTrack()==45) TrackIdreco=tracks.at(j).linkedTrackID();
if(tracks.at(j).processIDofLinkedTrack()==45 and tracks.size()>=3 and tracks.at(j).sector()==1) //and tracks.at(0).processIDofLinkedTrack()==45 and tracks.at(0).linkedTrackID()!=tracks.at(1).linkedTrackID()){
{yes2++;

        TVector3 out(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0);

if(tracks.at(j).linkedTrackID()==code_e) {the_rec_v.push_back(in.Angle(out)); chi_min_e.push_back(tracks.at(j).chi2perDegreeOfFreedom());}
if(tracks.at(j).linkedTrackID()==code_mu) {thmu_rec_v.push_back(in.Angle(out)); chi_min_mu.push_back(tracks.at(j).chi2perDegreeOfFreedom());}

if(tracks.size()>3 and tracks.at(j).linkedTrackID()==code_e) cout << "linkedID " << tracks.at(j).linkedTrackID() << " VS code_e " << code_e << endl;

std::vector<MUonERecoOutputTrackHit> hits_=tracks.at(j).hits();

for(int h=0;h<hits_.size();h++){
position.at(hits_.at(h).moduleID()).push_back(hits_.at(h).position());
		}
	}
}
if(the_rec_v.size()!=0){ auto it = min_element(chi_min_e.begin(),chi_min_e.end()); the_rec = the_rec_v.at(std::distance(chi_min_e.begin(), it));}
if(thmu_rec_v.size()!=0){ auto it = min_element(chi_min_mu.begin(),chi_min_mu.end()); thmu_rec = thmu_rec_v.at(std::distance(chi_min_mu.begin(), it));}

h_angle_MC1->Fill(the_gen,thmu_gen,MesmerEvent->wgt_full);

if(yes2>=2){reco+=MesmerEvent->wgt_full;
        h_anglee0->Fill(the_rec,MesmerEvent->wgt_full);
        h_anE_el->Fill(the_rec,Ee,MesmerEvent->wgt_full);
	h_angle_MC_rec->Fill(the_rec,the_gen,MesmerEvent->wgt_full);
	h_En_e->Fill(Ee,MesmerEvent->wgt_full);

if(the_rec<0.010) {h_E1->Fill(Ee,MesmerEvent->wgt_full);
		   h_angle_rec1->Fill(the_rec,thmu_rec,MesmerEvent->wgt_full);}
else if(the_rec>=0.010 and the_rec<0.030) {h_E2->Fill(Ee,MesmerEvent->wgt_full);
					   h_angle_rec2->Fill(the_rec,thmu_rec,MesmerEvent->wgt_full);}
else if(the_rec>=0.030 ) {h_E3->Fill(Ee,MesmerEvent->wgt_full);
			  h_angle_rec3->Fill(the_rec,thmu_rec,MesmerEvent->wgt_full);}

if(the_sdr<0.010) {h_E1MC->Fill(Ee,MesmerEvent->wgt_full);}//h_angle_MC1->Fill(the_gen,thmu_gen,MesmerEvent->wgt_full);}
else if(the_sdr>=0.010 and the_sdr<0.030) {h_E2MC->Fill(Ee,MesmerEvent->wgt_full);h_angle_MC2->Fill(the_sdr,thmu_sdr,MesmerEvent->wgt_full);}
else if(the_sdr>=0.030 ) {h_E3MC->Fill(Ee,MesmerEvent->wgt_full);h_angle_MC3->Fill(the_sdr,thmu_sdr,MesmerEvent->wgt_full);}
	}

if(yes2>2){more_reco+=MesmerEvent->wgt_full;}

if(yes2<2){reco1+=MesmerEvent->wgt_full;
cout <<"NOT RECONSTRUCTED event "<< i <<endl;
       if(TrackIdreco!=code_e){h_anglee1->Fill(the_gen,MesmerEvent->wgt_full);
				h_E4->Fill(Ee,MesmerEvent->wgt_full);}
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
cout << "Su " << signal << " eventi di segnale con piu tracce (stesso id), " << more_reco << " sono ricostruiti, con un rapporto del " << ratioM*100 << "%"<< endl;
cout << "Su " << signal << " eventi di segnale con 0 tracce di segnale reco, " << reco0 << ", con un rapporto del " << ratio0*100 << "%"<< endl;
cout << "Su " << signal << " eventi di con 1 sola traccia di segnale reco, " << reco1 << ", con un rapporto del " << ratio1*100 << "%"<< endl;


Int_t nxTh0 = h_anE_el->GetNbinsX();
Int_t nyTh0 = h_anE_el->GetNbinsY();
for (Int_t i=1; i<nxTh0+1; i++) {
for (Int_t j=1; j<nyTh0+1; j++) {
if (h_anE_el->GetBinContent(i,j)<1) h_anE_el->SetBinContent(i,j,0);}}

Int_t nxTh1a = h_angle_MC_rec->GetNbinsX();
Int_t nyTh1a = h_angle_MC_rec->GetNbinsY();
for (Int_t i=1; i<nxTh1a+1; i++) {
for (Int_t j=1; j<nyTh1a+1; j++) {
if (h_angle_MC_rec->GetBinContent(i,j)<1) h_angle_MC_rec->SetBinContent(i,j,0);}}


TCanvas a2("a2","a2",1400,1400);
a2.Divide(2,4);
a2.cd(1);
h_anglee0->Draw("hist");
a2.cd(2);
h_anglee1->Draw("hist");
a2.cd(3);
h_En_e->Draw("hist");
a2.cd(4);
h_angle_MC_rec->Draw("COLZ");
a2.cd(5);
h_anE_el->Draw("COLZ");
a2.SaveAs("pdf_exercise.pdf");

TCanvas a1("p","p",700,700);
a1.Divide(2,4);
a1.cd(1);
h_E1->Draw("hist");
gPad->SetLogy();
a1.cd(2);
h_E1MC->Draw("hist");
gPad->SetLogy();
a1.cd(3);
h_E2->Draw("hist");
gPad->SetLogy();
a1.cd(4);
h_E2MC->Draw("hist");
gPad->SetLogy();
a1.cd(5);
h_E3->Draw("hist");
gPad->SetLogy();
a1.cd(6);
h_E3MC->Draw("hist");
gPad->SetLogy();
a1.cd(7);
h_E4->Draw("hist");
gPad->SetLogy();
a1.SaveAs("different_E.pdf");


Int_t nxTh = h_angle_MC1->GetNbinsX();
Int_t nyTh = h_angle_MC1->GetNbinsY();
for (Int_t i=1; i<nxTh+1; i++) {
for (Int_t j=1; j<nyTh+1; j++) {
if (h_angle_MC1->GetBinContent(i,j)<1) h_angle_MC1->SetBinContent(i,j,0);}}

Int_t nxTh1 = h_angle_rec1->GetNbinsX();
Int_t nyTh1 = h_angle_rec1->GetNbinsY();
for (Int_t i=1; i<nxTh1+1; i++) {
for (Int_t j=1; j<nyTh1+1; j++) {
if (h_angle_rec1->GetBinContent(i,j)<1) h_angle_rec1->SetBinContent(i,j,0);}}

Int_t nxTh11 = h_angle_MC2->GetNbinsX();
Int_t nyTh11 = h_angle_MC2->GetNbinsY();
for (Int_t i=1; i<nxTh11+1; i++) {
for (Int_t j=1; j<nyTh11+1; j++) {
if (h_angle_MC2->GetBinContent(i,j)<1) h_angle_MC2->SetBinContent(i,j,0);}}

Int_t nxTh2 = h_angle_rec3->GetNbinsX();
Int_t nyTh2 = h_angle_rec3->GetNbinsY();
for (Int_t i=1; i<nxTh2+1; i++) {
for (Int_t j=1; j<nyTh2+1; j++) {
if (h_angle_rec3->GetBinContent(i,j)<1) h_angle_rec3->SetBinContent(i,j,0);}}

Int_t nxTh22 = h_angle_MC3->GetNbinsX();
Int_t nyTh22 = h_angle_MC3->GetNbinsY();
for (Int_t i=1; i<nxTh22+1; i++) {
for (Int_t j=1; j<nyTh22+1; j++) {
if (h_angle_MC3->GetBinContent(i,j)<1) h_angle_MC3->SetBinContent(i,j,0);}}

Int_t nxTh23 = h_angle_rec3->GetNbinsX();
Int_t nyTh23 = h_angle_rec3->GetNbinsY();
for (Int_t i=1; i<nxTh23+1; i++) {
for (Int_t j=1; j<nyTh23+1; j++) {
if (h_angle_rec3->GetBinContent(i,j)<1) h_angle_rec3->SetBinContent(i,j,0);}}


TCanvas a3("p","p",700,700);
a3.Divide(2,3);
a3.cd(1);
h_angle_MC1->Draw("COLZ");
a3.cd(2);
h_angle_rec1->Draw("COLZ");
a3.cd(3);
h_angle_MC2->Draw("COLZ");
a3.cd(4);
h_angle_rec2->Draw("COLZ");
a3.cd(5);
h_angle_MC3->Draw("COLZ");
a3.cd(6);
h_angle_rec3->Draw("COLZ");
a3.SaveAs("angle2D.pdf");

/*Int_t nxTh3 = h_angle_noreco->GetNbinsX();
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


