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

 	TFile *inputfile = new TFile("TRMesmer_box_offset_100k_2GeV_2.root");//TRMesmer_boxdiv_offset_2.root");
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

TH1D *vrtx_chi=new TH1D("chie","Chi2 per DOF of the kinematic vrtx for signal",100,0,1500);

double the_gen=0; double thmu_gen=0; double thmu=0; double the=0;
double the_rec=0;
double signal=0.; double reco=0.; double reco1=0.; double more_reco=0.; double reco0=0.;double reco3=0.;
double reco_v=0.; double more_reco_v=0.; double reco0_v=0.;
int yes_e_g=0;int yes_mu_g=0; int yes2=0; int yes_v=0;
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
vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
vector<MUonERecoOutputVertex> vrtx = ReconstructionOutput->reconstructedVertices();

int yes_mu=0;
int yes_e=0;

for(int j=0; j<tracks.size();j++)
{

if(tracks.at(j).processIDofLinkedTrack()==45 and tracks.at(j).sector()==1) TrackIdreco=tracks.at(j).linkedTrackID();
if(tracks.at(j).processIDofLinkedTrack()==45 and tracks.size()>=3 and tracks.at(j).sector()==1 and  tracks.at(j).percentageOfHitsSharedWithLinkedTrack()>=0) //and tracks.at(0).processIDofLinkedTrack()==45 and tracks.at(0).linkedTrackID()!=tracks.at(1).linkedTrackID()){
{yes2++; cout << "tracks.size " << tracks.size() << endl;

 int sum = tracks.at(j).numberOfXProjectionHits() + tracks.at(j).numberOfYProjectionHits() + tracks.at(j).numberOfStereoHits();

                 if(code_e==tracks.at(j).linkedTrackID()) { yes_e++;
                                        }

                 if(code_mu==tracks.at(j).linkedTrackID()) { yes_mu++;
                                        }
	}
}

if(yes_e>=1 and yes_mu>=1){reco+=MesmerEvent->wgt_full;
for(int j=0; j<vrtx.size();j++)
{
if(vrtx.at(j).stationIndex()==1) //and tracks.at(0).processIDofLinkedTrack()==45 and tracks.at(0).linkedTrackID()!=tracks.at(1).linkedTrackID()){
{
  if(vrtx.at(j).chi2perDegreeOfFreedom()<150)yes_v++; if(j==0)vrtx_chi->Fill(vrtx.at(j).chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);
 }
}
}

if(yes2<2 and TrackIdreco==-99){reco0+=MesmerEvent->wgt_full;
cout <<"NOT RECONSTRUCTED"<<endl;
	}

if(yes2<2 and TrackIdreco!=-99){reco1+=MesmerEvent->wgt_full;
	}

if(yes_v>=1 and yes_e>=1 and yes_mu>=1){reco_v+=MesmerEvent->wgt_full;
        }

if(yes_v>1 and yes_e>=1 and yes_mu>=1){more_reco_v+=MesmerEvent->wgt_full;
        }

if(yes_v<1 and yes_e>=1 and yes_mu>=1){reco0_v+=MesmerEvent->wgt_full;
cout <<"NOT RECONSTRUCTED vertex"<<endl;
        }



}
cout << "---------------------"<<endl;
yes_e_g=0;yes_mu_g=0;
code_e=-99;code_mu=-99;
point_el=0;point_mu=0;
TrackIdreco=-99;
yes2=0;yes_v=0;
} //end of general for

double ratio_v =reco_v/reco;
double ratio0_v =reco0_v/reco;
double ratioM_v =more_reco_v/reco;

cout << "Su " << signal << " eventi di segnale, " << reco_v << " sono ricostruiti con almeno 1 vertice, con un rapporto del " << ratio_v*100 << "%"<< endl;
cout << "Su " << signal << " eventi di segnale, " << more_reco_v << " sono ricostruiti con piu' di unun vertice, con un rapporto del " << ratioM_v*100 << "%"<< endl;
cout << "Su " << signal << " eventi di segnale, " << reco0_v << " hanno 0 vertici, con un rapporto del " << ratio0_v*100 << "%"<< endl;

TCanvas o1("o1","o1",700,700);

/*h_quality->SetLineWidth(5);
h_quality->Draw("hist");
h_quality_more->SetLineWidth(5);
h_quality_more->SetLineColor(kOrange);
h_quality_more->Draw("hist same");*/
vrtx_chi->Draw("hist");
//vrtx_chi->SetMinimum(0);
gPad->SetLogy();
//vrtx_chi->SetMinimum(0);
o1.SaveAs("vrtx_chi_sig.pdf");
/*TCanvas o("o","o",700,700);
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
a2.SaveAs("pdf_notReco/angle.pdf");

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


