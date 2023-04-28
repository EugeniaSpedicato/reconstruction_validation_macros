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


//	TFile *inputfile = new TFile("TRPP_2shared_veryloose.root");//narrowBP_20GeV_1shared.root");//trPoints,20GeV.root");//exampleProductionJob.root");narrowBP.root
//        TTree* cbmsim = (TTree*) inputfile->Get("cbmsim");

TChain * cbmsim = new TChain("cbmsim");
cbmsim->Add("TRPP_1k.root");
cbmsim->Add("TRMesmer_1k_unw.root");//TRMesmer_2.root");

        TClonesArray *MCTrack = 0;
        TClonesArray *SignalTracks = 0;
        TClonesArray *TrackerStripDigis = 0;
        TClonesArray *TrackerPoints = 0;
        TClonesArray *TrackerStubs = 0;
        MuE::Event *MesmerEvent = 0;
        MUonERecoOutput *ReconstructionOutput = 0;

        cbmsim->SetBranchAddress("MCTrack", &MCTrack);
        cbmsim->SetBranchAddress("SignalTracks", &SignalTracks);
        cbmsim->SetBranchAddress("MesmerEvent", &MesmerEvent);
        cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);

TH1F *Z_PP=new TH1F("zpp","Z position from adaptive fit mmuon+el/pos from PP",1000,1000,1100);
TH1F *Z_el=new TH1F("zel","Z position from adaptive fit mmuon+el from elastic",1000,1000,1100);
TH1F *Z_PP_gen=new TH1F("zpp","Z position from adaptive fit mmuon+el/pos from PP generated",1000,1000,1100);

TH1D *h_anglee0 = new TH1D("e0","polar reconstructed electron angle when event is reco",700,0,0.07);
TH2D *h_anE_el = new TH2D("eth0","reconstructed electron angle-energy when event is reco",700,0,0.07,300,0,60);
TH2D *h_angle_MC_rec = new TH2D("thrmc0","polar electron angle reco(X) VS generated(Y)",700,0,0.07,50,0,0.005);
TH1D *h_En_e = new TH1D("ene0"," generated electron energy when event is reco" ,300,0,60);
TH1D *h_anglee1 = new TH1D("e1","electron angle when event isn't reco",700,0,0.07);

TH2D *h_angle_MC1 = new TH2D("thmc3","polar angle electron vs positron generated",700,0,0.07,50,0,0.005);
TH2D *h_angle_MC2 = new TH2D("thmc2","polar angle electron vs muon generated",700,0,0.07,50,0,0.005);
TH2D *h_angle_MC3 = new TH2D("thmc4","polar angle electron vs muon generated",700,0,0.07,50,0,0.005);

TH1D *h_size=new TH1D("s","track multiplicity intID==5",10,0,10);
TH1D *h_part=new TH1D("p","reconstructed particle", 30,-15,15);

TH1D *opening = new TH1D("op0","opening angle when not reco in event noreco and 0 signal particle reco",400,0,0.04);
TH1D *h_ID=new TH1D("ID","reconstructed particle ID when no PPp PPm", 50,0,50);

double sig=0.; double bkg=0.; double three_trks=0.;
double ada_PP=0.;
double ada_sig=0.;


for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
		cbmsim->GetEntry(i);
		if(i%1000 == 0) cout<<"Entry "<<i<<endl;

double thep_gen=0; double them_gen=0;
double thel_gen=0; double thmu_gen=0;
int code_ep,code_em,code_el,code_mu; double Zem=0.; double Zep=0.;

	TVector3 pep,pem, pel,pmu, pmu_in, pmu_44;

	for(int n = 0; n < MCTrack->GetEntries(); n++) {
	 const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
         if(MCTr->interactionID()==0 and MCTr->pdgCode()==13) pmu_in.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz());

         if(MCTr->interactionID()==44) pmu_44.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz());

	 if(MCTr->interactionID()==5){
	 if(MCTr->pdgCode()==11) {pem.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); code_em=n; Zep=MCTr->startZ();}
         if(MCTr->pdgCode()==-11) {pep.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); code_ep=n;Zem=MCTr->startZ();}
	 }
	if(MCTr->interactionID()==45){
	if(MCTr->pdgCode()==11){pel.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); code_el=n;}
        if(MCTr->pdgCode()==13){pmu.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); code_mu=n;}
	}

	}

/*           int last_modXmu=0; int last_modXel=0; //= qx+mx*(189.9218);
           int last_modYmu=0; int last_modYel=0; // = qy+my*(193.7693);
           int last_modXem=0; int last_modXep=0; //= qx+mx*(189.9218);
           int last_modYem=0; int last_modYep=0; // = qy+my*(193.7693);

           int stereo_mu=0; int stereo_el=0;
           int stereo_em=0; int stereo_ep=0;


                         for(int s=0; s<TrackerPoints->GetEntries(); s++)
                         {const MUonETrackerPoint *TrackerPt = static_cast<const MUonETrackerPoint*>(TrackerPoints->At(s));
                          if(TrackerPt->trackPDGCode()==11 and TrackerPt->trackID()==code_el and TrackerPt->stationID()==1){
                                                                                                if(TrackerPt->moduleID()==4) last_modXel++; 
                                                                                                 if(TrackerPt->moduleID()==5) last_modYel++;
                                                                                                 if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_el++;}
                          if(TrackerPt->trackPDGCode()==13 and TrackerPt->trackID()==code_mu and TrackerPt->stationID()==1){
                                                                                                 if(TrackerPt->moduleID()==4) last_modXmu++;
                                                                                                 if(TrackerPt->moduleID()==5) last_modYmu++;
                                                                                                 if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_mu++;}
                          if(TrackerPt->trackPDGCode()==-11 and TrackerPt->trackID()==code_ep and TrackerPt->stationID()==1){
                                                                                                 if(TrackerPt->moduleID()==4) last_modXep++; 
                                                                                                 if(TrackerPt->moduleID()==5) last_modYep++;
                                                                                                 if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_ep++;}
                          if(TrackerPt->trackPDGCode()==11 and TrackerPt->trackID()==code_em and TrackerPt->stationID()==1){
                                                                                                 if(TrackerPt->moduleID()==4) last_modXem++; 
                                                                                                 if(TrackerPt->moduleID()==5) last_modYem++;
                                                                                                 if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_em++;}
                         }

double ok_mu=0.; double ok_el=0.; double ok_em=0.; double ok_ep=0.;

           if(last_modXmu==2 and last_modYmu==2 and stereo_mu>1) ok_mu=1;
           if(last_modXel==2 and last_modYel==2 and stereo_el>1) ok_el=1; 
           if(last_modXem==2 and last_modYem==2 and stereo_em>1) ok_em=1; 
           if(last_modXep==2 and last_modYep==2 and stereo_ep>1) ok_ep=1; 
*/

	TVector3 pem_dir,pep_dir, pel_dir,pmu_dir, pmu_in_dir, pmu_44_dir;
pep_dir=pep.Unit();
pem_dir=pem.Unit();
pel_dir=pel.Unit();
pmu_dir=pmu.Unit();
pmu_in_dir=pmu_in.Unit();
pmu_44_dir=pmu_44.Unit();

double them_sdr,thep_sdr, thel_sdr,thmu_sdr;
		them_gen=acos(pmu_in_dir.Dot(pem_dir));
                thep_gen=acos(pmu_in_dir.Dot(pep_dir));
		thel_gen=acos(pmu_in_dir.Dot(pel_dir));
		thmu_gen=acos(pmu_in_dir.Dot(pmu_dir));

//double opening = acos(pep_dir.Dot(pem_dir));

int yes_mu=0; int yes_em=0; int yes_ep=0; int yes_el=0; int yes_mu_in=0; int yes_mu_44=0;

TVector3 pem_rec,pep_rec, pel_rec,pmu_rec, pmu_in_rec, in,pmu_44_rec;
double them_rec,thep_rec, thel_rec,thmu_rec, thmu_in_rec,thmu_44_rec;
vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();

for(int j=0; j<tracks.size();j++)
{
double th_inx,th_iny;
if(tracks.at(j).processIDofLinkedTrack()==0 and tracks.at(j).sector()==0){
th_inx=tracks.at(j).xSlope();
th_iny=tracks.at(j).ySlope();
in.SetXYZ(th_inx,th_iny,1.0);
//th_in_recTR=in.Theta();
}
}


 for(int j=0; j<tracks.size();j++)
 {
  if(tracks.at(j).processIDofLinkedTrack()==5) {h_size->Fill(tracks.size());
						if(tracks.at(j).linkedTrackID()==code_ep){h_part->Fill(-11);}
						else if(tracks.at(j).linkedTrackID()==code_em){h_part->Fill(11);}
						else h_part->Fill(0);}


if(tracks.at(j).linkedTrackID()==code_mu and tracks.at(j).processIDofLinkedTrack()==45) {yes_mu=1; TVector3 p; p.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.); pmu_rec=p.Unit(); thmu_rec=acos(in.Dot(pmu_rec));}
if(tracks.at(j).linkedTrackID()==0 and tracks.at(j).sector()==1) {yes_mu_in=1; TVector3 p; p.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.); pmu_in_rec=p.Unit(); thmu_in_rec=acos(in.Dot(pmu_in_rec));}
if(tracks.at(j).processIDofLinkedTrack()==44 and tracks.at(j).sector()==1) {yes_mu_44=1; TVector3 p; p.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.); pmu_44_rec=p.Unit(); thmu_44_rec=acos(in.Dot(pmu_in_rec));}
if(tracks.at(j).linkedTrackID()==code_em) {yes_em=1; TVector3 p; p.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.); pem_rec=p.Unit(); them_rec=acos(in.Dot(pem_rec));}
if(tracks.at(j).linkedTrackID()==code_ep) {yes_ep=1; TVector3 p; p.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.); pep_rec=p.Unit(); thep_rec=acos(in.Dot(pep_rec));}
if(tracks.at(j).linkedTrackID()==code_el and tracks.at(j).processIDofLinkedTrack()==45) {yes_el=1; TVector3 p; p.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.); pel_rec=p.Unit(); thel_rec=acos(in.Dot(pel_rec));}

}

double ppp=0;
double ppm=0;
double s=0;
vector<MUonERecoOutputAdaptiveFitterVertex> vrtx = ReconstructionOutput->adaptiveFitterVertices();
vector<MUonERecoOutputVertex> vrtx_kin = ReconstructionOutput->reconstructedVertices();


if(tracks.size()==3 and yes_mu_in==1 and yes_em==1) {bkg++; ppm=1; h_angle_MC1->Fill(them_rec,thmu_in_rec); if(vrtx.size()!=0) ada_PP++; opening->Fill(acos(pem_rec.Dot(pmu_in_rec)));}
if(tracks.size()==3 and yes_mu_in==1 and yes_ep==1) {bkg++; ppp=1; h_angle_MC1->Fill(thep_rec,thmu_in_rec); if(vrtx.size()!=0) ada_PP++; opening->Fill(acos(pep_rec.Dot(pmu_in_rec)));}
if(tracks.size()==3 and yes_mu_44==1 and yes_em==1) {bkg++; ppm=1; h_angle_MC1->Fill(them_rec,thmu_44_rec); if(vrtx.size()!=0) ada_PP++; opening->Fill(acos(pem_rec.Dot(pmu_44_rec)));}
if(tracks.size()==3 and yes_mu_44==1 and yes_ep==1) {bkg++; ppp=1; h_angle_MC1->Fill(thep_rec,thmu_44_rec); if(vrtx.size()!=0) ada_PP++; opening->Fill(acos(pep_rec.Dot(pmu_44_rec)));}
if(tracks.size()==3 and yes_mu==1 and yes_el==1) {sig++; s=1; h_angle_MC3->Fill(thel_rec,thmu_rec); if(vrtx.size()!=0) ada_sig++; opening->Fill(acos(pel_rec.Dot(pmu_rec)));}

if(tracks.size()==3 and ppm==0 and ppp==0)
{

 for(int j=0; j<tracks.size();j++)
 {
  h_ID->Fill(tracks.at(j).processIDofLinkedTrack());
  cout <<j<< ") ID " << tracks.at(j).processIDofLinkedTrack() << " and trID " << tracks.at(j).linkedTrackID() <<  endl;
 }
}

//if(yes_mu_in==1 and yes_em==1) h_angle_MC2->Fill(them_gen,thmu_gen);
//if(yes_mu_in==1 and yes_ep==1) h_angle_MC2->Fill(thep_gen,thmu_gen);
//if(yes_mu==1 and yes_el==1) h_angle_MC3->Fill(thel_gen,thmu_gen);


if(tracks.size()==3) {three_trks++; }

for(int j=0; j<vrtx.size();j++)
{
	if(ppp==1 | ppm==1){Z_PP->Fill(vrtx.at(j).z()); if(ppm==1)Z_PP_gen->Fill(Zem); if(ppp==1)Z_PP_gen->Fill(Zep); }
	if(s==1)Z_el->Fill(vrtx.at(j).z(),MesmerEvent->wgt_full);
}


} //end of general for

cout << "Su " << three_trks << " eventi con tre tracce, " << bkg << " sono muone+elettrone/positrone da PP: " << (bkg/three_trks)*100 << endl;
cout << "Su " << three_trks << " eventi con tre tracce, " << sig << " sono muone+elettrone da elastico: " << (sig/three_trks)*100 << endl;
cout << "Segnale: " << (ada_sig/sig)*100 << " delle volte ho adaptive fit vertex" << endl;
cout << "PP: " << (ada_PP/bkg)*100 << " delle volte ho adaptive fit vertex" << endl;


Int_t nxTh = h_angle_MC1->GetNbinsX();
Int_t nyTh = h_angle_MC1->GetNbinsY();
for (Int_t i=1; i<nxTh+1; i++) {
for (Int_t j=1; j<nyTh+1; j++) {
if (h_angle_MC1->GetBinContent(i,j)<1) h_angle_MC1->SetBinContent(i,j,0);}}

Int_t nxTh2 = h_angle_MC2->GetNbinsX();
Int_t nyTh2 = h_angle_MC2->GetNbinsY();
for (Int_t i=1; i<nxTh2+1; i++) {
for (Int_t j=1; j<nyTh2+1; j++) {
if (h_angle_MC2->GetBinContent(i,j)<1) h_angle_MC2->SetBinContent(i,j,0);}}

TCanvas a3("p","p",700,700);
a3.Divide(1,3);
a3.cd(1);
h_angle_MC3->SetMarkerColor(kBlue);
h_angle_MC3->Draw();
h_angle_MC2->SetMarkerColor(kGreen);
h_angle_MC2->Draw("same");
h_angle_MC1->SetMarkerColor(kGreen);
h_angle_MC1->Draw("same");
a3.cd(2);
h_size->Draw("hist");
a3.cd(3);
h_part->Draw("hist");
a3.SaveAs("anglePP.pdf");

TCanvas a4("p4","p4",700,700);
Z_el->SetLineColor(kRed);
Z_el->Draw();
Z_PP->Draw("same");
Z_PP_gen->SetLineColor(kGreen);
Z_PP_gen->Draw("same");
a4.SaveAs("Z.pdf");

TCanvas a5("a5","a5",700,700);
a5.Divide(1,2);
a5.cd(1);
opening->Draw("hist");
a5.cd(2);
h_ID->Draw("hist");
a5.SaveAs("opening.pdf");

}

