#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TMath.h"
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

void macro_PPsig(){


//	TFile *inputfile = new TFile("TRPP_2shared_veryloose.root");//narrowBP_20GeV_1shared.root");//trPoints,20GeV.root");//exampleProductionJob.root");narrowBP.root
//        TTree* cbmsim = (TTree*) inputfile->Get("cbmsim");

TChain * cbmsim = new TChain("cbmsim");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/gen_digi_reco/WiP_v0140_commit_258faf6b_0_infmrad_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/bkg/mesmer_bkg_2hit_1M_1.root");//TRMesmer_2.root");

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

TH1D *opening = new TH1D("op0","opening angle when not reco in event noreco and 0 signal particle reco",400,0.,0.04);
TH1D *h_ID=new TH1D("ID","reconstructed particle ID when no PPp PPm", 50,0,50);

double sig=0.; double bkg=0.; double three_trks=0.;
double ada_PP=0.;
double ada_sig=0.;
int truth=99;

for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
		cbmsim->GetEntry(i);
		if(i%1000 == 0) cout<<"Entry "<<i<<endl;

double thep_gen=0; double them_gen=0;
double thel_gen=0; double thmu_gen=0;
int code_ep=99;
int code_el=99;
int code_mu=99; double Zem=0.; double Zep=0.;

	TVector3 pep,pel,pmu, pmu_in, pmu_44;

	for(int n = 0; n < MCTrack->GetEntries(); n++) {
	 const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
         if(MCTr->interactionID()==0 and MCTr->pdgCode()==-13) pmu_in.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz());

	 if(MCTr->interactionID()==45){
         if(MCTr->pdgCode()==-11) {pep.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); code_ep=n;Zem=MCTr->startZ();}
	 if(MCTr->pdgCode()==11){pel.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); code_el=n;}
         if(MCTr->pdgCode()==-13){pmu.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); code_mu=n;}
	}

	}
if(code_el!=99 and code_mu!=99) truth=1;
if(code_el!=99 and code_ep!=99) truth=0;


	TVector3 pep_dir, pel_dir,pmu_dir, pmu_in_dir;

pep_dir=pep.Unit();
pel_dir=pel.Unit();
pmu_dir=pmu.Unit();
pmu_in_dir=pmu_in.Unit();

double thep_sdr, thel_sdr,thmu_sdr;
                thep_gen=pmu_in_dir.Angle(pep_dir);
		thel_gen=pmu_in_dir.Angle(pel_dir);
		thmu_gen=pmu_in_dir.Angle(pmu_dir);

//double opening = acos(pep_dir.Dot(pel_dir));

int yes_mu=0; int yes_ep=0; int yes_el=0; int yes_mu_in=0; int yes_mu_44=0;

TVector3 pep_rec, pel_rec,pmu_rec, pmu_in_rec, in;
double thep_rec, thel_rec,thmu_rec, thmu_in_rec;
vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();

for(int j=0; j<tracks.size();j++)
{
double th_inx,th_iny;
if(tracks.at(j).processIDofLinkedTrack()==0 and tracks.at(j).sector()==0){
th_inx=tracks.at(j).xSlope();
th_iny=tracks.at(j).ySlope();
in.SetXYZ(th_inx,th_iny,1.0);
}
}


 for(int j=0; j<tracks.size();j++)
 {

if(tracks.at(j).linkedTrackID()==code_mu and tracks.at(j).processIDofLinkedTrack()==45) {yes_mu=1; TVector3 p; p.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.); pmu_rec=p.Unit(); thmu_rec=in.Angle(pmu_rec);}
if(tracks.at(j).linkedTrackID()==code_ep and tracks.at(j).processIDofLinkedTrack()==45) {yes_ep=1; TVector3 p; p.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.); pep_rec=p.Unit(); thep_rec=in.Angle(pep_rec);}
if(tracks.at(j).linkedTrackID()==code_el and tracks.at(j).processIDofLinkedTrack()==45) {yes_el=1; TVector3 p; p.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.); pel_rec=p.Unit(); thel_rec=in.Angle(pel_rec);}

}

double ppp=0;
double ppm=0;
double s=0;

MUonERecoOutputVertex vrtx = ReconstructionOutput->bestVertex();
double chi=vrtx.chi2perDegreeOfFreedom();

if(chi!=0){
 if(tracks.size()==3 and yes_mu==1 and yes_ep==1 and truth==0) {bkg++; ppp=1; h_angle_MC1->Fill(thep_rec,thmu_rec,MesmerEvent->wgt_full);  ada_PP++; opening->Fill(pep_rec.Angle(pmu_rec),MesmerEvent->wgt_full);}
 if(tracks.size()==3 and yes_mu==1 and yes_el==1 and truth==0) {bkg++; ppm=1; h_angle_MC1->Fill(thel_rec,thmu_rec,MesmerEvent->wgt_full);  ada_PP++; opening->Fill(pel_rec.Angle(pmu_rec),MesmerEvent->wgt_full);}
 if(tracks.size()==3 and yes_mu==1 and yes_el==1 and truth==1) {sig++; s=1; h_angle_MC3->Fill(thel_rec,thmu_rec,MesmerEvent->wgt_full);  ada_sig++; opening->Fill(pel_rec.Angle(pmu_rec),MesmerEvent->wgt_full);}

 if(tracks.size()==3) {three_trks++; }

	if(ppp==1 | ppm==1){Z_PP->Fill(vrtx.zPositionFit()); if(ppm==1)Z_PP_gen->Fill(Zem); if(ppp==1)Z_PP_gen->Fill(Zep); }
	if(s==1)Z_el->Fill(vrtx.zPositionFit(),MesmerEvent->wgt_full);
 }//enf of chi

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
a3.SaveAs("bkg_mesmer/anglePP.pdf");

TCanvas a4("p4","p4",700,700);
Z_el->SetLineColor(kRed);
Z_el->Draw();
Z_PP->Draw("same");
Z_PP_gen->SetLineColor(kGreen);
Z_PP_gen->Draw("same");
a4.SaveAs("bkg_mesmer/Z.pdf");

TCanvas a5("a5","a5",700,700);
a5.Divide(1,2);
a5.cd(1);
opening->Draw("hist");
a5.cd(2);
h_ID->Draw("hist");
a5.SaveAs("bkg_mesmer/opening.pdf");

}

