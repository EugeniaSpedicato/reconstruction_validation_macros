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


TChain * cbmsim = new TChain("cbmsim");
cbmsim->Add("TRMesmer_box_100k.root");

        TClonesArray *MCTrack = 0;
        TClonesArray *SignalTracks = 0;
        TClonesArray *TrackerStripDigis = 0;
        TClonesArray *TrackerPoints = 0;
        TClonesArray *TrackerStubs = 0;
        MUonERecoOutput *ReconstructionOutput = 0;

        cbmsim->SetBranchAddress("MCTrack", &MCTrack);
        cbmsim->SetBranchAddress("SignalTracks", &SignalTracks);
        cbmsim->SetBranchAddress("TrackerPoints", &TrackerPoints);
        cbmsim->SetBranchAddress("TrackerStubs", &TrackerStubs);
        cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);


TH1D *h_size=new TH1D("s","track multiplicity in second station",15,0,15);
TH1D *h_part=new TH1D("p","reconstructed particle ID in second station", 50,0,50);

TH1D *h_part1=new TH1D("p1","reconstructed particle ID multiplicity=1", 50,0,50);
TH1D *h_part2=new TH1D("p2","reconstructed particle ID multiplicity=2", 50,0,50);
TH1D *h_part3=new TH1D("p3","reconstructed particle ID multiplicity=3", 50,0,50);
TH1D *h_part_more=new TH1D("pm","reconstructed particle ID multiplicity>3", 50,0,50);

TH1D *h_quality1=new TH1D("quality1","reconstructed particle quality multiplicity=1", 101,0,101);
TH1D *h_quality2=new TH1D("quality2","reconstructed particle quality multiplicity=2", 101,0,101);
TH1D *h_quality3=new TH1D("quality3","reconstructed particle quality multiplicity=3", 101,0,101);
TH1D *h_quality_more=new TH1D("qualitym","reconstructed particle ID multiplicity>3", 101,0,101);

TH1D *h_nstubs=new TH1D("stubs1","nstubs multiplicity=1", 20,0,20);
TH1D *h_nstubs1=new TH1D("stubs1","nstubs multiplicity=2", 20,0,20);


TH1D *h_quality2_PP=new TH1D("quality2pp","reconstructed particle quality multiplicity=2 PP", 101,0,101);
TH1D *h_quality3_PP=new TH1D("quality3pp","reconstructed particle quality multiplicity=3 PP", 101,0,101);

TH1D *h_quality2_no=new TH1D("quality2no","reconstructed particle quality multiplicity=2 no correct PP", 101,0,101);
TH1D *h_quality3_no=new TH1D("quality3no","reconstructed particle quality multiplicity=3 no correct PP", 101,0,101);

TH1D *h_part2_bool=new TH1D("p3b","reconstructed particle ID multiplicity=2 are always mu+e?", 4,-2,2);
TH1D *h_part2_mue=new TH1D("p3mue","reconstructed particle ID multiplicity=2 when not mu+e?", 50,0,50);
TH1D *h_part3_mue=new TH1D("p3mue","reconstructed particle ID multiplicity=3 when not mu+e+e?", 50,0,50);

TH1D *h_chi_em_3=new TH1D("chi0","chi2/ndof reconstructed electron track when multiplicity==2(yellow),3(blue)",160,0,80);
TH1D *h_chi_em_4=new TH1D("chi1","chi2/ndof reconstructed electron track when multiplicity==2(yellow),3(blue)",160,0,80);

TH1D *h_chi_mu_3=new TH1D("chi2","chi2/ndof reconstructed positron track when multiplicity==2(yellow),3(blue)",160,0,80);
TH1D *h_chi_mu_4=new TH1D("chi3","chi2/ndof reconstructed positron track when multiplicity==2(yellow),3(blue)",160,0,80);

TH1D *h_chi_in_3=new TH1D("chi4","chi2/ndof reconstructed mu_in=0 track when multiplicity==2(yellow),3(blue)",160,0,80);
TH1D *h_chi_in_4=new TH1D("chi5","chi2/ndof reconstructed mu_in=0 track when multiplicity==2(yellow),3(blue)",160,0,80);

TH1D *h_chi_44_3=new TH1D("chi6","chi2/ndof reconstructed mu_in=44 track when multiplicity==2(yellow),3(blue)",160,0,80);
TH1D *h_chi_44_4=new TH1D("chi7","chi2/ndof reconstructed mu_in=44 track when multiplicity==2(yellow),3(blue)",160,0,80);

TH1D *opening = new TH1D("op0","opening angle when not reco in event noreco and 0 signal particle reco",400,0,0.04);
TH1D *h_ID=new TH1D("ID","reconstructed particle ID when no PPp PPm", 50,0,50);

TH1D *h_Em=new TH1D("1","Energy of MC muon 44 when there is 55 interaction (GeV)",180,135,190);
TH1D *h_Em_no=new TH1D("1no","Energy of MC muon 44 when no 55 interaction (GeV)",180,135,190);
TH1D *h_thm=new TH1D("2","Angle of the outgoing mu from PP wrt incoming mu (rad)",200,0,0.0002);
TH1D *h_thm_no=new TH1D("2","Angle of the outgoing mu no PP wrt incoming mu (rad)",200,0,0.0002);

TH1D *h_Ee=new TH1D("e0","Energy of electron from PP (GeV)",100,0,10);
TH1D *h_Ep=new TH1D("e1","Energy of positron from PP (GeV)",100,0,10);
TH1D *h_thp=new TH1D("t2","Angle of the electron from PP wrt incoming muon (rad)",100,0,0.01);
TH1D *h_the=new TH1D("t3","Angle of the positron from PP wrt incoming muon (rad)",100,0,0.01);

TH1D *th_mu44g = new TH1D("th_mu44g"," angle mu 44 wrt incoming muon (generated)",200,0,0.005);
TH1D *th_mu44 = new TH1D("th_mu44"," angle mu 44 wrt incoming muon (reconstructed)",200,0,0.005);
TH1D *th_e44g = new TH1D("th_e44e"," angle e 44 wrt incoming muon (generated)",200,0,0.05);
TH1D *th_e44 = new TH1D("th_mu44_r"," angle e 44 wrt incoming muon (reconstructed)",200,0,0.05);

TH1D *h_E2=new TH1D("E44in","Energy e_e+e_p+e44-ein (GeV)",100,-10,10);

double sig=0.; double bkg3=0.; double bkg4=0.; double more_tracks=0.; double bkg=0.;double bkg3no=0.; double bkg4no=0.;
double bkg2_g=0.; double bkg3_g=0.; double bkg4_g=0.; double more_tracks_g=0.;
double ada_PP=0.;
double ada_sig=0.;
double gen = 0.;

double c_44=0.;
double c_5_44=0.;
double c_5=0.;

double n_st0=0.; double n_st1=0.;

for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
		cbmsim->GetEntry(i);
		if(i%1000 == 0) cout<<"Entry "<<i<<endl;



double e_epp,e44,einpp,e_p;
int yes_44=0;
int yespp=0;
        for(int n = 0; n < MCTrack->GetEntries(); n++) {
         const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));


         if(MCTr->interactionID()==0 and MCTr->pdgCode()==-13) {einpp=MCTr->energy();}
         if(MCTr->interactionID()==44 and MCTr->pdgCode()==-13) { yes_44++;
                                                                 e44=MCTr->energy();}

         if(MCTr->interactionID()==5){
         if(MCTr->pdgCode()==11) {yespp++;e_epp=MCTr->energy(); cout << "motherID " << MCTr->motherID() << endl;}

         if(MCTr->pdgCode()==-11) {yespp++;e_p=MCTr->energy();}

         }

	}


if(yespp==2 and yes_44==1)
{
// h_E1->Fill(ein-e44);
double Eres=einpp-e44-e_epp-e_p;
double Pres=(pmu_in.Mag()-pmu_44.Mag())*(pmu_in.Mag()-pmu_44.Mag())-(pem.Mag()+pep.Mag())*(pem.Mag()+pep.Mag());
// h_E2->Fill(Eres);
cout << "Eres PP" << Eres << endl;
}
















double thmu_gen=0; double them_gen=0; double the44_gen=0;
int code_mu=0; int code_em=0; int code_44=0; int code_44_2=0;
double yes=0;
double e_e,e_mu,ein;
vector<double> en_ph;
en_ph.reserve(2);
double mass;
	TVector3 pmu,pem, pmu_in, pmu_44;
                cout << "Entries: " << MCTrack->GetEntries() <<endl;
	for(int n = 0; n < MCTrack->GetEntries(); n++) {
	 const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
                cout << n <<  ") IntID : " <<MCTr->interactionID() << " and pdf " << MCTr->pdgCode() << endl;

         if(MCTr->interactionID()==0 and MCTr->pdgCode()==-13){pmu_in.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); ein=MCTr->energy();
								double mass1=sqrt(ein*ein-pmu_in.Mag()*pmu_in.Mag()); cout << "mass mu "<< mass1 << endl;}
	 if(MCTr->interactionID()==45){
	 if(MCTr->pdgCode()==11) {yes++; pem.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); code_em=n; e_e=MCTr->energy();
                                                                mass=sqrt(e_e*e_e-pem.Mag()*pem.Mag()); cout << "mass e "<< mass << endl;}

         if(MCTr->pdgCode()==-13) {yes++;pmu.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); code_mu=n; e_mu=MCTr->energy();
                                                                double mass2=sqrt(e_mu*e_mu-pmu.Mag()*pmu.Mag()); cout << "mass mu out "<< mass2 << endl;}

         if(MCTr->pdgCode()==22) {yes++; en_ph.push_back(MCTr->energy());}

	 }

	}

if(yes>=2)
{
double e_ph=0;
for(int p=0;p<en_ph.size();p++){e_ph+=en_ph.at(p);}
double Eres=e_ph+e_e+e_mu-ein-mass;
 h_E2->Fill(Eres);
cout << "Eres " << Eres << endl;
}

           int last_modXmu=0;
           int last_modYmu=0;
           int last_modXem=0;
           int last_modYem=0;
           int stereo_em=0; int stereo_mu=0;

                         for(int s=0; s<TrackerPoints->GetEntries(); s++)
                         {const MUonETrackerPoint *TrackerPt = static_cast<const MUonETrackerPoint*>(TrackerPoints->At(s));

                          if(TrackerPt->trackPDGCode()==-13 and TrackerPt->trackID()==code_mu and TrackerPt->stationID()==1){
                                                                                                 if(TrackerPt->moduleID()==4) last_modXmu++;
                                                                                                 if(TrackerPt->moduleID()==5) last_modYmu++;
                                                                                                 if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_mu++;}
                          if(TrackerPt->trackPDGCode()==11 and TrackerPt->trackID()==code_em and TrackerPt->stationID()==1){
                                                                                                 if(TrackerPt->moduleID()==4) last_modXem++;
                                                                                                 if(TrackerPt->moduleID()==5) last_modYem++;
                                                                                                 if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_em++;}

                         }




double ok_em=0.; double ok_mu=0.;

           if(last_modXem==2 and last_modYem==2 and stereo_em>1) ok_em=1;
           if(last_modXmu==2 and last_modYmu==2 and stereo_mu>1) ok_mu=1;


	TVector3 pem_dir,pmu_dir, pmu_in_dir;
	pmu_dir=pmu.Unit();
	pem_dir=pem.Unit();
	pmu_in_dir=pmu_in.Unit();

double them_sdr,thep_sdr;
		them_gen=acos(pmu_in_dir.Dot(pem_dir));
                thmu_gen=acos(pmu_in_dir.Dot(pmu_dir));

cout << "------------------" << endl;

if(ok_em==1 and ok_mu==1){

bkg++;

 int yes_em=0; int yes_mu=0;  int yes_mu_in=0; 
TVector3 pem_rec,pmu_rec, pmu_in_rec, in;
double them_rec,thmu_rec, thmu_in_rec;
double chi_em,chi_mu,chi_in;

vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
for(int j=0; j<tracks.size();j++)
{
double th_inx,th_iny;
if(tracks.at(j).processIDofLinkedTrack()==0 and tracks.at(j).sector()==0){
th_inx=tracks.at(j).xSlope();
th_iny=tracks.at(j).ySlope();
in.SetXYZ(th_inx,th_iny,1.0);
in.Unit();
 }
}

double s=0;


int stubs0=0;
int stubs1=0;

for(int s=0; s<TrackerStubs->GetEntries(); s++){
const MUonETrackerStub *stub = static_cast<const MUonETrackerStub*>(TrackerStubs->At(s));
if(stub->stationID()==0)stubs0++;
if(stub->stationID()==1)stubs1++;
	}

vector<int> perc;
perc.reserve(tracks.size());

for(int j=0; j<tracks.size();j++)
 {
cout << "numero tracce " << tracks.size() << " ID processo " << tracks.at(j).processIDofLinkedTrack() << " sec " << tracks.at(j).sector() <<  endl;
if(tracks.at(j).percentageOfHitsSharedWithLinkedTrack()>=0){

perc.push_back(tracks.at(j).percentageOfHitsSharedWithLinkedTrack());

if(tracks.at(j).sector()==1) {h_part->Fill(tracks.at(j).processIDofLinkedTrack()); }
if(tracks.at(j).sector()==1 and tracks.size()==2) {h_part1->Fill(tracks.at(j).processIDofLinkedTrack()); h_quality1->Fill(tracks.at(j).percentageOfHitsSharedWithLinkedTrack());
							h_nstubs->Fill(stubs1);}
if(tracks.at(j).sector()==1 and tracks.size()==3) { h_part2->Fill(tracks.at(j).processIDofLinkedTrack()); h_quality2->Fill(tracks.at(j).percentageOfHitsSharedWithLinkedTrack());
							h_nstubs1->Fill(stubs1);}
if(tracks.at(j).sector()==1 and tracks.size()==4) { h_part3->Fill(tracks.at(j).processIDofLinkedTrack()); h_quality3->Fill(tracks.at(j).percentageOfHitsSharedWithLinkedTrack());}
if(tracks.at(j).sector()==1 and tracks.size()>4) { h_part_more->Fill(tracks.at(j).processIDofLinkedTrack()); h_quality_more->Fill(tracks.at(j).percentageOfHitsSharedWithLinkedTrack());}


if(tracks.at(j).linkedTrackID()==code_em and tracks.at(j).sector()==1 and tracks.at(j).processIDofLinkedTrack()==45) {yes_em=1; TVector3 p; p.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.); pem_rec=p.Unit();
									   them_rec=acos(in.Dot(pem_rec)); chi_em=tracks.at(j).chi2perDegreeOfFreedom();}
if(tracks.at(j).linkedTrackID()==code_mu and tracks.at(j).sector()==1 and tracks.at(j).processIDofLinkedTrack()==45) {yes_mu=1; TVector3 p; p.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.); pmu_rec=p.Unit();
									   thmu_rec=acos(in.Dot(pmu_rec)); chi_mu=tracks.at(j).chi2perDegreeOfFreedom(); }
	}
}

h_size->Fill(tracks.size()-1);
if(tracks.size()==2 and perc.size()==2) {bkg2_g++;}
if(tracks.size()==3 and perc.size()==3) {bkg3_g++;}
if(tracks.size()==4 and perc.size()==4) {bkg4_g++;}
if(tracks.size()>4) {more_tracks_g++;}



double ppm=0;

if(yes_mu==1 and yes_em==1 and tracks.size()>3) {bkg4++;}

if(tracks.size()==3){
		if(yes_mu==1 and yes_em==1) {bkg3++; ppm=1; h_part2_bool->Fill(1);
                                             h_chi_em_3->Fill(chi_em);  h_chi_em_3->Fill(chi_mu);
						 th_e44->Fill(them_rec);th_e44g->Fill(them_gen);
						th_mu44->Fill(thmu_rec);th_mu44g->Fill(thmu_gen);
			for(int j=0; j<tracks.size();j++){ h_quality2_PP->Fill(tracks.at(j).percentageOfHitsSharedWithLinkedTrack());}
			}
		else{h_part2_bool->Fill(0); bkg3no++;
                                for(int j=0; j<tracks.size();j++)
                        {if(tracks.at(j).sector()==1){h_part2_mue->Fill(tracks.at(j).processIDofLinkedTrack());
			h_quality2_no->Fill(tracks.at(j).percentageOfHitsSharedWithLinkedTrack()); }
                        }
    }
 }


	}
} //end of general for
cout << "Numero di eventi generati PP su N eventi input: " << bkg << " ---> " << (bkg/cbmsim->GetEntries())*100 << endl;
cout << "Su " << bkg << " eventi PP, " << bkg2_g << " molteplicita' = 1: " << (bkg2_g/bkg)*100 << endl;
cout << "Su " << bkg << " eventi PP, " << bkg3_g << " molteplicita' = 2: " << (bkg3_g/bkg)*100 << endl;
cout << "Su " << bkg << " eventi PP, " << bkg4_g << " molteplicita' = 3: " << (bkg4_g/bkg)*100 << endl;
cout << "Su " << bkg << " eventi PP, " << more_tracks_g << "molteplicita' > 3: " << (more_tracks_g/bkg)*100 << endl;


cout << "Su " << bkg << " eventi PP, " << bkg3 << " hanno molteplicita 2 e sono muone+elettrone di segnale: " << (bkg3/bkg)*100 << endl;
cout << "Su " << bkg << " eventi PP, " << bkg4 << " sono eventi con alta molteplicita' ma con tracce di segnale: " << (bkg4/bkg)*100 << endl;
cout << "Su " << bkg << " eventi PP, " << bkg3no << " hanno molteplicita 2 ma NON sono muone+elettrone di segnale: " << (bkg3no/bkg)*100 << endl;

//cout << "Segnale: " << (ada_sig/sig)*100 << " delle volte ho adaptive fit vertex" << endl;
//cout << "PP: " << (ada_PP/bkg)*100 << " delle volte ho adaptive fit vertex" << endl;
cout << "Generazione solo 44 e no 5: " << (c_44/cbmsim->GetEntries())*100 << endl;
cout << "Generazione 44 e 5: " << (c_5_44/cbmsim->GetEntries())*100 << endl;
cout << "Generazione no 44 e solo 5: " << (c_5/cbmsim->GetEntries())*100 << endl;

cout << "Volte in cui muone ID=0 ha punti solo stazione 0 " << (n_st0/bkg)*100 << endl;
cout << "Volte in cui muone ID=0 ha punti anche stazione 1 " << (n_st1/bkg)*100 << endl;


/*
TCanvas a5("a5","a5",700,700);
a5.Divide(3,3);
a5.cd(1);
h_part->Draw("hist");
gPad->SetLogy();
a5.cd(2);
h_size->Draw("hist");
gPad->SetLogy();
a5.cd(3);
h_part1->Draw("hist");
gPad->SetLogy();
a5.cd(4);
h_part2->Draw("hist");
gPad->SetLogy();
a5.cd(5);
h_part3->Draw("hist");
gPad->SetLogy();
a5.cd(6);
h_part_more->Draw("hist");
gPad->SetLogy();
a5.cd(7);
h_part2_bool->Draw("hist");
a5.cd(8);
h_part2_mue->Draw("hist");
gPad->SetLogy();
a5.cd(9);
h_part3_mue->Draw("hist");
gPad->SetLogy();
a5.SaveAs("idPerMultiplicity_PPsignal.pdf");


TCanvas a6("a6","a6",700,700);
a6.Divide(2,3);
a6.cd(1);
h_chi_em_3->SetLineColor(kOrange);
h_chi_em_3->Draw("hist");
h_chi_em_4->Draw("hist same");
gPad->SetLogy();
a6.cd(2);
h_chi_mu_3->SetLineColor(kOrange);
h_chi_mu_3->Draw("hist");
h_chi_mu_4->Draw("hist same");
gPad->SetLogy();
a6.cd(3);
h_chi_in_3->SetLineColor(kOrange);
h_chi_in_3->Draw("hist");
h_chi_in_4->Draw("hist same");
gPad->SetLogy();
a6.cd(4);
h_chi_44_3->SetLineColor(kOrange);
h_chi_44_3->Draw("hist");
h_chi_44_4->Draw("hist same");
gPad->SetLogy();
a6.SaveAs("chi2_PPsignal.pdf");


TCanvas a7("a7","a7",700,700);
a7.Divide(2,3);
a7.cd(1);
h_Em->Draw("hist");
h_Em_no->SetLineColor(kOrange);
h_Em_no->Draw("hist same");
a7.cd(3);
h_Ee->Draw("hist");
gPad->SetLogy();
a7.cd(5);
h_Ep->Draw("hist");
gPad->SetLogy();
a7.cd(2);
h_thm->Draw("hist");
h_thm_no->SetLineColor(kOrange);
h_thm_no->Draw("hist same");
a7.cd(4);
h_the->Draw("hist");
a7.cd(6);
h_thp->Draw("hist");
a7.SaveAs("caratt_44signal.pdf");


TCanvas a8("a8","a8",700,700);
a8.Divide(2,2);
a8.cd(1);
h_quality1->Draw("hist");
gPad->SetLogy();
a8.cd(2);
h_quality2->Draw("hist");
h_quality2_PP->SetLineColor(kPink);
h_quality2_PP->Draw("hist same");
h_quality2_no->SetLineColor(kGreen);
h_quality2_no->Draw("hist same");
gPad->SetLogy();
a8.cd(3);
h_quality3->Draw("hist");
h_quality3_PP->SetLineColor(kPink);
h_quality3_PP->Draw("hist same");
h_quality3_no->SetLineColor(kGreen);
h_quality3_no->Draw("hist same");
gPad->SetLogy();
a8.cd(4);
h_quality_more->Draw("hist");
gPad->SetLogy();
a8.SaveAs("qualitysignal.pdf");

TCanvas a9("a9","a9",700,700);
a9.Divide(1,2);
a9.cd(1);
h_nstubs->Draw("hist");
a9.cd(2);
h_nstubs1->Draw("hist");
a9.SaveAs("nstubssignal.pdf");
*/
TCanvas b1("b1","b1",700,700);
b1.Divide(1,2);
b1.cd(1);
th_mu44g->Draw("hist");
th_mu44->SetLineColor(kRed);
th_mu44->Draw("hist same");
gPad->SetLogy();
b1.cd(2);
th_e44g->Draw("hist");
th_e44->SetLineColor(kRed);
th_e44->Draw("hist same");
b1.SaveAs("anglesignal.pdf");

TCanvas b2("b2","b2",700,700);
//h_E1->Draw("hist");
//h_E2->SetLineColor(kGreen);
h_E2->Draw("hist");
//gPad->SetLogy();
b2.SaveAs("energySig.pdf");
}
