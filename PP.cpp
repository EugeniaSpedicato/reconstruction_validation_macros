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
cbmsim->Add("TRPP_100k.root");

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

TH1D *h_chi_ep_3=new TH1D("chi2","chi2/ndof reconstructed positron track when multiplicity==2(yellow),3(blue)",160,0,80);
TH1D *h_chi_ep_4=new TH1D("chi3","chi2/ndof reconstructed positron track when multiplicity==2(yellow),3(blue)",160,0,80);

TH1D *h_chi_in_3=new TH1D("chi4","chi2/ndof reconstructed mu_in=0 track when multiplicity==2(yellow),3(blue)",160,0,80);
TH1D *h_chi_in_4=new TH1D("chi5","chi2/ndof reconstructed mu_in=0 track when multiplicity==2(yellow),3(blue)",160,0,80);

TH1D *h_chi_44_3=new TH1D("chi6","chi2/ndof reconstructed mu_in=44 track when multiplicity==2(yellow),3(blue)",160,0,80);
TH1D *h_chi_44_4=new TH1D("chi7","chi2/ndof reconstructed mu_in=44 track when multiplicity==2(yellow),3(blue)",160,0,80);

TH1D *th_in44 = new TH1D("th_in44"," angle mu 44 wrt incoming muon (generated)",200,0,0.005);
TH1D *th_in44_r = new TH1D("th_in44_r"," angle mu 44 wrt incoming muon (reconstructed)",200,0,0.005);
TH1D *th_e = new TH1D("th_e"," angle e wrt incoming mu (generated)",100,0,0.01);
TH1D *th_e_r = new TH1D("th_e_r"," angle e wrt incoming mu (reconstructed)",100,0,0.01);


TH1D *h_Em=new TH1D("1","Energy of MC muon 44 when there is 55 interaction (GeV)",180,135,190);
TH1D *h_Em_no=new TH1D("1no","Energy of MC muon 44 when no 55 interaction (GeV)",180,135,190);
TH1D *h_thm=new TH1D("2","Angle of the outgoing mu from PP wrt incoming mu (rad)",200,0,0.0002);
TH1D *h_thm_no=new TH1D("2","Angle of the outgoing mu no PP wrt incoming mu (rad)",200,0,0.0002);

TH1D *h_Ee=new TH1D("e0","Energy of electron from PP (GeV)",100,0,10);
TH1D *h_Ep=new TH1D("e1","Energy of positron from PP (GeV)",100,0,10);
TH1D *h_thp=new TH1D("t2","Angle of the electron from PP wrt incoming muon (rad)",100,0,0.01);
TH1D *h_the=new TH1D("t3","Angle of the positron from PP wrt incoming muon (rad)",100,0,0.01);

TH1D *h_E1=new TH1D("E1","Energy of muon_in-muon_44 from PP (GeV)",300,-30,30);
TH1D *h_E2=new TH1D("E2","Energy of electron+positron from PP (GeV)",300,-30,30);
TH1D *h_E3=new TH1D("E3","Energy E_e+E_p+E_44-E_in (GeV)",100,-10,10);

double sig=0.; double bkg3=0.; double bkg4=0.; double more_tracks=0.; double bkg=0.;double bkg3no=0.; double bkg4no=0.;
double bkg3_g=0.; double bkg4_g=0.; double more_tracks_g=0.;
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

double thep_gen=0; double them_gen=0; double the44_gen=0;
int code_ep=0; int code_em=0; int code_44=0; int code_44_2=0;
double Zem=0.; double Zep=0.;
double yes=0;
int yes_44=0;

double e44,z44,x44,y44;
double e_e,e_p,ein;

	TVector3 pep,pem, pmu_in, pmu_44;
                cout << "Entries: " << MCTrack->GetEntries() <<endl;
	for(int n = 0; n < MCTrack->GetEntries(); n++) {
	 const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
                cout << n <<  ") IntID : " <<MCTr->interactionID() << " and pdf " << MCTr->pdgCode() << endl;
        cout << "mum " << MCTr->motherID() << endl;


         if(MCTr->interactionID()==0 and MCTr->pdgCode()==-13) {pmu_in.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); ein=MCTr->energy()-MCTr->mass();
								cout << "Z 44: " << MCTr->startZ() << endl;}
					//double mass=sqrt(ein*ein-pmu_in.Mag()*pmu_in.Mag()); cout << "mass mu "<< mass << endl;}
         if(MCTr->interactionID()==44 and MCTr->pdgCode()==-13) { yes_44++;pmu_44.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); code_44=n;
								 e44=MCTr->energy()-MCTr->mass(); cout << "Z 44: " << MCTr->startZ() << endl;}
                                                                //double mass=sqrt(e44*e44-pmu_44.Mag()*pmu_44.Mag()); cout << "mass mu 44 "<< mass << endl;}

	 if(MCTr->interactionID()==5){
	 if(MCTr->pdgCode()==11) {yes++; pem.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); code_em=n; Zep=MCTr->startZ(); e_e=MCTr->energy()-MCTr->mass();
						cout << "Z em: " << MCTr->startZ() << endl;}
                                                                //double mass=sqrt(e_e*e_e-pem.Mag()*pem.Mag()); cout << "mass e "<< mass << endl;}

         if(MCTr->pdgCode()==-11) {yes++;pep.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); code_ep=n; Zem=MCTr->startZ(); e_p=MCTr->energy()-MCTr->mass();
						cout << "Z ep: " << MCTr->startZ() << endl;}
                                                                //double mass=sqrt(e_p*e_p-pep.Mag()*pep.Mag()); cout << "mass e "<< mass << endl;}

	 }

	}


/*if(yes==2){
{
double Eres=e44+e_e+e_p-ein;
h_E1->Fill(ein-e44);
h_E2->Fill(e_e+e_p);
h_E3->Fill(Eres);
cout << "Eres " << Eres << endl;
}*/

           int last_modXmu=0; int last_modXel=0;
           int last_modYmu=0; int last_modYel=0;
           int last_modXem=0; int last_modXep=0;
           int last_modYem=0; int last_modYep=0;

           int stereo_mu=0; int stereo_el=0;
           int stereo_em=0; int stereo_ep=0;

int st0=0;
int st1=0;
                         for(int s=0; s<TrackerPoints->GetEntries(); s++)
                         {const MUonETrackerPoint *TrackerPt = static_cast<const MUonETrackerPoint*>(TrackerPoints->At(s));

                          if(TrackerPt->trackPDGCode()==-11 and TrackerPt->trackID()==code_ep and TrackerPt->stationID()==1){
                                                                                                 if(TrackerPt->moduleID()==4) last_modXep++;
                                                                                                 if(TrackerPt->moduleID()==5) last_modYep++;
                                                                                                 if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_ep++;}
                          if(TrackerPt->trackPDGCode()==11 and TrackerPt->trackID()==code_em and TrackerPt->stationID()==1){
                                                                                                 if(TrackerPt->moduleID()==4) last_modXem++;
                                                                                                 if(TrackerPt->moduleID()==5) last_modYem++;
                                                                                                 if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_em++;}

                          if(TrackerPt->trackPDGCode()==-13 and TrackerPt->trackID()==0) {if(TrackerPt->stationID()==0) st0=1;
											   if(TrackerPt->stationID()==1) st1=1;}

                         }




double ok_em=0.; double ok_ep=0.;

           if(last_modXem==2 and last_modYem==2 and stereo_em>1) ok_em=1;
           if(last_modXep==2 and last_modYep==2 and stereo_ep>1) ok_ep=1;


	TVector3 pem_dir,pep_dir, pmu_in_dir, pmu_44_dir;
	pep_dir=pep.Unit();
	pem_dir=pem.Unit();
	pmu_in_dir=pmu_in.Unit();
	pmu_44_dir=pmu_44.Unit();

double them_sdr,thep_sdr;
		them_gen=acos(pmu_in_dir.Dot(pem_dir));
                thep_gen=acos(pmu_in_dir.Dot(pep_dir));
                //the44_gen=acos(pmu_in_dir.Dot(pmu_44_dir));



if(yes==2 and yes_44==1) {c_5_44++;}
if(yes==1 and yes_44==1) c_5_44++;
if(yes==2 and yes_44==0) c_5++;

if(yes==0 and yes_44==1) {c_44++; h_Em_no->Fill(e44);h_thm_no->Fill(the44_gen);	}

cout << "------------------" << endl;

if(yes==2){// and yes_44==1){//if(ok_em==1 or ok_ep==1){

if(st0==1) n_st0++;
if(st1==1) n_st1++;
bkg++;
gen++;

//h_Em->Fill(e44);
h_Ee->Fill(e_e);
h_Ep->Fill(e_p);
h_the->Fill(them_gen);
h_thp->Fill(thep_gen);
//h_thm->Fill(the44_gen);



 int yes_em=0; int yes_ep=0;  int yes_mu_in=0; int yes_mu_44=0;

TVector3 pem_rec,pep_rec, pmu_in_rec, in,pmu_44_rec;
double them_rec,thep_rec, thmu_in_rec,thmu_44_rec;
double chi_em,chi_ep, chi_in,chi_44;

vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
for(int j=0; j<tracks.size();j++)
{
double th_inx,th_iny;
if(tracks.at(j).processIDofLinkedTrack()==0 and tracks.at(j).sector()==0){
th_inx=tracks.at(j).xSlope();
th_iny=tracks.at(j).ySlope();
in.SetXYZ(th_inx,th_iny,1.0);
in=in.Unit();
 }
}

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

if(tracks.at(j).processIDofLinkedTrack()==0 and tracks.at(j).sector()==1) {yes_mu_44=1; TVector3 p; p.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0);
									    pmu_44_rec=p.Unit(); thmu_44_rec=acos(pmu_44_rec.Dot(in)); //th_in44_r->Fill(p.Theta()); th_in44->Fill(the44_gen);
									    chi_44=tracks.at(j).chi2perDegreeOfFreedom();}

if(tracks.at(j).linkedTrackID()==code_em and tracks.at(j).sector()==1) {yes_em=1; TVector3 p; p.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.); pem_rec=p.Unit();
									   them_rec=acos(in.Dot(pem_rec)); chi_em=tracks.at(j).chi2perDegreeOfFreedom();}
if(tracks.at(j).linkedTrackID()==code_ep and tracks.at(j).sector()==1) {yes_ep=1; TVector3 p; p.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.); pep_rec=p.Unit();
									   thep_rec=acos(in.Dot(pep_rec)); chi_ep=tracks.at(j).chi2perDegreeOfFreedom();}

	}
}

//if(yes_mu_in==1 and yes_mu_44==1) cout << "contemporaneamente ID0 e ID44 nella seconda stazione" << endl;


double ppp=0;
double ppm=0;
vector<MUonERecoOutputAdaptiveFitterVertex> vrtx = ReconstructionOutput->adaptiveFitterVertices();
vector<MUonERecoOutputVertex> vrtx_kin = ReconstructionOutput->reconstructedVertices();

h_size->Fill(tracks.size()-1);
if(tracks.size()==3 and perc.size()==3) {bkg3_g++;}
if(tracks.size()==4 and perc.size()==4) {bkg4_g++;}
if(tracks.size()>4) {more_tracks_g++;}




if(tracks.size()==3){
if(yes_mu_44==1 and yes_em==1) {bkg3++; ppm=1; h_part2_bool->Fill(1); if(vrtx.size()!=0) ada_PP++; //opening->Fill(acos(pem_rec.Dot(pmu_44_rec)));
				th_e->Fill(them_gen); th_e_r->Fill(them_rec); 
				 h_chi_em_3->Fill(chi_em); h_chi_44_3->Fill(chi_44); 
				th_in44_r->Fill(thmu_44_rec); th_in44->Fill(the44_gen);
for(int j=0; j<tracks.size();j++){ h_quality2_PP->Fill(tracks.at(j).percentageOfHitsSharedWithLinkedTrack());}
}
else if(yes_mu_44==1 and yes_ep==1) {bkg3++; ppm=1; h_part2_bool->Fill(1);if(vrtx.size()!=0) ada_PP++; //opening->Fill(acos(pep_rec.Dot(pmu_44_rec)));
				 h_chi_ep_3->Fill(chi_ep); h_chi_44_3->Fill(chi_44);
				th_e->Fill(thep_gen); th_e_r->Fill(thep_rec);
				th_in44_r->Fill(thmu_44_rec); th_in44->Fill(the44_gen);
for(int j=0; j<tracks.size();j++){ h_quality2_PP->Fill(tracks.at(j).percentageOfHitsSharedWithLinkedTrack());}
}
else{h_part2_bool->Fill(0); bkg3no++;
				for(int j=0; j<tracks.size();j++)
                                {if(tracks.at(j).sector()==1){h_part2_mue->Fill(tracks.at(j).processIDofLinkedTrack()); h_quality2_no->Fill(tracks.at(j).percentageOfHitsSharedWithLinkedTrack()); }
				}
    }
 }

if(tracks.size()==4){
if(yes_mu_44==1 and yes_em==1 and yes_ep==1) {bkg4++; ppp=1; if(vrtx.size()!=0) ada_PP++; //opening->Fill(acos(pem_rec.Dot(pmu_44_rec)));
		                                                   h_chi_em_4->Fill(chi_em); h_chi_ep_4->Fill(chi_ep); h_chi_44_4->Fill(chi_44);
							th_e->Fill(thep_gen); th_e_r->Fill(thep_rec);
							th_e_r->Fill(them_rec); th_e->Fill(them_gen);
							th_in44_r->Fill(thmu_44_rec); th_in44->Fill(the44_gen);
for(int j=0; j<tracks.size();j++){ h_quality3_PP->Fill(tracks.at(j).percentageOfHitsSharedWithLinkedTrack());}
}
else{bkg4no++;
				for(int j=0; j<tracks.size();j++)
                                {if(tracks.at(j).sector()==1){h_part3_mue->Fill(tracks.at(j).processIDofLinkedTrack());h_quality3_no->Fill(tracks.at(j).percentageOfHitsSharedWithLinkedTrack());  }
 	} 			}
}


if(ppm==1 or ppp==1){
for(int j=0; j<vrtx_kin.size();j++)
{
if(vrtx_kin.at(j).stationIndex()==1)
{
 MUonERecoOutputTrack mu_in = vrtx_kin.at(j).incomingMuon();
 MUonERecoOutputTrack mu_out = vrtx_kin.at(j).outgoingMuon();
 MUonERecoOutputTrack e_out = vrtx_kin.at(j).outgoingElectron();

TVector3 in1(mu_in.xSlope(),mu_in.ySlope(),1.0);
TVector3 mu1(mu_out.xSlope(),mu_out.ySlope(),1.0);
TVector3 e1(e_out.xSlope(),e_out.ySlope(),1.0);
in1=in1.Unit();
mu1=mu1.Unit();
e1=e1.Unit();

}

}


}


if(ppm==0 and ppp==0){more_tracks++;}

	}
} //end of general for
cout << "Numero di eventi generati PP su N eventi input: " << gen << " ---> " << (gen/cbmsim->GetEntries())*100 << endl;
cout << "Su " << bkg << " eventi PP, " << bkg3_g << " molteplicita' = 2: " << (bkg3_g/bkg)*100 << endl;
cout << "Su " << bkg << " eventi PP, " << bkg4_g << " molteplicita' = 3: " << (bkg4_g/bkg)*100 << endl;
cout << "Su " << bkg << " eventi PP, " << more_tracks_g << "molteplicita' > 3: " << (more_tracks_g/bkg)*100 << endl;


cout << "Su " << bkg << " eventi PP, " << bkg3 << " sono muone+elettrone/positrone da PP: " << (bkg3/bkg)*100 << endl;
cout << "Su " << bkg << " eventi PP, " << bkg4 << " sono muone+elettrone+positrone da PP: " << (bkg4/bkg)*100 << endl;
cout << "Su " << bkg << " eventi PP, " << bkg3no << " hanno molteplicita 2 ma NON sono muone+elettrone/positrone da PP: " << (bkg3no/bkg)*100 << endl;
cout << "Su " << bkg << " eventi PP, " << bkg4no << " hanno molteplicita 3 ma NON sono muone+elettrone+positrone da PP: " << (bkg4no/bkg)*100 << endl;
cout << "Su " << bkg << " eventi PP, " << more_tracks << " sono senza requisiti minimi: " << (more_tracks/bkg)*100 << endl;
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
a5.SaveAs("idPerMultiplicity_PP.pdf");

TCanvas a6("a6","a6",700,700);
a6.Divide(2,3);
a6.cd(1);
h_chi_em_3->SetLineColor(kOrange);
h_chi_em_3->Draw("hist");
h_chi_em_4->Draw("hist same");
gPad->SetLogy();
a6.cd(2);
h_chi_ep_3->SetLineColor(kOrange);
h_chi_ep_3->Draw("hist");
h_chi_ep_4->Draw("hist same");
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
a6.SaveAs("chi2_PP.pdf");

--------
*/
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
a7.SaveAs("caratt_44.pdf");

/*
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
a8.SaveAs("quality.pdf");

TCanvas a9("a9","a9",700,700);
a9.Divide(1,2);
a9.cd(1);
h_nstubs->Draw("hist");
a9.cd(2);
h_nstubs1->Draw("hist");
a9.SaveAs("nstubs.pdf");


TCanvas b1("b1","b1",700,700);
b1.Divide(1,2);
b1.cd(1);
th_in44->Draw("hist");
th_in44_r->SetLineColor(kViolet);
th_in44_r->Draw("hist same");
//gPad->SetLogy();
b1.cd(2);
th_e->Draw("hist");
th_e_r->SetLineColor(kViolet);
th_e_r->Draw("hist same");
b1.SaveAs("angle44.pdf");

TCanvas b2("b2","b2",700,700);
b2.Divide(1,3);
b2.cd(1);
h_E1->Draw("hist");
gPad->SetLogy();
b2.cd(2);
//h_E2->SetLineColor(kGreen);
h_E2->Draw("hist");
gPad->SetLogy();
b2.cd(3);
h_E3->Draw("hist");
b2.SaveAs("energy.pdf");
*/

}
