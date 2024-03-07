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

void minbias(){

TChain * cbmsim = new TChain("cbmsim");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/minbias_1M.root");
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/minbias_1M_2.root");


        TClonesArray *MCTrack = 0;
        TClonesArray *TrackerStripDigis = 0;
        TClonesArray *TrackerPoints = 0;
        TClonesArray *TrackerStubs = 0;
        MUonERecoOutput *ReconstructionOutput = 0;

        cbmsim->SetBranchAddress("MCTrack", &MCTrack);
        cbmsim->SetBranchAddress("TrackerPoints", &TrackerPoints);
        cbmsim->SetBranchAddress("TrackerStripDigis", &TrackerStripDigis);
        cbmsim->SetBranchAddress("TrackerStubs", &TrackerStubs);
        cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);


                auto pos_on_track = [](double q, double m, double z){return (q + m*z);};

TH1D::SetDefaultSumw2(kTRUE);

TH1D *multiplicity_st0_g = new TH1D("multiplicity_st0_g", "number of station0 gen tracks", 20,0,20);
TH1D *multiplicity_st0_r = new TH1D("multiplicity_st0_r", "number of station0 reco tracks", 20,0,20);

TH1D *multiplicity_st1_g = new TH1D("multiplicity_st1_g", "number of station1 gen tracks", 20,0,20);
TH1D *multiplicity_st1_r = new TH1D("multiplicity_st1_r", "number of station1 reco tracks", 20,0,20);

double n=0.;
double n_0=0.;
double n_1=0.;
double n_2=0.;
double n_3=0.;
double n_4=0.;

double n_g=0.;
double n_0_g=0.;
double n_1_g=0.;
double n_2_g=0.;
double n_3_g=0.;
double n_4_g=0.;

for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
                cbmsim->GetEntry(i);
                if(i%100 == 0) cout<<"Entry "<<i<<endl;

int mu_gen=0;

        for(int n = 0; n < MCTrack->GetEntries(); n++) {
         const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
         if(MCTr->interactionID()==0 and MCTr->pdgCode()==13)mu_gen++;
	}

if(mu_gen==0) n_0_g++;
if(mu_gen==1) n_1_g++;
if(mu_gen==2) n_2_g++;
if(mu_gen==3) n_3_g++;
if(mu_gen>3) n_4_g++;
/*std::array<int,mu_gen> hit_modX(0);
std::array<int,mu_gen> hit_modY(0);
std::array<int,mu_gen> stereo(0);

int index=0;

         for(int n = 0; n < MCTrack->GetEntries(); n++) {
         const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
         if(MCTr->interactionID()==0 and MCTr->pdgCode()==13){
			for(int s=0; s<TrackerPoints->GetEntries(); s++)
                         {const MUonETrackerPoint *TrackerPt = static_cast<const MUonETrackerPoint*>(TrackerPoints->At(s));
                          if(TrackerPt->trackPDGCode()==13 and TrackerPt->trackID()==n and TrackerPt->stationID()==0){
                                                                                                 if(TrackerPt->moduleID()==0 or TrackerPt->moduleID()==4) hit_modX.at(index)++;
                                                                                                 if(TrackerPt->moduleID()==1 or TrackerPt->moduleID()==5) hit_modY.at(index)++;
                                                                                                 if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo.at(index)++;}
                 }
		index++;
	    }
	}*/

	if(mu_gen!=0){// and std::all_of(v.cbegin(), v.cend(), [](int i) { return i  == 4; }) hit_modX%4==0 and hit_modY%4==0 and stereo>=2){

multiplicity_st0_g->Fill(mu_gen);
n++;
int stub0=0; int stub1=0;

        for(int t=0; t<TrackerStubs->GetEntries(); t++)
                         {const MUonETrackerStub *stubs = static_cast<const MUonETrackerStub*>(TrackerStubs->At(t));
			   if(stubs->stationID()==0){stub0++;}
                           if(stubs->stationID()==1){stub1++;}
			}

//cout << "stub0 " << stub0 << " and stub1 " << stub1 << endl;

vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
int mu_in=0;

         for (auto&& track : tracks) {
        if(track.processIDofLinkedTrack()==0 and track.sector()==0)mu_in++;
         }//for

if(mu_in==0) n_0++;
if(mu_in==1) n_1++;
if(mu_in==2) n_2++;
if(mu_in==3) n_3++;
if(mu_in>3) n_4++;

multiplicity_st0_r->Fill(mu_in);

        }//if generated

}//for

cout <<"N. of events with 0 muons generated in st_0 " << n_0_g << " -> " << n_0_g/cbmsim->GetEntries()*100 << "% of events" << endl;
cout <<"N. of events with 1 muons generated in st_0 " << n_1_g << " -> " << n_1_g/cbmsim->GetEntries()*100 << "% of events" << endl;
cout <<"N. of events with 2 muons generated in st_0 " << n_2_g << " -> " << n_2_g/cbmsim->GetEntries()*100 << "% of events" << endl;
cout <<"N. of events with ==3 muons generated in st_0 " << n_3_g << " -> " << n_3_g/cbmsim->GetEntries()*100 << "% of events" << endl;
cout <<"N. of events with >3 muons generated in st_0 " << n_4_g << " -> " << n_4_g/cbmsim->GetEntries()*100 << "% of events" << endl;


cout <<"N. of events with 0 muons reconstructed in st_0 " << n_0 << " -> " << n_0/n*100 << "% of events with at least 1 generated muons" << endl;
cout <<"N. of events with 1 muons reconstructed in st_0 " << n_1 << " -> " << n_1/n*100 << "% of events with at least 1 generated muons" << endl;
cout <<"N. of events with 2 muons reconstructed in st_0 " << n_2 << " -> " << n_2/n*100 << "% of events with at least 1 generated muons" << endl;
cout <<"N. of events with ==3 muons reconstructed in st_0 " << n_3 << " -> " << n_3/n*100 << "% of events with at least 1 generated muons" << endl;
cout <<"N. of events with >3 muons reconstructed in st_0 " << n_4 << " -> " << n_4/n*100 << "% of events with at least 1 generated muons" << endl;



 auto legend = new TLegend(0.7,0.7,0.9,0.9);
   legend->AddEntry(multiplicity_st0_g,"multiplicity_st0_g","LEP");
   legend->AddEntry(multiplicity_st0_r,"multiplicity_st0_r","LEP");

TCanvas b("b","b",700,700);
multiplicity_st0_g->SetLineColor(kOrange);
multiplicity_st0_r->SetLineColor(kRed);

multiplicity_st0_g->Draw("hist");
multiplicity_st0_r->Draw("hist same");
legend->Draw();
b.SaveAs("minbias_pdf/mult.pdf");


}
