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

void efficiency_MC(int nhits, int index){

TChain * cbmsim = new TChain("cbmsim");


if(nhits==0){
  if(index==0)cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_0-5mrad_%dhit_NOoutchi2_1M.root",static_cast<char>(nhits)));
  else if(index==1)cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_5-10mrad_%dhit_NOoutchi2_1M.root",static_cast<char>(nhits)));
  else if(index==2)cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_10-15mrad_%dhit_NOoutchi2_1M.root",static_cast<char>(nhits)));
  else if(index==3)cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_15-20mrad_%dhit_NOoutchi2_1M.root",static_cast<char>(nhits)));
  else if(index==4)cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_20-25mrad_%dhit_NOoutchi2_1M.root",static_cast<char>(nhits)));
  else if(index==5)cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_25-32mrad_%dhit_NOoutchi2_1M.root",static_cast<char>(nhits)));
  else cout << "WRONG NUMBER" << endl;
}
else{
  if(index==0)cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_0-5mrad_1M_%dhitFirstModules_NOoutchi2_1M.root",static_cast<char>(nhits)));
  else if(index==1)cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_5-10mrad_1M_%dhitFirstModules_NOoutchi2_1M.root",static_cast<char>(nhits)));
  else if(index==2)cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_10-15mrad_1M_%dhitFirstModules_NOoutchi2_1M.root",static_cast<char>(nhits)));
  else if(index==3)cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_15-20mrad_1M_%dhitFirstModules_NOoutchi2_1M.root",static_cast<char>(nhits)));
  else if(index==4)cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_20-25mrad_1M_%dhitFirstModules_NOoutchi2_1M.root",static_cast<char>(nhits)));
  else if(index==5)cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_25-32mrad_1M_%dhitFirstModules_NOoutchi2_1M.root",static_cast<char>(nhits)));
  else cout << "WRONG NUMBER" << endl;
}


        TClonesArray *MCTrack = 0;
        TClonesArray *TrackerStripDigis = 0;
        TClonesArray *TrackerPoints = 0;
        TClonesArray *TrackerStubs = 0;
        MuE::Event *MesmerEvent = 0;
        MUonERecoOutput *ReconstructionOutput = 0;

        cbmsim->SetBranchAddress("MCTrack", &MCTrack);
        cbmsim->SetBranchAddress("TrackerPoints", &TrackerPoints);
        cbmsim->SetBranchAddress("TrackerStripDigis", &TrackerStripDigis);
        cbmsim->SetBranchAddress("TrackerStubs", &TrackerStubs);
        cbmsim->SetBranchAddress("MesmerEvent", &MesmerEvent);
        cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);


                auto pos_on_track = [](double q, double m, double z){return (q + m*z);};
   const Int_t NBINS = 14;
   const Int_t NBINS_mu = 14;
   Double_t edges_el[NBINS + 1] = {0.,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.010, 0.015, 0.020, 0.025, 0.032};
   Double_t edges_mu[NBINS_mu + 1] = {0.,0.0001,0.0002,0.0003,0.0004,0.0005,0.0006,0.0007,0.0008,0.0009,0.001,0.002,0.003,0.004,0.005};

TH1D::SetDefaultSumw2(kTRUE);

TH1D *theta_e=new TH1D("theta_e", "Reco elastic event: Electron scattering reco angles from MESMER",NBINS,edges_el);
TH1D *theta_mu=new TH1D("theta_mu", "Reco elastic event: Muon scattering reco angles from MESMER",NBINS_mu,edges_mu);
TH1D *h_opening=new TH1D("h_opening", "Opening angle reco events from MESMER",35,0.,0.035);

TH1D *theta_e_clone=new TH1D("theta_e_clone", "Reco elastic event: Electron scattering reco angles from MESMER with clones",NBINS,edges_el);
TH1D *theta_mu_clone=new TH1D("theta_mu_clone", "Reco elastic event: Muon scattering reco angles from MESMER with clones",NBINS_mu,edges_mu);
TH1D *h_opening_clone=new TH1D("h_opening_clone", "Opening angle reco events from MESMER with clones",35,0.,0.035);


TH1D *theta_e_single=new TH1D("theta_e_single", "Reco single electron: Electron scattering reco angles from MESMER",NBINS,edges_el);
TH1D *theta_mu_single=new TH1D("theta_mu_single", "Reco single muon: Muon scattering reco angles from MESMER",NBINS_mu,edges_mu);

TH1D *theta_e_single_clone=new TH1D("theta_e_single_clone", "Reco single electron: Electron scattering reco angles from MESMER clone added",NBINS,edges_el);
TH1D *theta_mu_single_clone=new TH1D("theta_mu_single_clone", "Reco single muon: Muon scattering reco angles from MESMER clone added",NBINS_mu,edges_mu); 


TH1D *theta_e_gen=new TH1D("theta_e_gen", "Electron scattering generated angles from MESMER",NBINS,edges_el);
TH1D *theta_mu_gen=new TH1D("theta_mu_gen", "Muon scattering generated angles from MESMER",NBINS_mu,edges_mu);
TH1D *h_opening_gen=new TH1D("h_opening_gen", "Opening angle generated events from MESMER",35,0.,0.035);

TH1D *theta_e_noreco=new TH1D("theta_e_noreco", "Electron scattering angles from MESMER when NOT reco",NBINS,edges_el);
TH1D *theta_mu_noreco=new TH1D("theta_mu_noreco", "Muon scattering angle from MESMER when NOT reco",NBINS_mu,edges_mu);


	double r_wnorm[6]={20.786765103274643,33.313091221576336,42.396733790329876,50.584815206143652,61.828110400735824,106.88513370134392};



double n_clones_el=0.;
double n_clones_mu=0.;
double total=0.;
double all=0.;
double n_el2=0.;
double n_el3=0.;
double n_cl=0.;
double n_one_el=0.;
double n_other=0.;
double n_both=0.;
double tot=0.;

for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
                cbmsim->GetEntry(i);
                if(i%100 == 0) cout<<"Entry "<<i<<endl;


Double_t code_mu_in=-99;
Double_t code_e=-99;
Double_t code_mu=-99;
        TVector3 p_muin_MC;
        TVector3 p_mu_MC;
        TVector3 p_e_MC;
        Double_t the_gen, thmu_gen,theX_gen,theY_gen,thmuX_gen,thmuY_gen;

           int hit_modXmu=0; int hit_modXe=0;
           int hit_modYmu=0; int hit_modYe=0;
           int stereo_mu=0; int stereo_e=0;
	   int hit_modXmuin=0;
	   int hit_modYmuin=0;
	   int stereo_muin=0;
	   double E_e=0.;

// Checking if in the MCTracks container there are elastic particles and if they are in acceptance (4 hits in X modules (== 2 stubs, as TrackerPoints give the hit per sensor) + 4 hits in Y modules + at least 2 hits in UV)

        for(int n = 0; n < MCTrack->GetEntries(); n++) {
         const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));

         if(MCTr->interactionID()==0 and MCTr->pdgCode()==-13) {code_mu_in=n; p_muin_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); p_muin_MC.Unit();E_e=MCTr->energy();
        for(int s=0; s<TrackerPoints->GetEntries(); s++)
                         {const MUonETrackerPoint *TrackerPt = static_cast<const MUonETrackerPoint*>(TrackerPoints->At(s));
         if(TrackerPt->trackPDGCode()==-13 and TrackerPt->trackID()==n and TrackerPt->stationID()==0){
                                                                                                 if(TrackerPt->moduleID()==0 or TrackerPt->moduleID()==4) hit_modXmuin++;
                                                                                                 else if(TrackerPt->moduleID()==1 or TrackerPt->moduleID()==5) hit_modYmuin++;
                                                                                                 else if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_muin++;}
                }
	}

         if(MCTr->interactionID()==45 and MCTr->pdgCode()==11) {code_e=n; p_e_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); p_e_MC.Unit(); the_gen=p_muin_MC.Angle(p_e_MC);
                                                                                 theX_gen=MCTr->ax();
                                                                                 theY_gen=MCTr->ay();

        for(int s=0; s<TrackerPoints->GetEntries(); s++)
                         {const MUonETrackerPoint *TrackerPt = static_cast<const MUonETrackerPoint*>(TrackerPoints->At(s));
         if(TrackerPt->trackPDGCode()==11 and TrackerPt->trackID()==code_e and TrackerPt->stationID()==1){
                                                                                                 if(TrackerPt->moduleID()==0 or TrackerPt->moduleID()==4) hit_modXe++;
                                                                                                 else if(TrackerPt->moduleID()==1 or TrackerPt->moduleID()==5) hit_modYe++;
                                                                                                 else if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_e++;}
        	}
	 }

         if(MCTr->interactionID()==45 and MCTr->pdgCode()==-13) {code_mu=n; p_mu_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); p_mu_MC.Unit(); thmu_gen=p_muin_MC.Angle(p_mu_MC);
                                                                                 thmuX_gen=MCTr->ax();
                                                                                 thmuY_gen=MCTr->ay();
        for(int s=0; s<TrackerPoints->GetEntries(); s++)
                         {const MUonETrackerPoint *TrackerPt = static_cast<const MUonETrackerPoint*>(TrackerPoints->At(s));
	 if(TrackerPt->trackPDGCode()==-13 and TrackerPt->trackID()==code_mu and TrackerPt->stationID()==1){
                                                                                                 if(TrackerPt->moduleID()==0 or TrackerPt->moduleID()==4) hit_modXmu++;
                                                                                                 else if(TrackerPt->moduleID()==1 or TrackerPt->moduleID()==5) hit_modYmu++;
                                                                                                 else if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_mu++;}
		}
	 }
	}


// Look at reconstruction if events are reconstructible (all three particles with necessary hits to be potentially reconstructed)

if(code_mu_in!=-99 and code_e!=-99 and code_mu!=-99 and hit_modXmuin==4 and hit_modYmuin==4 and stereo_muin>1 and hit_modXe==4 and hit_modYe==4 and stereo_e>1 and hit_modXmu==4 and hit_modYmu==4 and stereo_mu>1){

// wnorm is a number needed for normalization when we use different mesmer sample together (like in this case, 6 sample in differen kinematic region of the electron)

double wnorm=99.;

if(index==0){wnorm=r_wnorm[0];}
else if(index==1){wnorm=r_wnorm[1];}
else if(index==2){wnorm=r_wnorm[2];}
else if(index==3){wnorm=r_wnorm[3];}
else if(index==4){wnorm=r_wnorm[4];}
else if(index==5){wnorm=r_wnorm[5];}

 all+=MesmerEvent->wgt_LO*wnorm;


double opening_angle=p_mu_MC.Angle(p_e_MC);


Int_t yes_mu=0;
Int_t yes_e=0;
Int_t yes2=0;
Double_t th_inx,th_iny,x0_in,y0_in;
Double_t chi2_muin=999.;
Int_t stubs_muin=0;
Int_t e=0;
Int_t mu=0;
TVector3 p_muin,p_e,p_mu;
Double_t the_rec,theX_rec,theY_rec,thmu_rec,thmuX_rec,thmuY_rec;

// Define best vertex and its chi2
MUonERecoOutputVertex vrtx = ReconstructionOutput->bestVertex();
Double_t chi=vrtx.chi2perDegreeOfFreedom();


// Look at reconstructed track container to see if outgoing mu and e are reconstructed
vector<double> quality_e; quality_e.reserve(5);
vector<double> quality_mu; quality_mu.reserve(5);

vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
int other=0;
int mu_in=0;
         for (auto&& track : tracks) {

        if(code_mu_in==track.linkedTrackID() and track.sector()==0){
	mu_in++;
         th_inx=track.xSlope();
         th_iny=track.ySlope();
         x0_in=track.x0();
         y0_in=track.y0();
         chi2_muin=track.chi2perDegreeOfFreedom();
         stubs_muin=track.hits().size();
        p_muin.SetXYZ(th_inx,th_iny,1.0);
        p_muin=p_muin.Unit();
                        }
        if(track.processIDofLinkedTrack()==45 and track.sector()==1)
                {
                 if(code_e==track.linkedTrackID()) {yes2++; e++;
						theX_rec=track.xSlope(); theY_rec=track.ySlope();
						p_e.SetXYZ(theX_rec,theY_rec,1.0); p_e=p_e.Unit();the_rec=p_e.Angle(p_muin);
						quality_e.push_back(track.fractionOfHitsSharedWithLinkedTrack());}
                 if(code_mu==track.linkedTrackID()) {yes2++; mu++;
						thmuX_rec=track.xSlope(); thmuY_rec=track.ySlope();
						p_mu.SetXYZ(thmuX_rec,thmuY_rec,1.0); p_mu=p_mu.Unit();thmu_rec=p_mu.Angle(p_muin);
						quality_mu.push_back(track.fractionOfHitsSharedWithLinkedTrack());}
                }
        if(track.processIDofLinkedTrack()!=45 and track.sector()==1){other++;}

         }//for


// Check that the incoming muon is well reconstructed
if(stubs_muin>=5 and chi2_muin<2){//and thmu_gen>0.0002){

tot+=MesmerEvent->wgt_LO*wnorm;

theta_mu_gen->Fill(thmu_gen,MesmerEvent->wgt_LO*wnorm);
theta_e_gen->Fill(the_gen,MesmerEvent->wgt_LO*wnorm);
h_opening_gen->Fill(opening_angle,MesmerEvent->wgt_LO*wnorm);

// Qulity checks for the evaluation of single mu and single el reconstrustion efficiency
if(find_if(quality_mu.begin(),quality_mu.end(),[](double i){return i>=0.65;})!=end(quality_mu)){
 if(mu==1){theta_mu_single->Fill(thmu_gen,MesmerEvent->wgt_LO*wnorm);}
 if(mu>=1){theta_mu_single_clone->Fill(thmu_gen,MesmerEvent->wgt_LO*wnorm);}
}

if(find_if(quality_e.begin(),quality_e.end(),[](double i){return i>=0.65;})!=end(quality_e)){
 if(e==1){theta_e_single->Fill(the_gen,MesmerEvent->wgt_LO*wnorm);}
 if(e>=1){theta_e_single_clone->Fill(the_gen,MesmerEvent->wgt_LO*wnorm);}
}



if(yes2>=0){total+=MesmerEvent->wgt_LO*wnorm;}


// Qulity checks for the evaluation of elastic events reconstrustion efficiency
if(find_if(quality_e.begin(),quality_e.end(),[](double i){return i>=0.65;})!=end(quality_e) and find_if(quality_mu.begin(),quality_mu.end(),[](double i){return i>=0.65;})!=end(quality_mu)){

 //some number for multiplicity and what is reconstructed other==particle different from elastic lepton
 if(e==1 and mu==1 and other==0) n_el2+=MesmerEvent->wgt_LO*wnorm;
 if(e==1 and mu==1 and other!=0) n_el3+=MesmerEvent->wgt_LO*wnorm;
 if( (e>1 or mu>1) and other==0) n_cl+=MesmerEvent->wgt_LO*wnorm;
 if(yes2==0 and other>0) n_other+=MesmerEvent->wgt_LO*wnorm;
 if(yes2==1 and other==0) n_one_el+=MesmerEvent->wgt_LO*wnorm;
 if( (e>1 or mu>1) and other>0) n_both+=MesmerEvent->wgt_LO*wnorm;
 if(chi!=0 and yes2>=2 and e>1)n_clones_el+=MesmerEvent->wgt_LO*wnorm;
 if(chi!=0 and yes2>=2 and mu>1)n_clones_mu+=MesmerEvent->wgt_LO*wnorm;


// reconstruction efficiency for elastic events without clones of elastic particles
 if(e==1 and mu==1 and chi!=0){
 theta_mu->Fill(thmu_gen,MesmerEvent->wgt_LO*wnorm);
 theta_e->Fill(the_gen,MesmerEvent->wgt_LO*wnorm);
 h_opening->Fill(opening_angle,MesmerEvent->wgt_LO*wnorm);
                 }//if chi!=0

// reconstruction efficiency for elastic events including clones of elastic particles
 if(e>=1 and mu>=1 and chi!=0){
 theta_mu_clone->Fill(thmu_gen,MesmerEvent->wgt_LO*wnorm);
 theta_e_clone->Fill(the_gen,MesmerEvent->wgt_LO*wnorm);
 h_opening_clone->Fill(opening_angle,MesmerEvent->wgt_LO*wnorm);
                 }
			}//quality
                }//if mu_in

        }//if generated

}//for

cout << "N. of elastic events with 1 muon and 1 electron and no other "<< n_el2 << " -> " << n_el2/tot*100 << "%"<< endl;
cout << "N. of elastic events with 1 muon and 1 electron and other "<< n_el3 << " -> " << n_el3/tot*100 << "%"<< endl;
cout << "N. of elastic events with mu/e clones "<< n_cl << " -> " << n_cl/tot*100 << "%"<< endl;
cout << "N. of events with just 1 elastic "<< n_one_el << " -> " << n_one_el/tot*100 << "%"<< endl;
cout << "N. of elastic events with just other (interaction ID !=45) particles "<< n_other << " -> " << n_other/tot*100 << "%"<< endl;
cout << "N. of events with other (interaction ID !=45) particles and clones "<< n_both << " -> " << n_both/tot*100 << "%"<< endl;


cout <<"N. of electron clones " << n_clones_el << " over events with >= 2 e+mu " << total << " -> " << n_clones_el/total*100 << "%" << endl;
cout <<"N. of muon clones " << n_clones_mu << " over events with >= 2 e+mu " << total << " -> " << n_clones_mu/total*100 << "%" << endl;



TH1D * h0 = (TH1D*) theta_e_single->Clone();
TH1D * h0gen = (TH1D*) theta_e_gen->Clone();
h0->Divide(h0,h0gen,1,1,"B");
TH1D * h1 = (TH1D*) theta_mu_single->Clone();
TH1D * h1gen = (TH1D*) theta_mu_gen->Clone();
h1->Divide(h1,h1gen,1,1,"B");

TH1D * h0_clone = (TH1D*) theta_e_single_clone->Clone();
h0_clone->Divide(h0_clone,h0gen,1,1,"B");
TH1D * h1_clone = (TH1D*) theta_mu_single_clone->Clone();
h1_clone->Divide(h1_clone,h1gen,1,1,"B");



TH1D * h2 = (TH1D*) theta_e->Clone();
TH1D * h2gen = (TH1D*) theta_e_gen->Clone();
h2->Divide(h2,h2gen,1,1,"B");
TH1D * h3 = (TH1D*) theta_mu->Clone();
TH1D * h3gen = (TH1D*) theta_mu_gen->Clone();
h3->Divide(h3,h3gen,1,1,"B");

TH1D * h4 = (TH1D*) h_opening->Clone();
TH1D * h4gen = (TH1D*) h_opening_gen->Clone();
h4->Divide(h4,h4gen,1,1,"B");



TH1D * h2_clone = (TH1D*) theta_e_clone->Clone();
h2_clone->Divide(h2_clone,h2gen,1,1,"B");
TH1D * h3_clone = (TH1D*) theta_mu_clone->Clone();
h3_clone->Divide(h3_clone,h3gen,1,1,"B");

TH1D * h4_clone = (TH1D*) h_opening_clone->Clone();
h4_clone->Divide(h4_clone,h4gen,1,1,"B");

cout << "Efficiency as function of muon angle" << endl;
for(int b=1; b<h3->GetNbinsX(); b++){
cout << b << " ) " << h3->GetBinContent(b) << " +- " << h3->GetBinError(b)<< " ----> " << h3->GetBinError(b)/h3->GetBinContent(b)*100 << "%" << endl;
}
cout << "Efficiency as function of electron angle" << endl;
for(int b=1; b<h2->GetNbinsX(); b++){
cout << b << " ) " << h2->GetBinContent(b) << " +- " << h2->GetBinError(b)<< " ---> " << h2->GetBinError(b)/h2->GetBinContent(b)*100 << "%" <<endl ;
}


cout << "Efficiency as function of muon angle including clones" << endl;
for(int b=1; b<h3_clone->GetNbinsX(); b++){
cout << b << " ) " << h3_clone->GetBinContent(b) << " +- " << h3_clone->GetBinError(b)<< " ----> " << h3_clone->GetBinError(b)/h3_clone->GetBinContent(b)*100 << "%" << endl;
}
cout << "Efficiency as function of electron angle including clones" << endl;
for(int b=1; b<h2_clone->GetNbinsX(); b++){
cout << b << " ) " << h2_clone->GetBinContent(b) << " +- " << h2_clone->GetBinError(b)<< " ---> " << h2_clone->GetBinError(b)/h2_clone->GetBinContent(b)*100 << "%" <<endl ;
}



TCanvas a1("a1","a1",700,700);
a1.Divide(2,3);
a1.cd(1);
h1->SetMinimum(0.92);
h1->SetMaximum(1.02);
h1->Draw("E");
h1_clone->Draw("E same");
gStyle->SetOptStat(0);
a1.cd(2);
h0->SetMinimum(0.65);
h0->SetMaximum(1.02);
h0->Draw("E");
h0_clone->Draw("E same");
a1.cd(3);
h3->SetMinimum(0.82);
h3->SetMaximum(1.02);
h3_clone->SetLineColor(kPink);
h3->Draw("E");
h3_clone->Draw("E same");
gStyle->SetOptStat(0);
a1.cd(4);
h2->SetMinimum(0.65);
h2->SetMaximum(1.02);
h2->Draw("E");
h2_clone->SetLineColor(kPink);
h2_clone->Draw("E same");
gStyle->SetOptStat(0);
a1.cd(5);
h4->SetMinimum(0.5);
h4->SetMaximum(1.02);
h4->Draw("E");
h4_clone->SetLineColor(kPink);
h4_clone->Draw("E same");
gStyle->SetOptStat(0);
a1.SaveAs(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/eff_LO_%dhit_%d.pdf",static_cast<char>(nhits),static_cast<char>(index)));


h0_clone->SaveAs(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/el_single_eff_LO_%dhit_%d.root",static_cast<char>(nhits),static_cast<char>(index)));
h1_clone->SaveAs(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/mu_single_eff_LO_%dhit_%d.root",static_cast<char>(nhits),static_cast<char>(index)));
h2_clone->SaveAs(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/el_eff_LO_%dhit_%d.root",static_cast<char>(nhits),static_cast<char>(index)));
h3_clone->SaveAs(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/mu_eff_LO_%dhit_%d.root",static_cast<char>(nhits),static_cast<char>(index)));
h4_clone->SaveAs(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/op_eff_LO_%dhit_%d.root",static_cast<char>(nhits),static_cast<char>(index)));

}

