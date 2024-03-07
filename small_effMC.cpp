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

void effMC(){

TChain * cbmsim = new TChain("cbmsim");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_0-5mrad_1M_1hit_300outchi2.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_5-10mrad_1M_1hit_300outchi2.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_10-15mrad_1M_1hit_300outchi2.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_15-20mrad_1M_1hit_300outchi2.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_20-25mrad_1M_1hit_300outchi2.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_25-32mrad_1M_1hit_300outchi2.root");

//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_0-32mrad_1M_1hit.root");
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_20-32mrad_300k_1hit.root");
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_20-32mrad_300k_1hit.root");

//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_20-32mrad_1M_1hit_1000outchi2.root");

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
   const Int_t NBINS = 6;
   Double_t edges[NBINS + 1] = {0.0, 0.005, 0.010, 0.015, 0.020, 0.025, 0.032};

TH1D::SetDefaultSumw2(kTRUE);

TH1D *theta_e=new TH1D("theta_e", "Reco elastic event: Electron scattering reco angles from MESMER",NBINS,edges);
TH1D *theta_mu=new TH1D("theta_mu", "Reco elastic event: Muon scattering reco angles from MESMER",10,0.,0.005);
TH1D *h_opening=new TH1D("h_opening", "Opening angle reco events from MESMER",35,0.,0.035);

TH1D *theta_e_single=new TH1D("theta_e_single", "Reco single electron: Electron scattering reco angles from MESMER",NBINS,edges);
TH1D *theta_mu_single=new TH1D("theta_mu_single", "Reco single muon: Muon scattering reco angles from MESMER",10,0.,0.005);

TH1D *theta_e_gen=new TH1D("theta_e_gen", "Electron scattering generated angles from MESMER",NBINS,edges);
TH1D *theta_mu_gen=new TH1D("theta_mu_gen", "Muon scattering generated angles from MESMER",10,0.,0.005);
TH1D *h_opening_gen=new TH1D("h_opening_gen", "Opening angle generated events from MESMER",35,0.,0.035);

TH1D *theta_e_noreco=new TH1D("theta_e_noreco", "Electron scattering angles from MESMER when NOT reco",NBINS,edges);
TH1D *theta_mu_noreco=new TH1D("theta_mu_noreco", "Muon scattering angle from MESMER when NOT reco",10,0.,0.005);

TH1D *multiplicity = new TH1D("multiplicity", "number of outgoing reco tracks", 20,0,20);
TH1D *multiplicity_el = new TH1D("multiplicity_el", "number of outgoing reco tracks when elastic event is reconstructed with 2 tracks and without other", 20,0,20);
TH1D *multiplicity_cl = new TH1D("multiplicity_cl", "number of outgoing reco tracks when elastic event is reconstructed with clones", 20,0,20);
TH1D *multiplicity_other = new TH1D("multiplicity_other", "number of outgoing reco tracks with just other interaction ID", 20,0,20);
TH1D *multiplicity_both = new TH1D("multiplicity_both", "number of outgoing reco tracks when elastic event is reconstructed with 2 tracks + other", 20,0,20);
TH1D *multiplicity_one_el = new TH1D("multiplicity_one_el", "number of outgoing reco tracks when just one elastic is reco", 20,0,20);
TH1D *multiplicity_el3  = new TH1D("multiplicity_el3", "number of outgoing reco tracks when elastic event is reconstructed with other particles", 20,0,20);
TH1D *clones_el = new TH1D("clones_el", "number of electron's clone tracks when elastic event is reconstructed",NBINS,edges);
TH1D *clones_mu = new TH1D("clones_mu", "number of muon's clone tracks when elastic event is reconstructed",10,0.,0.005);

//old wnorm
//  double r_wnorm[6]={20.383577565024513,23.500259910776567,36.594795515579023,52.764785364090010,61.449970536155391,104.67651677602215};
//   double r_wnorm[6]={306.11916749125714,306.11916749125714,306.11916749125714,306.11916749125714,306.11916749125714,306.11916749125714};
//   double r_wnorm[6]={307.91166977973597,307.91166977973597,307.91166977973597,307.91166977973597,307.91166977973597,307.91166977973597};

//wnorm 1M NLO
// double r_wnorm[6]={21.522863161260670,38.514215499428630,41.710050452754516,49.388519620534474,61.328951873440843,105.50015224118995};
// double r_wnorm[6]={311.68129471616066,311.68129471616066,311.68129471616066,311.68129471616066,311.68129471616066,311.68129471616066};

//wnorm 1M NNLO
// double r_wnorm[6]={20.461867203505570,24.217332290303446,36.677249299608334,53.056125493170988,61.449970536155391,104.98858886225770};
//   double r_wnorm[6]={307.91166977973597,307.91166977973597,307.91166977973597,307.91166977973597,307.91166977973597,307.91166977973597};

//   double r_wnorm[6]={1.,1.,1.,1.,1.,1.};

//wnorm 1M/1M NLO 300outlier
	double r_wnorm[6]={20.786765103274643,33.313091221576336,42.396733790329876,50.584815206143652,61.828110400735824,106.88513370134392};

//1000chi2 1M single range
//      double r_wnorm[6]={169.10075957722117,169.10075957722117,169.10075957722117,169.10075957722117,169.10075957722117,169.10075957722117};


double sum_wgt=0.;
double sumW2[6]={0.};

double n_clones_el=0.;
double n_clones_mu=0.;
double total=0.;
double no_muin=0.;
double all=0.;
double n_el2=0.;
double n_el3=0.;
double n_cl=0.;
double n_one_el=0.;
double n_other=0.;
double n_zero=0.;
double n_both=0.;
double tot=0.;
double n_electron=0.;
double n_muon=0.;

for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
                cbmsim->GetEntry(i);
                if(i%100 == 0) cout<<"Entry "<<i<<endl;

int index=7;

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

        for(int n = 0; n < MCTrack->GetEntries(); n++) {
         const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
         if(MCTr->interactionID()==0 and MCTr->pdgCode()==13) {code_mu_in=n; p_muin_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); p_muin_MC.Unit();}
         if(MCTr->interactionID()==45 and MCTr->pdgCode()==11) {code_e=n; p_e_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); p_e_MC.Unit(); the_gen=p_muin_MC.Angle(p_e_MC);
                                                                                 theX_gen=MCTr->ax();
                                                                                 theY_gen=MCTr->ay();

        for(int s=0; s<TrackerPoints->GetEntries(); s++)
                         {const MUonETrackerPoint *TrackerPt = static_cast<const MUonETrackerPoint*>(TrackerPoints->At(s));
         if(TrackerPt->trackPDGCode()==11 and TrackerPt->trackID()==code_e and TrackerPt->stationID()==1){
                                                                                                 if(TrackerPt->moduleID()==0 or TrackerPt->moduleID()==4) hit_modXe++;
                                                                                                 if(TrackerPt->moduleID()==1 or TrackerPt->moduleID()==5) hit_modYe++;
                                                                                                 if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_e++;}
        	}
	 }
         if(MCTr->interactionID()==45 and MCTr->pdgCode()==-13) {code_mu=n; p_mu_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); p_mu_MC.Unit(); thmu_gen=p_muin_MC.Angle(p_mu_MC);
                                                                                 thmuX_gen=MCTr->ax();
                                                                                 thmuY_gen=MCTr->ay();
        for(int s=0; s<TrackerPoints->GetEntries(); s++)
                         {const MUonETrackerPoint *TrackerPt = static_cast<const MUonETrackerPoint*>(TrackerPoints->At(s));
	 if(TrackerPt->trackPDGCode()==-13 and TrackerPt->trackID()==code_mu and TrackerPt->stationID()==1){
                                                                                                 if(TrackerPt->moduleID()==0 or TrackerPt->moduleID()==4) hit_modXmu++;
                                                                                                 if(TrackerPt->moduleID()==1 or TrackerPt->moduleID()==5) hit_modYmu++;
                                                                                                 if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_mu++;}
		}
	 }
         i++;
	}


	if(code_mu_in!=-99 and code_e!=-99 and code_mu!=-99 and hit_modXe==4 and hit_modYe==4 and stereo_e>0 and hit_modXmu==4 and hit_modYmu==4 and stereo_mu>0){

double wnorm=0.;


if(the_gen>=0 and the_gen<0.005){index=0; wnorm=r_wnorm[0];}
if(the_gen>=0.005 and the_gen<0.010){index=1; wnorm=r_wnorm[1];}
if(the_gen>=0.010 and the_gen<0.015){index=2; wnorm=r_wnorm[2];}
if(the_gen>=0.015 and the_gen<0.020){index=3; wnorm=r_wnorm[3];}
if(the_gen>=0.020 and the_gen<0.025){index=4; wnorm=r_wnorm[4];}
if(the_gen>=0.025 and the_gen<0.032){index=5; wnorm=r_wnorm[5];}
 all+=MesmerEvent->wgt_full*wnorm;


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
MUonERecoOutputVertex vrtx = ReconstructionOutput->bestVertex();

Double_t chi=vrtx.chi2perDegreeOfFreedom();

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
                 if(code_e==track.linkedTrackID()) {yes2++; e++; theX_rec=track.xSlope(); theY_rec=track.ySlope();p_e.SetXYZ(theX_rec,theY_rec,1.0); p_e=p_e.Unit();the_rec=p_e.Angle(p_muin);}
                 if(code_mu==track.linkedTrackID()) {yes2++; mu++;thmuX_rec=track.xSlope(); thmuY_rec=track.ySlope();p_mu.SetXYZ(thmuX_rec,thmuY_rec,1.0); p_mu=p_mu.Unit();thmu_rec=p_mu.Angle(p_muin);}
                }
        if(track.processIDofLinkedTrack()!=45 and track.sector()==1){other++;}

         }//for

if(mu_in==0){no_muin+=MesmerEvent->wgt_full*wnorm;}

Double_t posxIN=pos_on_track(x0_in,th_inx,912.7);
Double_t posyIN=pos_on_track(y0_in,th_iny,912.7);


//if(stubs_muin==6 and abs(posxIN)<=1.5 and abs(posyIN)<=1.5 and chi2_muin<2 and thmu_gen>0.0002){
if(stubs_muin==6 and chi2_muin<2 and thmu_gen>0.0002){//and the_gen<0.01){


tot+=MesmerEvent->wgt_full*wnorm;

multiplicity->Fill(tracks.size()-1,MesmerEvent->wgt_full*wnorm);
theta_mu_gen->Fill(thmu_gen,MesmerEvent->wgt_full*wnorm);
theta_e_gen->Fill(the_gen,MesmerEvent->wgt_full*wnorm);

sumW2[index]+=(MesmerEvent->wgt_full*wnorm*MesmerEvent->wgt_full*wnorm);
sum_wgt+=MesmerEvent->wgt_full*wnorm;

h_opening_gen->Fill(opening_angle,MesmerEvent->wgt_full*wnorm);

if(e==0){theta_e_noreco->Fill(the_gen,MesmerEvent->wgt_full*wnorm);}

if(mu==0){theta_mu_noreco->Fill(thmu_gen,MesmerEvent->wgt_full*wnorm);}

if(e==1){theta_e_single->Fill(the_gen,MesmerEvent->wgt_full*wnorm);}

if(mu==1){theta_mu_single->Fill(thmu_gen,MesmerEvent->wgt_full*wnorm);}

if(yes2>=0){total+=MesmerEvent->wgt_full*wnorm;}

if(e==1 and mu==1 and other==0){n_el2+=MesmerEvent->wgt_full*wnorm; multiplicity_el->Fill(tracks.size()-1,MesmerEvent->wgt_full*wnorm);}
if(e==1 and mu==1 and other!=0){n_el3+=MesmerEvent->wgt_full*wnorm; multiplicity_el3->Fill(tracks.size()-1,MesmerEvent->wgt_full*wnorm);}
if( (e>1 or mu>1) and other==0){ n_cl+=MesmerEvent->wgt_full*wnorm;  multiplicity_cl->Fill(tracks.size()-1,MesmerEvent->wgt_full*wnorm);}
if(yes2==0 and other>0){ n_other+=MesmerEvent->wgt_full*wnorm; multiplicity_other->Fill(tracks.size()-1,MesmerEvent->wgt_full*wnorm);}
if(yes2==1 and other==0){ if(e==1)n_electron+=MesmerEvent->wgt_full*wnorm; if(mu==1)n_muon+=MesmerEvent->wgt_full*wnorm; 
			  n_one_el+=MesmerEvent->wgt_full*wnorm;  multiplicity_one_el->Fill(tracks.size()-1,MesmerEvent->wgt_full*wnorm);}
//if( (e>1 or mu>1) and other>0){ n_both+=MesmerEvent->wgt_full*wnorm;  multiplicity_both->Fill(tracks.size()-1,MesmerEvent->wgt_full*wnorm);}
if( (e>1 or mu>1) and other>0){ n_both+=MesmerEvent->wgt_full*wnorm;  multiplicity_both->Fill(tracks.size()-1,MesmerEvent->wgt_full*wnorm);}
//if(e>1 or mu>1 or (other>0 and e==1 and mu==1)){ n_both+=MesmerEvent->wgt_full*wnorm;  multiplicity_both->Fill(tracks.size()-1,MesmerEvent->wgt_full*wnorm);}

if(tracks.size()==1){ n_zero+=MesmerEvent->wgt_full*wnorm;}



if(chi!=0 and yes2>=2 and e>1){n_clones_el+=MesmerEvent->wgt_full*wnorm; clones_el->Fill(the_gen,MesmerEvent->wgt_full*wnorm);}
if(chi!=0 and yes2>=2 and mu>1){n_clones_mu+=MesmerEvent->wgt_full*wnorm; clones_mu->Fill(thmu_gen,MesmerEvent->wgt_full*wnorm);}

if(chi!=0 and e==1 and mu==1){
theta_mu->Fill(thmu_gen,MesmerEvent->wgt_full*wnorm);
theta_e->Fill(the_gen,MesmerEvent->wgt_full*wnorm);
h_opening->Fill(opening_angle,MesmerEvent->wgt_full*wnorm);

//                }             // if out particles
//else{h_2d_gen->Fill(the_gen,thmu_gen,MesmerEvent->wgt_full);}
                 }//if chi!=0

                }//if mu_in

        }//if generated

}//for

cout << "Sum of weights wgt_full " << sum_wgt << endl;
//cout << "Electron: " << endl;
//for(int n=0; n<6;n++){cout << n <<") SumW2 " << sqrt(sumW2[n]) << endl;}
cout << "N. of events where muin not reco " << no_muin << " over " << all<< " -> " << no_muin/all*100 << "%"<< endl;

cout << "N. of elastic events with 1 muon and 1 electron and no other "<< n_el2 << " -> " << n_el2/tot*100 << "%"<< endl;
cout << "N. of elastic events with 1 muon and 1 electron and other "<< n_el3 << " -> " << n_el3/tot*100 << "%"<< endl;
cout << "N. of elastic events with mu/e clones "<< n_cl << " -> " << n_cl/tot*100 << "%"<< endl;
cout << "N. of events with just 1 elastic "<< n_one_el << " -> " << n_one_el/tot*100 << "%"<< endl;
cout << "Of which muon: " << n_muon << " -> " << n_muon/n_one_el*100 << " and electron : " << n_electron << " -> " << n_electron/n_one_el*100 << endl;
cout << "N. of elastic events with just other (interaction ID !=45) particles "<< n_other << " -> " << n_other/tot*100 << "%"<< endl;
cout << "N. of events with other (interaction ID !=45) particles and clones "<< n_both << " -> " << n_both/tot*100 << "%"<< endl;
cout << "N. of events with no particles "<< n_zero << " -> " << n_zero/tot*100 << "%"<< endl;


cout <<"N. of electron clones " << n_clones_el << " over events with >= 2 e+mu " << total << " -> " << n_clones_el/total*100 << "%" << endl;
cout <<"N. of muon clones " << n_clones_mu << " over events with >= 2 e+mu " << total << " -> " << n_clones_mu/total*100 << "%" << endl;

for(int b=1; b<7; b++){
cout << b << " gen) " << theta_e_gen->GetBinContent(b) << " +- " << theta_e_gen->GetBinError(b)<< " ----> " << theta_e_gen->GetBinError(b)/theta_e_gen->GetBinContent(b)*100 << "%" << endl;
}

 auto legend = new TLegend(0.7,0.7,0.9,0.9);
   legend->AddEntry(multiplicity,"multiplicity","LEP");
   legend->AddEntry(multiplicity_el,"multiplicity_el","LEP");
   legend->AddEntry(multiplicity_cl,"multiplicity_cl","LEP");
   legend->AddEntry(multiplicity_other,"multiplicity_other","LEP");
   legend->AddEntry(multiplicity_both,"multiplicity_both","LEP");
   legend->AddEntry(multiplicity_one_el,"multiplicity_one_el","LEP");
   legend->AddEntry(multiplicity_el3,"multiplicity_el3","LEP");

TCanvas b("b","b",700,700);
b.Divide(2,3);
b.cd(1);
theta_mu_noreco->Draw("E");
gPad->SetLogy();
b.cd(2);
theta_e_noreco->Draw("E");
b.cd(3);
clones_mu->SetMinimum(1);
clones_mu->Draw("E");
b.cd(4);
clones_el->Draw("E");
clones_el->SetMinimum(1);
gPad->SetLogy();
b.cd(5);
multiplicity->Draw("E");
multiplicity_el->SetLineColor(kRed-3); 
multiplicity_el->Draw("E same");
multiplicity_cl->SetLineColor(kTeal-3); 
multiplicity_cl->Draw("E same");
multiplicity_other->SetLineColor(kViolet-3);
multiplicity_other->Draw("E same");
multiplicity_both->SetLineColor(kOrange); 
multiplicity_both->Draw("E same");
multiplicity_one_el->SetLineColor(kAzure+10); 
multiplicity_one_el->Draw("E same");
multiplicity_el3->SetLineColor(kPink-4); 
multiplicity_el3->Draw("E same");
gPad->SetLogy();
gStyle->SetOptStat(0);
legend->Draw("same");
b.SaveAs("new_mc_NLOvsLO/300chi2_not_reco_angles_full_0.2cut_reconstructibility_smallAngle.pdf");


TCanvas c("c","c",700,700);
c.Divide(2,3);
c.cd(1);
theta_mu_gen->SetLineColor(kPink);
theta_mu_gen->Draw("E");
theta_mu_single->Draw("E same");
gPad->SetLogy();
gStyle->SetOptStat(0);
c.cd(2);
theta_e_gen->SetLineColor(kPink);
theta_e_gen->Draw("E");
theta_e_single->Draw("E same");
gStyle->SetOptStat(0);
c.cd(3);
theta_mu_gen->SetLineColor(kPink);
theta_mu_gen->Draw("E");
theta_mu->Draw("E same");
gPad->SetLogy();
gStyle->SetOptStat(0);
c.cd(4);
theta_e_gen->SetLineColor(kPink);
theta_e_gen->Draw("E");
theta_e->Draw("E same");
gStyle->SetOptStat(0);
c.cd(5);
h_opening_gen->SetLineColor(kPink);
h_opening_gen->Draw("E");
h_opening->Draw("E same");
gStyle->SetOptStat(0);
gPad->SetLogy();
c.SaveAs("new_mc_NLOvsLO/300chi2_nonp_angle_full_MC_0.2cut_reconstructibility_smallAngle.pdf");
/*theta_e_gen->SaveAs("new_mc_NLOvsLO/300chi2_nonp_theta_e_gen_full_0.2cut.root");
theta_mu_gen->SaveAs("new_mc_NLOvsLO/300chi2_nonp_theta_mu_gen_full_0.2cut.root");
*/


TH1D * h0 = (TH1D*) theta_e_single->Clone();
TH1D * h0gen = (TH1D*) theta_e_gen->Clone();
h0->Divide(h0,h0gen,1,1,"B");
TH1D * h1 = (TH1D*) theta_mu_single->Clone();
TH1D * h1gen = (TH1D*) theta_mu_gen->Clone();
h1->Divide(h1,h1gen,1,1,"B");

TH1D * h2 = (TH1D*) theta_e->Clone();
TH1D * h2gen = (TH1D*) theta_e_gen->Clone();
h2->Divide(h2,h2gen,1,1,"B");
TH1D * h3 = (TH1D*) theta_mu->Clone();
TH1D * h3gen = (TH1D*) theta_mu_gen->Clone();
h3->Divide(h3,h3gen,1,1,"B");

TH1D * h4 = (TH1D*) h_opening->Clone();
TH1D * h4gen = (TH1D*) h_opening_gen->Clone();
h4->Divide(h4,h4gen,1,1,"B");

cout << "muon " << endl;
for(int b=1; b<7; b++){
cout << b << " ) " << h3->GetBinContent(b) << " +- " << h3->GetBinError(b)<< " ----> " << h3->GetBinError(b)/h3->GetBinContent(b)*100 << "%" << endl;
}
cout << "electron " << endl;
for(int b=1; b<7; b++){
cout << b << " ) " << h2->GetBinContent(b) << " +- " << h2->GetBinError(b)<< " ---> " << h2->GetBinError(b)/h2->GetBinContent(b)*100 << "%" <<endl ;
}

TCanvas a1("a1","a1",700,700);
a1.Divide(2,3);
a1.cd(1);
h1->SetMinimum(0.92);
h1->SetMaximum(1.02);
h1->Draw("E");
gStyle->SetOptStat(0);
a1.cd(2);
h0->SetMinimum(0.65);
h0->SetMaximum(1.02);
h0->Draw("E");
a1.cd(3);
h3->SetMinimum(0.82);
h3->SetMaximum(1.02);
h3->Draw("E");
gStyle->SetOptStat(0);
a1.cd(4);
h2->SetMinimum(0.65);
h2->SetMaximum(1.02);
h2->Draw("E");
gStyle->SetOptStat(0);
a1.cd(5);
h4->SetMinimum(0.5);
h4->SetMaximum(1.02);
h4->Draw("E");
gStyle->SetOptStat(0);
a1.SaveAs("new_mc_NLOvsLO/300chi2_nonp_eff_full_MC_0.2cut_reconstructibility_smallAngle.pdf");

/*h0->SaveAs("new_mc_NLOvsLO/300chi2_nonp_el_single_eff_full_MC_0.2cut.root");
h1->SaveAs("new_mc_NLOvsLO/300chi2_nonp_mu_single_eff_full_MC_0.2cut.root");
h2->SaveAs("new_mc_NLOvsLO/300chi2_nonp_el_eff_full_MC_0.2cut.root");
h3->SaveAs("new_mc_NLOvsLO/300chi2_nonp_mu_eff_full_MC_0.2cut.root");
h4->SaveAs("new_mc_NLOvsLO/300chi2_nonp_op_eff_full_MC_0.2cut.root");*/

}
