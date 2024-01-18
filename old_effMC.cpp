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
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_0-5mrad_100k_1hit.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_5-10mrad_100k_1hit.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_10-15mrad_100k_1hit.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_15-20mrad_100k_1hit.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_20-25mrad_100k_1hit.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_25-32mrad_100k_1hit.root");

//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/theta_0-32mrad_100k_1hit_NLO.root");

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

//old wnorm
//  double r_wnorm[6]={20.383577565024513,23.500259910776567,36.594795515579023,52.764785364090010,61.449970536155391,104.67651677602215};
//   double r_wnorm[6]={306.11916749125714,306.11916749125714,306.11916749125714,306.11916749125714,306.11916749125714,306.11916749125714};

//new wnorm
   double r_wnorm[6]={21.522863161260670,38.514215499428630,41.710050452754516,49.388519620534474,61.328951873440843,105.50015224118995};
//   double r_wnorm[6]={311.68129471616066,311.68129471616066,311.68129471616066,311.68129471616066,311.68129471616066,311.68129471616066};

//   double r_wnorm[6]={1.,1.,1.,1.,1.,1.};

double sum_wgt=0.;
double sum_wgt_nownorm=0.;

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

        for(int n = 0; n < MCTrack->GetEntries(); n++) {
         const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
         if(MCTr->interactionID()==0 and MCTr->pdgCode()==13) {code_mu_in=n; p_muin_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); p_muin_MC.Unit();}
         if(MCTr->interactionID()==45 and MCTr->pdgCode()==11) {code_e=n; p_e_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); p_e_MC.Unit(); the_gen=p_muin_MC.Angle(p_e_MC);
                                                                                 theX_gen=MCTr->ax();
                                                                                 theY_gen=MCTr->ay();}
         if(MCTr->interactionID()==45 and MCTr->pdgCode()==13) {code_mu=n; p_mu_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); p_mu_MC.Unit(); thmu_gen=p_muin_MC.Angle(p_mu_MC);
                                                                                 thmuX_gen=MCTr->ax();
                                                                                 thmuY_gen=MCTr->ay();}
         i++;
	}

	if(code_mu_in!=-99 and code_e!=-99 and code_mu!=-99){

double wnorm=1.;

if(the_gen>=0 and the_gen<0.005)wnorm=r_wnorm[0];
if(the_gen>=0.005 and the_gen<0.010)wnorm=r_wnorm[1];
if(the_gen>=0.010 and the_gen<0.015)wnorm=r_wnorm[2];
if(the_gen>=0.015 and the_gen<0.020)wnorm=r_wnorm[3];
if(the_gen>=0.020 and the_gen<0.025)wnorm=r_wnorm[4];
if(the_gen>=0.025 and the_gen<0.032)wnorm=r_wnorm[5];

cout << "the_gen " << the_gen << ", wnorm " << wnorm << endl;

double opening_angle=p_mu_MC.Angle(p_e_MC);


Int_t yes_mu=0;
Int_t yes_e=0;
Int_t yes2=0;
Double_t th_inx,th_iny,x0_in,y0_in;
Double_t chi2_muin;
Double_t stubs_muin;
Int_t e=0;
Int_t mu=0;
TVector3 p_muin,p_e,p_mu;
Double_t the_rec,theX_rec,theY_rec,thmu_rec,thmuX_rec,thmuY_rec;
MUonERecoOutputVertex vrtx = ReconstructionOutput->bestVertex();

Double_t chi=vrtx.chi2perDegreeOfFreedom();

vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();

         for (auto&& track : tracks) {

        if(code_mu_in==track.linkedTrackID() and track.sector()==0){
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
                {yes2++;
                 if(code_e==track.linkedTrackID()) {e=1; theX_rec=track.xSlope(); theY_rec=track.ySlope();p_e.SetXYZ(theX_rec,theY_rec,1.0); p_e=p_e.Unit();the_rec=p_e.Angle(p_muin);}
                 if(code_mu==track.linkedTrackID()) {mu=1;thmuX_rec=track.xSlope(); thmuY_rec=track.ySlope();p_mu.SetXYZ(thmuX_rec,thmuY_rec,1.0); p_mu=p_mu.Unit();thmu_rec=p_mu.Angle(p_muin);}
                }

         }//for


Double_t posxIN=pos_on_track(x0_in,th_inx,912.7);
Double_t posyIN=pos_on_track(y0_in,th_iny,912.7);


if(stubs_muin==6 and abs(posxIN)<=1.5 and abs(posyIN)<=1.5 and chi2_muin<2){//and thmu_gen>0.0002){

theta_mu_gen->Fill(thmu_gen,MesmerEvent->wgt_full*wnorm);
theta_e_gen->Fill(the_gen,MesmerEvent->wgt_full*wnorm);

sum_wgt_nownorm+=MesmerEvent->wgt_full;
sum_wgt+=MesmerEvent->wgt_full*wnorm;

h_opening_gen->Fill(opening_angle,MesmerEvent->wgt_full*wnorm);

if(e==0){theta_e_noreco->Fill(the_gen,MesmerEvent->wgt_full*wnorm);}

if(mu==0){theta_mu_noreco->Fill(thmu_gen,MesmerEvent->wgt_full*wnorm);}

if(e==1){theta_e_single->Fill(the_gen,MesmerEvent->wgt_full*wnorm);}

if(mu==1){theta_mu_single->Fill(thmu_gen,MesmerEvent->wgt_full*wnorm);}

if(chi!=0 and yes2==2){
theta_mu->Fill(thmu_gen,MesmerEvent->wgt_full*wnorm);
theta_e->Fill(the_gen,MesmerEvent->wgt_full*wnorm);
h_opening->Fill(opening_angle,MesmerEvent->wgt_full*wnorm);

//                }             // if out particles
//else{h_2d_gen->Fill(the_gen,thmu_gen,MesmerEvent->wgt_full);}
                 }//if chi!=0

                }//if mu_in

        }//if generated

}//for

cout << "Sum of weights " << sum_wgt << endl;
//cout << "Electron: " << endl;
//for(int b=1; b<7; b++){ cout << b+1 << ") " << theta_e_gen->GetBinContent(b)*306.11916749125714/sum_wgt_nownorm << " +- " << theta_e_gen->GetBinError(b)*306.11916749125714/sum_wgt_nownorm <<endl;} 
//<< " -> considering sum_wgt: " << theta_e_gen->GetBinContent(b)*306.11916749125714/sum_wgt_nownorm << endl;}

TCanvas b("b","b",700,700);
b.Divide(1,2);
b.cd(1);
theta_mu_noreco->Draw("E");
b.cd(2);
theta_e_noreco->Draw("E");
b.SaveAs("new_mc_NLOvsLO/not_reco_angles_full.pdf");

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
c.SaveAs("new_mc_NLOvsLO/nonp_angle_full_MC.pdf");
theta_e_gen->SaveAs("new_mc_NLOvsLO/nonp_theta_e_gen_full.root");
theta_mu_gen->SaveAs("new_mc_NLOvsLO/nonp_theta_mu_gen_full.root");


TH1D * h0 = (TH1D*) theta_e_single->Clone();
TH1D * h0gen = (TH1D*) theta_e_gen->Clone();
h0->Divide(h0gen);
TH1D * h1 = (TH1D*) theta_mu_single->Clone();
TH1D * h1gen = (TH1D*) theta_mu_gen->Clone();
h1->Divide(h1gen);

TH1D * h2 = (TH1D*) theta_e->Clone();
TH1D * h2gen = (TH1D*) theta_e_gen->Clone();
h2->Divide(h2gen);
TH1D * h3 = (TH1D*) theta_mu->Clone();
TH1D * h3gen = (TH1D*) theta_mu_gen->Clone();
h3->Divide(h3gen);

TH1D * h4 = (TH1D*) h_opening->Clone();
TH1D * h4gen = (TH1D*) h_opening_gen->Clone();
h4->Divide(h4gen);

TCanvas a1("a1","a1",700,700);
a1.Divide(2,3);
a1.cd(1);
h0->SetMinimum(0.);
h0->Draw("E");
gStyle->SetOptStat(0);
a1.cd(2);
h1->SetMinimum(0.);
h1->Draw("E");
a1.cd(3);
h2->SetMinimum(0.);
h2->Draw("E");
gStyle->SetOptStat(0);
a1.cd(4);
h3->SetMinimum(0.);
h3->Draw("E");
gStyle->SetOptStat(0);
a1.cd(5);
h4->SetMinimum(0.);
h4->Draw("E");
gStyle->SetOptStat(0);
a1.SaveAs("new_mc_NLOvsLO/nonp_eff_full_MC.pdf");
h0->SaveAs("new_mc_NLOvsLO/nonp_el_single_eff_full_MC.root");
h1->SaveAs("new_mc_NLOvsLO/nonp_mu_single_eff_full_MC.root");
h2->SaveAs("new_mc_NLOvsLO/nonp_el_eff_full_MC.root");
h3->SaveAs("new_mc_NLOvsLO/nonp_mu_eff_full_MC.root");
h4->SaveAs("new_mc_NLOvsLO/nonp_op_eff_full_MC.root");
}
