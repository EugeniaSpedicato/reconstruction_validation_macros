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

void gen_eff(){

TChain * cbmsim = new TChain("cbmsim");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency/theta_0-5mrad_100k_1hit.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency/theta_5-10mrad_100k_1hit.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency/theta_10-15mrad_100k_1hit.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency/theta_15-20mrad_100k_1hit.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency/theta_20-25mrad_100k_1hit.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency/theta_25-32mrad_100k_1hit.root");

//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency/theta_0-32mrad_100k_1hit.root");

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

//new wnorm 100k
//   double r_wnorm[6]={21.522863161260670,38.514215499428630,41.710050452754516,49.388519620534474,61.328951873440843,105.50015224118995};
//   double r_wnorm[6]={311.68129471616066,311.68129471616066,311.68129471616066,311.68129471616066,311.68129471616066,311.68129471616066};

//new wnorm 1M
//   double r_wnorm[6]={21.522863161260670,38.514215499428630,41.710050452754516,49.388519620534474,61.328951873440843,105.50015224118995};
//   double r_wnorm[6]={311.68129471616066,311.68129471616066,311.68129471616066,311.68129471616066,311.68129471616066,311.68129471616066};

//wnorm 100k NNLO
 double r_wnorm[6]={20.461867203505570,24.217332290303446,36.677249299608334,53.056125493170988,61.449970536155391,104.98858886225770};
//   double r_wnorm[6]={307.91166977973597,307.91166977973597,307.91166977973597,307.91166977973597,307.91166977973597,307.91166977973597};


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
         if(MCTr->interactionID()==45 and MCTr->pdgCode()==-13) {code_mu=n; p_mu_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); p_mu_MC.Unit(); thmu_gen=p_muin_MC.Angle(p_mu_MC);
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


theta_mu_gen->Fill(thmu_gen,MesmerEvent->wgt_full*wnorm);
theta_e_gen->Fill(the_gen,MesmerEvent->wgt_full*wnorm);
sum_wgt+=MesmerEvent->wgt_full*wnorm;

        }//if generated

}//for

cout << "Sum of weights " << sum_wgt << endl;

TCanvas c("c","c",700,700);
c.Divide(1,2);
c.cd(1);
theta_mu_gen->SetLineColor(kPink);
theta_mu_gen->Draw("E");
gPad->SetLogy();
gStyle->SetOptStat(0);
c.cd(2);
theta_e_gen->SetLineColor(kPink);
theta_e_gen->Draw("E");
gStyle->SetOptStat(0);
c.SaveAs("mc_NNLOvsLO/100k_gen_angle_full_MC.pdf");
theta_e_gen->SaveAs("mc_NNLOvsLO/100k_gen_theta_e_gen_full.root");
theta_mu_gen->SaveAs("mc_NNLOvsLO/100k_gen_theta_mu_gen_full.root");


}
