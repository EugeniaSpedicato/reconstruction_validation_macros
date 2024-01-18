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
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency/theta_0-5mrad_100k_1hit.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency/theta_5-10mrad_100k_1hit.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency/theta_10-15mrad_100k_1hit.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency/theta_15-20mrad_100k_1hit.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency/theta_20-25mrad_100k_1hit.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency/theta_25-32mrad_100k_1hit.root");

//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/theta_0-32mrad_600k_1hit.root");

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

std::array<TH1D*,6> theta_e;
for(int m=0; m<6; m++){
        string name="theta_e"+to_string(m);
        theta_e.at(m)=new TH1D(name.c_str(),"Reco elastic event: Electron scattering reco angles from MESMER",NBINS,edges);}
std::array<TH1D*,6> theta_mu;
for(int m=0; m<6; m++){
        string name="theta_mu"+to_string(m);
        theta_mu.at(m)=new TH1D(name.c_str(),"Reco elastic event: Muon scattering reco angles from MESMER",NBINS,edges);}
std::array<TH1D*,6> h_opening;
for(int m=0; m<6; m++){
        string name="h_opening"+to_string(m);
        h_opening.at(m)=new TH1D(name.c_str(),"Opening angle reco events from MESMER",35,0.,0.035);
        }
std::array<TH1D*,6> theta_e_single;
for(int m=0; m<6; m++){
        string name="theta_e_single"+to_string(m);
        theta_e_single.at(m)=new TH1D(name.c_str(),"Reco single electron: Electron scattering reco angles from MESMER",NBINS,edges);}
std::array<TH1D*,6> theta_mu_single;
for(int m=0; m<6; m++){
        string name="theta_mu_single"+to_string(m);
        theta_mu_single.at(m)=new TH1D(name.c_str(),"Reco single muon: Muon scattering reco angles from MESMER",NBINS,edges);}
std::array<TH1D*,6> theta_e_gen;
for(int m=0; m<6; m++){
        string name="theta_e_gen"+to_string(m);
        theta_e_gen.at(m)=new TH1D(name.c_str(),"Reco elastic event: Electron scattering generated angles from MESMER",NBINS,edges);}
std::array<TH1D*,6> theta_mu_gen;
for(int m=0; m<6; m++){
        string name="theta_mu_gen"+to_string(m);
        theta_mu_gen.at(m)=new TH1D(name.c_str(),"Reco elastic event: Muon scattering generated angles from MESMER",NBINS,edges);}
std::array<TH1D*,6> h_opening_gen;
for(int m=0; m<6; m++){
        string name="h_opening_gen"+to_string(m);
        h_opening_gen.at(m)=new TH1D(name.c_str(),"Opening angle generated events from MESMER",35,0.,0.035);
        }
std::array<TH1D*,6> theta_e_noreco;
for(int m=0; m<6; m++){
        string name="theta_e_noreco"+to_string(m);
        theta_e_noreco.at(m)=new TH1D(name.c_str(),"Electron scattering angles from MESMER when NOT reco",NBINS,edges);}

std::array<TH1D*,6> theta_mu_noreco;
for(int m=0; m<6; m++){
        string name="theta_mu_noreco"+to_string(m);
        theta_mu_noreco.at(m)=new TH1D(name.c_str(),"Muon scattering angle from MESMER when NOT reco",10,0.,0.005);}

   double r_wnorm[6]={20.383577565024513,23.500259910776567,36.594795515579023,52.764785364090010,61.449970536155391,104.67651677602215};
   int index=7;
//   double r_wnorm[6]={1.,1.,1.,1.,1.,1.};
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


if(the_gen>=0. and the_gen<5.){wnorm=r_wnorm[0]; index=0;}
if(the_gen>=5. and the_gen<10.){wnorm=r_wnorm[1]; index=1;}
if(the_gen>=10. and the_gen<15.){wnorm=r_wnorm[2]; index=2;}
if(the_gen>=15. and the_gen<20.){wnorm=r_wnorm[3]; index=3;}
if(the_gen>=20. and the_gen<25.){wnorm=r_wnorm[4]; index=4;}
if(the_gen>=25. and the_gen<32.){wnorm=r_wnorm[5]; index=5;}

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

theta_mu_gen[index]->Fill(thmu_gen,MesmerEvent->wgt_full);
theta_e_gen[index]->Fill(the_gen,MesmerEvent->wgt_full);

h_opening_gen[index]->Fill(opening_angle,MesmerEvent->wgt_full);

if(e==0){theta_e_noreco[index]->Fill(the_gen,MesmerEvent->wgt_full);}

if(mu==0){theta_mu_noreco[index]->Fill(thmu_gen,MesmerEvent->wgt_full);}

if(e==1){theta_e_single[index]->Fill(the_gen,MesmerEvent->wgt_full);}

if(mu==1){theta_mu_single[index]->Fill(thmu_gen,MesmerEvent->wgt_full);}

if(chi!=0 and yes2==2){
theta_mu[index]->Fill(thmu_gen,MesmerEvent->wgt_full);
theta_e[index]->Fill(the_gen,MesmerEvent->wgt_full);
h_opening[index]->Fill(opening_angle,MesmerEvent->wgt_full);

//                }             // if out particles
//else{h_2d_gen->Fill(the_gen,thmu_gen,MesmerEvent->wgt_full);}
                 }//if chi!=0

                }//if mu_in

        }//if generated

}//for

TCanvas b("b","b",700,700);
b.Divide(1,2);
for(int h=0; h<6; h++){
b.cd(1);
theta_mu_noreco[h]->Scale(r_wnorm[h]/20874.9);
if(h==0)theta_mu_noreco[h]->Draw("E");
else theta_mu_noreco[h]->Draw("E same");
b.cd(2);
theta_e_noreco[h]->Scale(r_wnorm[h]/20874.9);
if(h==0)theta_e_noreco[h]->Draw("E");
else theta_e_noreco[h]->Draw("E same");
}
b.SaveAs("mc_NLOvsLO/try_not_reco_angles_full.pdf");

TCanvas c("c","c",700,700);
c.Divide(2,3);
for(int h=0; h<6; h++){
c.cd(1);
theta_mu_gen[h]->Scale(r_wnorm[h]/20874.9);
theta_mu_gen[h]->SetLineColor(kPink);
if(h==0)theta_mu_gen[h]->Draw("E");
else theta_mu_gen[h]->Draw("E same");
theta_mu_single[h]->Scale(r_wnorm[h]/20874.9);
theta_mu_single[h]->Draw("E same");
gPad->SetLogy();
gStyle->SetOptStat(0);
c.cd(2);
theta_e_gen[h]->Scale(r_wnorm[h]/20874.9);
theta_e_gen[h]->SetLineColor(kPink);
if(h==0)theta_e_gen[h]->Draw("E");
else theta_e_gen[h]->Draw("E same");
theta_e_single[h]->Scale(r_wnorm[h]/20874.9);
theta_e_single[h]->Draw("E same");
gStyle->SetOptStat(0);
c.cd(3);
theta_mu_gen[h]->Scale(r_wnorm[h]/20874.9);
theta_mu_gen[h]->SetLineColor(kPink);
if(h==0)theta_mu_gen[h]->Draw("E");
else theta_mu_gen[h]->Draw("E same");
theta_mu[h]->Scale(r_wnorm[h]/20874.9);
theta_mu[h]->Draw("E same");
gPad->SetLogy();
gStyle->SetOptStat(0);
c.cd(4);
theta_e_gen[h]->Scale(r_wnorm[h]/20874.9);
theta_e_gen[h]->SetLineColor(kPink);
if(h==0)theta_e_gen[h]->Draw("E");
else theta_e_gen[h]->Draw("E same");
theta_e[h]->Scale(r_wnorm[h]/20874.9);
theta_e[h]->Draw("E same");
gStyle->SetOptStat(0);
c.cd(5);
h_opening_gen[h]->Scale(r_wnorm[h]/20874.9);
h_opening_gen[h]->SetLineColor(kPink);
if(h==0)h_opening_gen[h]->Draw("E");
else h_opening_gen[h]->Draw("E same");
h_opening[h]->Scale(r_wnorm[h]/20874.9);
h_opening[h]->Draw("E same");
gStyle->SetOptStat(0);
gPad->SetLogy();
}
c.SaveAs("mc_NLOvsLO/try_nonp_angle_full_MC.pdf");

/*TH1D * h0 = (TH1D*) theta_e_single->Clone();
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
a1.SaveAs("mc_NLOvsLO/nonp_eff_full_MC.pdf");
h0->SaveAs("mc_NLOvsLO/nonp_el_single_eff_full_MC.root");
h1->SaveAs("mc_NLOvsLO/nonp_mu_single_eff_full_MC.root");
h2->SaveAs("mc_NLOvsLO/nonp_el_eff_full_MC.root");
h3->SaveAs("mc_NLOvsLO/nonp_mu_eff_full_MC.root");
h4->SaveAs("mc_NLOvsLO/nonp_op_eff_full_MC.root");
*/

}

