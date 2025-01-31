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

                auto pos_on_track = [](double q, double m, double z){return (q + m*z);};

int minbias_p_pre_dati_noReco(string version, int nhits,string bend,string mc){

  int nthreads = 6;
  ROOT::EnableImplicitMT(nthreads);


TChain * cbmsim = new TChain("cbmsim");
TChain * cbmsim_g = new TChain("cbmsim");

if(version=="my_modifica_minbias" and nhits==0 and bend=="Bend"){
          cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/minbias_skim/asMC_minbias_10Mevents_0hit.root");
          cbmsim_g->Add("/mnt/raid10/DATA/espedica/fairmu/gen_digi/minbias_SIM-DIGI_10Mevents.root");

          cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/minbias_skim/asMC_minbias_11Mevents_0hit.root");
          cbmsim_g->Add("/mnt/raid10/DATA/espedica/fairmu/gen_digi/minbias_SIM-DIGI_11Mevents.root");

          cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/minbias_skim/asMC_minbias_8Mevents_0hit.root");
          cbmsim_g->Add("/mnt/raid10/DATA/espedica/fairmu/gen_digi/minbias_SIM-DIGI_8Mevents.root");
cbmsim->AddFriend(cbmsim_g);
}
else if(version=="my_modifica_minbias_preskim" and nhits==0 and bend=="Bend"){
          cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/minbias_skim/MCMinimumBias_SimDigiReco_8MeventsGNN2.root");
          cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/minbias_skim/MCMinimumBias_SimDigiReco_master10Mevents.root");
	  cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/minbias_skim/MCMinimumBias_SimDigiReco_11MeventsGNN.root");
}

/*if(version=="my_modifica_minbias" and nhits==0 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/minbias_skim/minbias_nuovo_algo_0hit.root");
}*/


ROOT::TTreeProcessorMT tp1(*cbmsim,nthreads);

      MUonERecoOutput *ReconstructionOutput = 0;



TH1::SetDefaultSumw2(kTRUE);
   const Int_t NBINS = 6;
   Double_t edges[NBINS + 1] = {0.0, 0.005, 0.010, 0.015, 0.020, 0.025, 0.032};
   ROOT::TThreadedObject<TH1D> theta_e("theta_e", "Electron scattering reco angles",32,0.,0.032);
   ROOT::TThreadedObject<TH1D> theta_mu("theta_mu", "Muon scattering reco angles",50,0.,0.005);
   ROOT::TThreadedObject<TH1D> h_opening("h_opening", "Opening angle reco events",32,0.,0.032);
   ROOT::TThreadedObject<TH1D> d_aco("d_aco_real", "Acoplanarity",600,-3.2,3.2);



   ROOT::TThreadedObject<TH1D> theta_e_pre("theta_e_pre", "Electron scattering reco angles pre-cuts",32,0.,0.032);
   ROOT::TThreadedObject<TH1D> theta_mu_pre("theta_mu_pre", "Muon scattering reco angles pre-cuts",50,0.,0.005);
   ROOT::TThreadedObject<TH1D> h_opening_pre("h_opening_pre", "Opening angle reco events pre-cuts",32,0.,0.032);
   ROOT::TThreadedObject<TH1D> d_aco_pre("d_aco_real_pre", "Acoplanarity pre-cuts",600,-3.2,3.2);

   ROOT::TThreadedObject<TH1D> theq("theq", "Equal angles",50,0.,0.005);

   ROOT::TThreadedObject<TH1D> h_z_pos_pre("h_z_pos_pre","Z selected events pos fit pre-selection",450, 890.,940.);
   ROOT::TThreadedObject<TH1D> h_z_pos("h_z_pos","Z selected events pos fit",450, 890.,940.);


   const Int_t NBINS2 = 12;
   Double_t edges2[NBINS2 + 1] = {0.,0.0002,0.0004,0.0006,0.0008,0.001,0.00125,0.0015,0.00175,0.002,0.003,0.004,0.005};


//   ROOT::TThreadedObject<TH2D> h_2d("h2D","theta mu vs theta E with all cuts", 16, 0.,0.032,NBINS2,edges2);
//   ROOT::TThreadedObject<TH2D> h_2d_pre("h2D_pre","theta mu vs theta E with all cuts", 16, 0.,0.032,NBINS2,edges2);

   ROOT::TThreadedObject<TH2D> h_2d("h2D","theta mu vs theta E with all cuts", 320, 0.,0.032,50,0.,0.005);
   ROOT::TThreadedObject<TH2D> h_2d_pre("h2D_pre","theta mu vs theta E with all cuts", 320, 0.,0.032,50,0.,0.005);


   ROOT::TThreadedObject<TH1D> h_chi2("h_chi2","chi2 mu in", 100,0,10);
   ROOT::TThreadedObject<TH1D> h_xslope("h_xslope","xSlope mu in ", 200,0.,0.02);

   ROOT::TThreadedObject<TH1D> h_e_pre("h_e_pre","h_e_pre",300,0.,0.035);
   ROOT::TThreadedObject<TH1D> h_e_post("h_e_post","h_e_post",300,0.,0.035);
   ROOT::TThreadedObject<TH1D> h_mu_pre("h_mu_pre","h_mu_pre", 100,0.,0.005);
   ROOT::TThreadedObject<TH1D> h_mu_post("h_mu_post","h_mu_post", 100,0.,0.005);
   ROOT::TThreadedObject<TH1D> h_inc_pre("h_inc_pre","h_inc_pre", 100,0.,0.005);
   ROOT::TThreadedObject<TH1D> h_inc_post("h_inc_post","h_inc_post", 100,0.,0.005);

/*   ROOT::TThreadedObject<TH1D> residualX_pre("residualX_pre","residual X_pre",400,-2.,2.);
   ROOT::TThreadedObject<TH1D> residualY_pre("residualY_pre","residual Y_pre",200,-1.,1.);

   ROOT::TThreadedObject<TH1D> residualX_post("residualX_post","residual X_post",200,-1.,1.);
   ROOT::TThreadedObject<TH1D> residualY_post("residualY_post","residual Y_post",200,-1.,1.);
*/

   ROOT::TThreadedObject<TH1D> residualU_critic("residualU_critic","residual U when track distance <0.2cm and track_e no Ustub",200,-1.,1.);
   ROOT::TThreadedObject<TH1D> residualU_no_critic("residualU_no_critic","residual U when track distance <0.2cm",400,-2.,2.);

   ROOT::TThreadedObject<TH1D> residualU("residualU","residual U when track distance >0.2cm",200,-1.,1.);
   ROOT::TThreadedObject<TH1D> residualU_no("residualU_no","residual U when track distance >0.2cm and track_e no Ustub",400,-2.,2.);


   ROOT::TThreadedObject<TH1D> h_dist_mod0_range0("h_dist_mod0_range0","h_dist mod0 el_mu in cm mod0 range 0-5mrad",80,-0.4,0.4);
   ROOT::TThreadedObject<TH1D> h_dist_mod0_range1("h_dist_mod0_range1","h_dist mod0 el_mu in cm mod0 range 5-7.5 mrad",80,-0.4,0.4);
   ROOT::TThreadedObject<TH1D> h_dist_mod0_range2("h_dist_mod0_range2","h_dist mod0 el_mu in cm mod0 range 7.5-10 mrad",80,-0.4,0.4);
   ROOT::TThreadedObject<TH1D> h_dist_mod0_range3("h_dist_mod0_range3","h_dist mod0 el_mu in cm mod0 range 10-15",160,-0.8,0.8);
   ROOT::TThreadedObject<TH1D> h_dist_mod0_range4("h_dist_mod0_range4","h_dist mod0 el_mu in cm mod0 range 15-25",160,-0.8,0.8);
   ROOT::TThreadedObject<TH1D> h_dist_mod0_range5("h_dist_mod0_range5","h_dist mod0 el_mu in cm mod0 range 25-35",160,-0.8,0.8);

   ROOT::TThreadedObject<TH1D> h_dist_mod2_range0("h_dist_mod2_range0","h_dist mod2 el_mu in cm mod0 range 0-5mrad",160,-0.8,0.8);
   ROOT::TThreadedObject<TH1D> h_dist_mod2_range1("h_dist_mod2_range1","h_dist mod2 el_mu in cm mod0 range 5-7.5 mrad",160,-0.8,0.8);
   ROOT::TThreadedObject<TH1D> h_dist_mod2_range2("h_dist_mod2_range2","h_dist mod2 el_mu in cm mod0 range 7.5-10 mrad",160,-0.8,0.8);
   ROOT::TThreadedObject<TH1D> h_dist_mod2_range3("h_dist_mod2_range3","h_dist mod2 el_mu in cm mod0 range 10-15",200,-1.,1.);
   ROOT::TThreadedObject<TH1D> h_dist_mod2_range4("h_dist_mod2_range4","h_dist mod2 el_mu in cm mod0 range 15-25",200,-1.,1.);
   ROOT::TThreadedObject<TH1D> h_dist_mod2_range5("h_dist_mod2_range5","h_dist mod2 el_mu in cm mod0 range 25-35",200,-1.,1.);


   ROOT::TThreadedObject<TH1D> h_dist0("h_dist0", "distance stub electron - stub muon in cm mod0",300,-1.5,1.5);
   ROOT::TThreadedObject<TH1D> h_dist1("h_dist1", "distance stub electron - stub muon in cm mod1",300,-1.5,1.5);
   ROOT::TThreadedObject<TH1D> h_dist2("h_dist2", "distance stub electron - stub muon in cm mod2",300,-1.5,1.5);
   ROOT::TThreadedObject<TH1D> h_dist3("h_dist3", "distance stub electron - stub muon in cm mod3",300,-1.5,1.5);
   ROOT::TThreadedObject<TH1D> h_dist4("h_dist4", "distance stub electron - stub muon in cm mod4",300,-1.5,1.5);
   ROOT::TThreadedObject<TH1D> h_dist5("h_dist5", "distance stub electron - stub muon in cm mod5",300,-1.5,1.5);


   ROOT::TThreadedObject<TH2D> h_pos_range0("h_pos_range0","h_pos el_mu in cm mod0 range 0-5mrad",400,-2,2,400,-2,2);
   ROOT::TThreadedObject<TH2D> h_pos_range1("h_pos_range1","h_pos el_mu in cm mod0 range 5-7.5 mrad",400,-2,2,400,-2,2);
   ROOT::TThreadedObject<TH2D> h_pos_range2("h_pos_range2","h_pos el_mu in cm mod0 range 7.5-10 mrad",400,-2,2,400,-2,2);
   ROOT::TThreadedObject<TH2D> h_pos_range3("h_pos_range3","h_pos el_mu in cm mod0 range 10-15",400,-2,2,400,-2,2);
   ROOT::TThreadedObject<TH2D> h_pos_range4("h_pos_range4","h_pos el_mu in cm mod0 range 15-25",400,-2,2,400,-2,2);
   ROOT::TThreadedObject<TH2D> h_pos_range5("h_pos_range5","h_pos el_mu in cm mod0 range 25-35",400,-2,2,400,-2,2);

   ROOT::TThreadedObject<TH1D> h_pos_el_range0("h_pos_el_range0","h_pos el in cm mod0 range 0-5mrad",400,-2,2);
   ROOT::TThreadedObject<TH1D> h_pos_el_range1("h_pos_el_range1","h_pos el in cm mod1 range 5-7.5mrad",400,-2,2);
   ROOT::TThreadedObject<TH1D> h_pos_el_range2("h_pos_el_range2","h_pos el in cm mod2 range 7.5-10mrad",400,-2,2);
   ROOT::TThreadedObject<TH1D> h_pos_el_range3("h_pos_el_range3","h_pos el in cm mod3 range 10-15mrad",400,-2,2);
   ROOT::TThreadedObject<TH1D> h_pos_el_range4("h_pos_el_range4","h_pos el in cm mod4 range 15-25mrad",400,-2,2);
   ROOT::TThreadedObject<TH1D> h_pos_el_range5("h_pos_el_range5","h_pos el in cm mod5 range 25-35mrad",400,-2,2);

   ROOT::TThreadedObject<TH1D> h_pos_mu_range0("h_pos_mu_range0","h_pos mu in cm mod0 range 0-5mrad",400,-2,2);
   ROOT::TThreadedObject<TH1D> h_pos_mu_range1("h_pos_mu_range1","h_pos mu in cm mod1 range 5-7.5mrad",400,-2,2);
   ROOT::TThreadedObject<TH1D> h_pos_mu_range2("h_pos_mu_range2","h_pos mu in cm mod2 range 7.5-10mrad",400,-2,2);
   ROOT::TThreadedObject<TH1D> h_pos_mu_range3("h_pos_mu_range3","h_pos mu in cm mod3 range 10-15mrad",400,-2,2);
   ROOT::TThreadedObject<TH1D> h_pos_mu_range4("h_pos_mu_range4","h_pos mu in cm mod4 range 15-25mrad",400,-2,2);
   ROOT::TThreadedObject<TH1D> h_pos_mu_range5("h_pos_mu_range5","h_pos mu in cm mod5 range 25-35mrad",400,-2,2);

   ROOT::TThreadedObject<TH2D> h_du_aco("h_du_aco","distance u VS aco",600,-3.2,3.2,200,-1.,1.);
   ROOT::TThreadedObject<TH2D> h_dv_aco("h_dv_aco","distance v VS aco",600,-3.2,3.2,200,-1.,1.);

   ROOT::TThreadedObject<TH1D> h_hist_distance0("h_hist_distance0","Distance between stub1 and stub0 in module0",200,-1.,1.);
   ROOT::TThreadedObject<TH1D> h_hist_distance_zoom0("h_hist_distance_zoom0","Distance between stub0 and stub1 in module0 zoomed",80,-0.4,0.4);
   ROOT::TThreadedObject<TH1D> h_hist_distance2("h_hist_distance2","Distance between stub1 and stub0 in module",200,-1.,1.);
   ROOT::TThreadedObject<TH1D> h_hist_distance_zoom2("h_hist_distance_zoom2","Distance between stub0 and stub1 in module2 zoomed",80,-0.4,0.4);

   ROOT::TThreadedObject<TH1D> h_nstubs_track_e("h_nstubs_track_e","nstubs track e",8,0,8);
   ROOT::TThreadedObject<TH1D> h_nstubs_track_mu("h_nstubs_track_mu","nstubs track mu",8,0,8);

   ROOT::TThreadedObject<TH1D> h_nstubs_s1("h_nstubs_s1","nstubs station1",18,0,18);
   ROOT::TThreadedObject<TH1D> h_nstubs_pre_s1("h_nstubs_pre_s1","nstubs station1 pre selection",45,0,45);

   ROOT::TThreadedObject<TH1D> count("count","count",2,0.,2.);

   ROOT::TThreadedObject<TH1D> vrtx_chi2("vrtx_chi2","vrtx_chi2",80,0,20);
   ROOT::TThreadedObject<TH1D> track_el("track_el","track el chi2",80,0,20);
   ROOT::TThreadedObject<TH1D> track_mu("track_mu","track mu chi2",80,0,20);

   ROOT::TThreadedObject<TH1D> sig("sig","sig",2,0,2);
   ROOT::TThreadedObject<TH1D> bkg_e("bkg_e","bkg_e",45,0,45);
   ROOT::TThreadedObject<TH1D> bkg_mu("bkg_mu","bkg_mu",45,0,45);
   ROOT::TThreadedObject<TH1D> sig_sel("sig_sel","sig_sel",2,0,2);
   ROOT::TThreadedObject<TH1D> count_lastSil("count_lastSil","count_lastSil",2,0,2);
   ROOT::TThreadedObject<TH1D> count_lastSil_sel("count_lastSil_sel","count_lastSil_sel",2,0,2);

auto trackatZ = [](double q, double m,double z) {return q + (z ) * m;};
double z_mod[6]={911.2+18.0218,911.2+21.8693,911.2+55.3635,911.2+56.6205,911.2+89.9218,911.2+93.7693};
double alpha[6]={0  *TMath::DegToRad(), 90 *TMath::DegToRad(),135*TMath::DegToRad(), 45 *TMath::DegToRad(),0  *TMath::DegToRad(), 90 *TMath::DegToRad()};

 auto myFunction = [&](TTreeReader &myReader) {


     TTreeReaderArray<MUonETrack> MCtracksRV(myReader, "MCTrack");
     TTreeReaderValue<std::vector<MUonERecoOutputTrack>> RVtracks(myReader, "ReconstructedTracks");
     TTreeReaderValue<MUonERecoOutputVertex> vrtx(myReader, "BestVertex");
     TTreeReaderValue<std::vector<MUonERecoOutputHit>> RVstubs(myReader, "ReconstructedHits");

     while (myReader.Next()) {

double reco_v=0.; double more_reco_v=0.; double reco0_v=0.;
int yes2=0; int yes_v=0;
int code_mu=-99; int code_e=-99; int code_mu_in=-99;
double z_fix=912.7;


 double chi=vrtx->chi2perDegreeOfFreedom();

int yes_mu=0;
int yes_e=0;
double th_inx,th_iny,x0_in,y0_in;
double chi2_muin;
int sec0=0;
int sec1=0;
int stubs_muin=0.;
double th_muin=0.;


double the_rec=-99.;
double thmu_rec=-99.;
TVector3 p_e,p_mu,p_muin;

        std::vector<MUonERecoOutputHit> hits_mu=vrtx->outgoingMuon().hits();
        std::vector<double> pos_mu; pos_mu.resize(6);
        for(int p=0;p<hits_mu.size();p++)pos_mu.push_back(hits_mu.at(p).position());
        std::vector<MUonERecoOutputHit> hits_e=vrtx->outgoingElectron().hits();
        std::vector<double> pos_e; pos_e.resize(6);
        for(int p=0;p<hits_e.size();p++)pos_e.push_back(hits_e.at(p).position());


         auto tracks = *RVtracks;
         for (auto&& track : tracks) {
        if(track.sector()==0) sec0++;
        if(track.sector()==1) sec1++;
        }

 MUonERecoOutputTrack mu_in_v = vrtx->incomingMuon();
 MUonERecoOutputTrack mu_out_v = vrtx->outgoingMuon();
 MUonERecoOutputTrack e_out_v = vrtx->outgoingElectron();
        TVector3 p_muin_v(mu_in_v.xSlope(),mu_in_v.ySlope(),1.0);
        TVector3 p_mu_v(mu_out_v.xSlope(),mu_out_v.ySlope(),1.0);
        TVector3 p_e_v(e_out_v.xSlope(),e_out_v.ySlope(),1.0);
        p_e_v.Unit();p_mu_v.Unit();p_muin_v.Unit();

MUonERecoOutputTrack t_mu;
MUonERecoOutputTrack t_e;

         for (auto&& track : tracks) {
        std::vector<MUonERecoOutputHit> hits_=track.hits();

        if(track.sector()==0 and sec0==1){
	stubs_muin=hits_.size();
        th_inx=track.xSlope();
        th_iny=track.ySlope();
        x0_in=track.x0();
        y0_in=track.y0();
        chi2_muin=track.chi2perDegreeOfFreedom();
        p_muin.SetXYZ(th_inx,th_iny,1.0);
        p_muin=p_muin.Unit();
        th_muin=p_muin.Theta();
	h_chi2->Fill(chi2_muin);
	h_xslope->Fill(th_inx);
                        }
	if(track.sector()==1)
	{
	std::vector<double> pos; pos.resize(6);

        for(int p=0;p<hits_.size();p++){pos.push_back(hits_.at(p).position());
      //if(hits_.at(p).moduleID()==0 or hits_.at(p).moduleID()==4) residualX_pre->Fill(hits_.at(p).position()-trackatZ(track.x0(),track.xSlope(),hits_.at(p).z()));
      //if(hits_.at(p).moduleID()==1 or hits_.at(p).moduleID()==5) residualY_pre->Fill(hits_.at(p).position()-trackatZ(track.y0(),track.ySlope(),hits_.at(p).z()));
	}

         if(std::equal(pos.begin(),pos.end(),pos_mu.begin())){p_mu.SetXYZ(track.xSlope(),track.ySlope(),1.0); p_mu.Unit(); thmu_rec=p_mu.Angle(p_muin);t_mu=track;}
        else if(std::equal(pos.begin(),pos.end(),pos_e.begin())){p_e.SetXYZ(track.xSlope(),track.ySlope(),1.0); p_e.Unit(); the_rec=p_e.Angle(p_muin);t_e=track;}



	yes2++;}
}

double posxIN=99.;//pos_on_track(x0_in,th_inx,z_fix);
double posyIN=99.;//pos_on_track(y0_in,th_iny,z_fix);

int stub0 = 0;
int stub1 = 0;

         auto stubs = *RVstubs;
         std::array<int,6> module_st1_2;module_st1_2.fill({0});
         std::array<int,6> module_st1;module_st1.fill({0});
        std::array<std::vector<float>,6> localX{{std::vector<float>(0.0)}};
         for (auto&& stub : stubs) {
                if(stub.stationID()==0){stub0++; if(stub.moduleID()==4){posxIN=stub.position(); } else if(stub.moduleID()==5){posyIN=stub.position();}   }
                if(stub.stationID()==1){stub1++; module_st1_2.at(stub.moduleID())+=1; module_st1.at(stub.moduleID())=1;int link= stub.moduleID(); localX.at(link).push_back(stub.position());}
                //if(stub.stationID()==1){stub1++; module_st1_2.at(stub.moduleID())+=1;int link= stub.moduleID(); localX.at(link).push_back(stub.position());}
        }


if(chi!=0){
        h_e_pre->Fill(the_rec);
        h_e_post->Fill(vrtx->electronTheta());
        h_mu_pre->Fill(thmu_rec);
        h_mu_post->Fill(vrtx->muonTheta());
        h_inc_pre->Fill(th_muin);
        h_inc_post->Fill(p_muin_v.Theta());
/*for(int p=0;p<hits_e.size();p++){
//	if(hits_e.at(p).moduleID()==0 or hits_e.at(p).moduleID()==4)residualX->Fill(hits_e.at(p).perpendicularResiduum());
//	if(hits_e.at(p).moduleID()==1 or hits_e.at(p).moduleID()==5)residualY->Fill(hits_e.at(p).perpendicularResiduum());

      if(hits_e.at(p).moduleID()==0 or hits_e.at(p).moduleID()==4) residualX_post->Fill(hits_e.at(p).position()-trackatZ(vrtx->outgoingElectron().x0(),vrtx->outgoingElectron().xSlope(),hits_e.at(p).z()));
      if(hits_e.at(p).moduleID()==1 or hits_e.at(p).moduleID()==5) residualY_post->Fill(hits_e.at(p).position()-trackatZ(vrtx->outgoingElectron().y0(),vrtx->outgoingElectron().ySlope(),hits_e.at(p).z()));


 }*/
}


 if(sec0==1 and stubs_muin==6 and abs(posxIN)<=1.5 and abs(posyIN)<=1.5 and chi2_muin<=2 and stub0==6 and th_muin<0.004){// and stub1<=15){
count->Fill(1.);


        for (auto&& mc_tracks : MCtracksRV) {
                if(mc_tracks.interactionID()==9 and mc_tracks.startZ()<(906.) and mc_tracks.startZ()>(902.)) { count_lastSil->Fill(1.); }
                else if(mc_tracks.interactionID()==9 and mc_tracks.startZ()<(914.) and mc_tracks.startZ()>(909.)) { sig->Fill(1.); }
        }



//bool allmod=std::all_of(std::begin(module_st1), std::end(module_st1), [](int i){return i==1;});
bool allmod2=std::all_of(std::begin(module_st1_2), std::end(module_st1_2), [](int i){return i==2;});
double u=99.;
double v=99.;


  if(allmod2 and stub1==12){
double tilt[6]={0.233,0.233,0.,0.,0.233,0.233};
for(int m=0;m<6;m++){
        double real_pos1= cos(tilt[m])*localX.at(m).at(1);
        double real_pos0= cos(tilt[m])*localX.at(m).at(0);
                            if(m==0){ h_hist_distance0->Fill(real_pos1-real_pos0);
                             h_hist_distance_zoom0->Fill(real_pos1-real_pos0);
                             }

                            if(m==2){ h_hist_distance2->Fill(real_pos1-real_pos0);
                             h_hist_distance_zoom2->Fill(real_pos1-real_pos0); u=real_pos1-real_pos0;
                             }

                            if(m==3){
				v=real_pos1-real_pos0;
                             }
                } //end for

}

if(chi!=0 and the_rec!=-99 and thmu_rec!=-99){

        h_z_pos_pre->Fill(vrtx->zPositionFit());


                                                double dotProduct_v = p_mu.Dot(p_e);
                                                TVector3 crossProduct_v = p_mu.Cross(p_e);
                                                double T_v = p_muin.Dot(crossProduct_v);
                                                TVector3 im_v= p_muin.Cross(p_mu);
                                                TVector3 ie_v= p_muin.Cross(p_e);
                                                T_v = T_v>0? 1:-1;
                                                double acoplanarity_v= T_v*(TMath::Pi()- acos( ((im_v).Dot(ie_v))/(im_v.Mag()*ie_v.Mag()) ));


 d_aco_pre->Fill(acoplanarity_v);
 h_2d_pre->Fill(the_rec,thmu_rec);
 theta_mu_pre->Fill(thmu_rec);
 theta_e_pre->Fill(the_rec);
 h_opening_pre->Fill(p_mu.Angle(p_e));
 h_nstubs_pre_s1->Fill(stub1);

double Elastic=0.5109989461*0.001*((1+(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*cos(the_rec)*cos(the_rec))/(1-(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*cos(the_rec)*cos(the_rec)));
double Elastic2=asin( (sin(the_rec)*sqrt(Elastic*Elastic-0.5109989461*0.001*0.5109989461*0.001))/sqrt( (160+0.5109989461*0.001-Elastic)*(160+0.5109989461*0.001-Elastic)-105.6583745 *0.001*105.6583745 *0.001 ) );

bool allmod=std::all_of(std::begin(module_st1), std::end(module_st1), [](int i){return i==1;});

int mod2=0;
double mod2_mu=99.;

//  if(allmod and abs(acoplanarity_v)<=0.4 and chi<20 and thmu_rec>0.0002 and stub1<=20 and thmu_rec<=Elastic2+0.0002 and thmu_rec>=Elastic2-0.0002 and the_rec<0.025 and the_rec>=0.003 and thmu_rec<=0.003){// and vrtx->zPositionFit()<917 and vrtx->zPositionFit()>907){
double zpos =vrtx->zPositionFit()-1.5;

//if(chi>0 and stub1<=15 and thmu_rec>=0.0003 and thmu_rec<=0.0011 and the_rec<0.020 and the_rec>=0.005 and thmu_rec<=Elastic2+0.0002 and thmu_rec>=Elastic2-0.0002 and p_mu.Angle(p_e)>0.005 and abs(acoplanarity_v)<=0.4 and zpos<(911.2+3) and zpos>(911.2-3)){
if(chi>0 and stub1<=15 and thmu_rec>=0.0002 and the_rec<0.020 and the_rec>=0.005 and thmu_rec<=Elastic2+0.0002 and thmu_rec>=Elastic2-0.0002 and p_mu.Angle(p_e)>0.005 and abs(acoplanarity_v)<=0.4 and zpos>906 ){

// if(allmod and chi<150 and stub1<=20 and thmu_rec>0.0002 and thmu_rec<=0.002 and the_rec<0.030 and the_rec>=0.003 and thmu_rec<=Elastic2+0.0002 and thmu_rec>=Elastic2-0.0002 and p_mu.Angle(p_e)>0.005 and abs(acoplanarity_v)<=0.4){
// if(allmod and t_e.chi2perDegreeOfFreedom()<30 and t_mu.chi2perDegreeOfFreedom()<10 and stub1<=150 and thmu_rec>0.0002 and thmu_rec<=0.002 and the_rec<0.030 and the_rec>=0.003 and thmu_rec<=Elastic2+0.0002 and thmu_rec>=Elastic2-0.0002 and p_mu.Angle(p_e)>0.005 and abs(acoplanarity_v)<=0.4){
// if(allmod and chi<20 and stub1<=20 and thmu_rec>0.0002 and thmu_rec<=0.002 and the_rec<0.030 and the_rec>=0.003 and thmu_rec<=Elastic2+0.0002 and thmu_rec>=Elastic2-0.0002 and p_mu.Angle(p_e)>0.005 and abs(acoplanarity_v)<=0.4){
// if(allmod and chi<20 and stub1<=20 and thmu_rec>0.0002 and thmu_rec<=0.002 and the_rec<0.015 and the_rec>=0.003 and thmu_rec<=Elastic2+0.0002 and thmu_rec>=Elastic2-0.0002 and p_mu.Angle(p_e)>0.005 and abs(acoplanarity_v)<=0.4){

// if(abs(u)>0.17 and abs(v)>0.17 and allmod2 and stub1==12 and thmu_rec>0.0005 and thmu_rec<=0.002 and the_rec<0.010 and the_rec>=0.003){
// if(allmod2 and stub1==12 and thmu_rec>0.0005 and thmu_rec<=0.002 and the_rec<0.010 and the_rec>=0.003 and thmu_rec<=Elastic2+0.0002 and thmu_rec>=Elastic2-0.0002 and p_mu.Angle(p_e)>0.005 and abs(acoplanarity_v)<=0.4){
//if(allmod and stub1==12 and abs(u)<0.2 ){//and thmu_rec>0.0005 and thmu_rec<=0.002 and the_rec<0.010 and the_rec>=0.003 and thmu_rec<=Elastic2+0.0002 and thmu_rec>=Elastic2-0.0002 and p_mu.Angle(p_e)>0.005){

//first
 if(yes2>=2){


vrtx_chi2->Fill(chi);
track_el->Fill(t_e.chi2perDegreeOfFreedom());
track_mu->Fill(t_mu.chi2perDegreeOfFreedom());


          double tracklocxy0 = (t_e.x0() + t_e.xSlope()*z_mod[2])*cos(alpha[2]) +
                               (t_e.y0() + t_e.ySlope()*z_mod[2])*sin(alpha[2]);



        h_nstubs_track_e->Fill(t_e.hits().size());
        h_nstubs_track_mu->Fill(t_mu.hits().size());

std::array<double,6> t_e_hits;t_e_hits.fill({-99.});
 std::array<double,6> t_mu_hits;t_mu_hits.fill({-99.});

 for(int h=0; h<t_mu.hits().size(); h++){
  if(t_mu.hits().at(h).stationID()==1) t_mu_hits.at(t_mu.hits().at(h).moduleID())=t_mu.hits().at(h).position();
 }

 for(int h=0; h<t_e.hits().size(); h++){
  if(t_e.hits().at(h).stationID()==1) t_e_hits.at(t_e.hits().at(h).moduleID())=t_e.hits().at(h).position();
 }



//analog
                             /*h_hist_distance0->Fill(t_e_hits.at(0)-t_mu_hits.at(0));
                             h_hist_distance_zoom0->Fill(t_e_hits.at(0)-t_mu_hits.at(0));
                             h_hist_distance2->Fill(t_e_hits.at(2)-t_mu_hits.at(2));
                             h_hist_distance_zoom2->Fill(t_e_hits.at(2)-t_mu_hits.at(2));*/


 if( t_e_hits.at(0)!=-99. and t_mu_hits.at(0)!=-99.) h_dist0->Fill(t_e_hits.at(0)-t_mu_hits.at(0));
 if( t_e_hits.at(1)!=-99. and t_mu_hits.at(1)!=-99.) h_dist1->Fill(t_e_hits.at(1)-t_mu_hits.at(1));
 if( t_e_hits.at(2)!=-99. and t_mu_hits.at(2)!=-99.) h_dist2->Fill(t_e_hits.at(2)-t_mu_hits.at(2));
 if( t_e_hits.at(3)!=-99. and t_mu_hits.at(3)!=-99.) h_dist3->Fill(t_e_hits.at(3)-t_mu_hits.at(3));
 if( t_e_hits.at(4)!=-99. and t_mu_hits.at(4)!=-99.) h_dist4->Fill(t_e_hits.at(4)-t_mu_hits.at(4));
 if( t_e_hits.at(5)!=-99. and t_mu_hits.at(5)!=-99.) h_dist5->Fill(t_e_hits.at(5)-t_mu_hits.at(5));


     if( t_e_hits.at(2)!=-99. and t_mu_hits.at(2)!=-99.)h_du_aco->Fill(acoplanarity_v,t_e_hits.at(2)-t_mu_hits.at(2));
     if( t_e_hits.at(3)!=-99. and t_mu_hits.at(3)!=-99.)h_dv_aco->Fill(acoplanarity_v,t_e_hits.at(3)-t_mu_hits.at(3));

/*
if( the_rec<=0.005 ){ if( t_e_hits.at(2)!=-99. and t_mu_hits.at(2)!=-99.){h_dist_mod2_range0->Fill(t_e_hits.at(2)-t_mu_hits.at(2)); h_pos_range0->Fill(t_e_hits.at(2),t_mu_hits.at(2));
                                                                         h_pos_el_range0->Fill(t_e_hits.at(2)); h_pos_mu_range0->Fill(t_mu_hits.at(2));}
                      if( t_e_hits.at(2)!=-99. and t_mu_hits.at(2)!=-99.){h_dist_mod0_range0->Fill(t_e_hits.at(0)-t_mu_hits.at(0));}}

if( the_rec>0.005 and the_rec<=0.0075 ){ if( t_e_hits.at(2)!=-99. and t_mu_hits.at(2)!=-99.){ h_dist_mod2_range1->Fill(t_e_hits.at(2)-t_mu_hits.at(2));h_pos_range1->Fill(t_e_hits.at(2),t_mu_hits.at(2));
                                                                         h_pos_el_range1->Fill(t_e_hits.at(2)); h_pos_mu_range1->Fill(t_mu_hits.at(2));}
                      if( t_e_hits.at(2)!=-99. and t_mu_hits.at(2)!=-99.){h_dist_mod0_range1->Fill(t_e_hits.at(0)-t_mu_hits.at(0));}}

if( the_rec>0.0075 and the_rec<=0.01 ){ if( t_e_hits.at(2)!=-99. and t_mu_hits.at(2)!=-99.){ h_dist_mod2_range2->Fill(t_e_hits.at(2)-t_mu_hits.at(2));h_pos_range2->Fill(t_e_hits.at(2),t_mu_hits.at(2));
                                                                         h_pos_el_range2->Fill(t_e_hits.at(2)); h_pos_mu_range2->Fill(t_mu_hits.at(2));}
                      if( t_e_hits.at(2)!=-99. and t_mu_hits.at(2)!=-99.){h_dist_mod0_range2->Fill(t_e_hits.at(0)-t_mu_hits.at(0));}}

if( the_rec>0.010 and the_rec<=0.015 ){ if( t_e_hits.at(2)!=-99. and t_mu_hits.at(2)!=-99.){ h_dist_mod2_range3->Fill(t_e_hits.at(2)-t_mu_hits.at(2));h_pos_range3->Fill(t_e_hits.at(2),t_mu_hits.at(2));
                                                                         h_pos_el_range3->Fill(t_e_hits.at(2)); h_pos_mu_range3->Fill(t_mu_hits.at(2));}
                      if( t_e_hits.at(2)!=-99. and t_mu_hits.at(2)!=-99.){h_dist_mod0_range3->Fill(t_e_hits.at(0)-t_mu_hits.at(0));}}

if( the_rec>0.015 and the_rec<=0.025 ){ if( t_e_hits.at(2)!=-99. and t_mu_hits.at(2)!=-99.){ h_dist_mod2_range4->Fill(t_e_hits.at(2)-t_mu_hits.at(2));h_pos_range4->Fill(t_e_hits.at(2),t_mu_hits.at(2));
                                                                         h_pos_el_range4->Fill(t_e_hits.at(2)); h_pos_mu_range4->Fill(t_mu_hits.at(2));}
                      if( t_e_hits.at(2)!=-99. and t_mu_hits.at(2)!=-99.){h_dist_mod0_range4->Fill(t_e_hits.at(0)-t_mu_hits.at(0));}}

if( the_rec>0.025 and the_rec<=0.035 ){ if( t_e_hits.at(2)!=-99. and t_mu_hits.at(2)!=-99.){ h_dist_mod2_range5->Fill(t_e_hits.at(2)-t_mu_hits.at(2));h_pos_range5->Fill(t_e_hits.at(2),t_mu_hits.at(2));
                                                                         h_pos_el_range5->Fill(t_e_hits.at(2)); h_pos_mu_range5->Fill(t_mu_hits.at(2));}
                      if( t_e_hits.at(2)!=-99. and t_mu_hits.at(2)!=-99.){h_dist_mod0_range5->Fill(t_e_hits.at(0)-t_mu_hits.at(0));}}
*/


 h_z_pos->Fill(vrtx->zPositionFit());
 d_aco->Fill(acoplanarity_v);
 h_2d->Fill(the_rec,thmu_rec);
 theta_mu->Fill(thmu_rec);
 theta_e->Fill(the_rec);
 h_opening->Fill(p_mu.Angle(p_e));
 h_nstubs_s1->Fill(stub1);

        for (auto&& mc_tracks : MCtracksRV) {
                if(mc_tracks.interactionID()==9 and mc_tracks.startZ()<(906.) and mc_tracks.startZ()>(902.)) { count_lastSil_sel->Fill(1.); }
        }


/*
if(zpos<(906.) and zpos>(902.))
{
double genz=0.;

        Int_t i=0;
        for (auto&& mc_tracks : MCtracksRV) {
                if(mc_tracks.interactionID()==9 and t_e.linkedTrackID()==i) { genz= mc_tracks.startZ(); cout <<"electron at Z " << mc_tracks.startZ() << endl;}
        i++;
	}

 h_z_pos_pre->Fill(genz);
 h_z_pos->Fill(vrtx->zPositionFit());
cout << "t_mu.processIDofLinkedTrack " << t_mu.processIDofLinkedTrack() << endl;
cout << "t_e.processIDofLinkedTrack " << t_e.processIDofLinkedTrack() << endl;
cout <<"----------"<<endl;
}*/

if( (t_mu.processIDofLinkedTrack()==0 and t_e.processIDofLinkedTrack()==9) or (t_mu.processIDofLinkedTrack()==9 and t_e.processIDofLinkedTrack()==0) ){sig_sel->Fill(1.);}
else{bkg_mu->Fill(t_mu.processIDofLinkedTrack()); bkg_e->Fill(t_e.processIDofLinkedTrack());}

if(thmu_rec>=the_rec-0.00005 and thmu_rec<=the_rec+0.00005)theq->Fill(thmu_rec);
				}//yes2>=2
			}//aco e chi
		}//chi!=0
	}//mu_in
yes2=0;
 } //end of general while
}; //end of ySlopefunction

  tp1.Process(myFunction);

auto count_lastSilM=count_lastSil.Merge();
auto count_lastSil_selM=count_lastSil_sel.Merge();

auto h_nstubs_pre_s1M=h_nstubs_pre_s1.Merge();
auto h_nstubs_s1M=h_nstubs_s1.Merge();


auto vrtx_chi2M=vrtx_chi2.Merge();
auto track_elM=track_el.Merge();
auto track_muM=track_mu.Merge();

auto countM = count.Merge();

auto h_hist_distance0M = h_hist_distance0.Merge();
auto h_hist_distance_zoom0M = h_hist_distance_zoom0.Merge();
auto h_hist_distance2M = h_hist_distance2.Merge();
auto h_hist_distance_zoom2M = h_hist_distance_zoom2.Merge();


auto h_nstubs_track_eM=h_nstubs_track_e.Merge();
auto h_nstubs_track_muM=h_nstubs_track_mu.Merge();

  auto h_dist_mod0_range0M= h_dist_mod0_range0.Merge();
  auto h_dist_mod0_range1M= h_dist_mod0_range1.Merge();
  auto h_dist_mod0_range2M= h_dist_mod0_range2.Merge();
  auto h_dist_mod0_range3M= h_dist_mod0_range3.Merge();
  auto h_dist_mod0_range4M= h_dist_mod0_range4.Merge();
  auto h_dist_mod0_range5M= h_dist_mod0_range5.Merge();

  auto h_du_acoM = h_du_aco.Merge();
  auto h_dv_acoM = h_dv_aco.Merge();


  auto h_dist_mod2_range0M= h_dist_mod2_range0.Merge();
  auto h_dist_mod2_range1M= h_dist_mod2_range1.Merge();
  auto h_dist_mod2_range2M= h_dist_mod2_range2.Merge();
  auto h_dist_mod2_range3M= h_dist_mod2_range3.Merge();
  auto h_dist_mod2_range4M= h_dist_mod2_range4.Merge();
  auto h_dist_mod2_range5M= h_dist_mod2_range5.Merge();

  auto h_dist0M=h_dist0.Merge();
  auto h_dist1M=h_dist1.Merge();
  auto h_dist2M=h_dist2.Merge();
  auto h_dist3M=h_dist3.Merge();
  auto h_dist4M=h_dist4.Merge();
  auto h_dist5M=h_dist5.Merge();

  auto h_pos_el_range0M=h_pos_el_range0.Merge();
  auto h_pos_el_range1M=h_pos_el_range1.Merge();
  auto h_pos_el_range2M=h_pos_el_range2.Merge();
/*  auto h_pos_el_range3M=h_pos_el_range3.Merge();
  auto h_pos_el_range4M=h_pos_el_range4.Merge();
  auto h_pos_el_range5M=h_pos_el_range5.Merge();*/

  auto h_pos_mu_range0M=h_pos_mu_range0.Merge();
  auto h_pos_mu_range1M=h_pos_mu_range1.Merge();
  auto h_pos_mu_range2M=h_pos_mu_range2.Merge();
/*  auto h_pos_mu_range3M=h_pos_mu_range3.Merge();
  auto h_pos_mu_range4M=h_pos_mu_range4.Merge();
  auto h_pos_mu_range5M=h_pos_mu_range5.Merge();*/

  auto h_pos_range0M=h_pos_range0.Merge();
  auto h_pos_range1M=h_pos_range1.Merge();
  auto h_pos_range2M=h_pos_range2.Merge();

  auto h_2dMerged = h_2d.Merge();
  auto theta_muMerged = theta_mu.Merge();
  auto theta_eMerged = theta_e.Merge();
  auto d_acoM = d_aco.Merge();
  auto h_openingM = h_opening.Merge();

  auto h_2d_preMerged = h_2d_pre.Merge();
  auto theta_mu_preMerged = theta_mu_pre.Merge();
  auto theta_e_preMerged = theta_e_pre.Merge();
  auto d_aco_preM = d_aco_pre.Merge();
  auto h_opening_preM = h_opening_pre.Merge();
  auto theqM = theq.Merge();
  auto h_z_pos_preM = h_z_pos_pre.Merge();
  auto h_z_posM = h_z_pos.Merge();
  auto h_chi2M = h_chi2.Merge();
  auto h_xslopeM = h_xslope.Merge();

  auto h_e_preM = h_e_pre.Merge();
  auto h_e_postM = h_e_post.Merge();
  auto h_mu_preM = h_mu_pre.Merge();
  auto h_mu_postM = h_mu_post.Merge();
  auto h_inc_preM = h_inc_pre.Merge();
  auto h_inc_postM = h_inc_post.Merge();
 


  auto residualU_criticM = residualU_critic.Merge();
  auto residualU_no_criticM =residualU_no_critic.Merge();
  auto residualUM_post = residualU.Merge();
  auto residualU_noM_post =residualU_no.Merge();

/*
TCanvas t("t","t",700,700);
t.Divide(2,2);
t.cd(1);
residualU_criticM->Draw("hist");
t.cd(2);
residualU_no_criticM->Draw("hist");
t.cd(3);
residualUM_post->Draw("hist");
t.cd(4);
residualU_noM_post->Draw("hist");
t.SaveAs(Form("comparison_RDMC/%s_RD_muin_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.pdf",version.c_str(),bend.c_str(),nhits,mc.c_str()));
*/

auto sig_selM=sig_sel.Merge();
auto sigM = sig.Merge();
auto bkg_eM = bkg_e.Merge();
auto bkg_muM = bkg_mu.Merge();

cout << "countM.Integral() " << countM->Integral() << endl;
cout << "bkg_eM.Integral() " << bkg_eM->Integral() << endl;
cout << "bkg_muM.Integral() " << bkg_muM->Integral() << endl;
cout << "sig_selM.Integral() " << sig_selM->Integral() << endl;
cout << "sigM.Integral() " << sigM->Integral() << endl;
cout << "count_lastSil.Integral() " << count_lastSilM->Integral() << endl;
cout << "count_lastSil_sel.Integral() " << count_lastSil_selM->Integral() << endl;

TCanvas bs("bs","bs",1400,700);
bs.Divide(2,1);
bs.cd(1);
bkg_muM->Draw("hist");
bs.cd(2);
bkg_eM->Draw("hist");
bs.SaveAs(Form("comparison_RDMC/%s_bkg_link_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));

vrtx_chi2M->SaveAs(Form("comparison_RDMC/%s_vrtx_chi2_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));
track_elM->SaveAs(Form("comparison_RDMC/%s_track_el_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));
track_muM->SaveAs(Form("comparison_RDMC/%s_track_mu_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));


h_2dMerged->SaveAs(Form("comparison_RDMC/%s_2D_RD_parallel_pre_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));

theta_muMerged->SaveAs(Form("comparison_RDMC/%s_theta_mu_RD_parallel_pre_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));

theta_eMerged->SaveAs(Form("comparison_RDMC/%s_theta_e_RD_parallel_pre_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));

d_acoM->SaveAs(Form("comparison_RDMC/%s_d_aco_RD_parallel_pre_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));

h_openingM->SaveAs(Form("comparison_RDMC/%s_opening_RD_parallel_pre_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));

h_z_pos_preM->SaveAs(Form("comparison_RDMC/%s_z_pre_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));

h_nstubs_pre_s1M->SaveAs(Form("comparison_RDMC/%s_nstubs_pre_s1_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));


h_2d_preMerged->SaveAs(Form("comparison_RDMC/%s_preCuts_2D_RD_parallel_pre_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));

theta_mu_preMerged->SaveAs(Form("comparison_RDMC/%s_preCuts_theta_mu_RD_parallel_pre_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));

theta_e_preMerged->SaveAs(Form("comparison_RDMC/%s_preCuts_theta_e_RD_parallel_pre_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));

d_aco_preM->SaveAs(Form("comparison_RDMC/%s_preCuts_d_aco_RD_parallel_pre_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));

h_opening_preM->SaveAs(Form("comparison_RDMC/%s_preCuts_opening_RD_parallel_pre_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));

h_z_posM->SaveAs(Form("comparison_RDMC/%s_z_post_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));

h_nstubs_s1M->SaveAs(Form("comparison_RDMC/%s_nstubs_s1_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));


TCanvas z("z","z",800,800);
h_z_posM->Draw("hist");
h_z_pos_preM->SetLineColor(kRed);
h_z_pos_preM->Draw("hist same");
z.SaveAs("minbias_z.pdf");

//theqM->SaveAs(Form("comparison_RDMC/%s_th_eqM_pre_Z.pdf",version.c_str(),bend.c_str(),nhits,mc.c_str()));

/*TCanvas a("a","a",700,700);
a.Divide(2,2);
a.cd(1);
h_z_pos_preM->SetLineColor(kOrange+10);
h_z_pos_preM->Draw("hist");
a.cd(2);
h_z_posM->SetLineColor(kOrange+10);
h_z_posM->Draw("hist");
a.SaveAs(Form("comparison_RDMC/%s_data_pos_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.pdf",version.c_str(),bend.c_str(),nhits,mc.c_str()));
*/


/*h_dist_mod2_range0M->SaveAs(Form("comparison_RDMC/%s_h_dist_mod2_range_0_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));
h_dist_mod2_range1M->SaveAs(Form("comparison_RDMC/%s_h_dist_mod2_range_1_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));
h_dist_mod2_range2M->SaveAs(Form("comparison_RDMC/%s_h_dist_mod2_range_2_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));
//h_dist_mod2_range3M->SaveAs(Form("comparison_RDMC/%s_h_dist_mod2_range_3_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));
//h_dist_mod2_range4M->SaveAs(Form("comparison_RDMC/%s_h_dist_mod2_range_4_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));
//h_dist_mod2_range5M->SaveAs(Form("comparison_RDMC/%s_h_dist_mod2_range_5_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));

h_dist_mod0_range0M->SaveAs(Form("comparison_RDMC/%s_h_dist_mod0_range_0_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));
h_dist_mod0_range1M->SaveAs(Form("comparison_RDMC/%s_h_dist_mod0_range_1_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));
h_dist_mod0_range2M->SaveAs(Form("comparison_RDMC/%s_h_dist_mod0_range_2_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));
//h_dist_mod0_range3M->SaveAs(Form("comparison_RDMC/%s_h_dist_mod0_range_3_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));
//h_dist_mod0_range4M->SaveAs(Form("comparison_RDMC/%s_h_dist_mod0_range_4_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));
//h_dist_mod0_range5M->SaveAs(Form("comparison_RDMC/%s_h_dist_mod0_range_5_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));
*/


 auto legend_e_mod0 = new TLegend(0.6,0.6,0.9,0.9);
legend_e_mod0->AddEntry("h_dist_mod0_range0M","0-5 mrad kBlue","L");
legend_e_mod0->AddEntry("h_dist_mod0_range1M","5-7.5 mrad kRed","L");
legend_e_mod0->AddEntry("h_dist_mod0_range2M","7.5-10 mrad kOrange","L");
legend_e_mod0->AddEntry("h_dist_mod0_range3M","10-15 mrad kGreen","L");
legend_e_mod0->AddEntry("h_dist_mod0_range4M","15-25 mrad kViolet","L");
legend_e_mod0->AddEntry("h_dist_mod0_range5M","25-35 mrad kBlack","L");

 auto legend_e_mod2 = new TLegend(0.6,0.6,0.9,0.9);
legend_e_mod2->AddEntry("h_dist_mod2_range0M","0-5 mrad kBlue","L");
legend_e_mod2->AddEntry("h_dist_mod2_range1M","5-7.5 mrad kRed","L");
legend_e_mod2->AddEntry("h_dist_mod2_range2M","7.5-10 mrad kOrange","L");
legend_e_mod2->AddEntry("h_dist_mod2_range3M","10-15 mrad kGreen","L");
legend_e_mod2->AddEntry("h_dist_mod2_range4M","15-25 mrad kViolet","L");
legend_e_mod2->AddEntry("h_dist_mod2_range5M","25-35 mrad kBlack","L");

/*TCanvas c("c","c",900,1800);
c.Divide(1,2);
c.cd(1);
h_dist_mod0_range0M->SetLineColor(kBlue);
h_dist_mod0_range1M->SetLineColor(kRed);
h_dist_mod0_range2M->SetLineColor(kOrange);
//h_dist_mod0_range3M->SetLineColor(kGreen);
//h_dist_mod0_range4M->SetLineColor(kViolet);
//h_dist_mod0_range5M->SetLineColor(kBlack);

//h_dist_mod0_range4M->Draw("hist");
h_dist_mod0_range2M->Draw("hist");
h_dist_mod0_range1M->Draw("hist same");
h_dist_mod0_range0M->Draw("hist same");
//h_dist_mod0_range3M->Draw("hist same");
//h_dist_mod0_range5M->Draw("hist same");

legend_e_mod0->Draw("same");
c.cd(2);
h_dist_mod2_range0M->SetLineColor(kBlue);
h_dist_mod2_range1M->SetLineColor(kRed);
h_dist_mod2_range2M->SetLineColor(kOrange);
//h_dist_mod2_range3M->SetLineColor(kGreen);
//h_dist_mod2_range4M->SetLineColor(kViolet);
//h_dist_mod2_range5M->SetLineColor(kBlack);

//h_dist_mod2_range4M->Draw("hist");
h_dist_mod2_range0M->Draw("hist");
h_dist_mod2_range1M->Draw("hist same");
h_dist_mod2_range2M->Draw("hist same");
//h_dist_mod2_range3M->Draw("hist same");
//h_dist_mod2_range5M->Draw("hist same");

legend_e_mod2->Draw("same");
c.SaveAs(Form("comparison_RDMC/%s_h_dist_ranges_el_mu_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.pdf",version.c_str(),bend.c_str(),nhits,mc.c_str()));

TCanvas p("p","p",2100,2100);
p.Divide(3,3);
p.cd(1);
h_pos_range0M->Draw("COLZ");
p.cd(2);
h_pos_el_range0M->Draw("hist");
p.cd(3);
h_pos_mu_range0M->Draw("hist");
p.cd(4);
h_pos_range1M->Draw("COLZ");
p.cd(5);
h_pos_el_range1M->Draw("hist");
p.cd(6);
h_pos_mu_range1M->Draw("hist");
p.cd(7);
h_pos_range2M->Draw("COLZ");
p.cd(8);
h_pos_el_range2M->Draw("hist");
p.cd(9);
h_pos_mu_range2M->Draw("hist");
p.SaveAs(Form("comparison_RDMC/%s_pos_MOD2_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.pdf",version.c_str(),bend.c_str(),nhits,mc.c_str()));
*/

TCanvas d("d","d",1400,2100);
d.Divide(2,3);
d.cd(1);
h_dist0M->Draw("hist");
h_dist0M->SaveAs(Form("comparison_RDMC/%s_distance_MOD0_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));
d.cd(2);
h_dist1M->Draw("hist");
h_dist1M->SaveAs(Form("comparison_RDMC/%s_distance_MOD1_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));
d.cd(3);
h_dist2M->Draw("hist");
h_dist2M->SaveAs(Form("comparison_RDMC/%s_distance_MOD2_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));
d.cd(4);
h_dist3M->Draw("hist");
h_dist3M->SaveAs(Form("comparison_RDMC/%s_distance_MOD3_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));
d.cd(5);
h_dist4M->Draw("hist");
h_dist4M->SaveAs(Form("comparison_RDMC/%s_distance_MOD4_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));
d.cd(6);
h_dist5M->Draw("hist");
h_dist5M->SaveAs(Form("comparison_RDMC/%s_distance_MOD5_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));

d.SaveAs(Form("comparison_RDMC/%s_distance_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.pdf",version.c_str(),bend.c_str(),nhits,mc.c_str()));

TCanvas ad("ad","ad",1400,700);
ad.Divide(2,1);
ad.cd(1);
h_du_acoM->Draw("COLZ");
ad.cd(2);
h_dv_acoM->Draw("COLZ");
ad.SaveAs(Form("comparison_RDMC/%s_dist_VS_aco_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.pdf",version.c_str(),bend.c_str(),nhits,mc.c_str()));


TCanvas n1("n1","n1",1000,1000);
n1.Divide(2,2);
n1.cd(1);
h_hist_distance0M->Draw("hist");
n1.cd(2);
h_hist_distance2M->Draw("hist");
n1.cd(3);
h_hist_distance_zoom0M->Draw("hist");
n1.cd(4);
h_hist_distance_zoom2M->Draw("hist");
n1.SaveAs(Form("comparison_RDMC/%s_distance_stubs_zoom_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.pdf",version.c_str(),bend.c_str(),nhits,mc.c_str()));

TCanvas ns("ns","ns",700,700);
ns.Divide(1,2);
ns.cd(1);
h_nstubs_track_eM->Draw("hist");
ns.cd(2);
h_nstubs_track_muM->Draw("hist");
ns.SaveAs(Form("comparison_RDMC/%s_nstubs_track_basic_chi150_flatArea_the20_the05_thmu02_nstubs15_nochivrtx_zPosNoPeak_el_%s_%dhit_%s.pdf",version.c_str(),bend.c_str(),nhits,mc.c_str()));

return 0;
}

