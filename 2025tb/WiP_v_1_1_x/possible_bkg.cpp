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


int possible_bkg(int run, string type, string filen, int nhits){

  int nthreads = 6;


  ROOT::EnableImplicitMT(nthreads);


TChain * cbmsim = new TChain("cbmsim");
TChain * cbmsim_g = new TChain("cbmsim");



cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/TB2025/run%i/%s/%s_%ihit_WiP_17_6.root",run,type.c_str(),filen.c_str(),nhits));
cout << Form("/mnt/raid10/DATA/espedica/fairmu/TB2025/run%i/%s/%s_%ihit_WiP_17_6.root",run,type.c_str(),filen.c_str(),nhits) << endl;

cout << "cbmsim->GetEntries() " << cbmsim->GetEntries() <<endl;

ROOT::TTreeProcessorMT tp1(*cbmsim,nthreads);

      MUonERecoOutputAnalysis *ReconstructionOutput = 0;

   ROOT::TThreadedObject<TH1D> count_tr_1("count_chi1","count_chi1",2,0,2);
   ROOT::TThreadedObject<TH1D> count_tr_2("count_chi2","count_chi2",2,0,2);
   ROOT::TThreadedObject<TH1D> count_tr_1_sig("count_chi1_sig","count_chi1_sig",2,0,2);
   ROOT::TThreadedObject<TH1D> count_tr_2_sig("count_chi2_sig","count_chi2_sig",2,0,2);
   ROOT::TThreadedObject<TH1D> count_mu_is_mu_0("count_mu_is_mu_0","count_mu_is_mu_0",2,0,2);
   ROOT::TThreadedObject<TH1D> count_e_is_mu_0("count_e_is_mu_0","count_e_is_mu_0",2,0,2);
   ROOT::TThreadedObject<TH1D> count_mu_is_mu_1("count_mu_is_mu_0","count_mu_is_mu_0",2,0,2);
   ROOT::TThreadedObject<TH1D> count_e_is_mu_1("count_e_is_mu_0","count_e_is_mu_0",2,0,2);
   ROOT::TThreadedObject<TH1D> count_mu_is_mu_0_sig("count_mu_is_mu_0_sig","count_mu_is_mu_0",2,0,2);
   ROOT::TThreadedObject<TH1D> count_e_is_mu_0_sig("count_e_is_mu_0_sig","count_e_is_mu_0",2,0,2);
   ROOT::TThreadedObject<TH1D> count_mu_is_mu_1_sig("count_mu_is_mu_0_sig","count_mu_is_mu_0",2,0,2);
   ROOT::TThreadedObject<TH1D> count_e_is_mu_1_sig("count_e_is_mu_0_sig","count_e_is_mu_0",2,0,2);

   ROOT::TThreadedObject<TH1D> h_IP_0("h_IP_0","Impact point outgoing mu - outgoing el tar 0",500,-5.,5.);
   ROOT::TThreadedObject<TH2D> h_IPop_0("h_IPop_0","Impact point outgoing mu - outgoing el VS opening angle tar 0",500,-5.,5.,160,0.,0.032);

   ROOT::TThreadedObject<TH1D> h_nhits_pre_0("h_nhits_pre_0","number hits pre tar 0",34,0,34);
   ROOT::TThreadedObject<TH1D> h_nhits_post_0("h_nhits_post_0","number hits post tar 0",34,0,34);
   ROOT::TThreadedObject<TH1D> h_ntracks_post_0("h_ntracks_post_0","number tracks post tar 0",34,0,34);
   ROOT::TThreadedObject<TH1D> h_chi2_vrtx_0("h_chi2_vrtx_0","chi2 vrtx tar 0",300,0.,30.);
   ROOT::TThreadedObject<TH1D> h_Dz_vrtx_0("h_Dz_vrtx_0","(z_vrtx-z_tar) vrtx tar 0",500,-2000.,2000.);
   ROOT::TThreadedObject<TH1D> h_aco_0("h_aco_0","acoplanarity tar 0",200,-1.,1.);
   ROOT::TThreadedObject<TH1D> h_th_min_0("h_th_min_0","theta min tar 0",250,0.,0.005);
   ROOT::TThreadedObject<TH1D> h_th_max_0("h_th_max_0","theta max tar 0",160,0.,0.032);
   ROOT::TThreadedObject<TH1D> h_mf_0("h_mf_0","Nstubs muon filter tar 0",30,0,30);
   ROOT::TThreadedObject<TH1D> h_mf_mod_0("h_mf_mod_0","Nstubs per mod muon filter tar 0",5,0,5);
   ROOT::TThreadedObject<TH1D> h_mu_mf_0("h_mu_mf_0","Nstubs muon track in MF tar 0",30,0,30);
   ROOT::TThreadedObject<TH1D> h_e_mf_0("h_e_mf_0","Nstubs electron track in MF tar 0",30,0,30);
   ROOT::TThreadedObject<TH1D> h_phi_e_0("h_phi_e_0","phi electron tar 0",180,-180,180);
   ROOT::TThreadedObject<TH1D> h_phi_mu_0("h_phi_mu_0","phi muon tar 0",180,-180,180);

   ROOT::TThreadedObject<TH1D> h_nhits_pre_0_sig("h_nhits_pre_0_sig","number hits pre tar 0",34,0,34);
   ROOT::TThreadedObject<TH1D> h_nhits_post_0_sig("h_nhits_post_0_sig","number hits post tar 0",34,0,34);
   ROOT::TThreadedObject<TH1D> h_ntracks_post_0_sig("h_ntracks_post_0_sig","number tracks post tar 0",34,0,34);
   ROOT::TThreadedObject<TH1D> h_chi2_vrtx_0_sig("h_chi2_vrtx_0_sig","chi2 vrtx tar 0",300,0.,30.);
   ROOT::TThreadedObject<TH1D> h_Dz_vrtx_0_sig("h_Dz_vrtx_0_sig","(z_vrtx-z_tar) vrtx tar 0",300,-30.,30.);
   ROOT::TThreadedObject<TH1D> h_aco_0_sig("h_aco_0_sig","acoplanarity tar 0",200,-1.,1.);
   ROOT::TThreadedObject<TH1D> h_th_min_0_sig("h_th_min_0_sig","theta min tar 0",250,0.,0.005);
   ROOT::TThreadedObject<TH1D> h_th_max_0_sig("h_th_max_0_sig","theta max tar 0",160,0.,0.032);
   ROOT::TThreadedObject<TH1D> h_mf_0_sig("h_mf_0_sig","Nstubs muon filter tar 0",30,0,30);
   ROOT::TThreadedObject<TH1D> h_mf_mod_0_sig("h_mf_mod_0_sig","Nstubs per mod muon filter tar 0",5,0,5);
   ROOT::TThreadedObject<TH1D> h_mu_mf_0_sig("h_mu_mf_0_sig","Nstubs muon track in MF tar 0",30,0,30);
   ROOT::TThreadedObject<TH1D> h_e_mf_0_sig("h_e_mf_0_sig","Nstubs electron track in MF tar 0",30,0,30);
   ROOT::TThreadedObject<TH1D> h_phi_e_0_sig("h_phi_e_0_sig","phi electron tar 0",180,-180,180);
   ROOT::TThreadedObject<TH1D> h_phi_mu_0_sig("h_phi_mu_0_sig","phi muon tar 0",180,-180,180);


   ROOT::TThreadedObject<TH1D> h_IP_1("h_IP_1","Impact point outgoing mu - outgoing el tar 1",500,-5.,5.);
   ROOT::TThreadedObject<TH2D> h_IPop_1("h_IPop_1","Impact point outgoing mu - outgoing el VS opening angle tar 1",500,-5.,5.,160,0.,0.032);

   ROOT::TThreadedObject<TH1D> h_nhits_pre_1("h_nhits_pre_1","number hits pre tar 1",34,0,34);
   ROOT::TThreadedObject<TH1D> h_nhits_post_1("h_nhits_post_1","number hits post tar 1",34,0,34);
   ROOT::TThreadedObject<TH1D> h_ntracks_post_1("h_ntracks_post_1","number tracks post tar 1",34,0,34);
   ROOT::TThreadedObject<TH1D> h_chi2_vrtx_1("h_chi2_vrtx_1","chi2 vrtx tar 1",300,0.,30.);
   ROOT::TThreadedObject<TH1D> h_Dz_vrtx_1("h_Dz_vrtx_1","(z_vrtx-z_tar) vrtx tar 1",500,-2000.,2000.);
   ROOT::TThreadedObject<TH1D> h_aco_1("h_aco_1","acoplanarity tar 1",200,-1.,1.);
   ROOT::TThreadedObject<TH1D> h_th_min_1("h_th_min_1","theta min tar 1",250,0.,0.005);
   ROOT::TThreadedObject<TH1D> h_th_max_1("h_th_max_1","theta max tar 1",160,0.,0.032);
   ROOT::TThreadedObject<TH1D> h_mf_1("h_mf_1","Nstubs muon filter tar 1",30,0,30);
   ROOT::TThreadedObject<TH1D> h_mf_mod_1("h_mf_mod_1","Nstubs per mod muon filter tar 1",5,0,5);
   ROOT::TThreadedObject<TH1D> h_mu_mf_1("h_mu_mf_1","Nstubs muon track in MF tar 1",30,0,30);
   ROOT::TThreadedObject<TH1D> h_e_mf_1("h_e_mf_1","Nstubs electron track in MF tar 1",30,0,30);
   ROOT::TThreadedObject<TH1D> h_phi_e_1("h_phi_e_1","phi electron tar 1",180,-180,180);
   ROOT::TThreadedObject<TH1D> h_phi_mu_1("h_phi_mu_1","phi muon tar 1",180,-180,180);

   ROOT::TThreadedObject<TH1D> h_nhits_pre_1_sig("h_nhits_pre_1_sig","number hits pre tar 1",34,0,34);
   ROOT::TThreadedObject<TH1D> h_nhits_post_1_sig("h_nhits_post_1_sig","number hits post tar 1",34,0,34);
   ROOT::TThreadedObject<TH1D> h_ntracks_post_1_sig("h_ntracks_post_1_sig","number tracks post tar 1",34,0,34);
   ROOT::TThreadedObject<TH1D> h_chi2_vrtx_1_sig("h_chi2_vrtx_1_sig","chi2 vrtx tar 1",300,0.,30.);
   ROOT::TThreadedObject<TH1D> h_Dz_vrtx_1_sig("h_Dz_vrtx_1_sig","(z_vrtx-z_tar) vrtx tar 1",300,-30.,30.);
   ROOT::TThreadedObject<TH1D> h_aco_1_sig("h_aco_1_sig","acoplanarity tar 1",200,-1.,1.);
   ROOT::TThreadedObject<TH1D> h_th_min_1_sig("h_th_min_1_sig","theta min tar 1",250,0.,0.005);
   ROOT::TThreadedObject<TH1D> h_th_max_1_sig("h_th_max_1_sig","theta max tar 1",160,0.,0.032);
   ROOT::TThreadedObject<TH1D> h_mf_1_sig("h_mf_1_sig","Nstubs muon filter tar 1",30,0,30);
   ROOT::TThreadedObject<TH1D> h_mf_mod_1_sig("h_mf_mod_1_sig","Nstubs per mod muon filter tar 1",5,0,5);
   ROOT::TThreadedObject<TH1D> h_mu_mf_1_sig("h_mu_mf_1_sig","Nstubs muon track in MF tar 1",30,0,30);
   ROOT::TThreadedObject<TH1D> h_e_mf_1_sig("h_e_mf_1_sig","Nstubs electron track in MF tar 1",30,0,30);
   ROOT::TThreadedObject<TH1D> h_phi_e_1_sig("h_phi_e_1_sig","phi electron tar 1",180,-180,180);
   ROOT::TThreadedObject<TH1D> h_phi_mu_1_sig("h_phi_mu_1_sig","phi muon tar 1",180,-180,180);

 auto myFunction = [&](TTreeReader &myReader) {


     TTreeReaderValue<MUonEEventHeader> RVeventHeader(myReader, "EventHeader");
     TTreeReaderValue<std::vector<MUonERecoOutputTrackAnalysis>> RVtracks(myReader, "ReconstructedTracks");
     TTreeReaderValue<std::vector<MUonERecoOutputVertexAnalysis>> vrtx(myReader, "BestVertex");
     TTreeReaderValue<std::vector<MUonERecoOutputHitAnalysis>> RVstubs(myReader, "ReconstructedHits");
//     TTreeReaderValue<std::vector<MUonETrackerStub>> tr_stubs(myReader,"TrackerStubs");

//     TTreeReaderValue<vector<short>> mf_modID(myReader, "MuonFilterModuleId");
//     TTreeReaderValue<vector<float>> mf_localX(myReader, "MuonFilterLocalX");

//     TTreeReaderValue<std::vector<bitset<8>>> RVtrig_bits(myReader, "TriggerBits");

//     TTreeReaderValue<Bool_t> trig_smi0(myReader, "TriggerSingleMuonInteraction0");
//     TTreeReaderValue<Bool_t> trig_smi1(myReader, "TriggerSingleMuonInteraction1");

count_tr_1->Fill(-99);
count_tr_2->Fill(-99);
count_tr_1_sig->Fill(-99);
count_tr_2_sig->Fill(-99);
count_mu_is_mu_0->Fill(-99);
count_e_is_mu_0->Fill(-99);
count_mu_is_mu_1->Fill(-99);
count_e_is_mu_1->Fill(-99);
count_mu_is_mu_0_sig->Fill(-99);
count_e_is_mu_0_sig->Fill(-99);
count_mu_is_mu_1_sig->Fill(-99);
count_e_is_mu_1_sig->Fill(-99);

     while (myReader.Next()) {

double chi=0.;
if(vrtx->size()!=0)chi=vrtx.at(0)->chi2();

double posxIN=99.;//pos_on_track(x0_in,th_inx,z_fix);
double posyIN=99.;//pos_on_track(y0_in,th_iny,z_fix);


 MUonERecoOutputTrackAnalysis mu_in_v = vrtx.at(0)->incomingMuon();
 MUonERecoOutputTrackAnalysis mu_out_v = vrtx.at(0)->outgoingMuon();
 MUonERecoOutputTrackAnalysis e_out_v = vrtx.at(0)->outgoingElectron();


std::array<int,3> nTr{0};
std::vector<MUonERecoOutputTrackAnalysis> track_0;
std::vector<MUonERecoOutputTrackAnalysis> track_1;
std::vector<MUonERecoOutputTrackAnalysis> track_vrtx1;
std::vector<MUonERecoOutputTrackAnalysis> track_vrtx2;

         auto& tracks = *RVtracks;
         for (auto&& track : tracks) {

       if(track.sector()==0){nTr.at(0)++; track_0.push_back(track);}
        else if(track.sector()==1){nTr.at(1)++; track_1.push_back(track);
				   if(track.index()==mu_out_v.index() or track.index()==e_out_v.index())track_vrtx1.push_back(track);
}
        else if(track.sector()==2){nTr.at(2)++;
				   if(track.index()==mu_out_v.index() or track.index()==e_out_v.index())track_vrtx2.push_back(track);}
        }


//if(chi==0)std::cout << " no vrtx " << myReader.GetCurrentEntry() << " sec0: " << nTr.at(0) << " sec1: " << nTr.at(1) << " sec2: " << nTr.at(2) << std::endl;

MUonERecoOutputHitAnalysis muin_4;
MUonERecoOutputHitAnalysis muin_5;

auto& hits = *RVstubs;

auto trig_bits = *RVeventHeader;
if(trig_bits.onlineTriggerSingleMuonInteraction0()==1 and nTr.at(0)==1){
 std::vector<Short_t> ids_=track_0.at(0).hitIds();
 for (auto&& stub : hits) {
        for(auto&& h : ids_){
	 if(h==stub.index()){
          if(stub.moduleID()==4){muin_4 = stub;}
          else if(stub.moduleID()==5){muin_5 = stub;}
	  }
         }
	}
}
else if(trig_bits.onlineTriggerSingleMuonInteraction1()==1 and nTr.at(1)==1){
 std::vector<Short_t> ids_=track_1.at(0).hitIds();
 for (auto&& stub : hits) {
        for(auto&& h : ids_){
         if(h==stub.index()){
          if(stub.moduleID()==4){muin_4 = stub;}
          else if(stub.moduleID()==5){muin_5 = stub;}
          }
         }
        }
}


int stub0=0;int stub1=0;int stub2=0; int mf_nStubs=0;

         for (auto&& hit : hits) {
//		std::cout << "hit muid size " << hit.muonIds().size() << std::endl;
        if(hit.stationID()==0){stub0++;}
        else if(hit.stationID()==1){stub1++;}
        else if(hit.stationID()==2){stub2++;}
        else if(hit.stationID()==3){mf_nStubs++;}
        }

//std::cout<<"ce sono" << std::endl;

if(trig_bits.onlineTriggerSingleMuonInteraction0()==1 and nTr.at(0)==1 and abs(muin_4.position())<=1.4985 and abs(muin_5.position())<=1.4985 and stub0<10){
/*
std::cout << "-----------" << std::endl;
if(track_vrtx1.at(0).index()==mu_out_v.index()){std::cout << "muon is muon? " << track_vrtx1.at(0).isMuon() << std::endl;}
else if(track_vrtx1.at(1).index()==mu_out_v.index()){std::cout << "muon is muon? " << track_vrtx1.at(1).isMuon() << std::endl;}
if(track_vrtx1.at(0).index()==e_out_v.index()){std::cout << "electron is muon? " << track_vrtx1.at(0).isMuon() << std::endl;}
else if(track_vrtx1.at(1).index()==e_out_v.index()){std::cout << "electron is muon? " << track_vrtx1.at(1).isMuon() << std::endl;}

*/


	if(chi!=0 and vrtx.at(0)->stationIndex()==1){
	double the_rec=vrtx.at(0)->electronTheta();
	double thmu_rec=vrtx.at(0)->muonTheta();
        double acoplanarity_v=vrtx.at(0)->modifiedAcoplanarity();
        double zpos =vrtx.at(0)->zPositionFit();

        double IP_out=-99.;


		double IPx=pos_on_track(track_vrtx1.at(1).x0(),track_vrtx1.at(1).xSlope(),(667.3-2.7)) - pos_on_track(track_vrtx1.at(0).x0(),track_vrtx1.at(0).xSlope(),(667.3-2.7));
                double IPy=pos_on_track(track_vrtx1.at(1).y0(),track_vrtx1.at(1).ySlope(),(667.3-2.7)) - pos_on_track(track_vrtx1.at(0).y0(),track_vrtx1.at(0).ySlope(),(667.3-2.7));

                IP_out=sqrt(IPx*IPx + IPy*IPy);

		TVector3 p_e(vrtx.at(0)->outgoingElectron().xSlope(),vrtx.at(0)->outgoingElectron().ySlope(),1.0); p_e.Unit();
		TVector3 p_mu(vrtx.at(0)->outgoingMuon().xSlope(),vrtx.at(0)->outgoingMuon().ySlope(),1.0); p_mu.Unit();

                int n_mu_mf=0; int n_e_mf=0;
                int muonid_mu=0;
                int muonid_e=0;
		int c_mu=0;
		int c_e=0;
		double phi_mu=1000.;
		double phi_e=1000.;


                for(auto&& track : track_vrtx1){
                        if(track.index()==mu_out_v.index()){
					TVector3 p(track.xSlope(),track.ySlope(),1.0); p.Unit(); phi_mu=p.Phi()*(180/M_PI);
					if(track.isMuon()){c_mu=1; muonid_mu=track.muonId();}
					}
                        else if(track.index()==e_out_v.index()){
					TVector3 p(track.xSlope(),track.ySlope(),1.0); p.Unit(); phi_e=p.Phi()*(180/M_PI);
					if(track.isMuon()){c_e=1; muonid_e=track.muonId();}
					}
                 }


		for (auto&& hit : hits) { if(hit.stationID()==3){
					if(std::find_if(hit.muonIds().begin(),hit.muonIds().end(),[&](int i){ return i==muonid_mu; })!=hit.muonIds().end())n_mu_mf++;
					else if(std::find_if(hit.muonIds().begin(),hit.muonIds().end(),[&](int i){ return i==muonid_e; })!=hit.muonIds().end())n_e_mf++;
								}
					}


//	 int mf_nStubs;
//	 auto& mf = *mf_modID;
//         for (auto&& mod : mf) {mf_nStubs++;}

		if(IP_out>0.2 && p_mu.Angle(p_e)<0.002){
		 count_tr_1->Fill(1.);

//std::cout << " TAR0 bkg " << myReader.GetCurrentEntry() << " sec0: " << nTr.at(0) << " sec1: " << nTr.at(1) << " sec2: " << nTr.at(2) << std::endl;

//cout << "n_mu_mf " << n_mu_mf << " muonid_mu " << muonid_mu << std::endl;
//cout << "n_e_mf " << n_e_mf << std::endl;

 		 h_phi_e_0->Fill(phi_e);
		 h_phi_mu_0->Fill(phi_mu);
		 h_nhits_pre_0->Fill(stub0);
		 h_nhits_post_0->Fill(stub1);
		 h_ntracks_post_0->Fill(nTr.at(1));
		 h_chi2_vrtx_0->Fill(chi);
		 h_Dz_vrtx_0->Fill(zpos-(667.3-2.7));
		 h_aco_0->Fill(acoplanarity_v);
		 h_th_min_0->Fill(thmu_rec);
		 h_th_max_0->Fill(the_rec);
                h_IP_0->Fill(IP_out);
		h_IPop_0->Fill(IP_out,p_mu.Angle(p_e));
		 h_mf_0->Fill(mf_nStubs);

		if(c_mu==1)count_mu_is_mu_0->Fill(1.);
		if(c_e==1)count_e_is_mu_0->Fill(1.);

                for (auto&& hit : hits) {h_mf_mod_0->Fill(hit.moduleID());}

		h_mu_mf_0->Fill(n_mu_mf);
		h_e_mf_0->Fill(n_e_mf);

		}
                else if(IP_out<0.2 && p_mu.Angle(p_e)>0.002){
                 count_tr_1_sig->Fill(1.);

//std::cout << " TAR0 sig " << myReader.GetCurrentEntry() << " sec0: " << nTr.at(0) << " sec1: " << nTr.at(1) << " sec2: " << nTr.at(2) << std::endl;
                 h_phi_e_0_sig->Fill(phi_e);
                 h_phi_mu_0_sig->Fill(phi_mu);
                 h_nhits_pre_0_sig->Fill(stub0);
                 h_nhits_post_0_sig->Fill(stub1);
                 h_ntracks_post_0_sig->Fill(nTr.at(1));
                 h_chi2_vrtx_0_sig->Fill(chi);
                 h_Dz_vrtx_0_sig->Fill(zpos-(667.3-2.7));
                 h_aco_0_sig->Fill(acoplanarity_v);
                 h_th_min_0_sig->Fill(thmu_rec);
                 h_th_max_0_sig->Fill(the_rec);
                 h_mf_0_sig->Fill(mf_nStubs);

                if(c_mu==1)count_mu_is_mu_0_sig->Fill(1.);
                if(c_e==1)count_e_is_mu_0_sig->Fill(1.);

                for (auto&& hit : hits) {h_mf_mod_0_sig->Fill(hit.moduleID());}

                h_mu_mf_0_sig->Fill(n_mu_mf);
                h_e_mf_0_sig->Fill(n_e_mf);

              }
	}
}
else if(trig_bits.onlineTriggerSingleMuonInteraction1()==1 and nTr.at(1)==1 and abs(muin_4.position())<=1.4985 and abs(muin_5.position())<=1.4985 and stub1<10){

//std::cout<<"ce sono 2" << std::endl;

	if(chi!=0 and vrtx.at(0)->stationIndex()==2){
	double the_rec=vrtx.at(0)->electronTheta();
	double thmu_rec=vrtx.at(0)->muonTheta();
	double acoplanarity_v=vrtx.at(0)->modifiedAcoplanarity();
	double zpos =vrtx.at(0)->zPositionFit();

        double IP_out=-99.;

		double IPx=pos_on_track(track_vrtx2.at(1).x0(),track_vrtx2.at(1).xSlope(),(784.6-3.9)) - pos_on_track(track_vrtx2.at(0).x0(),track_vrtx2.at(0).xSlope(),(784.6-3.9));
                double IPy=pos_on_track(track_vrtx2.at(1).y0(),track_vrtx2.at(1).ySlope(),(784.6-3.9)) - pos_on_track(track_vrtx2.at(0).y0(),track_vrtx2.at(0).ySlope(),(784.6-3.9));

                IP_out=sqrt(IPx*IPx + IPy*IPy);

                TVector3 p_e(vrtx.at(0)->outgoingElectron().xSlope(),vrtx.at(0)->outgoingElectron().ySlope(),1.0); p_e.Unit();
                TVector3 p_mu(vrtx.at(0)->outgoingMuon().xSlope(),vrtx.at(0)->outgoingMuon().ySlope(),1.0); p_mu.Unit();

                int n_mu_mf=0; int n_e_mf=0;
                int muonid_mu=0;
                int muonid_e=0;
                int c_mu=0;
                int c_e=0;

//std::cout<<"ce sono 3" << std::endl;

                for(auto&& track : track_vrtx2){
                        if(track.index()==mu_out_v.index() and track.isMuon()){c_mu=1; muonid_mu=track.muonId();}
                        else if(track.index()==e_out_v.index() and track.isMuon()){c_e=1; muonid_e=track.muonId();}
                 }

                for (auto&& hit : hits) { if(hit.stationID()==3){
                                                        if(std::find_if(hit.muonIds().begin(),hit.muonIds().end(),[&](int i){ return i==muonid_mu; })!=hit.muonIds().end())n_mu_mf++;
                                                        else if(std::find_if(hit.muonIds().begin(),hit.muonIds().end(),[&](int i){ return i==muonid_e; })!=hit.muonIds().end())n_e_mf++;
                                                                }
                                        }

//	 int mf_nStubs;
//         auto mf = *mf_modID;
//         for (auto&& mod : mf) {mf_nStubs++;}

                if(IP_out>0.2 && p_mu.Angle(p_e)<0.002){
//std::cout << " TAR1 bkg " << myReader.GetCurrentEntry()<< " sec0: " << nTr.at(0) << " sec1: " << nTr.at(1) << " sec2: " << nTr.at(2) << std::endl;
		 count_tr_2->Fill(1.);
		 h_nhits_pre_0->Fill(stub0);
		 h_nhits_pre_1->Fill(stub1);
		 h_nhits_post_1->Fill(stub2);
		 h_ntracks_post_1->Fill(nTr.at(2));
		 h_chi2_vrtx_1->Fill(chi);
		 h_Dz_vrtx_1->Fill(zpos-(784.6-3.9));
		 h_aco_1->Fill(acoplanarity_v);
		 h_th_min_1->Fill(thmu_rec);
		 h_th_max_1->Fill(the_rec);
                h_IP_1->Fill(IP_out);
                h_IPop_1->Fill(IP_out,p_mu.Angle(p_e));
                 h_mf_1->Fill(mf_nStubs);

                h_mu_mf_1->Fill(n_mu_mf);
                h_e_mf_1->Fill(n_e_mf);

		for (auto&& hit : hits) {h_mf_mod_1->Fill(hit.moduleID());}

                if(c_mu==1)count_mu_is_mu_1->Fill(1.);
                if(c_e==1)count_e_is_mu_1->Fill(1.);

		}
                else if(IP_out<0.2 && p_mu.Angle(p_e)>0.002){

//std::cout << " TAR1 sig " << myReader.GetCurrentEntry()<< " sec0: " << nTr.at(0) << " sec1: " << nTr.at(1) << " sec2: " << nTr.at(2) << std::endl;
                 count_tr_2_sig->Fill(1.);
                 h_nhits_pre_0_sig->Fill(stub0);
                 h_nhits_pre_1_sig->Fill(stub1);
                 h_nhits_post_1_sig->Fill(stub2);
                 h_ntracks_post_1_sig->Fill(nTr.at(2));
                 h_chi2_vrtx_1_sig->Fill(chi);
                 h_Dz_vrtx_1_sig->Fill(zpos-(784.6-3.9));
                 h_aco_1_sig->Fill(acoplanarity_v);
                 h_th_min_1_sig->Fill(thmu_rec);
                 h_th_max_1_sig->Fill(the_rec);
                 h_mf_1_sig->Fill(mf_nStubs);

                h_mu_mf_1_sig->Fill(n_mu_mf);
                h_e_mf_1_sig->Fill(n_e_mf);

                for (auto&& hit : hits) {h_mf_mod_1_sig->Fill(hit.moduleID());}

                if(c_mu==1)count_mu_is_mu_1_sig->Fill(1.);
                if(c_e==1)count_e_is_mu_1_sig->Fill(1.);

                }
	}
}





 } //end of general while
}; //end of ySlopefunction

  tp1.Process(myFunction);


  auto count_tr_1M=count_tr_1.Merge();
  auto count_tr_2M=count_tr_2.Merge();
  auto count_tr_1_sigM=count_tr_1_sig.Merge();
  auto count_tr_2_sigM=count_tr_2_sig.Merge();

  auto h_nhits_pre_0M=h_nhits_pre_0.Merge();
  auto h_nhits_post_0M=h_nhits_post_0.Merge();
  auto h_ntracks_post_0M=h_ntracks_post_0.Merge();
  auto h_chi2_vrtx_0M=h_chi2_vrtx_0.Merge();
  auto h_Dz_vrtx_0M=h_Dz_vrtx_0.Merge();
  auto h_aco_0M=h_aco_0.Merge();
  auto h_th_min_0M=h_th_min_0.Merge();
  auto h_th_max_0M=h_th_max_0.Merge();
  auto h_IP_0M=h_IP_0.Merge();
  auto h_IPop_0M=h_IPop_0.Merge();
  auto h_mf_0M=h_mf_0.Merge();
  auto h_mf_mod_0M=h_mf_mod_0.Merge();
  auto h_mu_mf_0M=h_mu_mf_0.Merge();
  auto h_e_mf_0M=h_e_mf_0.Merge();
  auto count_mu_is_mu_0M=count_mu_is_mu_0.Merge();
  auto count_e_is_mu_0M=count_e_is_mu_0.Merge();
  auto h_phi_e_0M=h_phi_e_0.Merge();
  auto h_phi_mu_0M=h_phi_mu_0.Merge();

  auto h_nhits_pre_1M=h_nhits_pre_1.Merge();
  auto h_nhits_post_1M=h_nhits_post_1.Merge();
  auto h_ntracks_post_1M=h_ntracks_post_1.Merge();
  auto h_chi2_vrtx_1M=h_chi2_vrtx_1.Merge();
  auto h_Dz_vrtx_1M=h_Dz_vrtx_1.Merge();
  auto h_aco_1M=h_aco_1.Merge();
  auto h_th_min_1M=h_th_min_1.Merge();
  auto h_th_max_1M=h_th_max_1.Merge();
  auto h_IP_1M=h_IP_1.Merge();
  auto h_IPop_1M=h_IPop_1.Merge();
  auto h_mf_1M=h_mf_1.Merge();
  auto h_mf_mod_1M=h_mf_mod_1.Merge();
  auto h_mu_mf_1M=h_mu_mf_1.Merge();
  auto h_e_mf_1M=h_e_mf_1.Merge();
  auto count_mu_is_mu_1M=count_mu_is_mu_1.Merge();
  auto count_e_is_mu_1M=count_e_is_mu_1.Merge();

  auto h_nhits_pre_0_sigM=h_nhits_pre_0_sig.Merge();
  auto h_nhits_post_0_sigM=h_nhits_post_0_sig.Merge();
  auto h_ntracks_post_0_sigM=h_ntracks_post_0_sig.Merge();
  auto h_chi2_vrtx_0_sigM=h_chi2_vrtx_0_sig.Merge();
  auto h_Dz_vrtx_0_sigM=h_Dz_vrtx_0_sig.Merge();
  auto h_aco_0_sigM=h_aco_0_sig.Merge();
  auto h_th_min_0_sigM=h_th_min_0_sig.Merge();
  auto h_th_max_0_sigM=h_th_max_0_sig.Merge();
  auto h_mf_0_sigM=h_mf_0_sig.Merge();
  auto h_mf_mod_0_sigM=h_mf_mod_0_sig.Merge();
  auto h_mu_mf_0_sigM=h_mu_mf_0_sig.Merge();
  auto h_e_mf_0_sigM=h_e_mf_0_sig.Merge();
  auto count_mu_is_mu_0_sigM=count_mu_is_mu_0_sig.Merge();
  auto count_e_is_mu_0_sigM=count_e_is_mu_0_sig.Merge();
  auto h_phi_e_0_sigM=h_phi_e_0_sig.Merge();
  auto h_phi_mu_0_sigM=h_phi_mu_0_sig.Merge();

  auto h_nhits_pre_1_sigM=h_nhits_pre_1_sig.Merge();
  auto h_nhits_post_1_sigM=h_nhits_post_1_sig.Merge();
  auto h_ntracks_post_1_sigM=h_ntracks_post_1_sig.Merge();
  auto h_chi2_vrtx_1_sigM=h_chi2_vrtx_1_sig.Merge();
  auto h_Dz_vrtx_1_sigM=h_Dz_vrtx_1_sig.Merge();
  auto h_aco_1_sigM=h_aco_1_sig.Merge();
  auto h_th_min_1_sigM=h_th_min_1_sig.Merge();
  auto h_th_max_1_sigM=h_th_max_1_sig.Merge();
  auto h_mf_1_sigM=h_mf_1_sig.Merge();
  auto h_mf_mod_1_sigM=h_mf_mod_1_sig.Merge();
  auto h_mu_mf_1_sigM=h_mu_mf_1_sig.Merge();
  auto h_e_mf_1_sigM=h_e_mf_1_sig.Merge();
  auto count_mu_is_mu_1_sigM=count_mu_is_mu_1_sig.Merge();
  auto count_e_is_mu_1_sigM=count_e_is_mu_1_sig.Merge();

//cout << "good_vrtx_1->Integral() " << good_vrtx_1->Integral() << endl;
//cout << "fiducial_1->Integral() " << fiducial_1->Integral() << endl;

    std::ofstream out(Form("txt/run%i/possible_bkg_report_%s_%s_nhits%i.txt",run,type.c_str(),filen.c_str(),nhits)); // apre (o crea) il file in scrittura
    if (!out) {
        std::cerr << "Errore nell'aprire il file!" << std::endl;
        return 1;
    }
    out << "Entries " << cbmsim->GetEntries() << std::endl;
    out << "Events with IP>0.2cm and op_angle<2mrad tar0 " << count_tr_1M->Integral()
	<< " with muons with isMuon " << count_mu_is_mu_0M->Integral() << " -> " << count_mu_is_mu_0M->Integral()/count_tr_1M->Integral()
	<< " with electron with isMuon " << count_e_is_mu_0M->Integral() << " -> " << count_e_is_mu_0M->Integral()/count_tr_1M->Integral() << std::endl;

    out << "Events with IP>0.2cm and op_angle<2mrad tar1 " << count_tr_2M->Integral()
	<< " with muons with isMuon " << count_mu_is_mu_1M->Integral() << " -> " << count_mu_is_mu_1M->Integral()/count_tr_2M->Integral()
	<< " with electron with isMuon " << count_e_is_mu_1M->Integral() << " -> " << count_e_is_mu_1M->Integral()/count_tr_2M->Integral() << std::endl;

    out << "Events with IP<0.2cm and op_angle>2mrad tar0 " << count_tr_1_sigM->Integral()
        << " with muons with isMuon " << count_mu_is_mu_0_sigM->Integral() << " -> " << count_mu_is_mu_0_sigM->Integral()/count_tr_1_sigM->Integral()
        << " with electron with isMuon " << count_e_is_mu_0_sigM->Integral() << " -> " << count_e_is_mu_0_sigM->Integral()/count_tr_1_sigM->Integral() << std::endl;

    out << "Events with IP<0.2cm and op_angle>2mrad tar1 " << count_tr_2_sigM->Integral()
        << " with muons with isMuon " << count_mu_is_mu_1_sigM->Integral() << " -> " << count_mu_is_mu_1_sigM->Integral()/count_tr_2_sigM->Integral()
        << " with electron with isMuon " << count_e_is_mu_1_sigM->Integral() << " -> " << count_e_is_mu_1_sigM->Integral()/count_tr_2_sigM->Integral() << std::endl;

    out << "Ratio IP>0.2cm and op_angle<2mrad tar0 over entries: " << count_tr_1M->Integral()/cbmsim->GetEntries() << std::endl;
    out << "Ratio IP>0.2cm and op_angle<2mrad tar1 over entries: " << count_tr_2M->Integral()/cbmsim->GetEntries() << std::endl;
    out << "Ratio IP<0.2cm and op_angle>2mrad tar0 over entries: " << count_tr_1_sigM->Integral()/cbmsim->GetEntries() << std::endl;
    out << "Ratio IP<0.2cm and op_angle>2mrad tar1 over entries: " << count_tr_2_sigM->Integral()/cbmsim->GetEntries() << std::endl;

    out.close();


    TFile out_root(Form("txt/run%i/possible_bkg_root_%s_%s_nhits%i.root",run,type.c_str(),filen.c_str(),nhits),"recreate");
      if(type=="single_muon_interaction_0"){
	h_nhits_pre_0M->Write();
	h_nhits_post_0M->Write();
	h_ntracks_post_0M->Write();
	h_chi2_vrtx_0M->Write();
	h_Dz_vrtx_0M->Write();
	h_aco_0M->Write();
	h_th_min_0M->Write();
	h_th_max_0M->Write();
	h_IP_0M->Write();
	h_IPop_0M->Write();
	h_mf_0M->Write();
	h_mf_mod_0M->Write();
        h_mu_mf_0M->Write();
        h_e_mf_0M->Write();
        h_phi_mu_0M->Write();
        h_phi_e_0M->Write();

	h_nhits_pre_0_sigM->Write();
	h_nhits_post_0_sigM->Write();
	h_ntracks_post_0_sigM->Write();
	h_chi2_vrtx_0_sigM->Write();
	h_Dz_vrtx_0_sigM->Write();
	h_aco_0_sigM->Write();
	h_th_min_0_sigM->Write();
	h_th_max_0_sigM->Write();
	h_mf_0_sigM->Write();
        h_mf_mod_0_sigM->Write();
	h_mu_mf_0_sigM->Write();
	h_e_mf_0_sigM->Write();
	h_phi_mu_0_sigM->Write();
	h_phi_e_0_sigM->Write();
      }
      else if(type=="single_muon_interaction_1"){

	h_nhits_pre_0M->Write();
        h_nhits_pre_1M->Write();
        h_nhits_post_1M->Write();
        h_ntracks_post_1M->Write();
        h_chi2_vrtx_1M->Write();
        h_Dz_vrtx_1M->Write();
        h_aco_1M->Write();
        h_th_min_1M->Write();
        h_th_max_1M->Write();
        h_IP_1M->Write();
        h_IPop_1M->Write();
        h_mf_1M->Write();
        h_mf_mod_1M->Write();
        h_mu_mf_1M->Write();
        h_e_mf_1M->Write();

        h_nhits_pre_0_sigM->Write();
        h_nhits_pre_1_sigM->Write();
        h_nhits_post_1_sigM->Write();
        h_ntracks_post_1_sigM->Write();
        h_chi2_vrtx_1_sigM->Write();
        h_Dz_vrtx_1_sigM->Write();
        h_aco_1_sigM->Write();
        h_th_min_1_sigM->Write();
        h_th_max_1_sigM->Write();
	h_mf_1_sigM->Write();
        h_mf_mod_1_sigM->Write();
        h_mu_mf_1_sigM->Write();
        h_e_mf_1_sigM->Write();
      }
out_root.Write();
out_root.Close();


/*TCanvas c("c","c",1000,500);
c.Divide(2,1);
c.cd(1);
h_el_1M->Draw("COLZ");
c.cd(2);
h_el_2M->Draw("COLZ");
c.SaveAs("2D.pdf");
*/
return 0;
}



/*
	double IP_mu=-99.;
	double IP_e=-99.;

	if(track_vrtx1.at(0).index()==mu_out_v.index() and track_vrtx1.at(1).index()==e_out_v.index()){
		double IPx_mu=pos_on_track(track_0.at(0).x0(),track_0.at(0).xSlope(),(667.3-2.7)) - pos_on_track(track_vrtx1.at(0).x0(),track_vrtx1.at(0).xSlope(),(667.3-2.7));
		double IPx_e=pos_on_track(track_0.at(0).x0(),track_0.at(0).xSlope(),(667.3-2.7)) - pos_on_track(track_vrtx1.at(1).x0(),track_vrtx1.at(1).xSlope(),(667.3-2.7));
                double IPy_mu=pos_on_track(track_0.at(0).y0(),track_0.at(0).ySlope(),(667.3-2.7)) - pos_on_track(track_vrtx1.at(0).y0(),track_vrtx1.at(0).ySlope(),(667.3-2.7));
                double IPy_e=pos_on_track(track_0.at(0).y0(),track_0.at(0).ySlope(),(667.3-2.7)) - pos_on_track(track_vrtx1.at(1).y0(),track_vrtx1.at(1).ySlope(),(667.3-2.7));

		IP_mu=sqrt(IPx_mu*IPx_mu + IPy_mu*IPy_mu);
		IP_e=sqrt(IPx_e*IPx_e + IPy_e*IPy_e);
                h_IP_mu0->Fill(IP_mu);
                h_IP_0->Fill(IP_e);
	}
	else if(track_vrtx1.at(1).index()==mu_out_v.index() and track_vrtx1.at(0).index()==e_out_v.index()){

                double IPx_e=pos_on_track(track_0.at(0).x0(),track_0.at(0).xSlope(),(667.3-2.7)) - pos_on_track(track_vrtx1.at(0).x0(),track_vrtx1.at(0).xSlope(),(667.3-2.7));
                double IPx_mu=pos_on_track(track_0.at(0).x0(),track_0.at(0).xSlope(),(667.3-2.7)) - pos_on_track(track_vrtx1.at(1).x0(),track_vrtx1.at(1).xSlope(),(667.3-2.7));
                double IPy_e=pos_on_track(track_0.at(0).y0(),track_0.at(0).ySlope(),(667.3-2.7)) - pos_on_track(track_vrtx1.at(0).y0(),track_vrtx1.at(0).ySlope(),(667.3-2.7));
                double IPy_mu=pos_on_track(track_0.at(0).y0(),track_0.at(0).ySlope(),(667.3-2.7)) - pos_on_track(track_vrtx1.at(1).y0(),track_vrtx1.at(1).ySlope(),(667.3-2.7));

                IP_mu=sqrt(IPx_mu*IPx_mu + IPy_mu*IPy_mu);
                IP_e=sqrt(IPx_e*IPx_e + IPy_e*IPy_e);
                h_IP_mu0->Fill(IP_mu);
                h_IP_0->Fill(IP_e);
	}
*/



/*
        if(track_vrtx2.at(0).index()==mu_out_v.index() and track_vrtx2.at(1).index()==e_out_v.index() and track_1.at(0).index()==mu_in_v.index()){

		double IPx_mu=pos_on_track(track_1.at(0).x0(),track_1.at(0).xSlope(),(784.6-3.9)) - pos_on_track(track_vrtx2.at(0).x0(),track_vrtx2.at(0).xSlope(),(784.6-3.9));
                double IPx_e=pos_on_track(track_1.at(0).x0(),track_1.at(0).xSlope(),(784.6-3.9)) - pos_on_track(track_vrtx2.at(1).x0(),track_vrtx2.at(1).xSlope(),(784.6-3.9));
                double IPy_mu=pos_on_track(track_1.at(0).y0(),track_1.at(0).ySlope(),(784.6-3.9)) - pos_on_track(track_vrtx2.at(0).y0(),track_vrtx2.at(0).ySlope(),(784.6-3.9));
                double IPy_e=pos_on_track(track_1.at(0).y0(),track_1.at(0).ySlope(),(784.6-3.9)) - pos_on_track(track_vrtx2.at(1).y0(),track_vrtx2.at(1).ySlope(),(784.6-3.9));

                IP_mu=sqrt(IPx_mu*IPx_mu + IPy_mu*IPy_mu);
                IP_e=sqrt(IPx_e*IPx_e + IPy_e*IPy_e);
                h_IP_mu1->Fill(IP_mu);
                h_IP_1->Fill(IP_e);

        }
	else if(track_vrtx2.at(1).index()==mu_out_v.index() and track_vrtx2.at(0).index()==e_out_v.index() and track_1.at(0).index()==mu_in_v.index()){

                double IPx_e=pos_on_track(track_1.at(0).x0(),track_1.at(0).xSlope(),(784.6-3.9)) - pos_on_track(track_vrtx2.at(0).x0(),track_vrtx2.at(0).xSlope(),(784.6-3.9));
                double IPx_mu=pos_on_track(track_1.at(0).x0(),track_1.at(0).xSlope(),(784.6-3.9)) - pos_on_track(track_vrtx2.at(1).x0(),track_vrtx2.at(1).xSlope(),(784.6-3.9));
                double IPy_e=pos_on_track(track_1.at(0).y0(),track_1.at(0).ySlope(),(784.6-3.9)) - pos_on_track(track_vrtx2.at(0).y0(),track_vrtx2.at(0).ySlope(),(784.6-3.9));
                double IPy_mu=pos_on_track(track_1.at(0).y0(),track_1.at(0).ySlope(),(784.6-3.9)) - pos_on_track(track_vrtx2.at(1).y0(),track_vrtx2.at(1).ySlope(),(784.6-3.9));

                IP_mu=sqrt(IPx_mu*IPx_mu + IPy_mu*IPy_mu);
                IP_e=sqrt(IPx_e*IPx_e + IPy_e*IPy_e);
                h_IP_mu1->Fill(IP_mu);
                h_IP_1->Fill(IP_e);
        }

*/

