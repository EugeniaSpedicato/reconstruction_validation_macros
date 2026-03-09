#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TSystemDirectory.h"
#include <TStyle.h>

using namespace std;

auto pos_on_track = [](double q, double m, double z){return (q + m*z);};


int MC_nocut5mrad_report_count_snakemake(string run, string type, int nhits){

  int nthreads = 6;

  ROOT::EnableImplicitMT(nthreads);


TChain * cbmsim = new TChain("cbmsim");
TChain * cbmsim_g = new TChain("cbmsim");


//if(type=="single_muon_interaction_0"){cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/MC_validation_tar0.root");}
//else if(type=="single_muon_interaction_1"){cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/MC_validation_tar1.root");}

if(type=="single_muon_interaction_0"){
cbmsim->Add("/eos/user/e/ehess/production_Mesmer/FairMUonE_WiP_v1.1.x/Condor_files/MCsignalJob_target0_1.root");
cbmsim->Add("/eos/user/e/ehess/production_Mesmer/FairMUonE_WiP_v1.1.x/Condor_files/MCsignalJob_target0_3.root");
cbmsim->Add("/eos/user/e/ehess/production_Mesmer/FairMUonE_WiP_v1.1.x/Condor_files/MCsignalJob_target0_6.root");
}
else if(type=="single_muon_interaction_1"){
cbmsim->Add("/eos/user/e/ehess/production_Mesmer/FairMUonE_WiP_v1.1.x/Condor_files/MCsignalJob_target1_2.root");
cbmsim->Add("/eos/user/e/ehess/production_Mesmer/FairMUonE_WiP_v1.1.x/Condor_files/MCsignalJob_target1_3.root");
cbmsim->Add("/eos/user/e/ehess/production_Mesmer/FairMUonE_WiP_v1.1.x/Condor_files/MCsignalJob_target1_4.root");
}

/*
cbmsim->Add(Form("/eos/experiment/mu-e/reco/2025/%i/%s/%s.root",run.c_str(),type.c_str(),filen.c_str()));
cout << Form("/eos/experiment/mu-e/reco/2025/%i/%s/%s.root",run.c_str(),type.c_str(),filen.c_str())<< endl;;
*/

cout << "cbmsim->GetEntries() " << cbmsim->GetEntries() <<endl;

ROOT::TTreeProcessorMT tp1(*cbmsim,nthreads);

      MUonERecoOutputAnalysis *ReconstructionOutput = 0;

   ROOT::TThreadedObject<TH1D> h_nvrtx_1("h_nvrtx_s2","nvrtxs station01",2,0,2);
   ROOT::TThreadedObject<TH1D> h_nvrtx_2("h_nvrtx_s2","nvrtxs station12",2,0,2);

   ROOT::TThreadedObject<TH1D> count_tr_1("count_chi1","count_chi1",2,0,2);
   ROOT::TThreadedObject<TH1D> count_tr_2("count_chi2","count_chi2",2,0,2);

   ROOT::TThreadedObject<TH1D> good_vrtx_1("good_vrtx_1","good_vrtx_1",2,0,2);
   ROOT::TThreadedObject<TH1D> good_vrtx_2("good_vrtx_2","good_vrtx_2",2,0,2);

   ROOT::TThreadedObject<TH1D> fiducial_1("fiducial_1","fiducial",2,0,2);
   ROOT::TThreadedObject<TH1D> fiducial_2("fiducial_2","fiducial",2,0,2);

   ROOT::TThreadedObject<TH1D> pre_elastic_1("pre_elastic_1","pre_elastic_1",2,0,2);
   ROOT::TThreadedObject<TH1D> pre_elastic_2("pre_elastic_2","pre_elastic_2",2,0,2);

   ROOT::TThreadedObject<TH1D> elastic_1("elastic_1","elastic_1",2,0,2);
   ROOT::TThreadedObject<TH1D> elastic_2("elastic_2","elastic_2",2,0,2);

   ROOT::TThreadedObject<TH2D> h_el_1("h_el_1","h_el_1",320, 0.,0.032,100,0.,0.005);
   ROOT::TThreadedObject<TH2D> h_el_2("h_el_2","h_el_2",320, 0.,0.032,100,0.,0.005);

   ROOT::TThreadedObject<TH1D> h_x_0("h_x_0","X1 fiducial tar 0",400,-2.,2.);
   ROOT::TThreadedObject<TH1D> h_y_0("h_y_0","Y1 fiducial tar 0",400,-2.,2.);
   ROOT::TThreadedObject<TH1D> h_thx_0("h_thx_0","theta X incoming fiducial tar 0",200,-0.002,0.002);
   ROOT::TThreadedObject<TH1D> h_thy_0("h_thy_0","theta Y incoming fiducial tar 0",200,-0.002,0.002);

   ROOT::TThreadedObject<TH1D> h_IP_0("h_IP_0","Impact point outgoing mu - outgoing el tar 0",500,-5.,5.);
   ROOT::TThreadedObject<TH2D> h_IPop_0("h_IPop_0","Impact point outgoing mu - outgoing el VS opening angle tar 0",500,-5.,5.,160,0.,0.032);

   ROOT::TThreadedObject<TH1D> h_nhits_pre_0("h_nhits_pre_0","number hits pre tar 0",34,0,34);
   ROOT::TThreadedObject<TH1D> h_nhits_post_0("h_nhits_post_0","number hits post tar 0",34,0,34);
   ROOT::TThreadedObject<TH1D> h_ntracks_post_0("h_ntracks_post_0","number tracks post tar 0",34,0,34);
   ROOT::TThreadedObject<TH1D> h_chi2_vrtx_0("h_chi2_vrtx_0","chi2 vrtx tar 0",300,0.,30.);
   ROOT::TThreadedObject<TH1D> h_Dz_vrtx_0("h_Dz_vrtx_0","(z_vrtx-z_tar) vrtx tar 0",300,-30.,30.);
   ROOT::TThreadedObject<TH1D> h_aco_0("h_aco_0","acoplanarity tar 0",200,-1.,1.);
   ROOT::TThreadedObject<TH1D> h_th_min_0("h_th_min_0","theta min tar 0",250,0.,0.005);
   ROOT::TThreadedObject<TH1D> h_th_max_0("h_th_max_0","theta max tar 0",160,0.,0.032);

   ROOT::TThreadedObject<TH1D> h_preE_th_min_0("h_preE_th_min_0","theta min tar 0 PreEl",250,0.,0.005);
   ROOT::TThreadedObject<TH1D> h_preE_th_max_0("h_preE_th_max_0","theta max tar 0 PreEl",160,0.,0.032);
   ROOT::TThreadedObject<TH1D> h_elastic_0("h_elastic_0","Elasticity tar0",200,-0.001,0.001);

   ROOT::TThreadedObject<TH1D> h_E_th_min_0("h_E_th_min_0","theta min tar 0 Elastic cut",250,0.,0.005);
   ROOT::TThreadedObject<TH1D> h_E_th_max_0("h_E_th_max_0","theta max tar 0 Elastic cu",160,0.,0.032);

   ROOT::TThreadedObject<TH1D> h_final_nhits_pre_0("h_final_nhits_pre_0","number hits pre tar 0 final",10,0,10);
   ROOT::TThreadedObject<TH1D> h_final_nhits_post_0("h_final_nhits_post_0","number hits post tar 0 final",25,0,25);
   ROOT::TThreadedObject<TH1D> h_final_ntracks_post_0("h_final_ntracks_post_0","number tracks post tar 0 final",10,0,10);
   ROOT::TThreadedObject<TH1D> h_final_chi2_vrtx_0("h_final_chi2_vrtx_0","chi2 vrtx tar 0 final",200,0.,20.);
   ROOT::TThreadedObject<TH1D> h_final_Dz_vrtx_0("h_final_Dz_vrtx_0","(z_vrtx-z_tar) vrtx tar 0 final",100,-10.,10.);
   ROOT::TThreadedObject<TH1D> h_final_aco_0("h_final_aco_0","acoplanarity tar 0 final",100,-0.3,0.3);
   ROOT::TThreadedObject<TH1D> h_final_chi2_tr0_0("h_final_chi2_tr0_0","chi2 track0 tar 0 final",200,0.,20.);
   ROOT::TThreadedObject<TH1D> h_final_chi2_trMax_0("h_final_chi2_trMax_0","chi2 track Max tar 0 final",200,0.,20.);
   ROOT::TThreadedObject<TH1D> h_final_chi2_trMin_0("h_final_chi2_trMin_0","chi2 track Min tar 0 final",200,0.,20.);
   ROOT::TThreadedObject<TH2D> h_final_rndm_0("h_final_rndm_0","h_final_rndm_0",200.,-0.01,0.01,100,0.002,0.012);
//   ROOT::TThreadedObject<TGraph> h_final_rndm_0;

   ROOT::TThreadedObject<TH1D> h_x_1("h_x_1","X1 fiducial tar 1",400,-2.,2.);
   ROOT::TThreadedObject<TH1D> h_y_1("h_y_1","Y1 fiducial tar 1",400,-2.,2.);
   ROOT::TThreadedObject<TH1D> h_thx_1("h_thx_1","theta X incoming fiducial tar 1",200,-0.002,0.002);
   ROOT::TThreadedObject<TH1D> h_thy_1("h_thy_1","theta Y incoming fiducial tar 1",200,-0.002,0.002);

   ROOT::TThreadedObject<TH1D> h_IP_1("h_IP_1","Impact point outgoing mu - outgoing el tar 1",500,-5.,5.);
   ROOT::TThreadedObject<TH2D> h_IPop_1("h_IPop_1","Impact point outgoing mu - outgoing el VS opening angle tar 1",500,-5.,5.,160,0.,0.032);

   ROOT::TThreadedObject<TH1D> h_nhits_pre_1("h_nhits_pre_1","number hits pre tar 1",34,0,34);
   ROOT::TThreadedObject<TH1D> h_nhits_post_1("h_nhits_post_1","number hits post tar 1",34,0,34);
   ROOT::TThreadedObject<TH1D> h_ntracks_post_1("h_ntracks_post_1","number tracks post tar 1",34,0,34);
   ROOT::TThreadedObject<TH1D> h_chi2_vrtx_1("h_chi2_vrtx_1","chi2 vrtx tar 1",300,0.,30.);
   ROOT::TThreadedObject<TH1D> h_Dz_vrtx_1("h_Dz_vrtx_1","(z_vrtx-z_tar) vrtx tar 1",300,-30.,30.);
   ROOT::TThreadedObject<TH1D> h_aco_1("h_aco_1","acoplanarity tar 1",200,-1.,1.);
   ROOT::TThreadedObject<TH1D> h_th_min_1("h_th_min_1","theta min tar 1",250,0.,0.005);
   ROOT::TThreadedObject<TH1D> h_th_max_1("h_th_max_1","theta max tar 1",160,0.,0.032);

   ROOT::TThreadedObject<TH1D> h_preEl_th_min_1("h_preE_th_min_1","theta min tar 1 PreEl",250,0.,0.005);
   ROOT::TThreadedObject<TH1D> h_preE_th_max_1("h_preE_th_max_1","theta max tar 1 PreEl",160,0.,0.032);
   ROOT::TThreadedObject<TH1D> h_elastic_1("h_elastic_1","Elasticity tar1",200,-0.001,0.001);

   ROOT::TThreadedObject<TH1D> h_E_th_min_1("h_E_th_min_1","theta min tar 1 Elastic cut",250,0.,0.005);
   ROOT::TThreadedObject<TH1D> h_E_th_max_1("h_E_th_max_1","theta max tar 1 Elastic cu",160,0.,0.032);

   ROOT::TThreadedObject<TH1D> h_final_nhits_pre_1("h_final_nhits_pre_1","number hits pre tar 1 final",10,0,10);
   ROOT::TThreadedObject<TH1D> h_final_nhits_post_1("h_final_nhits_post_1","number hits post tar 1 final",25,0,25);
   ROOT::TThreadedObject<TH1D> h_final_ntracks_post_1("h_final_ntracks_post_1","number tracks post tar 1 final",10,0,10);
   ROOT::TThreadedObject<TH1D> h_final_chi2_vrtx_1("h_final_chi2_vrtx_1","chi2 vrtx tar 1 final",200,0.,20.);
   ROOT::TThreadedObject<TH1D> h_final_Dz_vrtx_1("h_final_Dz_vrtx_1","(z_vrtx-z_tar) vrtx tar 1 final",100,-10.,10.);
   ROOT::TThreadedObject<TH1D> h_final_aco_1("h_final_aco_1","acoplanarity tar 1 final",100,-0.3,0.3);
   ROOT::TThreadedObject<TH1D> h_final_chi2_tr0_1("h_final_chi2_tr0_1","chi2 track0 tar 1 final",200,0.,20.);
   ROOT::TThreadedObject<TH1D> h_final_chi2_trMax_1("h_final_chi2_trMax_1","chi2 track Max tar 1 final",200,0.,20.);
   ROOT::TThreadedObject<TH1D> h_final_chi2_trMin_1("h_final_chi2_trMin_1","chi2 track Min tar 1 final",200,0.,20.);
   ROOT::TThreadedObject<TH2D> h_final_rndm_1("h_final_rndm_1","h_final_rndm_1",200.,-0.01,0.01,100,0.002,0.012);
//   ROOT::TThreadedObject<TGraph> h_final_rndm_1;

 auto myFunction = [&](TTreeReader &myReader) {

     TTreeReaderValue<MUonEEventHeader> RVeventHeader(myReader, "EventHeader");
     TTreeReaderValue<std::vector<MUonERecoOutputTrackAnalysis>> RVtracks(myReader, "ReconstructedTracks");
     TTreeReaderValue<std::vector<MUonERecoOutputVertexAnalysis>> RVvrtx(myReader, "BestVertex");
     TTreeReaderValue<std::vector<MUonERecoOutputHitAnalysis>> RVstubs(myReader, "ReconstructedHits");
//     TTreeReaderValue<std::vector<MUonETrackerStub>> tr_stubs(myReader,"TrackerStubs");
     TTreeReaderValue<bool> isReconstructed(myReader, "IsReconstructed");
        const char* wf="wgt_full";
        TTreeReaderValue<Double_t> wgt_f(myReader,wf);

//     TTreeReaderValue<Bool_t> trig_smi0(myReader, "trigger_single_muon_interaction_0");
//   TTreeReaderValue<Bool_t> trig_smi1(myReader, "trigger_single_muon_interaction_1");

count_tr_1->Fill(-99); h_nvrtx_1->Fill(-99);
count_tr_2->Fill(-99); h_nvrtx_2->Fill(-99);
good_vrtx_1->Fill(-99);
good_vrtx_2->Fill(-99);
fiducial_1->Fill(-99);
fiducial_2->Fill(-99);
pre_elastic_1->Fill(-99);
pre_elastic_2->Fill(-99);
elastic_1->Fill(-99);
elastic_2->Fill(-99);
h_el_1->Fill(-99,-99);
h_el_2->Fill(-99,-99);

     while (myReader.Next()) {

Long64_t entry = myReader.GetCurrentEntry();
auto wgt = *wgt_f;

 auto vrtx = *RVvrtx;
 double chi=0.;
if(vrtx.size()!=0)chi=vrtx.at(0).chi2();

auto& trig_bits = *RVeventHeader;

if(type=="single_muon_interaction_0")count_tr_1->Fill(1,wgt);
if(type=="single_muon_interaction_1")count_tr_2->Fill(1,wgt);


// if(chi!=0 and vrtx.at(0).stationIndex()==1 and *trig_smi0==1){ h_nvrtx_1->Fill(1); if(chi<20) {good_vrtx_1->Fill(1);} }
// if(chi!=0 and vrtx.at(0).stationIndex()==2 and *trig_smi1==1){ h_nvrtx_2->Fill(1); if(chi<20) {good_vrtx_2->Fill(1);} }

//if(*trig_smi0==1){count_tr_1->Fill(1);}
//else if(*trig_smi1==1){count_tr_2->Fill(1);}

double posxIN=99.;//pos_on_track(x0_in,th_inx,z_fix);
double posyIN=99.;//pos_on_track(y0_in,th_iny,z_fix);




std::array<int,3> nTr{0};
std::vector<MUonERecoOutputTrackAnalysis> track_0;
std::vector<MUonERecoOutputTrackAnalysis> track_1;
std::vector<MUonERecoOutputTrackAnalysis> track_vrtx1;
std::vector<MUonERecoOutputTrackAnalysis> track_vrtx2;

/*
         auto tracks = *RVtracks;
         for (auto&& track : tracks) {
        if(track.sector()==0){nTr.at(0)++; track_0.push_back(track);}
        else if(track.sector()==1){nTr.at(1)++; track_1.push_back(track);
				   if(track.index()==mu_out_v.index() or track.index()==e_out_v.index())track_vrtx1.push_back(track);}
        else if(track.sector()==2){nTr.at(2)++;
				   if(track.index()==mu_out_v.index() or track.index()==e_out_v.index())track_vrtx2.push_back(track);}
        }
*/

         auto tracks = *RVtracks;
         for (auto&& track : tracks) {
        if(track.sector()==0){nTr.at(0)++; track_0.push_back(track);}
        else if(track.sector()==1){nTr.at(1)++; track_1.push_back(track);track_vrtx1.push_back(track);
                }
        else if(track.sector()==2){nTr.at(2)++;track_vrtx2.push_back(track);
                }
        }

MUonERecoOutputHitAnalysis muin_4;
MUonERecoOutputHitAnalysis muin_5;

         auto hits = *RVstubs;

if(type=="single_muon_interaction_0" and nTr.at(0)==1){
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
else if(type=="single_muon_interaction_1" and nTr.at(1)==1){
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

int stub0=0;int stub1=0;int stub2=0;

         for (auto&& hit : hits) {
        if(hit.stationID()==0){stub0++;}
        else if(hit.stationID()==1){stub1++;}
        else if(hit.stationID()==2){stub2++;}
        }


if(type=="single_muon_interaction_0" and nTr.at(0)==1 and abs(muin_4.position())<=1.4985 and abs(muin_5.position())<=1.4985 and stub0<10){


	fiducial_1->Fill(1.,wgt);
	h_x_0->Fill(muin_4.position(),wgt);
	h_y_0->Fill(muin_5.position(),wgt);
	h_thx_0->Fill(track_0.at(0).xSlope(),wgt);
	h_thy_0->Fill(track_0.at(0).ySlope(),wgt);

auto& isReco = *isReconstructed;

//std::cout<<"isR entry --? " << entry << " ?--" << std::endl;


	if(chi!=0 and vrtx.at(0).stationIndex()==1){
//	if(chi!=0 and nTr.at(1)>1){
	//h_nvrtx_1->Fill(1);
 MUonERecoOutputTrackAnalysis mu_in_v = vrtx.at(0).incomingMuon();
 MUonERecoOutputTrackAnalysis mu_out_v = vrtx.at(0).outgoingMuon();
 MUonERecoOutputTrackAnalysis e_out_v = vrtx.at(0).outgoingElectron();

	if(chi<20) {good_vrtx_1->Fill(1,wgt);}

	double the_rec=vrtx.at(0).electronTheta();
	double thmu_rec=vrtx.at(0).muonTheta();
        double acoplanarity_v=vrtx.at(0).modifiedAcoplanarity();
        double zpos =vrtx.at(0).zPositionFit();

        double th1=0.;
        double th2=0.;
        if(entry%2==0){th1=the_rec;th2=thmu_rec;}
        else{th2=the_rec;th1=thmu_rec;}

	        double IP_out=vrtx.at(0).distanceBetweenOutgoingTracksAtTargetZ();
//		double IPx=pos_on_track(track_vrtx1.at(1).x0(),track_vrtx1.at(1).xSlope(),(667.3-2.7)) - pos_on_track(track_vrtx1.at(0).x0(),track_vrtx1.at(0).xSlope(),(667.3-2.7));
//              double IPy=pos_on_track(track_vrtx1.at(1).y0(),track_vrtx1.at(1).ySlope(),(667.3-2.7)) - pos_on_track(track_vrtx1.at(0).y0(),track_vrtx1.at(0).ySlope(),(667.3-2.7));
//              IP_out=sqrt(IPx*IPx + IPy*IPy);

                h_IP_0->Fill(IP_out,wgt);

		TVector3 p_e(vrtx.at(0).outgoingElectron().xSlope(),vrtx.at(0).outgoingElectron().ySlope(),1.0); p_e.Unit();
		TVector3 p_mu(vrtx.at(0).outgoingMuon().xSlope(),vrtx.at(0).outgoingMuon().ySlope(),1.0); p_mu.Unit();

		h_IPop_0->Fill(IP_out,p_mu.Angle(p_e),wgt);

		if(IP_out<0.2){


		 h_nvrtx_1->Fill(1,wgt);

		 h_nhits_pre_0->Fill(stub0,wgt);
		 h_nhits_post_0->Fill(stub1,wgt);
		 h_ntracks_post_0->Fill(nTr.at(1),wgt);
		 h_chi2_vrtx_0->Fill(chi,wgt);
		 h_Dz_vrtx_0->Fill(zpos-(667.3-2.7),wgt);
		 h_aco_0->Fill(acoplanarity_v,wgt);
		 h_th_min_0->Fill(thmu_rec,wgt);
		 h_th_max_0->Fill(the_rec,wgt);

		 if(chi<20 and stub1<=15 and abs(acoplanarity_v)<=0.3 and abs(zpos-(667.3-2.7))<3. and the_rec<0.02 and thmu_rec>0.0002){
	          pre_elastic_1->Fill(1,wgt);
		  double Elastic=0.5109989461*0.001*((1+(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*cos(the_rec)*cos(the_rec))/(1-(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*cos(the_rec)*cos(the_rec)));
		  double Elastic2=asin( (sin(the_rec)*sqrt(Elastic*Elastic-0.5109989461*0.001*0.5109989461*0.001))/sqrt( (160+0.5109989461*0.001-Elastic)*(160+0.5109989461*0.001-Elastic)-105.6583745 *0.001*105.6583745 *0.001 ) );

		  h_preE_th_min_0->Fill(thmu_rec,wgt);
		  h_preE_th_max_0->Fill(the_rec,wgt);
		  h_elastic_0->Fill(Elastic2-thmu_rec,wgt);
		  h_el_1->Fill(the_rec,thmu_rec,wgt);

			if(thmu_rec<=Elastic2+0.0002 and thmu_rec>=Elastic2-0.0002){ //if(thmu_rec>0.0002 and thmu_rec<=Elastic2+0.0002 and thmu_rec>=Elastic2-0.0002){
			        elastic_1->Fill(1,wgt);
	        	        h_E_th_min_0->Fill(thmu_rec,wgt);
		                h_E_th_max_0->Fill(the_rec,wgt);
                                h_final_nhits_pre_0->Fill(stub0,wgt);
                                h_final_nhits_post_0->Fill(stub1,wgt);
                                h_final_ntracks_post_0->Fill(nTr.at(1),wgt);
                                h_final_chi2_vrtx_0->Fill(chi,wgt);
                                h_final_Dz_vrtx_0->Fill(zpos-(667.3-2.7),wgt);
				h_final_aco_0->Fill(acoplanarity_v,wgt);
                                h_final_chi2_tr0_0->Fill(mu_in_v.chi2(),wgt);
                                h_final_chi2_trMax_0->Fill(e_out_v.chi2(),wgt);
				h_final_chi2_trMin_0->Fill(mu_out_v.chi2(),wgt);
				h_final_rndm_0->Fill(th1-th2,th1+th2,wgt);
				//h_final_rndm_0.Get()->SetPoint(h_final_rndm_0.Get()->GetN(),th1-th2,th1+th2);
				}

			}
		}
	}
}
else if(type=="single_muon_interaction_1" and nTr.at(1)==1 and abs(muin_4.position())<=1.4985 and abs(muin_5.position())<=1.4985 and stub1<10){


	fiducial_2->Fill(1.,wgt);
        h_x_1->Fill(muin_4.position(),wgt);
        h_y_1->Fill(muin_5.position(),wgt);
        h_thx_1->Fill(track_1.at(0).xSlope(),wgt);
        h_thy_1->Fill(track_1.at(0).ySlope(),wgt);

auto& isReco = *isReconstructed;


	if(chi!=0 and vrtx.at(0).stationIndex()==2){
//	if(chi!=0 and nTr.at(2)>1){
//h_nvrtx_2->Fill(1);

 MUonERecoOutputTrackAnalysis mu_in_v = vrtx.at(0).incomingMuon();
 MUonERecoOutputTrackAnalysis mu_out_v = vrtx.at(0).outgoingMuon();
 MUonERecoOutputTrackAnalysis e_out_v = vrtx.at(0).outgoingElectron();

	if(chi<20) {good_vrtx_2->Fill(1,wgt); }
	double the_rec=vrtx.at(0).electronTheta();
	double thmu_rec=vrtx.at(0).muonTheta();
	double acoplanarity_v=vrtx.at(0).modifiedAcoplanarity();
	double zpos =vrtx.at(0).zPositionFit();

	double th1=0.;
	double th2=0.;
	if(entry%2==0){th1=the_rec;th2=thmu_rec;}
	else{th2=the_rec;th1=thmu_rec;}

        double IP_out=vrtx.at(0).distanceBetweenOutgoingTracksAtTargetZ();
//		double IPx=pos_on_track(track_vrtx2.at(1).x0(),track_vrtx2.at(1).xSlope(),(784.6-3.9)) - pos_on_track(track_vrtx2.at(0).x0(),track_vrtx2.at(0).xSlope(),(784.6-3.9));
//              double IPy=pos_on_track(track_vrtx2.at(1).y0(),track_vrtx2.at(1).ySlope(),(784.6-3.9)) - pos_on_track(track_vrtx2.at(0).y0(),track_vrtx2.at(0).ySlope(),(784.6-3.9));
//              IP_out=sqrt(IPx*IPx + IPy*IPy);

                h_IP_1->Fill(IP_out,wgt);

                TVector3 p_e(vrtx.at(0).outgoingElectron().xSlope(),vrtx.at(0).outgoingElectron().ySlope(),1.0); p_e.Unit();
                TVector3 p_mu(vrtx.at(0).outgoingMuon().xSlope(),vrtx.at(0).outgoingMuon().ySlope(),1.0); p_mu.Unit();

                h_IPop_1->Fill(IP_out,p_mu.Angle(p_e),wgt);

                if(IP_out<0.2){
		 h_nvrtx_2->Fill(1,wgt);
		 h_nhits_pre_1->Fill(stub1,wgt);
		 h_nhits_post_1->Fill(stub2,wgt);
		 h_ntracks_post_1->Fill(nTr.at(2),wgt);
		 h_chi2_vrtx_1->Fill(chi,wgt);
		 h_Dz_vrtx_1->Fill(zpos-(784.6-3.9),wgt);
		 h_aco_1->Fill(acoplanarity_v,wgt);
		 h_th_min_1->Fill(thmu_rec,wgt);
		 h_th_max_1->Fill(the_rec,wgt);

		 if(chi<20 and stub2<=15 and abs(acoplanarity_v)<=0.3 and abs(zpos-(784.6-3.9))<3. and the_rec<0.02 and thmu_rec>0.0002){
		  pre_elastic_2->Fill(1,wgt);
		  double Elastic=0.5109989461*0.001*((1+(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*cos(the_rec)*cos(the_rec))/(1-(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*cos(the_rec)*cos(the_rec)));
		  double Elastic2=asin( (sin(the_rec)*sqrt(Elastic*Elastic-0.5109989461*0.001*0.5109989461*0.001))/sqrt( (160+0.5109989461*0.001-Elastic)*(160+0.5109989461*0.001-Elastic)-105.6583745 *0.001*105.6583745 *0.001 ) );

		  h_preEl_th_min_1->Fill(thmu_rec,wgt);
		  h_preE_th_max_1->Fill(the_rec,wgt);
		  h_elastic_1->Fill(Elastic2-thmu_rec,wgt);
				h_el_2->Fill(the_rec,thmu_rec,wgt);

                        if(thmu_rec<=Elastic2+0.0002 and thmu_rec>=Elastic2-0.0002){ //if(thmu_rec>0.0002 and thmu_rec<=Elastic2+0.0002 and thmu_rec>=Elastic2-0.0002){
                                elastic_2->Fill(1,wgt);
                                h_E_th_min_1->Fill(thmu_rec,wgt);
                                h_E_th_max_1->Fill(the_rec,wgt);
                                h_E_th_min_0->Fill(thmu_rec,wgt);
                                h_E_th_max_0->Fill(the_rec,wgt);
                                h_final_nhits_pre_1->Fill(stub1,wgt);
                                h_final_nhits_post_1->Fill(stub2,wgt);
                                h_final_ntracks_post_1->Fill(nTr.at(2),wgt);
                                h_final_chi2_vrtx_1->Fill(chi,wgt);
                                h_final_Dz_vrtx_1->Fill(zpos-(784.6-3.9),wgt);
                                h_final_aco_1->Fill(acoplanarity_v,wgt);
                                h_final_chi2_tr0_1->Fill(mu_in_v.chi2(),wgt);
                                h_final_chi2_trMax_1->Fill(e_out_v.chi2(),wgt);
                                h_final_chi2_trMin_1->Fill(mu_out_v.chi2(),wgt);
                                h_final_rndm_1->Fill(th1-th2,th1+th2,wgt);
				//h_final_rndm_1.Get()->SetPoint(h_final_rndm_1.Get()->GetN(),th1-th2,th1+th2);
				//h_el_2->Fill(the_rec,thmu_rec);
				}
                        }
		}
	}
}





 } //end of general while
}; //end of ySlopefunction

  tp1.Process(myFunction);


  auto h_nvrtx_1M=h_nvrtx_1.Merge();
  auto h_nvrtx_2M=h_nvrtx_2.Merge();
  auto count_tr_1M=count_tr_1.Merge();
  auto count_tr_2M=count_tr_2.Merge();
  auto good_vrtx_1M=good_vrtx_1.Merge();
  auto good_vrtx_2M=good_vrtx_2.Merge();
  auto fiducial_1M=fiducial_1.Merge();
  auto fiducial_2M=fiducial_2.Merge();
  auto pre_elastic_1M=pre_elastic_1.Merge();
  auto pre_elastic_2M=pre_elastic_2.Merge();
  auto elastic_1M=elastic_1.Merge();
  auto elastic_2M=elastic_2.Merge();
  auto h_el_1M=h_el_1.Merge();
  auto h_el_2M=h_el_2.Merge();

  auto h_x_0M=h_x_0.Merge();
  auto h_y_0M=h_y_0.Merge();
  auto h_thx_0M=h_thx_0.Merge();
  auto h_thy_0M=h_thy_0.Merge();
  auto h_nhits_pre_0M=h_nhits_pre_0.Merge();
  auto h_nhits_post_0M=h_nhits_post_0.Merge();
  auto h_ntracks_post_0M=h_ntracks_post_0.Merge();
  auto h_chi2_vrtx_0M=h_chi2_vrtx_0.Merge();
  auto h_Dz_vrtx_0M=h_Dz_vrtx_0.Merge();
  auto h_aco_0M=h_aco_0.Merge();
  auto h_th_min_0M=h_th_min_0.Merge();
  auto h_th_max_0M=h_th_max_0.Merge();
  auto h_preE_th_min_0M=h_preE_th_min_0.Merge();
  auto h_preE_th_max_0M=h_preE_th_max_0.Merge();
  auto h_elastic_0M=h_elastic_0.Merge();
  auto h_IP_0M=h_IP_0.Merge();
  auto h_IPop_0M=h_IPop_0.Merge();
  auto h_E_th_min_0M=h_E_th_min_0.Merge();
  auto h_E_th_max_0M=h_E_th_max_0.Merge();
  auto h_final_nhits_pre_0M=h_final_nhits_pre_0.Merge();
  auto h_final_nhits_post_0M=h_final_nhits_post_0.Merge();
  auto h_final_ntracks_post_0M=h_final_ntracks_post_0.Merge();
  auto h_final_chi2_vrtx_0M=h_final_chi2_vrtx_0.Merge();
  auto h_final_Dz_vrtx_0M=h_final_Dz_vrtx_0.Merge();
  auto h_final_aco_0M=h_final_aco_0.Merge();
  auto h_final_chi2_tr0_0M=h_final_chi2_tr0_0.Merge();
  auto h_final_chi2_trMax_0M=h_final_chi2_trMax_0.Merge();
  auto h_final_chi2_trMin_0M=h_final_chi2_trMin_0.Merge();
  auto h_final_rndm_0M=h_final_rndm_0.Merge();


  auto h_x_1M=h_x_1.Merge();
  auto h_y_1M=h_y_1.Merge();
  auto h_thx_1M=h_thx_1.Merge();
  auto h_thy_1M=h_thy_1.Merge();
  auto h_nhits_pre_1M=h_nhits_pre_1.Merge();
  auto h_nhits_post_1M=h_nhits_post_1.Merge();
  auto h_ntracks_post_1M=h_ntracks_post_1.Merge();
  auto h_chi2_vrtx_1M=h_chi2_vrtx_1.Merge();
  auto h_Dz_vrtx_1M=h_Dz_vrtx_1.Merge();
  auto h_aco_1M=h_aco_1.Merge();
  auto h_th_min_1M=h_th_min_1.Merge();
  auto h_th_max_1M=h_th_max_1.Merge();
  auto h_preEl_th_min_1M=h_preEl_th_min_1.Merge();
  auto h_preE_th_max_1M=h_preE_th_max_1.Merge();
  auto h_elastic_1M=h_elastic_1.Merge();
  auto h_IP_1M=h_IP_1.Merge();
  auto h_IPop_1M=h_IPop_1.Merge();
  auto h_E_th_min_1M=h_E_th_min_1.Merge();
  auto h_E_th_max_1M=h_E_th_max_1.Merge();
  auto h_final_nhits_pre_1M=h_final_nhits_pre_1.Merge();
  auto h_final_nhits_post_1M=h_final_nhits_post_1.Merge();
  auto h_final_ntracks_post_1M=h_final_ntracks_post_1.Merge();
  auto h_final_chi2_vrtx_1M=h_final_chi2_vrtx_1.Merge();
  auto h_final_Dz_vrtx_1M=h_final_Dz_vrtx_1.Merge();
  auto h_final_aco_1M=h_final_aco_1.Merge();
  auto h_final_chi2_tr0_1M=h_final_chi2_tr0_1.Merge();
  auto h_final_chi2_trMax_1M=h_final_chi2_trMax_1.Merge();
  auto h_final_chi2_trMin_1M=h_final_chi2_trMin_1.Merge();
  auto h_final_rndm_1M=h_final_rndm_1.Merge();

//cout << "good_vrtx_1->Integral() " << good_vrtx_1->Integral() << endl;
//cout << "fiducial_1->Integral() " << fiducial_1->Integral() << endl;

    std::ofstream out(Form("txt/run%s/nocut5mrad_report_%s_nhits%i.txt",run.c_str(),type.c_str(),nhits)); // apre (o crea) il file in scrittura
    if (!out) {
        std::cerr << "Errore nell'aprire il file!" << std::endl;
        return 1;
    }
    out << "Entries " << cbmsim->GetEntries() << std::endl;
    out << "Events with trigger 0 " << count_tr_1M->Integral() << std::endl;
    out << "Events with trigger 1 " << count_tr_2M->Integral() << std::endl;
    out << "Events with fiducial 0 " << fiducial_1M->Integral() << std::endl;
    out << "Events with fiducial 1 " << fiducial_2M->Integral() << std::endl;
    out << "Events with reco vrtx in st01 " << h_nvrtx_1M->Integral()
        << " with chi<20 " << good_vrtx_1M->Integral() << std::endl;
    out << "Events with reco vrtx in st12 " << h_nvrtx_2M->Integral()
        << " with chi<20 " << good_vrtx_2M->Integral() << std::endl;
    out << "Events pre_elastic 0 " << pre_elastic_1M->Integral() << std::endl;
    out << "Events pre_elastic 1 " << pre_elastic_2M->Integral() << std::endl;
    out << "Events elastic 0 " << elastic_1M->Integral() << std::endl;
    out << "Events elastic 1 " << elastic_2M->Integral() << std::endl;


    out << "Ratio fiducial tr0 over entries: " << fiducial_1M->Integral()/cbmsim->GetEntries() << std::endl;
    out << "Ratio fiducial tr1 over entries: " << fiducial_2M->Integral()/cbmsim->GetEntries() << std::endl;
    out << "Ratio any vrtx tr0 over fiducial: " << h_nvrtx_1M->Integral()/fiducial_1M->Integral() << std::endl;
    out << "Ratio any vrtx tr1 over fiducial: " << h_nvrtx_2M->Integral()/fiducial_2M->Integral() << std::endl;
    out << "Ratio pre_elastic tr0 over best vrtx: " << pre_elastic_1M->Integral()/h_nvrtx_1M->Integral() << std::endl;
    out << "Ratio pre_elastic tr1 over best vrtx: " << pre_elastic_2M->Integral()/h_nvrtx_2M->Integral() << std::endl;
    out << "Ratio elastic tr0 over pre_elastic: " << elastic_1M->Integral()/pre_elastic_1M->Integral() << std::endl;
    out << "Ratio elastic tr1 over pre_elastic: " << elastic_2M->Integral()/pre_elastic_2M->Integral() << std::endl;


    out << "Ratio any vrtx tr0 over fiducial: " << h_nvrtx_1M->Integral()/fiducial_1M->Integral() << std::endl;
    out << "Ratio any vrtx tr1 over fiducial: " << h_nvrtx_2M->Integral()/fiducial_2M->Integral() << std::endl;
    out << "Ratio pre_elastic tr0 over fiducial: " << pre_elastic_1M->Integral()/fiducial_1M->Integral() << std::endl;
    out << "Ratio pre_elastic tr1 over fiducial: " << pre_elastic_2M->Integral()/fiducial_2M->Integral() << std::endl;
    out << "Ratio elastic tr0 over fiducial: " << elastic_1M->Integral()/fiducial_1M->Integral() << std::endl;
    out << "Ratio elastic tr1 over fiducial: " << elastic_2M->Integral()/fiducial_2M->Integral() << std::endl;

    out.close();


    TFile out_root(Form("txt/run%s/nocut5mrad_root_%s_nhits%i.root",run.c_str(),type.c_str(),nhits),"recreate");
      if(type=="single_muon_interaction_0"){
	h_x_0M->Write();
	h_y_0M->Write();
	h_thx_0M->Write();
	h_thy_0M->Write();
	h_nhits_pre_0M->Write();
	h_nhits_post_0M->Write();
	h_ntracks_post_0M->Write();
	h_chi2_vrtx_0M->Write();
	h_Dz_vrtx_0M->Write();
	h_aco_0M->Write();
	h_th_min_0M->Write();
	h_th_max_0M->Write();
	h_preE_th_min_0M->Write();
	h_preE_th_max_0M->Write();
	h_elastic_0M->Write();
	h_IP_0M->Write();
	h_IPop_0M->Write();
	h_el_1M->Write();
        h_E_th_min_0M->Write();
        h_E_th_max_0M->Write();
        h_final_nhits_pre_0M->Write();
        h_final_nhits_post_0M->Write();
        h_final_ntracks_post_0M->Write();
        h_final_chi2_vrtx_0M->Write();
        h_final_Dz_vrtx_0M->Write();
        h_final_aco_0M->Write();
	h_final_chi2_tr0_0M->Write();
	h_final_chi2_trMax_0M->Write();
	h_final_chi2_trMin_0M->Write();
//	h_final_rndm_0M->Fit("pol2");
	auto f1 = new TF1("f1","[0]+[1]*x*x",-0.002,0.002);
	h_final_rndm_0M->Fit("f1","R"," ",-0.002,0.002);
	h_final_rndm_0M->Write();
      }
      else if(type=="single_muon_interaction_1"){
        h_x_1M->Write();
        h_y_1M->Write();
        h_thx_1M->Write();
        h_thy_1M->Write();
        h_nhits_pre_1M->Write();
        h_nhits_post_1M->Write();
        h_ntracks_post_1M->Write();
        h_chi2_vrtx_1M->Write();
        h_Dz_vrtx_1M->Write();
        h_aco_1M->Write();
        h_th_min_1M->Write();
        h_th_max_1M->Write();
        h_preEl_th_min_1M->Write();
        h_preE_th_max_1M->Write();
        h_elastic_1M->Write();
        h_IP_1M->Write();
        h_IPop_1M->Write();
	h_el_2M->Write();
        h_E_th_min_1M->Write();
        h_E_th_max_1M->Write();
        h_final_nhits_pre_1M->Write();
        h_final_nhits_post_1M->Write();
        h_final_ntracks_post_1M->Write();
        h_final_chi2_vrtx_1M->Write();
        h_final_Dz_vrtx_1M->Write();
        h_final_aco_1M->Write();
        h_final_chi2_tr0_1M->Write();
        h_final_chi2_trMax_1M->Write();
        h_final_chi2_trMin_1M->Write();
//        h_final_rndm_1M->Fit("pol2");
        auto f1 = new TF1("f1","[0]+[1]*x*x",-0.002,0.002);
        h_final_rndm_1M->Fit("f1","R"," ",-0.002,0.002);
        h_final_rndm_1M->Write();
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

