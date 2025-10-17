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


int report_count_snakemake(int run, string type, string filen, int nhits){

  int nthreads = 7;


  ROOT::EnableImplicitMT(nthreads);


TChain * cbmsim = new TChain("cbmsim");
TChain * cbmsim_g = new TChain("cbmsim");



cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/TB2025/run%i/%s/%s_%ihit.root",run,type.c_str(),filen.c_str(),nhits));
cout << Form("/mnt/raid10/DATA/espedica/fairmu/TB2025/run%i/%s/%s_%ihit.root",run,type.c_str(),filen.c_str(),nhits) << endl;

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

   ROOT::TThreadedObject<TH1D> h_nhits_pre_0("h_nhits_pre_0","number hits pre tar 0",34,0,34);
   ROOT::TThreadedObject<TH1D> h_nhits_post_0("h_nhits_post_0","number hits post tar 0",34,0,34);
   ROOT::TThreadedObject<TH1D> h_ntracks_post_0("h_ntracks_post_0","number tracks post tar 0",34,0,34);
   ROOT::TThreadedObject<TH1D> h_chi2_vrtx_0("h_chi2_vrtx_0","chi2 vrtx tar 0",300,0.,30.);
   ROOT::TThreadedObject<TH1D> h_Dz_vrtx_0("h_Dz_vrtx_0","(z_vrtx-z_tar) vrtx tar 0",300,-30.,30.);
   ROOT::TThreadedObject<TH1D> h_aco_0("h_aco_0","acoplanarity tar 0",200,-1.,1.);
   ROOT::TThreadedObject<TH1D> h_th_min_0("h_th_min_0","theta min tar 0",250,0.,0.005);
   ROOT::TThreadedObject<TH1D> h_th_max_0("h_th_max_0","theta max tar 0",160,0.,0.032);

   ROOT::TThreadedObject<TH1D> h_preEl_th_min_0("h_preE_th_min_0","theta min tar 0 PreEl",250,0.,0.005);
   ROOT::TThreadedObject<TH1D> h_preE_th_max_0("h_preE_th_max_0","theta max tar 0 PreEl",160,0.,0.032);
   ROOT::TThreadedObject<TH1D> h_elastic_0("h_elastic_0","Elasticity tar0",200,-0.001,0.001);

   ROOT::TThreadedObject<TH1D> h_x_1("h_x_1","X1 fiducial tar 1",400,-2.,2.);
   ROOT::TThreadedObject<TH1D> h_y_1("h_y_1","Y1 fiducial tar 1",400,-2.,2.);
   ROOT::TThreadedObject<TH1D> h_thx_1("h_thx_1","theta X incoming fiducial tar 1",200,-0.002,0.002);
   ROOT::TThreadedObject<TH1D> h_thy_1("h_thy_1","theta Y incoming fiducial tar 1",200,-0.002,0.002);

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


 auto myFunction = [&](TTreeReader &myReader) {


     TTreeReaderValue<std::vector<MUonERecoOutputTrackAnalysis>> RVtracks(myReader, "ReconstructedTracks");
     TTreeReaderValue<MUonERecoOutputVertexAnalysis> vrtx(myReader, "BestVertex");
     TTreeReaderValue<std::vector<MUonERecoOutputVertexAnalysis>> RVall_vrtx(myReader, "ReconstructedVertices");
     TTreeReaderValue<std::vector<MUonERecoOutputHitAnalysis>> RVstubs(myReader, "ReconstructedHits");
//     TTreeReaderValue<std::vector<MUonETrackerStub>> tr_stubs(myReader,"TrackerStubs");

     TTreeReaderValue<Bool_t> trig_smi0(myReader, "TriggerSingleMuonInteraction0");
     TTreeReaderValue<Bool_t> trig_smi1(myReader, "TriggerSingleMuonInteraction1");

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



 double chi=vrtx->chi2perDegreeOfFreedom();

if(*trig_smi0==1)count_tr_1->Fill(1);
if(*trig_smi1==1)count_tr_2->Fill(1);


// if(chi!=0 and vrtx->stationIndex()==1 and *trig_smi0==1){ h_nvrtx_1->Fill(1); if(chi<20) {good_vrtx_1->Fill(1);} }
// if(chi!=0 and vrtx->stationIndex()==2 and *trig_smi1==1){ h_nvrtx_2->Fill(1); if(chi<20) {good_vrtx_2->Fill(1);} }

//if(chi!=0 and vrtx->stationIndex()==1 and *trig_smi0==1){ h_nvrtx_1->Fill(1);}
//else if(chi!=0 and vrtx->stationIndex()==2 and *trig_smi1==1){ h_nvrtx_2->Fill(1);}

double posxIN=99.;//pos_on_track(x0_in,th_inx,z_fix);
double posyIN=99.;//pos_on_track(y0_in,th_iny,z_fix);


 MUonERecoOutputTrackAnalysis mu_in_v = vrtx->incomingMuon();
 MUonERecoOutputTrackAnalysis mu_out_v = vrtx->outgoingMuon();
 MUonERecoOutputTrackAnalysis e_out_v = vrtx->outgoingElectron();


std::array<int,3> nTr{0};
std::vector<MUonERecoOutputTrackAnalysis> track_0;
std::vector<MUonERecoOutputTrackAnalysis> track_1;
std::vector<MUonERecoOutputTrackAnalysis> track_2;

         auto tracks = *RVtracks;
         for (auto&& track : tracks) {
        if(track.sector()==0){nTr.at(0)++; track_0.push_back(track);}
        else if(track.sector()==1){nTr.at(1)++; track_1.push_back(track);}
        else if(track.sector()==2){nTr.at(2)++; track_2.push_back(track);}
//        if(track.sector()==2){nTr.at(2)++;}
        }

MUonERecoOutputHitAnalysis muin_4;
MUonERecoOutputHitAnalysis muin_5;

if(*trig_smi0==1 and nTr.at(0)==1){
	for(auto&& h : track_0.at(0).hits()){
	 if(h.moduleID()==4){muin_4 = h;}
	 else if(h.moduleID()==5){muin_5 = h;}
	}
}
else if(*trig_smi1==1 and nTr.at(1)==1){
	for(auto&& h : track_1.at(0).hits()){
         if(h.moduleID()==4){muin_4 = h;}
         else if(h.moduleID()==5){muin_5 = h;}
        }
}


int stub0=0;int stub1=0;int stub2=0;

         auto hits = *RVstubs;
         for (auto&& hit : hits) {
        if(hit.stationID()==0){stub0++;}
        else if(hit.stationID()==1){stub1++;}
        else if(hit.stationID()==2){stub2++;}
        }

if(*trig_smi0==1 and nTr.at(0)==1 and abs(muin_4.position())<=1.4985 and abs(muin_5.position())<=1.4985 and stub0<10){
	fiducial_1->Fill(1.);
	h_x_0->Fill(muin_4.position());
	h_y_0->Fill(muin_5.position());
	h_thx_0->Fill(track_0.at(0).xSlope());
	h_thy_0->Fill(track_0.at(0).ySlope());


	if(chi!=0 and vrtx->stationIndex()==1){ h_nvrtx_1->Fill(1); if(chi<20) {good_vrtx_1->Fill(1);}
	double the_rec=vrtx->electronTheta();
	double thmu_rec=vrtx->muonTheta();
        double acoplanarity_v=vrtx->modifiedAcoplanarity();
        double zpos =vrtx->zPositionFit();

h_nhits_pre_0->Fill(stub0);
h_nhits_post_0->Fill(stub1);
h_ntracks_post_0->Fill(nTr.at(1));
h_chi2_vrtx_0->Fill(chi);
h_Dz_vrtx_0->Fill(zpos-(667.3-2.7));
h_aco_0->Fill(acoplanarity_v);
h_th_min_0->Fill(thmu_rec);
h_th_max_0->Fill(the_rec);

		if(chi<20 and stub1<=15 and abs(acoplanarity_v)<=0.3 and abs(zpos-(667.3-2.7))<2. and the_rec<0.02 and the_rec>0.005){
	        pre_elastic_1->Fill(1);
		double Elastic=0.5109989461*0.001*((1+(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*cos(the_rec)*cos(the_rec))/(1-(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*cos(the_rec)*cos(the_rec)));
		double Elastic2=asin( (sin(the_rec)*sqrt(Elastic*Elastic-0.5109989461*0.001*0.5109989461*0.001))/sqrt( (160+0.5109989461*0.001-Elastic)*(160+0.5109989461*0.001-Elastic)-105.6583745 *0.001*105.6583745 *0.001 ) );

h_preEl_th_min_0->Fill(thmu_rec);
h_preE_th_max_0->Fill(the_rec);
h_elastic_0->Fill(Elastic2-thmu_rec);

			if(thmu_rec<=Elastic2+0.0002 and thmu_rec>=Elastic2-0.0002){ //if(thmu_rec>0.0002 and thmu_rec<=Elastic2+0.0002 and thmu_rec>=Elastic2-0.0002){
			        elastic_1->Fill(1);
				h_el_1->Fill(the_rec,thmu_rec);
			}
		}
	}
}
else if(*trig_smi1==1 and nTr.at(1)==1 and abs(muin_4.position())<=1.4985 and abs(muin_5.position())<=1.4985 and stub1<10){
	fiducial_2->Fill(1.);
        h_x_1->Fill(muin_4.position());
        h_y_1->Fill(muin_5.position());
        h_thx_1->Fill(track_1.at(0).xSlope());
        h_thy_1->Fill(track_1.at(0).ySlope());

	if(chi!=0 and vrtx->stationIndex()==2){ h_nvrtx_2->Fill(1); if(chi<20) {good_vrtx_2->Fill(1);}
	double the_rec=vrtx->electronTheta();
	double thmu_rec=vrtx->muonTheta();
	double acoplanarity_v=vrtx->modifiedAcoplanarity();
	double zpos =vrtx->zPositionFit();

h_nhits_pre_1->Fill(stub1);
h_nhits_post_1->Fill(stub2);
h_ntracks_post_1->Fill(nTr.at(2));
h_chi2_vrtx_1->Fill(chi);
h_Dz_vrtx_1->Fill(zpos-(784.6-3.9));
h_aco_1->Fill(acoplanarity_v);
h_th_min_1->Fill(thmu_rec);
h_th_max_1->Fill(the_rec);

		if(chi<20 and stub2<=15 and abs(acoplanarity_v)<=0.3 and abs(zpos-(784.6-3.9))<2. and the_rec<0.02 and the_rec>0.005){
		pre_elastic_2->Fill(1);
		double Elastic=0.5109989461*0.001*((1+(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*cos(the_rec)*cos(the_rec))/(1-(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*cos(the_rec)*cos(the_rec)));
		double Elastic2=asin( (sin(the_rec)*sqrt(Elastic*Elastic-0.5109989461*0.001*0.5109989461*0.001))/sqrt( (160+0.5109989461*0.001-Elastic)*(160+0.5109989461*0.001-Elastic)-105.6583745 *0.001*105.6583745 *0.001 ) );

h_preEl_th_min_1->Fill(thmu_rec);
h_preE_th_max_1->Fill(the_rec);
h_elastic_1->Fill(Elastic2-thmu_rec);

                        if(thmu_rec<=Elastic2+0.0002 and thmu_rec>=Elastic2-0.0002){ //if(thmu_rec>0.0002 and thmu_rec<=Elastic2+0.0002 and thmu_rec>=Elastic2-0.0002){
                                elastic_2->Fill(1);
				h_el_2->Fill(the_rec,thmu_rec);
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
  auto h_preEl_th_min_0M=h_preEl_th_min_0.Merge();
  auto h_preE_th_max_0M=h_preE_th_max_0.Merge();
  auto h_elastic_0M=h_elastic_0.Merge();

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
//cout << "good_vrtx_1->Integral() " << good_vrtx_1->Integral() << endl;
//cout << "fiducial_1->Integral() " << fiducial_1->Integral() << endl;

    std::ofstream out(Form("txt/run%i/report_%s_%s_nhits%i.txt",run,type.c_str(),filen.c_str(),nhits)); // apre (o crea) il file in scrittura
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


    TFile out_root(Form("txt/run%i/root_%s_%s_nhits%i.root",run,type.c_str(),filen.c_str(),nhits),"recreate");
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
	h_preEl_th_min_0M->Write();
	h_preE_th_max_0M->Write();
	h_elastic_0M->Write();
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
