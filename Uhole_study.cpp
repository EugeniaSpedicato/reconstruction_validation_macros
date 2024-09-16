
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

int Uhole_study(string version,int nhits,string bend){

  int nthreads = 6;
  ROOT::EnableImplicitMT(nthreads);


TChain * cbmsim = new TChain("cbmsim");

if(version=="default" and nhits==2){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig.root");
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_1.root");
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_2.root");
}
else if(version=="default" and nhits==0 and bend=="Bend"){
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_1.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_0hit.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_0hit_1.root");
}
else if(version=="default" and nhits==0 and bend=="noBend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_noBend_0hit.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_noBend_0hit_1.root");
}
else if(version=="ricMis"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_commit_b2ed7c3b.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_commit_b2ed7c3b_1.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_commit_b2ed7c3b_2.root");
}
else if(version=="chi2out50" and nhits==0 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_0hit_chi2out50.root");
}
else if(version=="wip11" and nhits==0 and bend=="Bend"){

cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_3234-3235_new_1hit_skimGA_inputs_0.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_3234-3235_new_1hit_skimGA_inputs_1.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_3234-3235_new_1hit_skimGA_inputs_2.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_3234-3235_new_1hit_skimGA_inputs_3.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_3234-3235_new_1hit_skimGA_inputs_4.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_3234-3235_new_1hit_skimGA_inputs_5.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_3234-3235_new_1hit_skimGA_inputs_6.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_3234-3235_new_1hit_skimGA_inputs_7.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_3234-3235_new_1hit_skimGA_inputs_8.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_3234-3235_new_1hit_skimGA_inputs_9.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_3234-3235_new_1hit_skimGA_inputs_10.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_3234-3235_new_1hit_skimGA_inputs_11.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_3234-3235_new_1hit_skimGA_inputs_12.root");
}
else if(version=="ricMis_noAlign" and nhits==0 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_commit_ricMis_0hit_noalignment_5M.root");
}
else if(version=="uOpposit" and nhits==0 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_0hit_uOpposit.root");
}
else if(version=="noalin" and nhits==0 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_0hit_noalin.root");
}
else if(version=="alin" and nhits==0 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_0hit_alin.root");
}
else if(version=="default" and nhits==-1 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_minus1hit.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_minus1hit_1.root");
}
else if(version=="fileMateusz" and nhits==0 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/15_eventsProblems_0hit.root");
}
else if(version=="converted" and nhits==0 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/commit_b2ed7c3b_MCsignal_RECO_convertedToRealDataFormat_0hit.root");
}
else if(version=="tiltAlin" and nhits==0 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_tiltAlin_0.root");
}
else if(version=="run17" and nhits==-1 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/commit_6c846c96_MCsignal_RECO_m1hit_run17_partial.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/commit_6c846c96_MCsignal_RECO_m1hit_run17_partial_1.root");
}
else if(version=="run8" and nhits==-1 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/commit_6c846c96_MCsignal_RECO_m1hit_run8_partial.root");
}
else if(version=="default_chi2out50" and nhits==2 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/commit_6c846c96_MCsignal_RECO_2hit_run6_chi2out50.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/commit_6c846c96_MCsignal_RECO_2hit_run6_chi2out50_1.root");
}
else if(version=="default_bb35b5de" and nhits==2 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/commit_bb35b5de_RECO_2hit_run6.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/commit_bb35b5de_RECO_2hit_run6_1.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/commit_bb35b5de_RECO_2hit_run6_2.root");
}
else if(version=="default_bb35b5de" and nhits==-1 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/commit_bb35b5de_RECO_m1hit_run6.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/commit_bb35b5de_RECO_m1hit_run6_1.root");
}
else if(version=="default_bb35b5de" and nhits==0 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/commit_bb35b5de_RECO_0hit_run6.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/commit_bb35b5de_RECO_0hit_run6_1.root");
}
else if(version=="my_modifica" and nhits==0 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/all_run6_nuovo_algo_0hit.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/all_run6_nuovo_algo_0hit_1.root");
}

ROOT::TTreeProcessorMT tp1(*cbmsim,nthreads);

      MUonERecoOutput *ReconstructionOutput = 0;


   ROOT::TThreadedObject<TH1D> residual0("residual0","residual 0",80,-0.2,0.2);
   ROOT::TThreadedObject<TH1D> residual1("residual1","residual 1",80,-0.2,0.2);
   ROOT::TThreadedObject<TH1D> residual2("residual2","residual 2",80,-0.2,0.2);
   ROOT::TThreadedObject<TH1D> residual3("residual3","residual 3",80,-0.2,0.2);
   ROOT::TThreadedObject<TH1D> residual4("residual4","residual 4",80,-0.2,0.2);
   ROOT::TThreadedObject<TH1D> residual5("residual5","residual 5",80,-0.2,0.2);

   ROOT::TThreadedObject<TH1D> residual0T("residual0T","residual 0 reco track",80,-0.2,0.2);
   ROOT::TThreadedObject<TH1D> residual1T("residual1T","residual 1 reco track",80,-0.2,0.2);
   ROOT::TThreadedObject<TH1D> residual2T("residual2T","residual 2 reco track",80,-0.2,0.2);
   ROOT::TThreadedObject<TH1D> residual3T("residual3T","residual 3 reco track",80,-0.2,0.2);
   ROOT::TThreadedObject<TH1D> residual4T("residual4T","residual 4 reco track",80,-0.2,0.2);
   ROOT::TThreadedObject<TH1D> residual5T("residual5T","residual 5 reco track",80,-0.2,0.2);

   ROOT::TThreadedObject<TH1D> residual0D("residual0D","residual distance 0",80,-0.2,0.2);
   ROOT::TThreadedObject<TH1D> residual1D("residual1D","residual distance 1",80,-0.2,0.2);
   ROOT::TThreadedObject<TH1D> residual2D("residual2D","residual distance 2",80,-0.2,0.2);
   ROOT::TThreadedObject<TH1D> residual3D("residual3D","residual distance 3",80,-0.2,0.2);
   ROOT::TThreadedObject<TH1D> residual4D("residual4D","residual distance 4",80,-0.2,0.2);
   ROOT::TThreadedObject<TH1D> residual5D("residual5D","residual distance 5",80,-0.2,0.2);



   ROOT::TThreadedObject<TH1D> h_hist_distance0("h_hist_distance0","Distance between stub1 and stub0 in module0",80,-0.2,0.2);
   ROOT::TThreadedObject<TH1D> h_hist_distance_zoom0("h_hist_distance_zoom0","Distance between stub0 and stub1 in module0 zoomed",80,-0.4,0.4);
   ROOT::TThreadedObject<TH1D> h_position0("h_position0","localX of the stub with smaller reisual in module0",1000,-5.,5.);


   ROOT::TThreadedObject<TH1D> h_hist_distance1("h_hist_distance1","Distance between stub1 and stub0 in module1",80,-0.2,0.2);
   ROOT::TThreadedObject<TH1D> h_hist_distance_zoom1("h_hist_distance_zoom1","Distance between stub0 and stub1 in module1 zoomed",80,-0.4,0.4);
   ROOT::TThreadedObject<TH1D> h_position1("h_position1","localX of the stub with smaller reisual in module1",1000,-5.,5.);



   ROOT::TThreadedObject<TH1D> h_hist_distance2("h_hist_distance2","Distance between stub1 and stub0 in module",80,-0.2,0.2);
   ROOT::TThreadedObject<TH1D> h_hist_distance_zoom2("h_hist_distance_zoom2","Distance between stub0 and stub1 in module2 zoomed",80,-0.4,0.4);
   ROOT::TThreadedObject<TH1D> h_position2("h_position2","localX of the stub with smaller reisual in module2",1000,-5.,5.);




   ROOT::TThreadedObject<TH1D> h_hist_distance3("h_hist_distance3","Distance between stub1 and stub0 in module3",80,-0.2,0.2);
   ROOT::TThreadedObject<TH1D> h_hist_distance_zoom3("h_hist_distance_zoom3","Distance between stub0 and stub1 in module3 zoomed",80,-0.4,0.4);
   ROOT::TThreadedObject<TH1D> h_position3("h_position3","localX of the stub with smaller reisual in module3",1000,-5.,5.);


   ROOT::TThreadedObject<TH1D> h_hist_distance4("h_hist_distance4","Distance between stub1 and stub0 in module4",80,-0.2,0.2);
   ROOT::TThreadedObject<TH1D> h_hist_distance_zoom4("h_hist_distance_zoom4","Distance between stub0 and stub1 in module4 zoomed",80,-0.4,0.4);
   ROOT::TThreadedObject<TH1D> h_position4("h_position4","localX of the stub with smaller reisual in module4",1000,-5.,5.);



   ROOT::TThreadedObject<TH1D> h_hist_distance5("h_hist_distance5","Distance between stub1 and stub0 in module5",80,-0.2,0.2);
   ROOT::TThreadedObject<TH1D> h_hist_distance_zoom5("h_hist_distance_zoom5","Distance between stub0 and stub1 in module5 zoomed",80,-0.4,0.4);
   ROOT::TThreadedObject<TH1D> h_position5("h_position5","localX of the stub with smaller reisual in module5",1000,-5.,5.);

	   ROOT::TThreadedObject<TH1D> perp_res0("perp_res0","perp_res mod0",80,-0.2,0.2);
	   ROOT::TThreadedObject<TH1D> perp_res1("perp_res1","perp_res mod1",80,-0.2,0.2);
           ROOT::TThreadedObject<TH1D> perp_res2("perp_res2","perp_res mod2",80,-0.2,0.2);
           ROOT::TThreadedObject<TH1D> perp_res3("perp_res3","perp_res mod3",80,-0.2,0.2);
           ROOT::TThreadedObject<TH1D> perp_res4("perp_res4","perp_res mod4",80,-0.2,0.2);
           ROOT::TThreadedObject<TH1D> perp_res5("perp_res5","perp_res mod3",80,-0.2,0.2);
   ROOT::TThreadedObject<TH2D> h_2d("h2D","theta mu vs theta E with all cuts", 320, 0.,0.032,50,0.,0.005);

           ROOT::TThreadedObject<TH1D> not_used_better("not_used_better","not_used_better",2,0,2);
           ROOT::TThreadedObject<TH1D> used_better("used_better","used_better",2,0,2);

           std::vector<uint64_t> nID_v;

auto trackatZ = [](double q, double m,double z) {return q + (z ) * m;};

//double tilt[6]={0.233,0.233,0.,0.,0.233,0.233};//,0.233,-0.233,0.,0.,0.233,-0.233};
double posZ[6]={18.0218,21.8693,55.3635,56.6205,89.9218,93.7693};//, 119.2218, 123.0693,156.4535, 157.8205,191.1218, 194.9693};
double rot[6]={0  *TMath::DegToRad(),90 *TMath::DegToRad(),135*TMath::DegToRad(),45 *TMath::DegToRad(),0  *TMath::DegToRad(),90 *TMath::DegToRad()};
double Zsens[6]={-0.09,+0.09,+0.2,-0.2,-0.09,+0.09};
double alpha[6]={0  *TMath::DegToRad(), 90 *TMath::DegToRad(), 135*TMath::DegToRad(), 45 *TMath::DegToRad(),0  *TMath::DegToRad(), 90 *TMath::DegToRad()};


double tilt[6]={0.233,0.233,0.,0.,0.233,0.233};
//double offsetX[6]={-0.041151187339928,0.,0.07044990681274964,0.07366436386828577,-0.02216311771211892,0., //prima staz
double offsetX[6]={0.05309770641837935,0.,0.01831731454809257,0.1400357237114471,-0.0564074004379715,0.};//seconda staz
//double offsetY[6]={0.,-0.01962944331210647,-0.07050890136986443, 0.07323347608588825,0.,-0.01011926023475611,//prima staz
double offsetY[6]={0.,0.2324750543793093,-0.01838517522060464,0.1391336460004981,0.,-0.04404067629216706};
//double offset_alpha[6]={-0.1164515010010689*TMath::DegToRad(),-0.03411354091431071*TMath::DegToRad(),0.1369298349690934*TMath::DegToRad(),-0.4577137342911528*TMath::DegToRad(),0.007725616313160143*TMath::DegToRad(),-0.1283674677074805*TMath::DegToRad(),//prima staz
double offset_alpha[6]={0.08712552021046044*TMath::DegToRad(),0.005738391422644699*TMath::DegToRad(),0.1284617932104747*TMath::DegToRad(),0.2515280983199862*TMath::DegToRad(),-0.19965187647738691*TMath::DegToRad(),-0.179018408807327*TMath::DegToRad()};//seconda staz

/*
//alignment tilt
double tilt[6]={0.233-0.0029263031563715987,0.233-0.01494246656307502,0.+0.99568515498159038,0.+1.7808494322921509,0.233-0.008612741424450638,0.233-0.011802765130691497};
double offsetX[6]={0.053097706418379353,0.,0.018317314548092569,0.1400357237114471,-0.056407400437971503,0.};//seconda staz
double offsetY[6]={0.,0.23247505437930929,-0.018385175220604639,0.1391336460004981,0.,-0.044040676292167062};
double offset_alpha[6]={0.087125520210460439*TMath::DegToRad(),0.0057383914226446994*TMath::DegToRad(),0.12846179321047471*TMath::DegToRad(),0.25152809831998618*TMath::DegToRad(),-0.19965187647738691*TMath::DegToRad(),-0.179018408807327*TMath::DegToRad()};//seconda staz
*/


auto trackFit = [&] (std::vector<MUonERecoOutputHit> stubs){

        double tmpA[4*(4+1)/2] = {0.};
        ROOT::Math::SVector<double, 4> B;
        for(int j = 0; j<stubs.size(); j++){

                int linkID= stubs.at(j).moduleID();
                if(stubs.at(j).stationID()==0 or linkID==2 or linkID==3) continue;
//                if(stubs.at(j).stationID()==0) continue;
                double cos_term = cos(rot[linkID]+offset_alpha[linkID]);
                double sin_term = sin(rot[linkID]+offset_alpha[linkID]);
//BEND
                double z_term   = stubs.at(j).z();
                double ui	= stubs.at(j).positionPerpendicular();
                double du2	= 0.002*0.002;

                tmpA[0] += cos_term*cos_term/du2;
                tmpA[1] += cos_term*sin_term/du2;
                tmpA[2] += sin_term*sin_term/du2;
                tmpA[3] += cos_term*cos_term*z_term/du2;
                tmpA[4] += cos_term*sin_term*z_term/du2;
                tmpA[5] += cos_term*cos_term*z_term*z_term/du2;
                tmpA[6] += cos_term*sin_term*z_term/du2;
                tmpA[7] += sin_term*sin_term*z_term/du2;
                tmpA[8] += cos_term*sin_term*z_term*z_term/du2;
                tmpA[9] += sin_term*sin_term*z_term*z_term/du2;

                B(0)    += (ui*cos_term + offsetX[linkID]*cos_term*cos_term + offsetY[linkID]*sin_term*cos_term)/du2;//x0
                B(1)    += (ui*sin_term + offsetX[linkID]*sin_term*cos_term + offsetY[linkID]*sin_term*sin_term)/du2;//y0
                B(2)    += (ui*z_term*cos_term + offsetX[linkID]*z_term*cos_term*cos_term + offsetY[linkID]*z_term*sin_term*cos_term)/du2;//xSlope
                B(3)    += (ui*z_term*sin_term + offsetX[linkID]*z_term*sin_term*cos_term + offsetY[linkID]*z_term*sin_term*sin_term)/du2;//ySlope
        }
       //note MatRepSym, hence not all terms are needed to define the matrix (it's symmetric)
        ROOT::Math::SMatrix<double, 4, 4, ROOT::Math::MatRepSym<double, 4>> A(tmpA, 4*(4+1)/2);

        //I'm inverting A, so that it becomes the covariance matrix:
        if(!A.Invert()) {
                cout<<"ERROR!! FAILED TO INVERT THE DERIVATIVES MATRIX!! CANNOT CREATE THE COVARIANCE MATRIX"<<endl; //failed to invert
        }

	ROOT::Math::SVector<Double_t, 4> c_svd = A*B;

        double *bestFitPars = new double[4];
        for(int i = 0; i < 4; i++) bestFitPars[i] = c_svd[i];

        return bestFitPars;};




auto trackPos = [&](int linkID,double globalz, double *bestFitPars) {
        double xtrack = bestFitPars[2]*globalz + bestFitPars[0];
        double ytrack = bestFitPars[3]*globalz + bestFitPars[1];
        return (xtrack - offsetX[linkID])*cos(rot[linkID]+offset_alpha[linkID]) +
               (ytrack - offsetY[linkID])*sin(rot[linkID]+offset_alpha[linkID]);
};






 auto myFunction = [&](TTreeReader &myReader) {

     TTreeReaderValue<std::vector<MUonERecoOutputTrack>> RVtracks(myReader, "ReconstructedTracks");
     TTreeReaderValue<std::vector<MUonERecoOutputHit>> RVstubs(myReader, "ReconstructedHits");
     TTreeReaderValue<MUonERecoOutputVertex> vrtx(myReader, "BestVertex");

     while (myReader.Next()) {

	 double chi=vrtx->chi2perDegreeOfFreedom();

        std::array<uint16_t,6> nstubs_per_mod_BX_st1{0};
        std::array<std::vector<double>,6> localX{{std::vector<double>(0.0)}};
        std::array<std::vector<double>,6> localZ{{std::vector<double>(0.0)}};
        int n_stubs=0.;
        auto stubs = *RVstubs;
        std::vector<MUonERecoOutputHit> all_hits;
        for(auto&& stub : stubs)
        {
	if(stub.stationID()==1){
                int link= stub.moduleID();
        localX.at(link).push_back(stub.positionPerpendicular());
	localZ.at(link).push_back(stub.z());
        nstubs_per_mod_BX_st1.at(link)+=1;
	                    n_stubs+=1;
	all_hits.push_back(stub);}
	}

int  size1=0;
        std::array<std::vector<double>,6> posT{{std::vector<double>(0.0)}};
        std::array<std::vector<double>,6> posTZ{{std::vector<double>(0.0)}};
        std::vector<MUonERecoOutputHit> hits1;
         auto tracks = *RVtracks;

int sec0=0;
int stubs_muin=0.;

         for (auto&& track : tracks) {
        std::vector<MUonERecoOutputHit> hits_=track.hits();

  if(track.sector()==0){
        stubs_muin=hits_.size();
	sec0++;
                        }
  if(track.sector()==1)
        {
          for(int p=0;p<hits_.size();p++){posT.at(hits_.at(p).moduleID()).push_back(hits_.at(p).positionPerpendicular());posTZ.at(hits_.at(p).moduleID()).push_back(hits_.at(p).z());}
        size1=hits_.size();
	hits1=hits_;
	}
 }

if(sec0==1 and stubs_muin==6 and n_stubs==12 and std::all_of(std::begin(nstubs_per_mod_BX_st1), std::end(nstubs_per_mod_BX_st1), [](int i){return i==2;}) ){


double u=99.;
double v=99.;
for(int m=0;m<6;m++){
      	double real_pos1= localX.at(m).at(1);
        double real_pos0= localX.at(m).at(0);
                            if(m==0){ h_hist_distance0->Fill(real_pos1-real_pos0);
                             h_hist_distance_zoom0->Fill(real_pos1-real_pos0);
                             h_position0->Fill(real_pos1);h_position0->Fill(real_pos0);}

                            if(m==1){ h_hist_distance1->Fill(real_pos1-real_pos0);
                             h_hist_distance_zoom1->Fill(real_pos1-real_pos0);
                             h_position1->Fill(real_pos1);h_position1->Fill(real_pos0);}

                            if(m==2){ //h_hist_distance2->Fill(real_pos1-real_pos0);
                             //h_hist_distance_zoom2->Fill(real_pos1-real_pos0);
                             h_position2->Fill(real_pos1);h_position2->Fill(real_pos0); u=real_pos1-real_pos0;}

                            if(m==3){ //h_hist_distance3->Fill(real_pos1-real_pos0);
                             //h_hist_distance_zoom3->Fill(real_pos1-real_pos0);
                             h_position3->Fill(real_pos1);h_position3->Fill(real_pos0); v=real_pos1-real_pos0;}

                            if(m==4){ h_hist_distance4->Fill(real_pos1-real_pos0);
                             h_hist_distance_zoom4->Fill(real_pos1-real_pos0);
                             h_position4->Fill(real_pos1);h_position4->Fill(real_pos0);}

                            if(m==5){ h_hist_distance5->Fill(real_pos1-real_pos0);
                             h_hist_distance_zoom5->Fill(real_pos1-real_pos0);
                             h_position5->Fill(real_pos1);h_position5->Fill(real_pos0);}

		} //end for


if(chi==0 and size1==6){
double tilt[6]={0.233,0.233,0.,0.,0.233,0.233};
        std::array<std::vector<double>,6> not_used{{std::vector<double>(0.0)}};
        std::array<std::vector<double>,6> not_usedZ{{std::vector<double>(0.0)}};


int aiut=0;
for(int p=0;p<6;p++){
if(localX.at(p).at(0)!=posT.at(p).at(0)){not_used.at(p).push_back(localX.at(p).at(0));not_usedZ.at(p).push_back(localZ.at(p).at(0));}
else if(localX.at(p).at(1)!=posT.at(p).at(0)){not_used.at(p).push_back(localX.at(p).at(1));not_usedZ.at(p).push_back(localZ.at(p).at(1));}
else{aiut=1;}
}

if(aiut==0 and (abs(u) <0.17 or abs(v)<0.17)){
        std::vector<MUonERecoOutputHit> not_used_hits;
         for (int s=0; s<all_hits.size(); s++) {

if(all_hits.at(s).positionPerpendicular()==not_used.at(all_hits.at(s).moduleID()).at(0) ){
not_used_hits.push_back(all_hits.at(s));
		}
        }
     auto paramsT = trackFit(hits1);
     auto params = trackFit(not_used_hits);


  double mx=params[2];
  double my=params[3];
  double mxT=paramsT[2];
  double myT=paramsT[3];



/*TVector3 pmuin;
TVector3 pmu;
TVector3 pe;

         for (auto&& track : tracks) {
        std::vector<MUonERecoOutputHit> hits_=track.hits();
  if(track.sector()==0)         {pmuin.SetXYZ(track.xSlope(),track.ySlope(),1.0);}
	}


pmu.SetXYZ(mxT,myT,1.0);
pe.SetXYZ(mx,my,1.0);

if(pe.Angle(pmuin)>pmu.Angle(pmuin)){h_2d->Fill(pe.Angle(pmuin),pmu.Angle(pmuin));}
else if(pe.Angle(pmuin)<pmu.Angle(pmuin) ){h_2d->Fill(pmu.Angle(pmuin),pe.Angle(pmuin));}


if( pe.Angle(pmuin)>pmu.Angle(pmuin) or pe.Angle(pmuin)<pmu.Angle(pmuin)){*/

TVector3 pmuin;
TVector3 pmu;
TVector3 pe;

         for (auto&& track : tracks) {
        std::vector<MUonERecoOutputHit> hits_=track.hits();
  if(track.sector()==0)         {pmuin.SetXYZ(track.xSlope(),track.ySlope(),1.0);}
        }


pmu.SetXYZ(mxT,myT,1.0);
pe.SetXYZ(mx,my,1.0);

double the_rec=0.;
double thmu_rec=0.;

if(pe.Angle(pmuin)>pmu.Angle(pmuin) and pmu.Angle(pmuin)>0.0003){the_rec=pe.Angle(pmuin); thmu_rec=pmu.Angle(pmuin);}
else if(pe.Angle(pmuin)<pmu.Angle(pmuin) and pe.Angle(pmuin)>0.0003){the_rec=pmu.Angle(pmuin); thmu_rec=pe.Angle(pmuin);}

double Elastic=0.5109989461*0.001*((1+(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*cos(the_rec)*cos(the_rec))/(1-(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*cos(the_rec)*cos(the_rec)));
double Elastic2=asin( (sin(the_rec)*sqrt(Elastic*Elastic-0.5109989461*0.001*0.5109989461*0.001))/sqrt( (160+0.5109989461*0.001-Elastic)*(160+0.5109989461*0.001-Elastic)-105.6583745 *0.001*105.6583745 *0.001 ) );


if( thmu_rec>0.0003 and thmu_rec<=Elastic2+0.0002 and thmu_rec>=Elastic2-0.0002){

//cout << "distance not_used " << not_used.at(2).at(0)- trackPos(2,not_usedZ.at(2).at(0),params) << endl;
//cout << "distance used " << posT.at(2).at(0)- trackPos(2,not_usedZ.at(2).at(0),params) <<endl;

if(not_used.at(2).at(0)- trackPos(2,not_usedZ.at(2).at(0),params) < posT.at(2).at(0)- trackPos(2,not_usedZ.at(2).at(0),params) ){not_used_better->Fill(1.);}
else{used_better->Fill(1.);}

h_2d->Fill(the_rec,thmu_rec);



 residual0T->Fill(posT.at(0).at(0)- trackPos(0,posTZ.at(0).at(0),paramsT) );
 residual1T->Fill(posT.at(1).at(0)- trackPos(1,posTZ.at(1).at(0),paramsT) );
 residual2T->Fill(posT.at(2).at(0)- trackPos(2,posTZ.at(2).at(0),paramsT) );
 residual3T->Fill(posT.at(3).at(0)- trackPos(3,posTZ.at(3).at(0),paramsT) );
 residual4T->Fill(posT.at(4).at(0)- trackPos(4,posTZ.at(4).at(0),paramsT) );
 residual5T->Fill(posT.at(5).at(0)- trackPos(5,posTZ.at(5).at(0),paramsT) );


/* residual0T->Fill(not_used.at(0).at(0)- trackPos(0,posTZ.at(0).at(0),paramsT) );
 residual1T->Fill(not_used.at(1).at(0)- trackPos(1,posTZ.at(1).at(0),paramsT) );
 residual2T->Fill(not_used.at(2).at(0)- trackPos(2,posTZ.at(2).at(0),paramsT) );
 residual3T->Fill(not_used.at(3).at(0)- trackPos(3,posTZ.at(3).at(0),paramsT) );
 residual4T->Fill(not_used.at(4).at(0)- trackPos(4,posTZ.at(4).at(0),paramsT) );
 residual5T->Fill(not_used.at(5).at(0)- trackPos(5,posTZ.at(5).at(0),paramsT) );
*/

 residual0->Fill(not_used.at(0).at(0)- trackPos(0,not_usedZ.at(0).at(0),params) );
 residual1->Fill(not_used.at(1).at(0)- trackPos(1,not_usedZ.at(1).at(0),params) );
 residual2->Fill(not_used.at(2).at(0)- trackPos(2,not_usedZ.at(2).at(0),params) );
 residual3->Fill(not_used.at(3).at(0)- trackPos(3,not_usedZ.at(3).at(0),params) );
 residual4->Fill(not_used.at(4).at(0)- trackPos(4,not_usedZ.at(4).at(0),params) );
 residual5->Fill(not_used.at(5).at(0)- trackPos(5,not_usedZ.at(5).at(0),params) );


/* residual0->Fill(posT.at(0).at(0)- trackPos(0,not_usedZ.at(0).at(0),params) );
 residual1->Fill(posT.at(1).at(0)- trackPos(1,not_usedZ.at(1).at(0),params) );
 residual2->Fill(posT.at(2).at(0)- trackPos(2,not_usedZ.at(2).at(0),params) );
 residual3->Fill(posT.at(3).at(0)- trackPos(3,not_usedZ.at(3).at(0),params) );
 residual4->Fill(posT.at(4).at(0)- trackPos(4,not_usedZ.at(4).at(0),params) );
 residual5->Fill(posT.at(5).at(0)- trackPos(5,not_usedZ.at(5).at(0),params) );*/

 residual0D->Fill(trackPos(0,posTZ.at(0).at(0),paramsT) - trackPos(0,not_usedZ.at(0).at(0),params) );
 residual1D->Fill(trackPos(1,posTZ.at(1).at(0),paramsT) - trackPos(1,not_usedZ.at(1).at(0),params) );
 residual2D->Fill(trackPos(2,posTZ.at(2).at(0),paramsT) - trackPos(2,not_usedZ.at(2).at(0),params) );
 residual3D->Fill(trackPos(3,posTZ.at(3).at(0),paramsT) - trackPos(3,not_usedZ.at(3).at(0),params) );
 residual4D->Fill(trackPos(4,posTZ.at(4).at(0),paramsT) - trackPos(4,not_usedZ.at(4).at(0),params) );
 residual5D->Fill(trackPos(5,posTZ.at(5).at(0),paramsT) - trackPos(5,not_usedZ.at(5).at(0),params) );

h_hist_distance_zoom2->Fill(posT.at(2).at(0)-not_used.at(2).at(0));
h_hist_distance_zoom3->Fill(posT.at(3).at(0)-not_used.at(3).at(0));
h_hist_distance2->Fill(u);
h_hist_distance3->Fill(v);



         for (auto&& track : tracks) {
        std::vector<MUonERecoOutputHit> hits_=track.hits();
 	if(track.sector()==1)
         {
	for(int h=0;h<hits_.size();h++){
                if(hits_.at(h).moduleID()==0)perp_res0->Fill(hits_.at(h).perpendicularResiduum());
                if(hits_.at(h).moduleID()==1)perp_res1->Fill(hits_.at(h).perpendicularResiduum());
                if(hits_.at(h).moduleID()==2)perp_res2->Fill(hits_.at(h).perpendicularResiduum());
                if(hits_.at(h).moduleID()==3)perp_res3->Fill(hits_.at(h).perpendicularResiduum());
		if(hits_.at(h).moduleID()==4)perp_res4->Fill(hits_.at(h).perpendicularResiduum());
		if(hits_.at(h).moduleID()==5)perp_res5->Fill(hits_.at(h).perpendicularResiduum());
if(h==0){uint64_t nID = static_cast<uint64_t>(hits_.at(h).superID())*3564 + static_cast<uint64_t>(hits_.at(h).bx()); nID_v.push_back(nID);}
			}
		 }
	 }


	}//0.2mrad
  }//aiut==0
}

	}//end ==12
}//end while

}; //end of myfunction

  tp1.Process(myFunction);


cout <<"nID_v size " << nID_v.size() << endl;
/*for( int f=0; f<nID_v.size(); f++){
cout << nID_v.at(f) <<", ";
}*/

   TFile *fout = TFile::Open("/home/espedica/tr2023_promptanalysis/RiccardoP/nID_v.root", "RECREATE");

   fout->WriteObject(&nID_v, "nID_v"); // I store the vector in the TFile
   fout->Close();

auto h_2dM=h_2d.Merge();
auto perp_res0M=perp_res0.Merge();
auto perp_res1M=perp_res1.Merge();
auto perp_res2M=perp_res2.Merge();
auto perp_res3M=perp_res3.Merge();
auto perp_res4M=perp_res4.Merge();
auto perp_res5M=perp_res5.Merge();


auto not_used_betterM=not_used_better.Merge();
auto used_betterM=used_better.Merge();

cout << "Not used hit closer to not tracked track " << not_used_betterM->Integral() << endl;
cout << "Used hit closer to not tracked track " << used_betterM->Integral() << endl;


auto h_hist_distance0M = h_hist_distance0.Merge();
auto h_hist_distance_zoom0M = h_hist_distance_zoom0.Merge();
auto h_position0M =h_position0.Merge();

auto h_hist_distance1M = h_hist_distance1.Merge();
auto h_hist_distance_zoom1M = h_hist_distance_zoom1.Merge();
auto h_position1M =h_position1.Merge();

auto h_hist_distance2M = h_hist_distance2.Merge();
auto h_hist_distance_zoom2M = h_hist_distance_zoom2.Merge();
auto h_position2M =h_position2.Merge();

auto h_hist_distance3M = h_hist_distance3.Merge();
auto h_hist_distance_zoom3M = h_hist_distance_zoom3.Merge();
auto h_position3M =h_position3.Merge();

auto h_hist_distance4M = h_hist_distance4.Merge();
auto h_hist_distance_zoom4M = h_hist_distance_zoom4.Merge();
auto h_position4M =h_position4.Merge();

auto h_hist_distance5M = h_hist_distance5.Merge();
auto h_hist_distance_zoom5M = h_hist_distance_zoom5.Merge();
auto h_position5M =h_position5.Merge();

TCanvas t2("t2","t2",700,700);
h_2dM->Draw("COLZ");
t2.SaveAs(Form("hole/%s_h2d_elasticCurve%s_%dhit.pdf",version.c_str(),bend.c_str(),nhits));

TCanvas n1("n1","n1",1000,1000);
n1.Divide(2,3);
n1.cd(1);
h_hist_distance0M->Draw("hist");
n1.cd(2);
h_hist_distance1M->Draw("hist");
n1.cd(3);
h_hist_distance2M->Draw("hist");
n1.cd(4);
h_hist_distance3M->Draw("hist");
n1.cd(5);
h_hist_distance4M->Draw("hist");
n1.cd(6);
h_hist_distance5M->Draw("hist");
n1.SaveAs(Form("hole/%s_distance_stubs_elasticCurve%s_%dhit.pdf",version.c_str(),bend.c_str(),nhits));

TCanvas n1_zoom("n1_zoom","n1_zoom",1000,1000);
n1_zoom.Divide(2,3);
n1_zoom.cd(1);
h_hist_distance2M->Draw("hist");//h_hist_distance_zoom0M->Draw("hist");
n1_zoom.cd(2);
h_hist_distance3M->Draw("hist");//h_hist_distance_zoom1M->Draw("hist");
n1_zoom.cd(3);
h_hist_distance_zoom2M->Draw("hist");
n1_zoom.cd(4);
h_hist_distance_zoom3M->Draw("hist");
n1_zoom.cd(5);
//h_hist_distance_zoom4M->Draw("hist");
n1_zoom.cd(6);
//h_hist_distance_zoom5M->Draw("hist");
n1_zoom.SaveAs(Form("hole/%s_distance_stubs_zoom_elasticCurve%s_%dhit.pdf",version.c_str(),bend.c_str(),nhits));


auto residual0M = residual0.Merge();
auto residual1M = residual1.Merge();
auto residual2M = residual2.Merge();
auto residual3M = residual3.Merge();
auto residual4M = residual4.Merge();
auto residual5M = residual5.Merge();

TCanvas r("r","r",1400,2100);
r.Divide(2,3);
r.cd(1);
residual0M->Draw("hist");
r.cd(2);
residual1M->Draw("hist");
r.cd(3);
residual2M->Draw("hist");
r.cd(4);
residual3M->Draw("hist");
r.cd(5);
residual4M->Draw("hist");
r.cd(6);
residual5M->Draw("hist");
r.SaveAs(Form("hole/%s_residuals_notrack_elasticCurve%s_%dhit.pdf",version.c_str(),bend.c_str(),nhits));


auto residual0TM = residual0T.Merge();
auto residual1TM = residual1T.Merge();
auto residual2TM = residual2T.Merge();
auto residual3TM = residual3T.Merge();
auto residual4TM = residual4T.Merge();
auto residual5TM = residual5T.Merge();

TCanvas rT("rT","rT",1400,2100);
rT.Divide(2,3);
rT.cd(1);
residual0TM->Draw("hist");
rT.cd(2);
residual1TM->Draw("hist");
rT.cd(3);
residual2TM->Draw("hist");
rT.cd(4);
residual3TM->Draw("hist");
rT.cd(5);
residual4TM->Draw("hist");
rT.cd(6);
residual5TM->Draw("hist");
rT.SaveAs(Form("hole/%s_residuals_track_elasticCurve%s_%dhit.pdf",version.c_str(),bend.c_str(),nhits));

auto residual0DM = residual0D.Merge();
auto residual1DM = residual1D.Merge();
auto residual2DM = residual2D.Merge();
auto residual3DM = residual3D.Merge();
auto residual4DM = residual4D.Merge();
auto residual5DM = residual5D.Merge();

TCanvas rD("rD","rD",1400,2100);
rD.Divide(2,3);
rD.cd(1);
residual0DM->Draw("hist");
rD.cd(2);
residual1DM->Draw("hist");
rD.cd(3);
residual2DM->Draw("hist");
rD.cd(4);
residual3DM->Draw("hist");
rD.cd(5);
residual4DM->Draw("hist");
rD.cd(6);
residual5DM->Draw("hist");
rD.SaveAs(Form("hole/%s_residuals_disttrack_elasticCurve%s_%dhit.pdf",version.c_str(),bend.c_str(),nhits));

TCanvas pr("pr","pr",1400,2100);
pr.Divide(2,3);
pr.cd(1);
perp_res0->Draw("hist");
pr.cd(2);
perp_res1->Draw("hist");
pr.cd(3);
perp_res2->Draw("hist");
pr.cd(4);
perp_res3->Draw("hist");
pr.cd(5);
perp_res4->Draw("hist");
pr.cd(6);
perp_res5->Draw("hist");
pr.SaveAs(Form("hole/%s_perp_res_elasticCurve%s_%dhit.pdf",version.c_str(),bend.c_str(),nhits));


return 1;
}

