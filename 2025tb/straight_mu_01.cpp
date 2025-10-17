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

int straight_mu_01(string sample, string type){

  int nthreads = 6;
if(sample=="data") nthreads = 16;

  ROOT::EnableImplicitMT(nthreads);

TChain * cbmsim = new TChain("cbmsim");

if(sample=="run8" and type=="passing_muon"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run8/passing_muon/reco/passing_muon_muedaq03-1750191365.root");
}
else if(sample=="data_flip"){ //run5
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/wip16_run6_goldenMu_newAlignment_0hit_Flip.root");
}

ROOT::TTreeProcessorMT tp1(*cbmsim,nthreads);

      MUonERecoOutputAnalysis *ReconstructionOutput = 0;


TH1::SetDefaultSumw2(kTRUE);

   ROOT::TThreadedObject<TH1D> theta_mu2("theta_mu2", "Muon scattering reco angles second station",100,0.,0.005);//50,0.,0.005);
   ROOT::TThreadedObject<TH1D> theta_mu("theta_mu", "Muon scattering reco angles first station",100,0.,0.005);//50,0.,0.005);
   ROOT::TThreadedObject<TH1D> diff_theta("diff_theta", "Muon scattering reco angles first station - second station",200,-0.00025,0.00025);//50,0.,0.005);
   ROOT::TThreadedObject<TH1D> diff_thetax("diff_thetax", "Muon scattering reco angles first station - second station THETA_X",200,-0.00025,0.00025);//50,0.,0.005);
   ROOT::TThreadedObject<TH1D> diff_thetay("diff_thetay", "Muon scattering reco angles first station - second station THETA_Y",200,-0.00025,0.00025);//50,0.,0.005);
   ROOT::TThreadedObject<TH1D> diff_theta_smeared("diff_theta_smeared", "Muon scattering reco angles first station - second station smeared",200,-0.00025,0.00025);//50,0.,0.005);
   ROOT::TThreadedObject<TH1D> ms_effect("ms_effect","Effect ms",200,0.018,0.038);


 auto myFunction = [&](TTreeReader &myReader) {

     TTreeReaderValue<std::vector<MUonERecoOutputTrackAnalysis>> RVtracks(myReader, "ReconstructedTracks");
     TTreeReaderValue<std::vector<MUonERecoOutputHitAnalysis>> RVstubs(myReader, "ReconstructedHits");


     while (myReader.Next()) {



double z_fix=912.7;


double th_inx,th_iny,x0_in,y0_in;
double th_inx2,th_iny2,x0_in2,y0_in2;

double th_inx_smeared,th_iny_smeared;
double th_inx2_smeared,th_iny2_smeared;


double chi2_muin;
double chi2_muin2;
int sec0=0;
int sec1=0;
int stubs_muin=0.;
double th_muin=0.;
int stubs_muin2=0.;
double th_muin2=0.;

double th_muin_smeared=0.;
double th_muin2_smeared=0.;

TVector3 p_muin,p_muin2;
TVector3 p_muin_smeared,p_muin2_smeared;

         auto tracks = *RVtracks;
         for (auto&& track : tracks) {
        if(track.sector()==0) sec0++;
        if(track.sector()==1) sec1++;
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
        }



MUonERecoOutputTrackAnalysis t_mu;
MUonERecoOutputTrackAnalysis t_mu2;

         for (auto&& track : tracks) {
        std::vector<MUonERecoOutputHitAnalysis> hits_=track.hits();

        if(track.sector()==0 and sec0==1){
	stubs_muin=hits_.size();
        th_inx=track.xSlope();
        th_iny=track.ySlope();
        th_inx_smeared=track.xSlope()+gRandom->Gaus(0.,1.7e-05);//gRandom->Gaus(0.,2.3e-05);
        th_iny_smeared=track.ySlope()+gRandom->Gaus(0.,1.7e-05);
        x0_in=track.x0();
        y0_in=track.y0();
        chi2_muin=track.chi2perDegreeOfFreedom();
        p_muin.SetXYZ(th_inx,th_iny,1.0);
        p_muin=p_muin.Unit();
        th_muin=p_muin.Theta();
        p_muin_smeared.SetXYZ(th_inx_smeared,th_iny_smeared,1.0);
        p_muin_smeared=p_muin_smeared.Unit();
        th_muin_smeared=p_muin_smeared.Theta();
	t_mu=track;
                        }
	if(track.sector()==1 and sec1==1)// and stub1<8)
	{
        stubs_muin2=hits_.size();
        th_inx2=track.xSlope();
        th_iny2=track.ySlope();
        th_inx2_smeared=track.xSlope()+gRandom->Gaus(0.,1.7e-05);//notar:1.6e-05, tar:2.3e-05
        th_iny2_smeared=track.ySlope()+gRandom->Gaus(0.,1.7e-05);
        x0_in2=track.x0();
        y0_in2=track.y0();
        chi2_muin2=track.chi2perDegreeOfFreedom();
        p_muin2.SetXYZ(th_inx2,th_iny2,1.0);
        p_muin2=p_muin2.Unit();
        th_muin2=p_muin2.Theta();
        p_muin2_smeared.SetXYZ(th_inx2_smeared,th_iny2_smeared,1.0);
        p_muin2_smeared=p_muin2_smeared.Unit();
        th_muin2_smeared=p_muin2_smeared.Theta();
       t_mu2=track;
	}
 }


if( abs(posxIN)<=3 and abs(posyIN)<=3 and stub0==6 and sec0==1 and sec1==1 and th_muin2!=0 and th_muin<0.004){

        theta_mu->Fill(th_muin);
        theta_mu2->Fill(th_muin2);
	diff_theta->Fill(th_muin-th_muin2);
	diff_theta_smeared->Fill(th_muin_smeared-th_muin2_smeared);
 if(sample=="data"){
	diff_thetax->Fill(th_inx-th_inx2);
	diff_thetay->Fill(th_iny-th_iny2);
	}
  else{
        diff_thetax->Fill(th_inx_smeared-th_inx2_smeared);
        diff_thetay->Fill(th_iny_smeared-th_iny2_smeared);
	}
double p=sqrt(160.*160. - 0.105658*0.105658);
double ms=pow(0.0136/p * sqrt(1.5/19.32) * (1. + 0.038 * log(1.5/19.32)),2)+pow(0.0136/p * sqrt(0.64*6/9.37) * (1. + 0.038 * log(0.64*6/9.37)),2);

	ms_effect->Fill(ms);


	}//mu_in

  } //end of general while
}; //end of ySlopefunction

  tp1.Process(myFunction);

  auto theta_muMerged = theta_mu.Merge();
  auto theta_mu2Merged = theta_mu2.Merge();
  auto diff_thetaMerged = diff_theta.Merge();
  auto diff_theta_smearedM=diff_theta_smeared.Merge();
  auto diff_thetaxM=diff_thetax.Merge();
  auto diff_thetayM=diff_thetay.Merge();


theta_muMerged->SaveAs(Form("dir_straight_mu/theta_mu_gm_%s.root",sample.c_str()));
theta_mu2Merged->SaveAs(Form("dir_straight_mu/theta_mu2_gm_%s.root",sample.c_str()));
diff_thetaMerged->SaveAs(Form("dir_straight_mu/diff_theta_mu_gm_%s.root",sample.c_str()));
if(sample=="mc" or sample=="mc_notar" or sample=="mc9000") {diff_theta_smearedM->SaveAs(Form("dir_straight_mu/diff_theta_smeared_mu_gm_%s.root",sample.c_str()));
							    diff_thetaxM->SaveAs(Form("dir_straight_mu/diff_thetaX_smeared_mu_gm_%s.root",sample.c_str()));
							    diff_thetayM->SaveAs(Form("dir_straight_mu/diff_thetaY_smeared_mu_gm_%s.root",sample.c_str()));}

TF1 *g1 = new TF1("g1", "gaus");
TF1 *g2 = new TF1("g2", "gaus");
diff_thetaMerged->Fit("g1","R","",-0.00015,0.00015);
if(sample=="mc" or sample=="mc_notar" or sample=="mc9000") diff_theta_smearedM->Fit("g2","R","",-0.00015,0.00015);

cout << "sigma: " << g1->GetParameter(2) << endl;
if(sample=="mc" or sample=="mc_notar" or sample=="mc9000") cout << "sigma smeared: " << g2->GetParameter(2) << endl;



TCanvas a("a","a",1500,1000);
a.Divide(2,3);
a.cd(1);
theta_muMerged->Draw("hist");
a.cd(2);
theta_mu2Merged->Draw("hist");
a.cd(3);
diff_thetaMerged->Draw("hist");
if(sample=="mc" or sample=="mc_notar" or sample=="mc9000"){diff_theta_smearedM->SetLineColor(kRed);
		 diff_theta_smearedM->Draw("hist sames");}
a.cd(4);
ms_effect->Draw("hist");
a.cd(5);
diff_thetaxM->Draw("hist");
a.cd(6);
diff_thetayM->Draw("hist");
a.SaveAs(Form("dir_straight_mu/all_theta_mu_gm_%s.pdf",sample.c_str()));

return 0;
}


