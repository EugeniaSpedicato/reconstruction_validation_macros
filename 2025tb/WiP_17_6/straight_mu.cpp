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

int straight_mu(string sample, string name){

  int nthreads = 6;
//if(sample=="data") nthreads = 16;

  ROOT::EnableImplicitMT(nthreads);

TChain * cbmsim = new TChain("cbmsim");

if(sample=="data"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run32/passing_muon/muedaq01-1754247289_0hit_WiP_17_6.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run32/passing_muon/muedaq01-1754265666_0hit_WiP_17_6.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run32/passing_muon/muedaq01-1754274334_0hit_WiP_17_6.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run32/passing_muon/muedaq01-1754278680_0hit_WiP_17_6.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run32/passing_muon/muedaq01-1754283088_0hit_WiP_17_6.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run32/passing_muon/muedaq02-1754252562_0hit_WiP_17_6.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run32/passing_muon/muedaq02-1754256880_0hit_WiP_17_6.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run32/passing_muon/muedaq02-1754270029_0hit_WiP_17_6.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run32/passing_muon/muedaq03-1754261271_0hit_WiP_17_6.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/run32/passing_muon/muedaq04-1754242824_0hit_WiP_17_6.root");
}
else if(sample=="MC_ideal"){
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/MC/minBias.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/MC/minimumbias_ideal_passingMuons_noECAL_WiP_v1_1_x.root");
}
else if(sample=="MC_realistic"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/TB2025/MC/minimumbias_realistic_passingMuons_noECAL_WiP_v1_1_x.root");
}

ROOT::TTreeProcessorMT tp1(*cbmsim,nthreads);

      MUonERecoOutputAnalysis *ReconstructionOutput = 0;


TH1::SetDefaultSumw2(kTRUE);

   ROOT::TThreadedObject<TH1D> theta_mu3("theta_mu3", "Muon scattering reco angles third station",100,0.,0.005);//50,0.,0.005);
   ROOT::TThreadedObject<TH1D> theta_mu2("theta_mu2", "Muon scattering reco angles second station",100,0.,0.005);//50,0.,0.005);
   ROOT::TThreadedObject<TH1D> theta_mu("theta_mu", "Muon scattering reco angles first station",100,0.,0.005);//50,0.,0.005);
   ROOT::TThreadedObject<TH1D> thetaX_mu3("thetaX_mu3", "Muon scattering reco X angles third station",300,-0.0015,0.0015);//50,0.,0.005);
   ROOT::TThreadedObject<TH1D> thetaX_mu2("thetaX_mu2", "Muon scattering reco X angles second station",300,-0.0015,0.0015);//50,0.,0.005);
   ROOT::TThreadedObject<TH1D> thetaX_mu("thetaX_mu", "Muon scattering reco X angles first station",300,-0.0015,0.0015);//50,0.,0.005);
   ROOT::TThreadedObject<TH1D> thetaY_mu3("thetaY_mu3", "Muon scattering reco Y angles third station",300,-0.0015,0.0015);//50,0.,0.005);
   ROOT::TThreadedObject<TH1D> thetaY_mu2("thetaY_mu2", "Muon scattering reco Y angles second station",300,-0.0015,0.0015);//50,0.,0.005);
   ROOT::TThreadedObject<TH1D> thetaY_mu("thetaY_mu", "Muon scattering reco Y angles first station",300,-0.0015,0.0015);//50,0.,0.005);
   ROOT::TThreadedObject<TH1D> diff_theta("diff_theta", "Muon scattering reco angles first station - second station",200,-0.00025,0.00025);//50,0.,0.005);
   ROOT::TThreadedObject<TH1D> diff_theta2("diff_theta2", "Muon scattering reco angles second station - third station",200,-0.00025,0.00025);//50,0.,0.005);
   ROOT::TThreadedObject<TH1D> diff_theta3("diff_theta3", "Muon scattering reco angles first station - third station",200,-0.00025,0.00025);//50,0.,0.005);
   ROOT::TThreadedObject<TH1D> diff_thetax("diff_thetax", "Muon scattering reco angles first station - second station THETA_X",200,-0.00025,0.00025);//50,0.,0.005);
   ROOT::TThreadedObject<TH1D> diff_thetay("diff_thetay", "Muon scattering reco angles first station - second station THETA_Y",200,-0.00025,0.00025);//50,0.,0.005);
   ROOT::TThreadedObject<TH1D> diff_thetax2("diff_thetax2", "Muon scattering reco angles second station - third station THETA_X",200,-0.00025,0.00025);//50,0.,0.005);
   ROOT::TThreadedObject<TH1D> diff_thetay2("diff_thetay2", "Muon scattering reco angles second station - third station THETA_Y",200,-0.00025,0.00025);//50,0.,0.005);
   ROOT::TThreadedObject<TH1D> diff_thetax3("diff_thetax3", "Muon scattering reco angles first station - third station THETA_X",200,-0.00025,0.00025);//50,0.,0.005);
   ROOT::TThreadedObject<TH1D> diff_thetay3("diff_thetay3", "Muon scattering reco angles first station - third station THETA_Y",200,-0.00025,0.00025);//50,0.,0.005);
   ROOT::TThreadedObject<TH1D> diff_theta_smeared("diff_theta_smeared", "Muon scattering reco angles first station - second station smeared",200,-0.00025,0.00025);//50,0.,0.005);
   ROOT::TThreadedObject<TH1D> diff_theta_smeared2("diff_theta_smeared2", "Muon scattering reco angles second station - third station smeared",200,-0.00025,0.00025);//50,0.,0.005);
   ROOT::TThreadedObject<TH1D> ms_effect("ms_effect","Effect ms",200,0.018,0.038);
   ROOT::TThreadedObject<TH1D> ms_effect2("ms_effect2","Effect ms",200,0.018,0.038);
   ROOT::TThreadedObject<TH2D> h_2D_thetaY_12("h_2D_thetaY_12", "Muon scattering reco Y angles 1st and 2nd station",300,-0.0015,0.0015,300,-0.0015,0.0015);//50,0.,0.005);
   ROOT::TThreadedObject<TH2D> h_2D_thetaY_23("h_2D_thetaY_23", "Muon scattering reco Y angles 2nd and 3rd station",300,-0.0015,0.0015,300,-0.0015,0.0015);//50,0.,0.005);
   ROOT::TThreadedObject<TH2D> h_2D_thetaY_13("h_2D_thetaY_13", "Muon scattering reco Y angles 1st and 3rd station",300,-0.0015,0.0015,300,-0.0015,0.0015);//50,0.,0.005);
   ROOT::TThreadedObject<TH1D> h_IPx_t0("h_IPx_t0", "Distance X tracks at tar 0 in m",500,-0.0005,0.0005);
   ROOT::TThreadedObject<TH1D> h_IPy_t0("h_IPy_t0", "Distance Y tracks at tar 0 in m",500,-0.0005,0.0005);
   ROOT::TThreadedObject<TH1D> h_IPx_t1("h_IPx_t1", "Distance X tracks at tar 1 in m",500,-0.0005,0.0005);
   ROOT::TThreadedObject<TH1D> h_IPy_t1("h_IPy_t1", "Distance Y tracks at tar 1 in m",500,-0.0005,0.0005);
   ROOT::TThreadedObject<TH1D> h_IPx_t0_02("h_IPx_t0_02", "Distance X tracks at tar 0 tracks 0 and 2 in m",500,-0.0005,0.0005);
   ROOT::TThreadedObject<TH1D> h_IPy_t0_02("h_IPy_t0_02", "Distance Y tracks at tar 0 tracks 0 and 2 in m",500,-0.0005,0.0005);
   ROOT::TThreadedObject<TH1D> h_energy("h_energy", "calo energy in passing muons sample",400,0.,200.);

std::cout << "entries " << cbmsim->GetEntries() << std::endl;


 auto myFunction = [&](TTreeReader &myReader) {

     TTreeReaderValue<std::vector<MUonERecoOutputTrackAnalysis>> RVtracks(myReader, "ReconstructedTracks");
     TTreeReaderValue<std::vector<MUonERecoOutputHitAnalysis>> RVstubs(myReader, "ReconstructedHits");
     TTreeReaderValue<Double_t> RVenergy(myReader, "ReconstructedCalorimeterCluster.ClusterEnergy");


     while (myReader.Next()) {

Long64_t entry = myReader.GetCurrentEntry();

double z_fix=912.7;


double th_inx,th_iny,x0_in,y0_in;
double th_inx2,th_iny2,x0_in2,y0_in2;
double th_inx3,th_iny3,x0_in3,y0_in3;

double th_inx_smeared,th_iny_smeared;
double th_inx2_smeared,th_iny2_smeared;
double th_inx3_smeared,th_iny3_smeared;


double chi2_muin;
double chi2_muin2;
double chi2_muin3;
int sec0=0;
int sec1=0;
int sec2=0;
int stubs_muin=0.;
double th_muin=0.;
int stubs_muin2=0.;
double th_muin2=0.;
int stubs_muin3=0.;
double th_muin3=0.;

double th_muin_smeared=0.;
double th_muin2_smeared=0.;
double th_muin3_smeared=0.;

TVector3 p_muin,p_muin2,p_muin3;
TVector3 p_muin_smeared,p_muin2_smeared,p_muin3_smeared;

         auto tracks = *RVtracks;
         for (auto&& track : tracks) {
        if(track.sector()==0) sec0++;
        else if(track.sector()==1) sec1++;
        else if(track.sector()==2) sec2++;
        }

double posxIN=99.;//pos_on_track(x0_in,th_inx,z_fix);
double posyIN=99.;//pos_on_track(y0_in,th_iny,z_fix);
double posxIN2=99.;//pos_on_track(x0_in,th_inx,z_fix);
double posyIN2=99.;//pos_on_track(y0_in,th_iny,z_fix);
double posxIN3=99.;//pos_on_track(x0_in,th_inx,z_fix);
double posyIN3=99.;//pos_on_track(y0_in,th_iny,z_fix);

int stub0 = 0;
int stub1 = 0;
int stub2 = 0;

         auto stubs = *RVstubs;
         std::array<int,6> module_st1_3;module_st1_3.fill({0});
         std::array<int,6> module_st1_2;module_st1_2.fill({0});
         std::array<int,6> module_st1;module_st1.fill({0});
         std::array<int,6> module_st2;module_st2.fill({0});
         std::array<std::vector<float>,6> localX{{std::vector<float>(0.0)}};
         std::array<std::vector<float>,6> localX1{{std::vector<float>(0.0)}};
         for (auto&& stub : stubs) {
                if(stub.stationID()==0){stub0++; if(stub.moduleID()==4){posxIN=stub.position(); } else if(stub.moduleID()==5){posyIN=stub.position();}   }
                else if(stub.stationID()==1){stub1++; if(stub.moduleID()==4){posxIN2=stub.position(); } else if(stub.moduleID()==5){posyIN2=stub.position();}
						      module_st1_2.at(stub.moduleID())+=1; module_st1.at(stub.moduleID())=1;int link= stub.moduleID(); localX.at(link).push_back(stub.position());}
                else if(stub.stationID()==2){stub2++; if(stub.moduleID()==4){posxIN3=stub.position(); } else if(stub.moduleID()==5){posyIN3=stub.position();}
						      module_st1_3.at(stub.moduleID())+=1; module_st2.at(stub.moduleID())=1;int link= stub.moduleID(); localX1.at(link).push_back(stub.position());}

        }



MUonERecoOutputTrackAnalysis t_mu;
MUonERecoOutputTrackAnalysis t_mu2;
MUonERecoOutputTrackAnalysis t_mu3;

double smearing = 0;//8.e-06;// 1.7e-05;
double smearing1 = 0.;// 1.7e-05;
double smearing2 = 0.;// 1.7e-05;


         for (auto&& track : tracks) {
        std::vector<Short_t> hits_=track.hitIds();
        if(track.sector()==0 and sec0==1){
	stubs_muin=hits_.size();
        th_inx=track.xSlope();
        th_iny=track.ySlope();
        th_inx_smeared=track.xSlope()+gRandom->Gaus(0.,smearing);//gRandom->Gaus(0.,2.3e-05);
        th_iny_smeared=track.ySlope()+gRandom->Gaus(0.,smearing);
        x0_in=track.x0();
        y0_in=track.y0();
        chi2_muin=track.chi2();
        p_muin.SetXYZ(th_inx,th_iny,1.0);
        p_muin=p_muin.Unit();
        th_muin=p_muin.Theta();
        p_muin_smeared.SetXYZ(th_inx_smeared,th_iny_smeared,1.0);
        p_muin_smeared=p_muin_smeared.Unit();
        th_muin_smeared=p_muin_smeared.Theta();
	t_mu=track;
                        }
	else if(track.sector()==1 and sec1==1)// and stub1<8)
	{
        stubs_muin2=hits_.size();
        th_inx2=track.xSlope();
        th_iny2=track.ySlope();//+gRandom->Gaus(0.,10.e-06);
        th_inx2_smeared=track.xSlope()+gRandom->Gaus(0.,smearing1);//notar:1.6e-05, tar:2.3e-05
        th_iny2_smeared=track.ySlope()+gRandom->Gaus(0.,smearing1);
        x0_in2=track.x0();
        y0_in2=track.y0();
        chi2_muin2=track.chi2();
        p_muin2.SetXYZ(th_inx2,th_iny2,1.0);
        p_muin2=p_muin2.Unit();
        th_muin2=p_muin2.Theta();
        p_muin2_smeared.SetXYZ(th_inx2_smeared,th_iny2_smeared,1.0);
        p_muin2_smeared=p_muin2_smeared.Unit();
        th_muin2_smeared=p_muin2_smeared.Theta();
       t_mu2=track;
	}
	else if(track.sector()==2 and sec2==1)// and stub1<8)
	{
        stubs_muin3=hits_.size();
        th_inx3=track.xSlope();
        th_iny3=track.ySlope();
        th_inx3_smeared=track.xSlope()+gRandom->Gaus(0.,smearing2);//notar:1.6e-05, tar:2.3e-05
        th_iny3_smeared=track.ySlope()+gRandom->Gaus(0.,smearing2);
        x0_in3=track.x0();
        y0_in3=track.y0();
        chi2_muin3=track.chi2();
        p_muin3.SetXYZ(th_inx3,th_iny3,1.0);
        p_muin3=p_muin3.Unit();
        th_muin3=p_muin3.Theta();
        p_muin3_smeared.SetXYZ(th_inx3_smeared,th_iny3_smeared,1.0);
        p_muin3_smeared=p_muin3_smeared.Unit();
        th_muin3_smeared=p_muin3_smeared.Theta();
       t_mu3=track;
	}
 }

/*std::cout << "sec0 " << sec0 << std::endl;
std::cout << "sec1 " << sec1 << std::endl;
std::cout << "sec2 " << sec2 << std::endl;
std::cout << "stub0 " << stub0 << std::endl;
std::cout << "stub1 " << stub1 << std::endl;
std::cout << "stub2 " << stub2 << std::endl;*/

	auto energy=*RVenergy;

//if( abs(posxIN)<=3 and abs(posyIN)<=3 and abs(posxIN2)<=3 and abs(posyIN2)<=3 and abs(posxIN3)<=3 and abs(posyIN3)<=3 and
//if( ( (abs(posxIN)>3 and (posxIN2)>3 and abs(posxIN3)>3) or (abs(posyIN)>3 and abs(posyIN2)>3 and abs(posyIN3)>3) ) and
if(
	stub0==6 and stub1==6 and stub2<=15 and
	sec0==1 and sec1==1 //and sec2==2
//	and abs(th_inx)>0.0002 and abs(th_inx2)>0.0002 and abs(th_inx3)>0.0002
//	and abs(th_iny)>0.0002 and abs(th_iny2)>0.0002 and abs(th_iny3)>0.0002
	//and abs(th_inx)<0.0001 and abs(th_inx2)<0.0001 and abs(th_inx3)<0.0001
        //and abs(th_iny)<0.0001 and abs(th_iny2)<0.0001 and abs(th_iny3)<0.0001
	){

	if(energy>30)h_energy->Fill(energy);


//	h_energy->Fill(energy);

                double IPx_t0=pos_on_track(t_mu.x0(),t_mu.xSlope(),(667.3-2.7))  - pos_on_track(t_mu2.x0(),t_mu2.xSlope(),(667.3-2.7));
                double IPy_t0=pos_on_track(t_mu.y0(),t_mu.ySlope(),(667.3-2.7))  - pos_on_track(t_mu2.y0(),t_mu2.ySlope(),(667.3-2.7));

                double IPx_t1=pos_on_track(t_mu2.x0(),t_mu2.xSlope(),(784.6-3.9))  - pos_on_track(t_mu3.x0(),t_mu3.xSlope(),(784.6-3.9));
                double IPy_t1=pos_on_track(t_mu2.y0(),t_mu2.ySlope(),(784.6-3.9))  - pos_on_track(t_mu3.y0(),t_mu3.ySlope(),(784.6-3.9));

		double IPx_t0_02=pos_on_track(t_mu.x0(),t_mu.xSlope(),(667.3-2.7))  - pos_on_track(t_mu3.x0(),t_mu3.xSlope(),(667.3-2.7));
                double IPy_t0_02=pos_on_track(t_mu.y0(),t_mu.ySlope(),(667.3-2.7))  - pos_on_track(t_mu3.y0(),t_mu3.ySlope(),(667.3-2.7));

//if( (th_iny-th_iny2)>-5.3e-5 and (th_iny-th_iny2)<-4.0e-5 )
//if( (th_iny2-th_iny3)>-4.0e-5 and (th_iny2-th_iny3)<-1.5e-5 )std::cout<<"entry --? " << entry << " ?--" << std::endl;

//	if(abs(posxIN)>3 and (posxIN2)>3 and abs(posxIN3)>3){
        thetaX_mu3->Fill(th_inx3);
        thetaX_mu2->Fill(th_inx2);
        thetaX_mu->Fill(th_inx);
	diff_thetax->Fill(th_inx-th_inx2);
	diff_thetax2->Fill(th_inx2-th_inx3);
	diff_thetax3->Fill(th_inx-th_inx3);
//	}
//	else if(abs(posyIN)>3 and abs(posyIN2)>3 and abs(posyIN3)>3){
        thetaY_mu3->Fill(th_iny3);
        thetaY_mu2->Fill(th_iny2);
        thetaY_mu->Fill(th_iny);
	diff_thetay->Fill(th_iny-th_iny2);
	diff_thetay2->Fill(th_iny2-th_iny3);
	diff_thetay3->Fill(th_iny-th_iny3);
//	}

        theta_mu->Fill(th_muin);
        theta_mu2->Fill(th_muin2);
        theta_mu3->Fill(th_muin3);
	diff_theta->Fill(th_muin-th_muin2);
	diff_theta_smeared->Fill(th_muin_smeared-th_muin2_smeared);

	diff_theta2->Fill(th_muin2-th_muin3);
	diff_theta_smeared2->Fill(th_muin2_smeared-th_muin3_smeared);

	diff_theta3->Fill(th_muin-th_muin3);

	h_2D_thetaY_12->Fill(th_iny,th_iny2);
	h_2D_thetaY_23->Fill(th_iny2,th_iny3);
	h_2D_thetaY_13->Fill(th_iny,th_iny3);

h_IPx_t0_02->Fill(IPx_t0_02*0.01);
h_IPy_t0_02->Fill(IPy_t0_02*0.01);
h_IPx_t0->Fill(IPx_t0*0.01);
h_IPy_t0->Fill(IPy_t0*0.01);
h_IPx_t1->Fill(IPx_t1*0.01);
h_IPy_t1->Fill(IPy_t1*0.01);

	}//mu_in


  } //end of general while
}; //end of ySlopefunction

  tp1.Process(myFunction);

  auto theta_muMerged = theta_mu.Merge();
  auto theta_mu2Merged = theta_mu2.Merge();
  auto diff_thetaM = diff_theta.Merge();
  auto diff_theta_smearedM=diff_theta_smeared.Merge();
  auto diff_thetaxM=diff_thetax.Merge();
  auto diff_thetayM=diff_thetay.Merge();

  auto theta_mu3Merged = theta_mu3.Merge();
  auto diff_theta_smeared2M=diff_theta_smeared2.Merge();
  auto diff_theta2Merged = diff_theta2.Merge();
  auto diff_thetax2M=diff_thetax2.Merge();
  auto diff_thetay2M=diff_thetay2.Merge();


  auto diff_theta3Merged = diff_theta3.Merge();
  auto diff_thetax3M=diff_thetax3.Merge();
  auto diff_thetay3M=diff_thetay3.Merge();

  auto thetaX_mu3Merged = thetaX_mu3.Merge();
  auto thetaY_mu3Merged = thetaY_mu3.Merge();
  auto thetaX_mu2Merged = thetaX_mu2.Merge();
  auto thetaY_mu2Merged = thetaY_mu2.Merge();
  auto thetaX_muMerged = thetaX_mu.Merge();
  auto thetaY_muMerged = thetaY_mu.Merge();

  auto h_2D_thetaY_12M=h_2D_thetaY_12.Merge();
  auto h_2D_thetaY_23M=h_2D_thetaY_23.Merge();
  auto h_2D_thetaY_13M=h_2D_thetaY_13.Merge();

  auto h_IPx_t0M=h_IPx_t0.Merge();
  auto h_IPy_t0M=h_IPy_t0.Merge();
  auto h_IPx_t1M=h_IPx_t1.Merge();
  auto h_IPy_t1M=h_IPy_t1.Merge();
  auto h_IPx_t0_02M=h_IPx_t0_02.Merge();
  auto h_IPy_t0_02M=h_IPy_t0_02.Merge();

theta_muMerged->SaveAs(Form("dir_straight_mu/theta_mu_gm_%s.root",sample.c_str()));
theta_mu2Merged->SaveAs(Form("dir_straight_mu/theta_mu2_gm_%s.root",sample.c_str()));
diff_thetaM->SaveAs(Form("dir_straight_mu/diff_theta_mu_gm_%s.root",sample.c_str()));

theta_mu3Merged->SaveAs(Form("dir_straight_mu/theta_mu3_gm_%s.root",sample.c_str()));
diff_theta2Merged->SaveAs(Form("dir_straight_mu/diff_theta_mu2_gm_%s.root",sample.c_str()));

diff_theta_smearedM->SaveAs(Form("dir_straight_mu/diff_theta_smeared_mu_gm_%s.root",sample.c_str()));
diff_theta_smeared2M->SaveAs(Form("dir_straight_mu/diff_theta_smeared_mu2_gm_%s.root",sample.c_str()));

TF1 *g1 = new TF1("g1", "gaus");
TF1 *g2 = new TF1("g2", "gaus");
TF1 *g3 = new TF1("g3", "gaus");
TF1 *g4 = new TF1("g4", "gaus");

diff_thetaM->Fit("g1","R","",-0.00015,0.00015);
diff_theta2Merged->Fit("g3","R","",-0.00015,0.00015);

diff_theta_smearedM->Fit("g2","R","",-0.00015,0.00015);
diff_theta_smeared2M->Fit("g4","R","",-0.00015,0.00015);

cout << "sigma: " << g1->GetParameter(2) << endl;
cout << "sigma2: " << g3->GetParameter(2) << endl;
cout << "sigma smeared: " << g2->GetParameter(2) << endl;
cout << "sigma2 smeared: " << g4->GetParameter(2) << endl;

TCanvas a("a","a",1500,1500);
a.Divide(2,3);
a.cd(1);
diff_thetax2M->SetLineColor(kRed);
diff_thetax3M->SetLineColor(kViolet);
diff_thetax2M->Draw("hist");
diff_thetax3M->Draw("hist sames");
diff_thetaxM->Draw("hist sames");
a.cd(2);
diff_thetay2M->SetLineColor(kRed);
diff_thetay3M->SetLineColor(kViolet);
diff_thetay2M->Draw("hist");
//diff_thetay3M->Draw("hist sames");
diff_thetayM->Draw("hist sames");
a.cd(3);
diff_thetaM->Draw("hist");
diff_theta2Merged->SetLineColor(kRed);
diff_theta2Merged->Draw("hist sames");
//diff_thetaM->SetLineColor(kGreen);
diff_thetaM->Draw("hist sames");
//diff_theta3Merged->SetLineColor(kViolet);
//diff_theta3Merged->Draw("hist sames");
a.cd(4);
thetaX_muMerged->Draw("hist");
thetaX_mu2Merged->SetLineColor(kRed);
thetaX_mu2Merged->Draw("hist sames");
thetaX_mu3Merged->SetLineColor(kViolet);
//thetaX_mu3Merged->Draw("hist sames");
a.cd(5);
thetaY_muMerged->Draw("hist");
thetaY_mu2Merged->SetLineColor(kRed);
thetaY_mu2Merged->Draw("hist sames");
thetaY_mu3Merged->SetLineColor(kViolet);
//thetaY_mu3Merged->Draw("hist sames");

a.SaveAs(Form("dir_straight_mu/all_stations_gm_%s_%s.root",sample.c_str(),name.c_str()));

TCanvas b("b","b",100,300);
b.Divide(1,3);
b.cd(1);
h_2D_thetaY_12M->Draw("COLZ");
b.cd(2);
h_2D_thetaY_23M->Draw("COLZ");
b.cd(3);
h_2D_thetaY_13M->Draw("COLZ");
b.SaveAs(Form("dir_straight_mu/2D_gm_%s_%s.root",sample.c_str(),name.c_str()));

TCanvas c("c","c",200,300);
c.Divide(2,3);
c.cd(1);
h_IPx_t0M->Draw("hist");
c.cd(2);
h_IPy_t0M->Draw("hist");
c.cd(3);
h_IPx_t1M->Draw("hist");
c.cd(4);
h_IPy_t1M->Draw("hist");
c.cd(5);
h_IPx_t0_02M->Draw("hist");
c.cd(6);
h_IPy_t0_02M->Draw("hist");

c.SaveAs(Form("dir_straight_mu/impactPoint_%s_%s.root",sample.c_str(),name.c_str()));

TCanvas e("e","e",500,500);
h_energy->Draw("hist");
e.SaveAs(Form("dir_straight_mu/energy_%s_%s.root",sample.c_str(),name.c_str()));


/*
TCanvas a("a","a",1500,1000);
a.Divide(2,3);
a.cd(1);
theta_muMerged->Draw("hist");
a.cd(2);
theta_mu2Merged->Draw("hist");
a.cd(3);
diff_thetaM->Draw("hist");
diff_theta_smearedM->SetLineColor(kGreen);
diff_theta_smearedM->Draw("hist sames");
diff_theta_smeared2M->SetLineColor(kRed);
diff_theta_smeared2M->Draw("hist sames");
a.cd(4);
ms_effect->Draw("hist");
a.cd(5);
diff_thetaxM->Draw("hist");
a.cd(6);
diff_thetayM->Draw("hist");
a.SaveAs(Form("dir_straight_mu/all_theta_mu_gm_%s.pdf",sample.c_str()));



TCanvas a2("a2","a2",1500,1000);
a2.Divide(2,3);
a2.cd(1);
theta_mu2Merged->Draw("hist");
a2.cd(2);
theta_mu3Merged->Draw("hist");
a2.cd(3);
diff_theta2Merged->Draw("hist");
diff_theta_smearedM->SetLineColor(kGreen);
diff_theta_smearedM->Draw("hist sames");
diff_theta_smeared2M->SetLineColor(kRed);
diff_theta_smeared2M->Draw("hist sames");
a2.cd(4);
ms_effect2->Draw("hist");
a2.cd(5);
diff_thetax2M->Draw("hist");
a2.cd(6);
diff_thetay2M->Draw("hist");
a2.SaveAs(Form("dir_straight_mu/all_theta_mu_gm2_%s.pdf",sample.c_str()));
*/
return 0;
}


