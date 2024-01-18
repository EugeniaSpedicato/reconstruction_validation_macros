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

int p_dati_noReco(){

  int nthreads = 8;
  ROOT::EnableImplicitMT(nthreads);

//  ROOT::TTreeProcessorMT tp1("/mnt/raid10/DATA/espedica/fairmu/dataReconstruction_3234-3235_new_1hit_skimGA.root","cbmsim",nthreads);
//  ROOT::TTreeProcessorMT tp2("/mnt/raid10/DATA/espedica/fairmu/dataReconstruction_3234-3235_new_1hit_skimGA_1.root","cbmsim",nthreads);
    ROOT::TTreeProcessorMT tp1("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_3234-3235_new_1hit_skimGA_inputs_12.root","cbmsim",nthreads);

      MUonERecoOutput *ReconstructionOutput = 0;


   ROOT::TThreadedObject<TH2D> h_2d("h2D","theta mu vs theta E with all cuts", 300,0.,0.035,150,0.,0.005);

TH1::SetDefaultSumw2(kTRUE);
   const Int_t NBINS = 6;
   Double_t edges[NBINS + 1] = {0.0, 0.005, 0.010, 0.015, 0.020, 0.025, 0.032};
   ROOT::TThreadedObject<TH1D> d_eff("d_eff_MC", "Efficiency as a function of the electron's angle",NBINS,edges);
   ROOT::TThreadedObject<TH1D> theta_e("theta_e", "Electron scattering reco angles from MESMER",10,0.,0.035);
   ROOT::TThreadedObject<TH1D> theta_mu("theta_mu", "Muon scattering reco angles from MESMER",20,0.,0.005);

   ROOT::TThreadedObject<TH1D> signal("signal", "sum of weights", 2,0,2);

 auto myFunction = [&](TTreeReader &myReader) {

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

         auto tracks = *RVtracks;
         for (auto&& track : tracks) {
        if(track.sector()==0) sec0++;
        if(track.sector()==1) sec1++;
        }

         for (auto&& track : tracks) {

        if(track.sector()==0 and sec0==1){
	std::vector<MUonERecoOutputHit> hits_=track.hits();
	stubs_muin=hits_.size();
        th_inx=track.xSlope();
        th_iny=track.ySlope();
        x0_in=track.x0();
        y0_in=track.y0();
        chi2_muin=track.chi2perDegreeOfFreedom();
                        }
	if(track.sector()==1 and sec1==2)
	{yes2++;}
}

double posxIN=pos_on_track(x0_in,th_inx,z_fix);
double posyIN=pos_on_track(y0_in,th_iny,z_fix);

int stub0 = 0;
         auto stubs = *RVstubs;
         for (auto&& stub : stubs) {if(stub.stationID()==0)stub0++;}

 if(sec0==1 and stubs_muin==6 and abs(posxIN)<=1.5 and abs(posyIN)<=1.5 and chi2_muin<2 and stub0==6){
signal->Fill(1);

if(chi!=0){
 MUonERecoOutputTrack mu_in = vrtx->incomingMuon();
 MUonERecoOutputTrack mu_out = vrtx->outgoingMuon();
 MUonERecoOutputTrack e_out = vrtx->outgoingElectron();
        TVector3 p_muin(mu_in.xSlope(),mu_in.ySlope(),1.0);
        TVector3 p_mu(mu_out.xSlope(),mu_out.ySlope(),1.0);
        TVector3 p_e(e_out.xSlope(),e_out.ySlope(),1.0);

                                                double dotProduct_v = p_mu.Dot(p_e);
                                                TVector3 crossProduct_v = p_mu.Cross(p_e);
                                                double T_v = p_muin.Dot(crossProduct_v);
                                                TVector3 im_v= p_muin.Cross(p_mu);
                                                TVector3 ie_v= p_muin.Cross(p_e);
                                                T_v = T_v>0? 1:-1;
                                                double acoplanarity_v= T_v*(TMath::Pi()- acos( ((im_v).Dot(ie_v))/(im_v.Mag()*ie_v.Mag()) ));

  if(abs(acoplanarity_v)<=1 and chi<20 and vrtx->muonTheta()>0.0002){
//first
 if(yes2>=2){d_eff->Fill(vrtx->electronTheta());

 h_2d->Fill(vrtx->electronTheta(),vrtx->muonTheta());
 theta_mu->Fill(vrtx->muonTheta());
 theta_e->Fill(vrtx->electronTheta());}

			}//aco e chi
		}//chi!=0
	}//mu_in
yes2=0;
 } //end of general while
}; //end of myfunction

  tp1.Process(myFunction);
  //tp2.Process(myFunction);

  auto d_effMerged = d_eff.Merge();
  auto h_2dMerged = h_2d.Merge();
  auto theta_muMerged = theta_mu.Merge();
  auto theta_eMerged = theta_e.Merge();
  auto signalM = signal.Merge();

//cout << "sum of weights " << signalM->Integral(0,3) << endl;
cout << "sum of weights " << signalM->Integral() << endl;


TCanvas a("a","a",700,700);
d_effMerged->Draw("E");
a.SaveAs("d_eff_RD.pdf");
d_effMerged->SaveAs("d_eff_RD_inputs_12.root");

TCanvas b("b","b",700,700);
h_2dMerged->Draw();
h_2dMerged->SaveAs("2D_RD_inputs_12.root");

TCanvas c("c","c",700,700);
theta_muMerged->Draw("E");
theta_muMerged->SaveAs("theta_mu_RD_inputs_12.root");

TCanvas d("d","d",700,700);
theta_eMerged->Draw("E");
theta_eMerged->SaveAs("theta_e_RD_inputs_12.root");

//signalM->SaveAs("signal_RD_input0.root");


return 0;
}
