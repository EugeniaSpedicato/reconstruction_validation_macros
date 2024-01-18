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

int p_notReco(){

  int nthreads = 8;
  ROOT::EnableImplicitMT(nthreads);

  //ROOT::TTreeProcessorMT tp("/mnt/raid10/DATA/espedica/fairmu/Mesmer_new_100k_1hit_bend.root","cbmsim",nthreads);
  ROOT::TTreeProcessorMT tp("/mnt/raid10/DATA/espedica/fairmu/Mesmer_new_1M_1hit_bend.root","cbmsim",nthreads);
        MUonERecoOutput *ReconstructionOutput = 0;



   ROOT::TThreadedObject<TH2D> h_2d("h2D","theta mu vs theta E with all cuts", 300,0.,0.035,150,0.,0.005);

TH1::SetDefaultSumw2(kTRUE);
   const Int_t NBINS = 6;
   Double_t edges[NBINS + 1] = {0.0, 0.005, 0.010, 0.015, 0.020, 0.025, 0.032};
   ROOT::TThreadedObject<TH1D> d_eff("d_eff_MC", "Efficiency as a function of the electron's angle",NBINS,edges);
   ROOT::TThreadedObject<TH1D> theta_e("theta_e", "Electron scattering reco angles from MESMER",10,0.,0.035);
   ROOT::TThreadedObject<TH1D> theta_mu("theta_mu", "Muon scattering reco angles from MESMER",20,0.,0.005);

   ROOT::TThreadedObject<TH1D> theta_e_gen("theta_e_gen", "Electron scattering generated angles from MESMER",10,0.,0.035);
   ROOT::TThreadedObject<TH1D> theta_mu_gen("theta_mu_gen", "Muon scattering generated angles from MESMER",20,0.,0.005);
   ROOT::TThreadedObject<TH1D> th_mu_ris("th_mu_ris", "muon resolution Emu75,85GeV",300,-0.0006,0.0006);

   ROOT::TThreadedObject<TH1D> h_res("res", "(the_rec-the_true) 0<theta_e<5 GeV",100,-0.01,0.01);
   ROOT::TThreadedObject<TH1D> h_res1("res1", "(the_rec-the_true) 5<theta_e<10 GeV",50,-0.0025,0.0025);
   ROOT::TThreadedObject<TH1D> h_res2("res2", "(the_rec-the_true) 10<theta_e<15 GeV",40,-0.01,0.01);
   ROOT::TThreadedObject<TH1D> h_res3("res3", "(the_rec-the_true) 15<theta_e<20 GeV",40,-0.01,0.01);
   ROOT::TThreadedObject<TH1D> h_res4("res4", "(the_rec-the_true) 20<theta_e<25 GeV",40,-0.01,0.01);
   ROOT::TThreadedObject<TH1D> h_res5("res5", "(the_rec-the_true) 25<theta_e<32 GeV",40,-0.01,0.01);

   ROOT::TThreadedObject<TH1D> signal("signal", "sum of weights", 2,0,2);

  auto myFunction = [&](TTreeReader &myReader) {
     //TTreeReaderValue<std::vector<MUonERecoOutputTrack>> tracksRV(myReader, "ReconstructionOutput");
     TTreeReaderValue<std::vector<MUonERecoOutputTrack>> tracksRV(myReader, "ReconstructedTracks");
     TTreeReaderArray<MUonETrack> MCtracksRV(myReader, "MCTrack");
     TTreeReaderValue<MUonERecoOutputVertex> vrtx(myReader, "BestVertex");
     TTreeReaderValue<MuE::Event> mesmerEvent(myReader, "MesmerEvent");

     //TTreeReaderValue<std::vector<Int_t>> MCTrack.pdgCode()(myReader, "MCTrack.PdgCode");

     // For performance reasons, a copy of the pointer associated to this thread on the
     // stack is used
     auto my_h_res = h_res.Get();


     while (myReader.Next()) {

int yes2=0; int yes_v=0;
int code_mu=-99; int code_e=-99; int code_mu_in=-99;
int TrackIdreco=-99;
double z_fix=912.7;

        TVector3 p_muin_MC;
        TVector3 p_mu_MC;
        TVector3 p_e_MC;
	double the_gen, thmu_gen;
	double emu=0.;

        Int_t i=0;
        for (auto&& MCTrack : MCtracksRV) {
	 if(MCTrack.interactionID()==0 and MCTrack.pdgCode()==13) {code_mu_in=i; p_muin_MC.SetXYZ(MCTrack.px(),MCTrack.py(),MCTrack.pz()); p_muin_MC.Unit();}
         if(MCTrack.interactionID()==45 and MCTrack.pdgCode()==11) {code_e=i; p_e_MC.SetXYZ(MCTrack.px(),MCTrack.py(),MCTrack.pz()); p_e_MC.Unit(); the_gen=p_muin_MC.Angle(p_e_MC);}
         if(MCTrack.interactionID()==45 and MCTrack.pdgCode()==13) {code_mu=i; p_mu_MC.SetXYZ(MCTrack.px(),MCTrack.py(),MCTrack.pz()); p_mu_MC.Unit(); thmu_gen=p_muin_MC.Angle(p_mu_MC); emu=MCTrack.energy();}
         i++;	}

 if(code_mu_in!=-99 and code_mu!=-99 and code_e!=-99){

double chi=vrtx->chi2perDegreeOfFreedom();

int yes_mu=0;
int yes_e=0;
double th_inx,th_iny,x0_in,y0_in;
double chi2_muin;
double stubs_muin;
int sec0=0; int sec1=0;

         auto tracks = *tracksRV;

         for (auto&& track : tracks) { if(track.sector()==0) sec0++; if(track.sector()==1) sec1++;}

         for (auto&& track : tracks) {

if(sec0==1 and track.sector()==0){
        th_inx=track.xSlope();
        th_iny=track.ySlope();
        x0_in=track.x0();
        y0_in=track.y0();
        chi2_muin=track.chi2perDegreeOfFreedom();
        std::vector<MUonERecoOutputHit> hits_=track.hits();
        stubs_muin=hits_.size();
                        }
if(track.processIDofLinkedTrack()==45 and tracks.size()==3 and track.sector()==1)
{yes2++;
                 if(code_e==track.linkedTrackID()) yes_e++;
                 if(code_mu==track.linkedTrackID()) yes_mu++;
	}
}

double posxIN=pos_on_track(x0_in,th_inx,z_fix);
double posyIN=pos_on_track(y0_in,th_iny,z_fix);

//h_xy->Fill(posxIN,posyIN);

if(stubs_muin==6 and abs(posxIN)<=1.5 and abs(posyIN)<=1.5 and chi2_muin<2){

signal->Fill(1,mesmerEvent->wgt_full);

                                                double dotProduct_MC = p_mu_MC.Dot(p_e_MC);
                                                TVector3 crossProduct_MC = p_mu_MC.Cross(p_e_MC);
                                                double T_MC = p_muin_MC.Dot(crossProduct_MC);
                                                TVector3 im_MC= p_muin_MC.Cross(p_mu_MC);
                                                TVector3 ie_MC= p_muin_MC.Cross(p_e_MC);
                                                T_MC = T_MC>0? 1:-1;
                                                double acoplanarity_MC= T_MC*(TMath::Pi()- acos( ((im_MC).Dot(ie_MC))/(im_MC.Mag()*ie_MC.Mag()) ));

if(abs(acoplanarity_MC)<=1 and thmu_gen>0.0002){
theta_mu_gen->Fill(thmu_gen,mesmerEvent->wgt_full);
theta_e_gen->Fill(the_gen,mesmerEvent->wgt_full);
 }

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

 d_eff->Fill(vrtx->electronTheta(),mesmerEvent->wgt_full);
 h_2d->Fill(vrtx->electronTheta(),vrtx->muonTheta(),mesmerEvent->wgt_full);
theta_mu->Fill(vrtx->muonTheta(),mesmerEvent->wgt_full);
theta_e->Fill(vrtx->electronTheta(),mesmerEvent->wgt_full);

//first
if(vrtx->electronTheta()>0.0 and vrtx->electronTheta()<=0.005){h_res->Fill(vrtx->electronTheta()-the_gen,mesmerEvent->wgt_full);}

//second
if(vrtx->electronTheta()>0.005 and vrtx->electronTheta()<=0.01){h_res1->Fill(vrtx->electronTheta()-the_gen,mesmerEvent->wgt_full);}

//third
if(vrtx->electronTheta()>0.01 and vrtx->electronTheta()<=0.015){h_res2->Fill(vrtx->electronTheta()-the_gen,mesmerEvent->wgt_full);}

//fourth
if(vrtx->electronTheta()>0.015 and vrtx->electronTheta()<=0.02){h_res3->Fill(vrtx->electronTheta()-the_gen,mesmerEvent->wgt_full);}

//fifth
if(vrtx->electronTheta()>0.02 and vrtx->electronTheta()<=0.025){h_res4->Fill(vrtx->electronTheta()-the_gen,mesmerEvent->wgt_full);}

//sixth
if(vrtx->electronTheta()>0.025 and vrtx->electronTheta()<=0.032){h_res5->Fill(vrtx->electronTheta()-the_gen,mesmerEvent->wgt_full);}

			}//aco e chi
		}//chi!=0
	}//chiusura mu_in
    }// chiusura if code_x!0-99
code_e=-99;code_mu=-99;code_mu_in=-99;
yes2=0;yes_v=0;
 }// end while
};//end my function

  tp.Process(myFunction);

  auto h_resMerged   = h_res.Merge();
  auto h_res1Merged   = h_res1.Merge();
  auto h_res2Merged   = h_res2.Merge();
  auto h_res3Merged   = h_res3.Merge();
  auto h_res4Merged   = h_res4.Merge();
  auto h_res5Merged   = h_res5.Merge();

  auto d_effMerged = d_eff.Merge();
  auto h_2dMerged = h_2d.Merge();
  auto theta_muMerged = theta_mu.Merge();
  auto theta_eMerged = theta_e.Merge();
  auto theta_mu_genMerged = theta_mu_gen.Merge();
  auto theta_e_genMerged = theta_e_gen.Merge();
  auto signalM = signal.Merge();

//cout << "sum of weights " << signalM->Integral(0,3) << endl;
cout << "sum of weights " << signalM->Integral() << endl;

TCanvas r("r","r",700,700);
r.Divide(2,3);
r.cd(1);
h_resMerged->Draw("hist");
r.cd(2);
h_res1Merged->Draw("hist");
r.cd(3);
h_res2Merged->Draw("hist");
r.cd(4);
h_res3Merged->Draw("hist");
r.cd(5);
h_res4Merged->Draw("hist");
r.cd(6);
h_res5Merged->Draw("hist");
r.SaveAs("res_bend.pdf");

TCanvas a("a","a",700,700);
d_effMerged->Draw("E");
a.SaveAs("d_eff_MC.pdf");
d_effMerged->SaveAs("d_eff_MC.root");

TCanvas b("b","b",700,700);
h_2dMerged->Draw();
h_2dMerged->SaveAs("2D_MC.root");

TCanvas c("c","c",700,700);
theta_muMerged->Draw("E");
theta_muMerged->SaveAs("theta_mu_MC.root");

TCanvas d("d","d",700,700);
theta_eMerged->Draw("E");
theta_eMerged->SaveAs("theta_e_MC.root");

TCanvas e("e","e",700,700);
theta_mu_genMerged->Draw("E");
theta_mu_genMerged->SaveAs("theta_mu_gen_MC.root");

TCanvas f("f","f",700,700);
theta_e_genMerged->Draw("E");
theta_e_genMerged->SaveAs("theta_e_gen_MC.root");

  return 0;

}
