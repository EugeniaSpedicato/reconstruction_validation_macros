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
#include <ROOT/TTreeProcessorMT.hxx>

using namespace std;

auto pos_on_track = [](double q, double m, double z){return (q + m*z);};

int IDparallel()
{
  // First enable implicit multi-threading globally, so that the implicit parallelisation is on.
  // The parameter of the call specifies the number of threads to use.
  int nthreads = 4;
  ROOT::EnableImplicitMT(nthreads);
// Create one ROOT::TThreadedObject per histogram to fill during the processing of the tree

  ROOT::TThreadedObject<TH1D> h_res("res", "(the_rec-the_true) 0<theta_e<5 GeV",100,-0.01,0.01);
  ROOT::TThreadedObject<TH1D> h_res1("res1", "(the_rec-the_true) 5<theta_e<10 GeV",100,-0.01,0.01);
  ROOT::TThreadedObject<TH1D> h_res2("res2", "(the_rec-the_true) 10<theta_e<15 GeV",250,-0.04,0.04);
  ROOT::TThreadedObject<TH1D> h_res3("res3", "(the_rec-the_true) 15<theta_e<20 GeV",250,-0.04,0.04);
  ROOT::TThreadedObject<TH1D> h_res4("res4", "(the_rec-the_true) 20<theta_e<25 GeV",250,-0.04,0.04);
  ROOT::TThreadedObject<TH1D> h_res5("res5", "(the_rec-the_true) 25<theta_e<32 GeV",250,-0.04,0.04);

ROOT::TThreadedObject<TH1D> thetaX("theta1" , "theta X of the muon for events in the tail of angular residuum theta<5mrad" , 100,-0.01,0.01);
ROOT::TThreadedObject<TH1D> thetaY("theta2" , "theta Y of the muon for events in the tail of angular residuum theta<5mrad" , 100,-0.01,0.01);
ROOT::TThreadedObject<TH1D> thetaXg("theta1g" , "theta X of the muon for events in the tail of angular residuum theta<5mrad" , 100,-0.01,0.01);
ROOT::TThreadedObject<TH1D> thetaYg("theta2g" , "theta Y of the muon for events in the tail of angular residuum theta<5mrad" , 100,-0.01,0.01);

ROOT::TThreadedObject<TH1D> thetaX2("theta1_large" , "theta X of the muon for events in the tail of angular residuum theta>5mrad" , 150,-0.1,0.1);
ROOT::TThreadedObject<TH1D> thetaY2("theta2_large" , "theta Y of the muon for events in the tail of angular residuum theta>5mrad" , 150,-0.1,0.1);
ROOT::TThreadedObject<TH1D> thetaXg2("theta1g_large" , "theta X of the muon for events in the tail of angular residuum theta>5mrad" , 150,-0.1,0.1);
ROOT::TThreadedObject<TH1D> thetaYg2("theta2g_large" , "theta Y of the muon for events in the tail of angular residuum theta>5mrad" , 150,-0.1,0.1);


  // Create a TTreeProcessor: specify the file and the tree in it
 ROOT::TTreeProcessorMT tp("/mnt/raid10/DATA/espedica/fairmu/Mesmer_sample_1M.root","cbmsim",nthreads);
  //ROOT::TTreeProcessorMT tp("/mnt/raid10/DATA/espedica/fairmu/Mesmer_sample_noBend_1hit_1M.root","cbmsim",nthreads);
        MUonERecoOutput *ReconstructionOutput = 0;

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

Double_t code_mu_in=-99;
Double_t code_e=-99;
Double_t code_mu=-99;
        TVector3 p_muin_MC;
        TVector3 p_mu_MC;
        TVector3 p_e_MC;
        Double_t the_gen, thmu_gen,theX_gen,theY_gen,thmuX_gen,thmuY_gen;

         //auto MCTracks = MCtracksRV;
	Int_t i=0;
        for (auto&& MCTrack : MCtracksRV) {
         if(MCTrack.interactionID()==0 and MCTrack.pdgCode()==13) {code_mu_in=i; p_muin_MC.SetXYZ(MCTrack.px(),MCTrack.py(),MCTrack.pz()); p_muin_MC.Unit();}
         if(MCTrack.interactionID()==45 and MCTrack.pdgCode()==11) {code_e=i; p_e_MC.SetXYZ(MCTrack.px(),MCTrack.py(),MCTrack.pz()); p_e_MC.Unit(); the_gen=p_muin_MC.Angle(p_e_MC);
									         theX_gen=MCTrack.ax();
									         theY_gen=MCTrack.ay();}
         if(MCTrack.interactionID()==45 and MCTrack.pdgCode()==-13) {code_mu=i; p_mu_MC.SetXYZ(MCTrack.px(),MCTrack.py(),MCTrack.pz()); p_mu_MC.Unit(); thmu_gen=p_muin_MC.Angle(p_mu_MC);
                                                                                 thmuX_gen=MCTrack.ax();
                                                                                 thmuY_gen=MCTrack.ay();}
	 i++;
        }

	if(code_mu_in!=-99 and code_e!=-99 and code_mu!=-99){

Int_t yes_mu=0;
Int_t yes_e=0;
Int_t yes2=0;
Double_t th_inx,th_iny,x0_in,y0_in;
Double_t chi2_muin;
Double_t stubs_muin;
TVector3 p_muin,p_e,p_mu;
Double_t the_rec,theX_rec,theY_rec,thmu_rec,thmuX_rec,thmuY_rec;

Double_t chi=vrtx->chi2perDegreeOfFreedom();

         auto tracks = *tracksRV;

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
	if(track.processIDofLinkedTrack()==45 and tracks.size()==3 and track.sector()==1)
		{yes2++;
                 if(code_e==track.linkedTrackID()) {theX_rec=track.xSlope(); theY_rec=track.ySlope();p_e.SetXYZ(theX_rec,theY_rec,1.0); p_e=p_e.Unit();the_rec=p_e.Angle(p_muin);}
                 if(code_mu==track.linkedTrackID()) {thmuX_rec=track.xSlope(); thmuY_rec=track.ySlope();p_mu.SetXYZ(thmuX_rec,thmuY_rec,1.0); p_mu=p_mu.Unit();thmu_rec=p_mu.Angle(p_muin);}
        	}
         }//for


Double_t posxIN=pos_on_track(x0_in,th_inx,912.7);
Double_t posyIN=pos_on_track(y0_in,th_iny,912.7);

if(stubs_muin==6 and abs(posxIN)<=1.5 and abs(posyIN)<=1.5 and chi2_muin<2){

 if(chi!=0){

                                                Double_t dotProduct_v = p_mu.Dot(p_e);
                                                TVector3 crossProduct_v = p_mu.Cross(p_e);
                                                Double_t T_v = p_muin.Dot(crossProduct_v);
                                                TVector3 im_v= p_muin.Cross(p_mu);
                                                TVector3 ie_v= p_muin.Cross(p_e);
                                                T_v = T_v>0? 1:-1;
                                                Double_t acoplanarity_v= T_v*(TMath::Pi()- acos( ((im_v).Dot(ie_v))/(im_v.Mag()*ie_v.Mag()) ));

 if(abs(acoplanarity_v)<=1 and chi<20 and thmu_rec>0.0002){
	if(the_rec>0.0 and the_rec<=0.005){h_res->Fill(the_rec-the_gen,mesmerEvent->wgt_full);
									if(the_rec-the_gen>0.001){thetaX->Fill(theX_rec,mesmerEvent->wgt_full);
									thetaY->Fill(theY_rec,mesmerEvent->wgt_full);
                                                                        thetaXg->Fill(theX_gen,mesmerEvent->wgt_full);
                                                                        thetaYg->Fill(theY_gen,mesmerEvent->wgt_full);}
	}
	else{                                                           if(the_rec-the_gen<-0.01 ){thetaX2->Fill(theX_rec,mesmerEvent->wgt_full);
                                                                        thetaY2->Fill(theY_rec,mesmerEvent->wgt_full);
                                                                        thetaXg2->Fill(theX_gen,mesmerEvent->wgt_full);
                                                                        thetaYg2->Fill(theY_gen,mesmerEvent->wgt_full);}
	}
	if(the_rec>0.005 and the_rec<=0.01){h_res1->Fill(the_rec-the_gen,mesmerEvent->wgt_full);}
	if(the_rec>0.01 and the_rec<=0.015){h_res2->Fill(the_rec-the_gen,mesmerEvent->wgt_full);}
	if(the_rec>0.015 and the_rec<=0.02){h_res3->Fill(the_rec-the_gen,mesmerEvent->wgt_full);}
	if(the_rec>0.02 and the_rec<=0.025){h_res4->Fill(the_rec-the_gen,mesmerEvent->wgt_full);}
	if(the_rec>0.025 and the_rec<=0.032){h_res5->Fill(the_rec-the_gen,mesmerEvent->wgt_full);}
		  }// if out particles
		 }//if chi!=0
		}//if mu_in

	}//if generated
     }//while


  };// function


  // Launch the parallel processing of the tree
  tp.Process(myFunction);
  // Use the ROOT::TThreadedObject::Merge method to merge the thread private histograms
  // into the final result
  auto h_resMerged   = h_res.Merge();
  auto h_res1Merged   = h_res1.Merge();
  auto h_res2Merged   = h_res2.Merge();
  auto h_res3Merged   = h_res3.Merge();
  auto h_res4Merged   = h_res4.Merge();
  auto h_res5Merged   = h_res5.Merge();

  auto thetaXMerged = thetaX.Merge();
  auto thetaYMerged = thetaY.Merge();
  auto thetaXgMerged = thetaXg.Merge();
  auto thetaYgMerged = thetaYg.Merge();

  auto thetaX2Merged = thetaX2.Merge();
  auto thetaY2Merged = thetaY2.Merge();
  auto thetaXg2Merged = thetaXg2.Merge();
  auto thetaYg2Merged = thetaYg2.Merge();

TCanvas a("a","a",700,700);
a.Divide(2,3);
a.cd(1);
h_resMerged->SetMinimum(1.);
h_resMerged->Draw("hist");
gPad->SetLogy();
a.cd(2);
h_res1Merged->SetMinimum(1.);
h_res1Merged->Draw("hist");
gPad->SetLogy();
a.cd(3);
h_res2Merged->SetMinimum(1.);
h_res2Merged->Draw("hist");
gPad->SetLogy();
a.cd(4);
h_res3Merged->SetMinimum(1.);
h_res3Merged->Draw("hist");
gPad->SetLogy();
a.cd(5);
h_res4Merged->SetMinimum(1.);
h_res4Merged->Draw("hist");
gPad->SetLogy();
a.cd(6);
h_res5Merged->SetMinimum(1.);
h_res5Merged->Draw("hist");
gPad->SetLogy();
a.SaveAs("res_bend_ID.pdf");

TCanvas b("b","b",1000,700);
b.Divide(2,2);
b.cd(1);
thetaXMerged->SetLineColor(kBlue);
thetaXgMerged->SetLineColor(kPink);
thetaXgMerged->Draw("hist");
thetaXMerged->Draw("hist same");
//gPad->SetLogy();

b.cd(2);
thetaYMerged->SetLineColor(kBlue);
thetaYgMerged->SetLineColor(kPink);
thetaYgMerged->Draw("hist");
thetaYMerged->Draw("hist same"); 
//gPad->SetLogy();

b.cd(3);
thetaX2Merged->SetLineColor(kBlue);
thetaXg2Merged->SetLineColor(kPink);
thetaX2Merged->Draw("hist");
thetaXg2Merged->Draw("hist same"); 
//gPad->SetLogy();

b.cd(4);
thetaY2Merged->SetLineColor(kBlue);
thetaYg2Merged->SetLineColor(kPink);
thetaY2Merged->Draw("hist");
thetaYg2Merged->Draw("hist same"); 
//gPad->SetLogy();

b.SaveAs("projections_bend_ID.pdf");

  return 0;
}

