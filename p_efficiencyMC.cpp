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

int p_efficiencyMC()
{
ROOT::EnableThreadSafety();
  // First enable implicit multi-threading globally, so that the implicit parallelisation is on.
  // The parameter of the call specifies the number of threads to use.
  int nthreads = 4;
  ROOT::EnableImplicitMT(nthreads);
// Create one ROOT::TThreadedObject per histogram to fill during the processing of the tree

TH1::SetDefaultSumw2(kTRUE);

   const Int_t NBINS = 6;
   Double_t edges[NBINS + 1] = {0.0, 0.005, 0.010, 0.015, 0.020, 0.025, 0.032};

   ROOT::TThreadedObject<TH1D> theta_e("theta_e", "Electron scattering reco angles from MESMER",NBINS,edges);
   ROOT::TThreadedObject<TH1D> theta_mu("theta_mu", "Muon scattering reco angles from MESMER",10,0.,0.005);

   ROOT::TThreadedObject<TH1D> theta_e_gen("theta_e_gen", "Electron scattering generated angles from MESMER",NBINS,edges);
   ROOT::TThreadedObject<TH1D> theta_mu_gen("theta_mu_gen", "Muon scattering generated angles from MESMER",10,0.,0.005);

   ROOT::TThreadedObject<TH2D> h_2d_gen("h2D","theta mu vs theta E with all cuts generated", 300,0.,0.035,150,0.,0.005);
//   ROOT::TThreadedObject<TH2D> h_2d_rec("h2D","theta mu vs theta E with all cuts reconstructed", 300,0.,0.035,150,0.,0.005);

  // Create a TTreeProcessor: specify the file and the tree in it
// ROOT::TTreeProcessorMT tp("/mnt/raid10/DATA/espedica/fairmu/Mesmer_new_1M_1hit_bend.root","cbmsim",nthreads);

  ROOT::TTreeProcessorMT tp1("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_0-5mrad_100k_1hit.root","cbmsim",nthreads);
  ROOT::TTreeProcessorMT tp2("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_5-10mrad_100k_1hit.root","cbmsim",nthreads);
  ROOT::TTreeProcessorMT tp3("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_10-15mrad_100k_1hit.root","cbmsim",nthreads);
  ROOT::TTreeProcessorMT tp4("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_15-20mrad_100k_1hit.root","cbmsim",nthreads);
  ROOT::TTreeProcessorMT tp5("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_20-25mrad_100k_1hit.root","cbmsim",nthreads);
  ROOT::TTreeProcessorMT tp6("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_25-32mrad_100k_1hit.root","cbmsim",nthreads);

        MUonERecoOutput *ReconstructionOutput = 0;

double r_wnorm[6]={21.522863161260670,38.514215499428630,41.710050452754516,49.388519620534474,61.328951873440843,105.50015224118995};

  auto myFunction = [&](TTreeReader &myReader) {
     //TTreeReaderValue<std::vector<MUonERecoOutputTrack>> tracksRV(myReader, "ReconstructionOutput");
     TTreeReaderValue<std::vector<MUonERecoOutputTrack>> tracksRV(myReader, "ReconstructedTracks");
     TTreeReaderArray<MUonETrack> MCtracksRV(myReader, "MCTrack");
     TTreeReaderValue<MUonERecoOutputVertex> vrtx(myReader, "BestVertex");
     TTreeReaderValue<MuE::Event> mesmerEvent(myReader, "MesmerEvent");

     //TTreeReaderValue<std::vector<Int_t>> MCTrack.pdgCode()(myReader, "MCTrack.PdgCode");

     // For performance reasons, a copy of the pointer associated to this thread on the
     // stack is used

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

double wnorm=1.;

if(the_gen>=0 and the_gen<0.005)wnorm=r_wnorm[0];
if(the_gen>=0.005 and the_gen<0.010)wnorm=r_wnorm[1];
if(the_gen>=0.010 and the_gen<0.015)wnorm=r_wnorm[2];
if(the_gen>=0.015 and the_gen<0.020)wnorm=r_wnorm[3];
if(the_gen>=0.020 and the_gen<0.025)wnorm=r_wnorm[4];
if(the_gen>=0.025 and the_gen<0.032)wnorm=r_wnorm[5];


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

if(stubs_muin==6 and abs(posxIN)<=1.5 and abs(posyIN)<=1.5 and chi2_muin<2){//and thmu_gen>0.0002){

theta_mu_gen->Fill(thmu_gen,mesmerEvent->wgt_full*wnorm);
theta_e_gen->Fill(the_gen,mesmerEvent->wgt_full*wnorm);

 if(chi!=0){

                                                Double_t dotProduct_v = p_mu.Dot(p_e);
                                                TVector3 crossProduct_v = p_mu.Cross(p_e);
                                                Double_t T_v = p_muin.Dot(crossProduct_v);
                                                TVector3 im_v= p_muin.Cross(p_mu);
                                                TVector3 ie_v= p_muin.Cross(p_e);
                                                T_v = T_v>0? 1:-1;
                                                Double_t acoplanarity_v= T_v*(TMath::Pi()- acos( ((im_v).Dot(ie_v))/(im_v.Mag()*ie_v.Mag()) ));

//if(abs(acoplanarity_v)<=1 and chi<20 and thmu_rec>0.0002){

theta_mu->Fill(thmu_gen,mesmerEvent->wgt_full*wnorm);
theta_e->Fill(the_gen,mesmerEvent->wgt_full*wnorm);

//		  }		// if out particles
//else{h_2d_gen->Fill(the_gen,thmu_gen,mesmerEvent->wgt_full);}
		 }//if chi!=0
		}//if mu_in

	}//if generated
     }//while


  };// function


  // Launch the parallel processing of the tree
  tp1.Process(myFunction);
  tp2.Process(myFunction);
  tp3.Process(myFunction);
  tp4.Process(myFunction);
  tp5.Process(myFunction);
  tp6.Process(myFunction);

  // Use the ROOT::TThreadedObject::Merge method to merge the thread private histograms
  // into the final result

auto theta_muMerged = theta_mu.Merge();
auto theta_eMerged = theta_e.Merge();
auto theta_mu_genMerged = theta_mu_gen.Merge();
auto theta_e_genMerged = theta_e_gen.Merge();


// auto h_2d_genMerged = h_2d_gen.Merge();

/*Int_t nxTh = h_2d_genMerged->GetNbinsX();
Int_t nyTh = h_2d_genMerged->GetNbinsY();
for (Int_t i=1; i<nxTh+1; i++) {
for (Int_t j=1; j<nyTh+1; j++) {
if (h_2d_genMerged->GetBinContent(i,j)<1) h_2d_genMerged->SetBinContent(i,j,0);}}

TCanvas b("b","b",700,700);
h_2d_gen->Draw();
b.SaveAs("2d_acocut.pdf");*/

TCanvas c("c","c",700,700);
c.Divide(1,2);
c.cd(1);
theta_mu_genMerged->SetLineColor(kPink);
theta_mu_genMerged->Draw("E");
theta_muMerged->Draw("E same");
gPad->SetLogy();
gStyle->SetOptStat(0);
c.cd(2);
theta_e_genMerged->SetLineColor(kPink);
theta_e_genMerged->Draw("E");
theta_eMerged->Draw("E same");
gStyle->SetOptStat(0);
c.SaveAs("angle_LO_MC_cut.pdf");

//theta_eMerged->SetBins(8,0.,0.032);
//theta_e_genMerged->SetBins(8,0.,0.032);

//theta_muMerged->SetBins(9,0.0002,0.0029);
//theta_mu_genMerged->SetBins(9,0.0002,0.0029);

TH1D * h2 = (TH1D*) theta_eMerged->Clone();
TH1D * h2gen = (TH1D*) theta_e_genMerged->Clone();
h2->Divide(h2gen);
TH1D * h3 = (TH1D*) theta_muMerged->Clone();
TH1D * h3gen = (TH1D*) theta_mu_genMerged->Clone();
h3->Divide(h3gen);
TCanvas a1("a1","a1",700,700);
a1.Divide(1,2);
a1.cd(1);
h3->SetMinimum(0.);
h3->Draw("E");
gStyle->SetOptStat(0);
a1.cd(2);
h2->SetMinimum(0.);
h2->Draw("E");
gStyle->SetOptStat(0);
a1.SaveAs("eff_LO_MC_cut.pdf");
h3->SaveAs("mu_eff_LO_MC_cut.root");
h2->SaveAs("el_eff_LO_MC_cut.root");

  return 0;
}

