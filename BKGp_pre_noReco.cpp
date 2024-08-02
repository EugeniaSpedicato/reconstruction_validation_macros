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

//int BKGp_pre_noReco(string path, int n){
int BKGp_pre_noReco(string path){

  int nthreads = 8;
  ROOT::EnableImplicitMT(nthreads);

/*   TFile *f1 = new TFile("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_commit_b892679b_0_infmrad_2hitFirstModules.root");
   TFile *f2 = new TFile("/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_b892679b_0_infmrad.root");//"/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_b892679b_25_32mrad.root");
   TTree *t1 = (TTree*)f1->Get("cbmsim");
   TTree *t2 = (TTree*)f2->Get("cbmsim");
   t1->AddFriend(t2);*/

TChain * cbmsim = new TChain("cbmsim");
TChain * cbmsim_g = new TChain("cbmsim");
TString gen_filename;
 TTree *t2;

//          cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_commit_2f4e96f4_MCbackground_bestConfig_%d.root",n));
//          cbmsim_g->Add("/mnt/raid10/DATA/espedica/fairmu/gen_digi/commit_2f4e96f4_MCbackground_SIM-DIGI_%d.root",n));

          cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_commit_2f4e96f4_MCbackground_bestConfig_wnorm14_0.root");
          cbmsim_g->Add("/mnt/raid10/DATA/espedica/fairmu/gen_digi/commit_2f4e96f4_MCbackground_SIM-DIGI_wnorm14_0.root");

          cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_commit_2f4e96f4_MCbackground_bestConfig_wnorm14_1.root");
          cbmsim_g->Add("/mnt/raid10/DATA/espedica/fairmu/gen_digi/commit_2f4e96f4_MCbackground_SIM-DIGI_wnorm14_1.root");

          cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_commit_2f4e96f4_MCbackground_bestConfig_wnorm14_2.root");
          cbmsim_g->Add("/mnt/raid10/DATA/espedica/fairmu/gen_digi/commit_2f4e96f4_MCbackground_SIM-DIGI_wnorm14_2.root");

          cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_commit_2f4e96f4_MCbackground_bestConfig_wnorm14_3.root");
          cbmsim_g->Add("/mnt/raid10/DATA/espedica/fairmu/gen_digi/commit_2f4e96f4_MCbackground_SIM-DIGI_wnorm14_3.root");

          cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_commit_2f4e96f4_MCbackground_bestConfig_wnorm14_4.root");
          cbmsim_g->Add("/mnt/raid10/DATA/espedica/fairmu/gen_digi/commit_2f4e96f4_MCbackground_SIM-DIGI_wnorm14_4.root");

          cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_commit_2f4e96f4_MCbackground_bestConfig_wnorm14_5.root");
          cbmsim_g->Add("/mnt/raid10/DATA/espedica/fairmu/gen_digi/commit_2f4e96f4_MCbackground_SIM-DIGI_wnorm14_5.root");

          cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_commit_2f4e96f4_MCbackground_bestConfig_wnorm14_6.root");
          cbmsim_g->Add("/mnt/raid10/DATA/espedica/fairmu/gen_digi/commit_2f4e96f4_MCbackground_SIM-DIGI_wnorm14_6.root");

          cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_commit_2f4e96f4_MCbackground_bestConfig_wnorm14_7.root");
          cbmsim_g->Add("/mnt/raid10/DATA/espedica/fairmu/gen_digi/commit_2f4e96f4_MCbackground_SIM-DIGI_wnorm14_7.root");


/*
          cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_5_10mrad_2hitFirstModules.root");
          cbmsim_g->Add(Form("/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_5_10mrad.root"); 

          cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_10_15mrad_2hitFirstModules.root");
          cbmsim_g->Add(Form("/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_10_15mrad.root"); 

          cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_15_20mrad_2hitFirstModules.root");
          cbmsim_g->Add(Form("/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_15_20mrad.root"); 

          cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_20_25mrad_2hitFirstModules.root");
          cbmsim_g->Add(Form("/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_20_25mrad.root"); 

          cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_25_32mrad_2hitFirstModules.root");
          cbmsim_g->Add(Form("/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_25_32mrad.root");

          cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_32_infmrad_2hitFirstModules.root");
          cbmsim_g->Add(Form("/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_32_infmrad.root");
*/
	  cbmsim->AddFriend(cbmsim_g);

ROOT::TTreeProcessorMT tp1(*cbmsim,nthreads);

// ROOT::TTreeProcessorMT tp1("/mnt/raid10/DATA/espedica/fairmu/gen_digi_reco/WiP_v0140_commit_258faf6b_0_infmrad_1M.root","cbmsim",nthreads);
// total:            1337.95331109692 +-         4.77058943221 μb

//ROOT::TTreeProcessorMT tp1(*t1,nthreads);
///mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/single_sample_2hitFirstModules_NOoutchi2_MCcorrections_v3_first_modules.root","cbmsim",nthreads);
//ROOT::TTreeProcessorMT tp2("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/single_sample_2hitFirstModules_NOoutchi2_MCcorrections_v3_1.root","cbmsim",nthreads);



/*ROOT::TTreeProcessorMT tp3("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/single_sample_2hitFirstModules_NOoutchi2_MCcorrections_v2_2.root","cbmsim",nthreads);
ROOT::TTreeProcessorMT tp4("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/single_sample_2hitFirstModules_NOoutchi2_MCcorrections_v2_3.root","cbmsim",nthreads);
ROOT::TTreeProcessorMT tp5("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/single_sample_2hitFirstModules_NOoutchi2_MCcorrections_v2_4.root","cbmsim",nthreads);
ROOT::TTreeProcessorMT tp6("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/single_sample_2hitFirstModules_NOoutchi2_MCcorrections_v2_5.root","cbmsim",nthreads);
ROOT::TTreeProcessorMT tp7("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/single_sample_2hitFirstModules_NOoutchi2_MCcorrections_v2_6.root","cbmsim",nthreads);*/

//ROOT::TTreeProcessorMT tp2("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_skim_run_6_v2_2hits_nochi2_MCcorr_gioskim_saskia_1.root","cbmsim",nthreads);


      MUonERecoOutput *ReconstructionOutput = 0;
        MuE::Event *MesmerEvent = 0;



TH1::SetDefaultSumw2(kTRUE);
   const Int_t NBINS = 6;
   Double_t edges[NBINS + 1] = {0.0, 0.005, 0.010, 0.015, 0.020, 0.025, 0.032};
   ROOT::TThreadedObject<TH1D> theta_e("theta_e", "Electron scattering reco angles",10,0.,0.035);
   ROOT::TThreadedObject<TH1D> theta_mu("theta_mu", "Muon scattering reco angles",20,0.,0.005);
   ROOT::TThreadedObject<TH1D> h_opening("h_opening", "Opening angle reco events",35,0.,0.035);
   ROOT::TThreadedObject<TH1D> d_aco("d_aco_real", "Acoplanarity",600,-3.2,3.2);



   ROOT::TThreadedObject<TH1D> theta_e_pre("theta_e_pre", "Electron scattering reco angles pre-cuts",10,0.,0.035);
   ROOT::TThreadedObject<TH1D> theta_mu_pre("theta_mu_pre", "Muon scattering reco angles pre-cuts",20,0.,0.005);
   ROOT::TThreadedObject<TH1D> h_opening_pre("h_opening_pre", "Opening angle reco events pre-cuts",35,0.,0.035);
   ROOT::TThreadedObject<TH1D> d_aco_pre("d_aco_real_pre", "Acoplanarity pre-cuts",600,-3.2,3.2);

   ROOT::TThreadedObject<TH1D> theq("theq", "Equal angles",20,0.,0.005);
   ROOT::TThreadedObject<TH1D> count("count","count",2,0.,2.);
   ROOT::TThreadedObject<TH1D> count_after("count_after","count_after",2,0.,2.);
   ROOT::TThreadedObject<TH1D> count_pre("count_pre","count_pre",2,0.,2.);

   ROOT::TThreadedObject<TH1D> h_chi2("h_chi2","chi2 mu in", 100,0,10);
   ROOT::TThreadedObject<TH1D> h_xslope("h_xslope","xSlope mu in ", 200,0.,0.02);


   const Int_t NBINS2 = 12;
   Double_t edges2[NBINS2 + 1] = {0.,0.0002,0.0004,0.0006,0.0008,0.001,0.00125,0.0015,0.00175,0.002,0.003,0.004,0.005};


  /* ROOT::TThreadedObject<TH2D> h_2d("h2D","theta mu vs theta E with all cuts", 16, 0.,0.032,NBINS2,edges2);
   ROOT::TThreadedObject<TH2D> h_2d_pre("h2D_pre","theta mu vs theta E with all cuts", 16, 0.,0.032,NBINS2,edges2);*/

   ROOT::TThreadedObject<TH2D> h_2d("h2D","theta mu vs theta E with all cuts", 320, 0.,0.032,50,0.,0.005);
   ROOT::TThreadedObject<TH2D> h_2d_pre("h2D_pre","theta mu vs theta E with all cuts", 320, 0.,0.032,50,0.,0.005);

/* double wnorm=0;
if(n==0) wnorm=31.20593;
if(n==1) wnorm=17.73709;
if(n==2) wnorm=9.6284348;
*/

 auto myFunction = [&](TTreeReader &myReader) {

     TTreeReaderValue<TClonesArray> RVMCTrack(myReader, "MCTrack");
     TTreeReaderValue<std::vector<MUonERecoOutputTrack>> RVtracks(myReader, "ReconstructedTracks");
     TTreeReaderValue<MUonERecoOutputVertex> vrtx(myReader, "BestVertex");
     TTreeReaderValue<std::vector<MUonERecoOutputHit>> RVstubs(myReader, "ReconstructedHits");
     TTreeReaderValue<Double_t> wgt_f(myReader,"wgt_full");

     while (myReader.Next()) {

double truth0=0;
double truth1=0;
double truth2=0;
auto MCTrack=*RVMCTrack;


        for(int n = 0; n < MCTrack.GetEntries(); n++) {
         const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack.At(n));
         if(MCTr->interactionID()==45){
         if(MCTr->pdgCode()==-11){truth0=1;}
         if(MCTr->pdgCode()==11){truth1=1;}
         if(MCTr->pdgCode()==-13){truth2=1;}
         }
        }
if( (truth0==1 and truth1==1) or (truth2==1 and truth1==1) ){

double reco_v=0.; double more_reco_v=0.; double reco0_v=0.;
int yes2=0; int yes_v=0;
int code_mu=-99; int code_e=-99; int code_mu_in=-99;
double z_fix=912.7;

auto wgt = *wgt_f;

 double chi=vrtx->chi2perDegreeOfFreedom();

int yes_mu=0;
int yes_e=0;
double th_inx,th_iny,x0_in,y0_in;
double chi2_muin;
int sec0=0;
int sec1=0;
int stubs_muin=0.;
double th_muin=0.;


double the_rec=-99.;
double thmu_rec=-99.;
TVector3 p_e,p_mu,p_muin;

        std::vector<MUonERecoOutputHit> hits_mu=vrtx->outgoingMuon().hits();
        std::vector<double> pos_mu; pos_mu.resize(6);
        for(int p=0;p<hits_mu.size();p++)pos_mu.push_back(hits_mu.at(p).position());
        std::vector<MUonERecoOutputHit> hits_e=vrtx->outgoingElectron().hits();
        std::vector<double> pos_e; pos_e.resize(6);
        for(int p=0;p<hits_e.size();p++)pos_e.push_back(hits_e.at(p).position());

count_pre->Fill(1.,wgt);

         auto tracks = *RVtracks;
         for (auto&& track : tracks) {
        if(track.sector()==0) sec0++;
        if(track.sector()==1) sec1++;
        }

         for (auto&& track : tracks) {
        std::vector<MUonERecoOutputHit> hits_=track.hits();

        if(track.sector()==0 and sec0==1){
	stubs_muin=hits_.size();
        th_inx=track.xSlope();
        th_iny=track.ySlope();
        x0_in=track.x0();
        y0_in=track.y0();
        chi2_muin=track.chi2perDegreeOfFreedom();
        p_muin.SetXYZ(th_inx,th_iny,1.0);
        p_muin=p_muin.Unit();
        th_muin=p_muin.Theta();
        h_chi2->Fill(chi2_muin);
        h_xslope->Fill(th_inx);

                        }
	if(track.sector()==1)
	{
	std::vector<double> pos; pos.resize(6);
        for(int p=0;p<hits_.size();p++)pos.push_back(hits_.at(p).position());

         if(std::equal(pos.begin(),pos.end(),pos_mu.begin())){p_mu.SetXYZ(track.xSlope(),track.ySlope(),1.0); p_mu.Unit(); thmu_rec=p_mu.Angle(p_muin);}
        else if(std::equal(pos.begin(),pos.end(),pos_e.begin())){p_e.SetXYZ(track.xSlope(),track.ySlope(),1.0); p_e.Unit(); the_rec=p_e.Angle(p_muin);}

	yes2++;}
}

double posxIN=99.;//pos_on_track(x0_in,th_inx,z_fix);
double posyIN=99.;//pos_on_track(y0_in,th_iny,z_fix);

int stub0 = 0;
int stub1 = 0;


         auto stubs = *RVstubs;
         for (auto&& stub : stubs) {
		if(stub.stationID()==0){stub0++; if(stub.moduleID()==4){posxIN=stub.position(); } else if(stub.moduleID()==5){posyIN=stub.position();}   }
		if(stub.stationID()==1)stub1++;
	}

 if(sec0==1 and stubs_muin==6 and abs(posxIN)<=1.5 and abs(posyIN)<=1.5 and chi2_muin<=2 and stub0==6 and th_muin<0.004){// and stub1<=15){

count->Fill(1.,wgt);

if(chi!=0 and the_rec!=-99 and thmu_rec!=-99){

                                                double dotProduct_v = p_mu.Dot(p_e);
                                                TVector3 crossProduct_v = p_mu.Cross(p_e);
                                                double T_v = p_muin.Dot(crossProduct_v);
                                                TVector3 im_v= p_muin.Cross(p_mu);
                                                TVector3 ie_v= p_muin.Cross(p_e);
                                                T_v = T_v>0? 1:-1;
                                                double acoplanarity_v= T_v*(TMath::Pi()- acos( ((im_v).Dot(ie_v))/(im_v.Mag()*ie_v.Mag()) ));


 d_aco_pre->Fill(acoplanarity_v,wgt);
 h_2d_pre->Fill(the_rec,thmu_rec,wgt);
 theta_mu_pre->Fill(thmu_rec,wgt);
 theta_e_pre->Fill(the_rec,wgt);
 h_opening_pre->Fill(p_mu.Angle(p_e),wgt);


  if(abs(acoplanarity_v)<=1 and chi<20 and thmu_rec>0.0002 and stub1<=15){// and vrtx->zPositionFit()<915 and vrtx->zPositionFit()>907){
//  if(thmu_rec>0.0002 and the_rec<0.032){// and vrtx->zPositionFit()<915 and vrtx->zPositionFit()>907){
//first
 if(yes2>=2){

count_after->Fill(1.,wgt);

 d_aco->Fill(acoplanarity_v,wgt);
 h_2d->Fill(the_rec,thmu_rec,wgt);
 theta_mu->Fill(thmu_rec,wgt);
 theta_e->Fill(the_rec,wgt);
 h_opening->Fill(p_mu.Angle(p_e),wgt);

if(thmu_rec>=the_rec-0.00005 and thmu_rec<=the_rec+0.00005)theq->Fill(thmu_rec);
				}//yes2>=2
			}//aco e chi
		}//chi!=0
	}//mu_in
yes2=0;
  }//end truth
 } //end of general while
}; //end of myfunction

  tp1.Process(myFunction);
  //tp2.Process(myFunction);

  auto h_2dMerged = h_2d.Merge();
  auto theta_muMerged = theta_mu.Merge();
  auto theta_eMerged = theta_e.Merge();
  auto d_acoM = d_aco.Merge();
  auto h_openingM = h_opening.Merge();

  auto h_2d_preMerged = h_2d_pre.Merge();
  auto theta_mu_preMerged = theta_mu_pre.Merge();
  auto theta_e_preMerged = theta_e_pre.Merge();
  auto d_aco_preM = d_aco_pre.Merge();
  auto h_opening_preM = h_opening_pre.Merge();
  auto theqM = theq.Merge();
  auto h_chi2M = h_chi2.Merge();
  auto h_xslopeM = h_xslope.Merge();
  auto countM = count.Merge();
  auto count_afterM = count_after.Merge();
  auto count_preM = count_pre.Merge();
/*
double sigma=0.;
if(n==0) sigma=11.782165; //0
if(n==1) sigma=11.248182; //1
if(n==2) sigma=11.412257; //2
*/

/*
double sigma=0.;
if(n==0) sigma=12.79992003470; //0
if(n==1) sigma=12.51073350106; //1
if(n==3) sigma=13.58726609942; //3
*/


double sigma=12.9221;

cout << "count_preM.Integral() " << count_preM->Integral() << endl;
cout << "count_preM.Integral() under+over " << count_preM->Integral(0,count_preM->GetNbinsX()+1) << endl;

cout << "count_preM.GetEffectiveEntries()() " << count_preM->GetEffectiveEntries() << endl;
//cout << "count_preM.Integral()/entries " << count_preM->Integral()/cbmsim->GetEntries() << endl;

cout << "countM.Integral() " << countM->Integral() << endl;
cout << "count_afterM.Integral() " << count_afterM->Integral() << endl;
cout << "Events that passes selection over the ones that pass the fiducial " << count_afterM->Integral()/countM->Integral()*100 << "%" << endl;
cout << "Events that passes selection over all " << count_afterM->Integral()/count_preM->Integral()*100 << "%" << endl;
cout << "Taking into account the CS = "<<sigma<<" micron : " << count_afterM->Integral()/count_preM->Integral()*sigma << endl;





TCanvas t1("t1","t1",700,700);
t1.Divide(1,2);
t1.cd(1);
h_chi2M->Draw("hist");
t1.cd(2);
h_xslopeM->Draw("hist");
t1.SaveAs("MC_muin.pdf");

cout << "pre " << h_2dMerged->Integral() << endl;
if(path=="bkg_mesmer"){
h_2dMerged->Scale(sigma/count_preM->Integral());
theta_muMerged->Scale(sigma/count_preM->Integral());
theta_eMerged->Scale(sigma/count_preM->Integral());
d_acoM->Scale(sigma/count_preM->Integral());
h_openingM->Scale(sigma/count_preM->Integral());
h_2d_preMerged->Scale(sigma/count_preM->Integral());
theta_mu_preMerged->Scale(sigma/count_preM->Integral());
theta_e_preMerged->Scale(sigma/count_preM->Integral());
d_aco_preM->Scale(sigma/count_preM->Integral());
h_opening_preM->Scale(sigma/count_preM->Integral());
}

cout << "post " << h_2dMerged->Integral() << endl;



Int_t nx = h_2dMerged->GetNbinsX();
Int_t ny = h_2dMerged->GetNbinsY();
for (Int_t i=1; i<nx+1; i++) {
for (Int_t j=1; j<ny+1; j++) {
if (h_2dMerged->GetBinContent(i,j)<1) h_2dMerged->SetBinContent(i,j,0);}}

Int_t nxpre = h_2d_preMerged->GetNbinsX();
Int_t nypre = h_2d_preMerged->GetNbinsY();
for (Int_t i=1; i<nxpre+1; i++) {
for (Int_t j=1; j<nypre+1; j++) {
if (h_2d_preMerged->GetBinContent(i,j)<1) h_2d_preMerged->SetBinContent(i,j,0);}}


TCanvas t("t","t",700,700);
t.Divide(2,3);
t.cd(1);
theta_mu_preMerged->SetLineColor(kRed);
theta_mu_preMerged->Draw("hist");
theta_muMerged->Draw("hist same");
theta_mu_preMerged->SetMinimum(0.);
t.cd(2);
theta_e_preMerged->SetLineColor(kRed);
theta_e_preMerged->Draw("hist");
theta_eMerged->Draw("hist same");
theta_e_preMerged->SetMinimum(0.);
t.cd(3);
d_aco_preM->SetLineColor(kRed);
d_aco_preM->Draw("hist");
d_acoM->Draw("hist same");
d_aco_preM->SetMinimum(0.);
t.cd(4);
h_opening_preM->SetLineColor(kRed);
h_opening_preM->Draw("hist");
h_openingM->Draw("hist same");
h_opening_preM->SetMinimum(0.);
t.cd(5);
h_2d_preMerged->SetMarkerColor(kRed);
h_2d_preMerged->Draw();
h_2dMerged->SetMarkerColor(kBlue);
h_2dMerged->Draw("same");
gStyle->SetOptStat(0);
t.SaveAs(Form("%s/info_bkg.pdf",path.c_str()));





h_2dMerged->SaveAs(Form("%s/2D_MC_parallel_bkg.root",path.c_str()));

theta_muMerged->SaveAs(Form("%s/theta_mu_MC_parallel_bkg.root",path.c_str()));

theta_eMerged->SaveAs(Form("%s/theta_e_MC_parallel_bkg.root",path.c_str()));

d_acoM->SaveAs(Form("%s/d_aco_MC_parallel_bkg.root",path.c_str()));

h_openingM->SaveAs(Form("%s/opening_MC_parallel_bkg.root",path.c_str()));


h_2d_preMerged->SaveAs(Form("%s/preCuts_2D_MC_parallel_bkg.root",path.c_str()));

theta_mu_preMerged->SaveAs(Form("%s/preCuts_theta_mu_MC_parallel_bkg.root",path.c_str()));

theta_e_preMerged->SaveAs(Form("%s/preCuts_theta_e_MC_parallel_bkg.root",path.c_str()));

d_aco_preM->SaveAs(Form("%s/preCuts_d_aco_MC_parallel_bkg.root",path.c_str()));

h_opening_preM->SaveAs(Form("%s/preCuts_opening_MC_parallel_bkg.root",path.c_str()));

//theqM->SaveAs(Form("%s/th_eqM_pre.pdf");

return 0;
}
