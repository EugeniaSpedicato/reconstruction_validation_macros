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

void resolutions(int nhits, int index){

TChain * cbmsim = new TChain("cbmsim");
TString gen_filename;



/*
if(nhits==0){
  if(index==0)cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_0-5mrad_%dhit_NOoutchi2_1M.root",static_cast<char>(nhits)));
  else if(index==1)cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_5-10mrad_%dhit_NOoutchi2_1M.root",static_cast<char>(nhits)));
  else if(index==2)cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_10-15mrad_%dhit_NOoutchi2_1M.root",static_cast<char>(nhits)));
  else if(index==3)cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_15-20mrad_%dhit_NOoutchi2_1M.root",static_cast<char>(nhits)));
  else if(index==4)cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_20-25mrad_%dhit_NOoutchi2_1M.root",static_cast<char>(nhits)));
  else if(index==5)cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_25-32mrad_%dhit_NOoutchi2_1M.root",static_cast<char>(nhits)));
  else cout << "WRONG NUMBER" << endl;
}
else{
  if(index==0)cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_0-5mrad_1M_%dhitFirstModules_NOoutchi2_1M.root",static_cast<char>(nhits)));
  else if(index==1)cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_5-10mrad_1M_%dhitFirstModules_NOoutchi2_1M.root",static_cast<char>(nhits)));
  else if(index==2)cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_10-15mrad_1M_%dhitFirstModules_NOoutchi2_1M.root",static_cast<char>(nhits)));
  else if(index==3)cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_15-20mrad_1M_%dhitFirstModules_NOoutchi2_1M.root",static_cast<char>(nhits)));
  else if(index==4)cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_20-25mrad_1M_%dhitFirstModules_NOoutchi2_1M.root",static_cast<char>(nhits)));
  else if(index==5)cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_25-32mrad_1M_%dhitFirstModules_NOoutchi2_1M.root",static_cast<char>(nhits)));
  else cout << "WRONG NUMBER" << endl;
}
*/


if(index==0){cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_commit_258faf6b_0_5mrad_%dhitFirstModules_allNewFeatures.root",static_cast<char>(nhits)));
          gen_filename="/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_258faf6b_0_5mrad.root";}

else if(index==1){cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_commit_258faf6b_5_10mrad_%dhitFirstModules_allNewFeatures.root",static_cast<char>(nhits)));
          gen_filename="/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_258faf6b_5_10mrad.root";}

else if(index==2){cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_commit_258faf6b_10_15mrad_%dhitFirstModules_allNewFeatures.root",static_cast<char>(nhits)));
          gen_filename="/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_258faf6b_10_15mrad.root";}

else if(index==3){cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_commit_258faf6b_15_20mrad_%dhitFirstModules_allNewFeatures.root",static_cast<char>(nhits)));
          gen_filename="/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_258faf6b_15_20mrad.root";}

else if(index==4){cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_commit_258faf6b_20_25mrad_%dhitFirstModules_allNewFeatures.root",static_cast<char>(nhits)));
          gen_filename="/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_258faf6b_20_25mrad.root";}

else if(index==5){cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_commit_258faf6b_25_32mrad_%dhitFirstModules_allNewFeatures.root",static_cast<char>(nhits)));
          gen_filename="/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_258faf6b_25_32mrad.root";}

/*
if(index==0){cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/reco/old_WiP/WiP_v0140_commit_0996edd4_0_5mrad_%dhitFirstModules.root",static_cast<char>(nhits)));
          gen_filename="/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_258faf6b_0_5mrad.root";}

else if(index==1){cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/reco/old_WiP/WiP_v0140_commit_0996edd4_5_10mrad_%dhitFirstModules.root",static_cast<char>(nhits)));
          gen_filename="/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_258faf6b_5_10mrad.root";}

else if(index==2){cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/reco/old_WiP/WiP_v0140_commit_0996edd4_10_15mrad_%dhitFirstModules.root",static_cast<char>(nhits)));
          gen_filename="/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_258faf6b_10_15mrad.root";}

else if(index==3){cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/reco/old_WiP/WiP_v0140_commit_0996edd4_15_20mrad_%dhitFirstModules.root",static_cast<char>(nhits)));
          gen_filename="/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_258faf6b_15_20mrad.root";}

else if(index==4){cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/reco/old_WiP/WiP_v0140_commit_0996edd4_20_25mrad_%dhitFirstModules.root",static_cast<char>(nhits)));
          gen_filename="/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_258faf6b_20_25mrad.root";}

else if(index==5){cbmsim->Add(Form("/mnt/raid10/DATA/espedica/fairmu/reco/old_WiP/WiP_v0140_commit_0996edd4_25_32mrad_%dhitFirstModules.root",static_cast<char>(nhits)));
          gen_filename="/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_258faf6b_25_32mrad.root";}
*/

else cout << "ERROR"<< endl;


   TFile *f2 = new TFile(gen_filename);
   TTree *t2 = (TTree*)f2->Get("cbmsim");

cout << "entries cbmsim " << cbmsim->GetEntries() << endl;
cout << "entries t2 " << t2->GetEntries() << endl;


   cbmsim->AddFriend(t2);
   cout << t2->GetEntries() << endl;
   cout << cbmsim->GetEntries() << endl;
   cbmsim->Print();



        TClonesArray *MCTrack = 0;
        TClonesArray *TrackerStripDigis = 0;
        TClonesArray *TrackerPoints = 0;
        TClonesArray *TrackerStubs = 0;
        MuE::Event *MesmerEvent = 0;
        MUonERecoOutput *ReconstructionOutput = 0;

        cbmsim->SetBranchAddress("MCTrack", &MCTrack);
        cbmsim->SetBranchAddress("TrackerPoints", &TrackerPoints);
        cbmsim->SetBranchAddress("TrackerStripDigis", &TrackerStripDigis);
        cbmsim->SetBranchAddress("TrackerStubs", &TrackerStubs);
        cbmsim->SetBranchAddress("MesmerEvent", &MesmerEvent);
        cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);


   const Int_t NBINS_mu = 14;
   Double_t edges_mu[NBINS_mu + 1] = {0.,0.0001,0.0002,0.0003,0.0004,0.0005,0.0006,0.0007,0.0008,0.0009,0.001,0.002,0.003,0.004,0.005};
   const Int_t NBINS = 14;
   Double_t edges_el[NBINS + 1] = {0.,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.010, 0.015, 0.020, 0.025, 0.032};

TH1D::SetDefaultSumw2(kTRUE);

        std::vector<TH1D*> h_res_vrtx_el(14);
        std::vector<TH1D*> h_res_vrtx_mu(14);
        std::vector<TH1D*> h_res_el_pre(14);
        std::vector<TH1D*> h_res_mu_pre(14);
        std::vector<TH1D*> h_res_vrtx_el_clones(14);
        std::vector<TH1D*> h_res_vrtx_mu_clones(14);
        std::vector<TH1D*> h_res_el_pre_clones(14);
        std::vector<TH1D*> h_res_mu_pre_clones(14);

h_res_vrtx_el.at(0)=new TH1D("res0_vrtx_el", "(the_rec_vrtx-the_true), theta_el[0,1]mrad",100,-0.001,0.001);
h_res_vrtx_el.at(1)=new TH1D("res1_vrtx_el", "(the_rec_vrtx-the_true), theta_el[1,2]mrad",100,-0.001,0.001);
h_res_vrtx_el.at(2)=new TH1D("res2_vrtx_el", "(the_rec_vrtx-the_true), theta_el[2,3]mrad",100,-0.001,0.001);
h_res_vrtx_el.at(3)=new TH1D("res3_vrtx_el", "(the_rec_vrtx-the_true), theta_el[3,4]mrad",100,-0.001,0.001);
h_res_vrtx_el.at(4)=new TH1D("res4_vrtx_el", "(the_rec_vrtx-the_true), theta_el[4,5]mrad",200,-0.002,0.002);
h_res_vrtx_el.at(5)=new TH1D("res5_vrtx_el", "(the_rec_vrtx-the_true), theta_el[5,6]mrad",200,-0.002,0.002);
h_res_vrtx_el.at(6)=new TH1D("res6_vrtx_el", "(the_rec_vrtx-the_true), theta_el[6,7]mrad",200,-0.002,0.002);
h_res_vrtx_el.at(7)=new TH1D("res7_vrtx_el", "(the_rec_vrtx-the_true), theta_el[7,8]mrad",200,-0.002,0.002);
h_res_vrtx_el.at(8)=new TH1D("res8_vrtx_el", "(the_rec_vrtx-the_true), theta_el[8,9]mrad",200,-0.002,0.002);
h_res_vrtx_el.at(9)=new TH1D("res9_vrtx_el", "(the_rec_vrtx-the_true), theta_el[9,10]mrad",200,-0.002,0.002);
h_res_vrtx_el.at(10)=new TH1D("res10_vrtx_el", "(the_rec_vrtx-the_true), theta_el[10,15]mrad",200,-0.01,0.01);
h_res_vrtx_el.at(11)=new TH1D("res11_vrtx_el", "(the_rec_vrtx-the_true), theta_el[15,20]mrad",200,-0.02,0.02);
h_res_vrtx_el.at(12)=new TH1D("res12_vrtx_el", "(the_rec_vrtx-the_true), theta_el[20,25]mrad",200,-0.02,0.02);
h_res_vrtx_el.at(13)=new TH1D("res13_vrtx_el", "(the_rec_vrtx-the_true), theta_el[25,32]mrad",200,-0.02,0.02);


h_res_vrtx_mu.at(0)=new TH1D("res0_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_mu[0.,0.1]mrad",120,-0.0004,0.0004);
h_res_vrtx_mu.at(1)=new TH1D("res1_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_mu[0.1,0.2]mrad",120,-0.0004,0.0004);
h_res_vrtx_mu.at(2)=new TH1D("res2_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_mu[0.2,0.3]mrad",120,-0.0004,0.0004);
h_res_vrtx_mu.at(3)=new TH1D("res3_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_mu[0.3,0.4]mrad",120,-0.0004,0.0004);
h_res_vrtx_mu.at(4)=new TH1D("res4_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_mu[0.4,0.5]mrad",120,-0.0004,0.0004);
h_res_vrtx_mu.at(5)=new TH1D("res5_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_mu[0.5,0.6]mrad",120,-0.0004,0.0004);
h_res_vrtx_mu.at(6)=new TH1D("res6_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_mu[0.6,0.7]mrad",120,-0.0004,0.0004);
h_res_vrtx_mu.at(7)=new TH1D("res7_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_mu[0.7,0.8]mrad",120,-0.0004,0.0004);
h_res_vrtx_mu.at(8)=new TH1D("res8_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_mu[0.8,0.9]mrad",120,-0.0004,0.0004);
h_res_vrtx_mu.at(9)=new TH1D("res9_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_mu[0.9,1]mrad",120,-0.0004,0.0004);
h_res_vrtx_mu.at(10)=new TH1D("res10_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_mu[0.3,0.4]mrad",120,-0.0004,0.0004);
h_res_vrtx_mu.at(11)=new TH1D("res11_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_mu[0.4,0.5]mrad",120,-0.0004,0.0004);
h_res_vrtx_mu.at(12)=new TH1D("res12_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_mu[3,4]mrad",120,-0.0004,0.0004);
h_res_vrtx_mu.at(13)=new TH1D("res13_vrtx_mu", "(thmu_rec_vrtx-thmu_true), theta_mu[4,5]mrad",120,-0.0004,0.0004);


h_res_el_pre.at(0)=new TH1D("res0_el_pre", "(the_rec-the_true), theta_el[0,1]mrad",100,-0.001,0.001);
h_res_el_pre.at(1)=new TH1D("res1_el_pre", "(the_rec-the_true), theta_el[1,2]mrad",100,-0.001,0.001);
h_res_el_pre.at(2)=new TH1D("res2_el_pre", "(the_rec-the_true), theta_el[2,3]mrad",100,-0.001,0.001);
h_res_el_pre.at(3)=new TH1D("res3_el_pre", "(the_rec-the_true), theta_el[3,4]mrad",100,-0.001,0.001);
h_res_el_pre.at(4)=new TH1D("res4_el_pre", "(the_rec-the_true), theta_el[4,5]mrad",200,-0.002,0.002);
h_res_el_pre.at(5)=new TH1D("res5_el_pre", "(the_rec-the_true), theta_el[5,6]mrad",200,-0.002,0.002);
h_res_el_pre.at(6)=new TH1D("res6_el_pre", "(the_rec-the_true), theta_el[6,7]mrad",200,-0.002,0.002);
h_res_el_pre.at(7)=new TH1D("res7_el_pre", "(the_rec-the_true), theta_el[7,8]mrad",200,-0.002,0.002);
h_res_el_pre.at(8)=new TH1D("res8_el_pre", "(the_rec-the_true), theta_el[8,9]mrad",200,-0.002,0.002);
h_res_el_pre.at(9)=new TH1D("res9_el_pre", "(the_rec-the_true), theta_el[9,10]mrad",200,-0.002,0.002);
h_res_el_pre.at(10)=new TH1D("res10_el_pre", "(the_rec-the_true), theta_el[10,15]mrad",200,-0.01,0.01);
h_res_el_pre.at(11)=new TH1D("res11_el_pre", "(the_rec-the_true), theta_el[15,20]mrad",200,-0.02,0.02);
h_res_el_pre.at(12)=new TH1D("res12_el_pre", "(the_rec-the_true), theta_el[20,25]mrad",200,-0.02,0.02);
h_res_el_pre.at(13)=new TH1D("res13_el_pre", "(the_rec-the_true), theta_el[25,32]mrad",200,-0.02,0.02);


h_res_mu_pre.at(0)=new TH1D("res0_mu_pre", "(thmu_rec-thmu_true), theta_mu[0.,0.1]mrad",120,-0.0004,0.0004);
h_res_mu_pre.at(1)=new TH1D("res1_mu_pre", "(thmu_rec-thmu_true), theta_mu[0.1,0.2]mrad",120,-0.0004,0.0004);
h_res_mu_pre.at(2)=new TH1D("res2_mu_pre", "(thmu_rec-thmu_true), theta_mu[0.2,0.3]mrad",120,-0.0004,0.0004);
h_res_mu_pre.at(3)=new TH1D("res3_mu_pre", "(thmu_rec-thmu_true), theta_mu[0.3,0.4]mrad",120,-0.0004,0.0004);
h_res_mu_pre.at(4)=new TH1D("res4_mu_pre", "(thmu_rec-thmu_true), theta_mu[0.4,0.5]mrad",120,-0.0004,0.0004);
h_res_mu_pre.at(5)=new TH1D("res5_mu_pre", "(thmu_rec-thmu_true), theta_mu[0.5,0.6]mrad",120,-0.0004,0.0004);
h_res_mu_pre.at(6)=new TH1D("res6_mu_pre", "(thmu_rec-thmu_true), theta_mu[0.6,0.7]mrad",120,-0.0004,0.0004);
h_res_mu_pre.at(7)=new TH1D("res7_mu_pre", "(thmu_rec-thmu_true), theta_mu[0.7,0.8]mrad",120,-0.0004,0.0004);
h_res_mu_pre.at(8)=new TH1D("res8_mu_pre", "(thmu_rec-thmu_true), theta_mu[0.8,0.9]mrad",120,-0.0004,0.0004);
h_res_mu_pre.at(9)=new TH1D("res9_mu_pre", "(thmu_rec-thmu_true), theta_mu[0.9,1]mrad",120,-0.0004,0.0004);
h_res_mu_pre.at(10)=new TH1D("res10_mu_pre", "(thmu_rec-thmu_true), theta_mu[1,2]mrad",120,-0.0004,0.0004);
h_res_mu_pre.at(11)=new TH1D("res11_mu_pre", "(thmu_rec-thmu_true), theta_mu[2,3]mrad",120,-0.0004,0.0004);
h_res_mu_pre.at(12)=new TH1D("res12_mu_pre", "(thmu_rec-thmu_true), theta_mu[3,4]mrad",120,-0.0004,0.0004);
h_res_mu_pre.at(13)=new TH1D("res13_mu_pre", "(thmu_rec-thmu_true), theta_mu[4,5]mrad",120,-0.0004,0.0004);


h_res_vrtx_el_clones.at(0)=new TH1D("res0_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[0,1]mrad with clones",100,-0.001,0.001);
h_res_vrtx_el_clones.at(1)=new TH1D("res1_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[1,2]mrad with clones",100,-0.001,0.001);
h_res_vrtx_el_clones.at(2)=new TH1D("res2_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[2,3]mrad with clones",100,-0.001,0.001);
h_res_vrtx_el_clones.at(3)=new TH1D("res3_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[3,4]mrad with clones",100,-0.001,0.001);
h_res_vrtx_el_clones.at(4)=new TH1D("res4_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[4,5]mrad with clones",200,-0.002,0.002);
h_res_vrtx_el_clones.at(5)=new TH1D("res5_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[5,6]mrad with clones",200,-0.002,0.002);
h_res_vrtx_el_clones.at(6)=new TH1D("res6_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[6,7]mrad with clones",200,-0.002,0.002);
h_res_vrtx_el_clones.at(7)=new TH1D("res7_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[7,8]mrad with clones",200,-0.002,0.002);
h_res_vrtx_el_clones.at(8)=new TH1D("res8_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[8,9]mrad with clones",200,-0.002,0.002);
h_res_vrtx_el_clones.at(9)=new TH1D("res9_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[9,10]mrad with clones",200,-0.002,0.002);
h_res_vrtx_el_clones.at(10)=new TH1D("res10_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[10,15]mrad with clones",200,-0.01,0.01);
h_res_vrtx_el_clones.at(11)=new TH1D("res11_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[15,20]mrad with clones",200,-0.02,0.02);
h_res_vrtx_el_clones.at(12)=new TH1D("res12_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[20,25]mrad with clones",200,-0.02,0.02);
h_res_vrtx_el_clones.at(13)=new TH1D("res13_vrtx_el_clones", "(the_rec_vrtx-the_true), theta_el[25,32]mrad with clones",200,-0.02,0.02);

h_res_vrtx_mu_clones.at(0)=new TH1D("res0_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_mu[0.,0.1]mrad with clones",120,-0.0004,0.0004);
h_res_vrtx_mu_clones.at(1)=new TH1D("res1_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_mu[0.1,0.2]mrad with clones",120,-0.0004,0.0004);
h_res_vrtx_mu_clones.at(2)=new TH1D("res2_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_mu[0.2,0.3]mrad with clones",120,-0.0004,0.0004);
h_res_vrtx_mu_clones.at(3)=new TH1D("res3_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_mu[0.3,0.4]mrad with clones",120,-0.0004,0.0004);
h_res_vrtx_mu_clones.at(4)=new TH1D("res4_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_mu[0.4,0.5]mrad with clones",120,-0.0004,0.0004);
h_res_vrtx_mu_clones.at(5)=new TH1D("res5_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_mu[0.5,0.6]mrad with clones",120,-0.0004,0.0004);
h_res_vrtx_mu_clones.at(6)=new TH1D("res6_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_mu[0.6,0.7]mrad with clones",120,-0.0004,0.0004);
h_res_vrtx_mu_clones.at(7)=new TH1D("res7_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_mu[0.7,0.8]mrad with clones",120,-0.0004,0.0004);
h_res_vrtx_mu_clones.at(8)=new TH1D("res8_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_mu[0.8,0.9]mrad with clones",120,-0.0004,0.0004);
h_res_vrtx_mu_clones.at(9)=new TH1D("res9_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_mu[0.9,1]mrad with clones",120,-0.0004,0.0004);
h_res_vrtx_mu_clones.at(10)=new TH1D("res10_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_mu[1,2]mrad with clones",120,-0.0004,0.0004);
h_res_vrtx_mu_clones.at(11)=new TH1D("res11_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_mu[2,3]mrad with clones",120,-0.0004,0.0004);
h_res_vrtx_mu_clones.at(12)=new TH1D("res12_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_mu[3,4]mrad with clones",120,-0.0004,0.0004);
h_res_vrtx_mu_clones.at(13)=new TH1D("res13_vrtx_mu_clones", "(thmu_rec_vrtx-thmu_true), theta_mu[4,5]mrad with clones",120,-0.0004,0.0004);

h_res_el_pre_clones.at(0)=new TH1D("res0_el_pre_clones", "(the_rec-the_true), theta_el[0,1]mrad with clones",100,-0.001,0.001);
h_res_el_pre_clones.at(1)=new TH1D("res1_el_pre_clones", "(the_rec-the_true), theta_el[1,2]mrad with clones",100,-0.001,0.001);
h_res_el_pre_clones.at(2)=new TH1D("res2_el_pre_clones", "(the_rec-the_true), theta_el[2,3]mrad with clones",100,-0.001,0.001);
h_res_el_pre_clones.at(3)=new TH1D("res3_el_pre_clones", "(the_rec-the_true), theta_el[3,4]mrad with clones",100,-0.001,0.001);
h_res_el_pre_clones.at(4)=new TH1D("res4_el_pre_clones", "(the_rec-the_true), theta_el[4,5]mrad with clones",200,-0.002,0.002);
h_res_el_pre_clones.at(5)=new TH1D("res5_el_pre_clones", "(the_rec-the_true), theta_el[5,6]mrad with clones",200,-0.002,0.002);
h_res_el_pre_clones.at(6)=new TH1D("res6_el_pre_clones", "(the_rec-the_true), theta_el[6,7]mrad with clones",200,-0.002,0.002);
h_res_el_pre_clones.at(7)=new TH1D("res7_el_pre_clones", "(the_rec-the_true), theta_el[7,8]mrad with clones",200,-0.002,0.002);
h_res_el_pre_clones.at(8)=new TH1D("res8_el_pre_clones", "(the_rec-the_true), theta_el[8,9]mrad with clones",200,-0.002,0.002);
h_res_el_pre_clones.at(9)=new TH1D("res9_el_pre_clones", "(the_rec-the_true), theta_el[9,10]mrad with clones",200,-0.002,0.002);
h_res_el_pre_clones.at(10)=new TH1D("res10_el_pre_clones", "(the_rec-the_true), theta_el[10,15]mrad with clones",200,-0.01,0.01);
h_res_el_pre_clones.at(11)=new TH1D("res11_el_pre_clones", "(the_rec-the_true), theta_el[15,20]mrad with clones",200,-0.02,0.02);
h_res_el_pre_clones.at(12)=new TH1D("res12_el_pre_clones", "(the_rec-the_true), theta_el[20,25]mrad with clones",200,-0.02,0.02);
h_res_el_pre_clones.at(13)=new TH1D("res13_el_pre_clones", "(the_rec-the_true), theta_el[25,32]mrad with clones",200,-0.02,0.02);

h_res_mu_pre_clones.at(0)=new TH1D("res0_mu_pre_clones", "(thmu_rec-thmu_true), theta_mu[0.,0.1]mrad with clones",120,-0.0004,0.0004);
h_res_mu_pre_clones.at(1)=new TH1D("res1_mu_pre_clones", "(thmu_rec-thmu_true), theta_mu[0.1,0.2]mrad with clones",120,-0.0004,0.0004);
h_res_mu_pre_clones.at(2)=new TH1D("res2_mu_pre_clones", "(thmu_rec-thmu_true), theta_mu[0.2,0.3]mrad with clones",120,-0.0004,0.0004);
h_res_mu_pre_clones.at(3)=new TH1D("res3_mu_pre_clones", "(thmu_rec-thmu_true), theta_mu[0.3,0.4]mrad with clones",120,-0.0004,0.0004);
h_res_mu_pre_clones.at(4)=new TH1D("res4_mu_pre_clones", "(thmu_rec-thmu_true), theta_mu[0.4,0.5]mrad with clones",120,-0.0004,0.0004);
h_res_mu_pre_clones.at(5)=new TH1D("res5_mu_pre_clones", "(thmu_rec-thmu_true), theta_mu[0.5,0.6]mrad with clones",120,-0.0004,0.0004);
h_res_mu_pre_clones.at(6)=new TH1D("res6_mu_pre_clones", "(thmu_rec-thmu_true), theta_mu[0.6,0.7]mrad with clones",120,-0.0004,0.0004);
h_res_mu_pre_clones.at(7)=new TH1D("res7_mu_pre_clones", "(thmu_rec-thmu_true), theta_mu[0.7,0.8]mrad with clones",120,-0.0004,0.0004);
h_res_mu_pre_clones.at(8)=new TH1D("res8_mu_pre_clones", "(thmu_rec-thmu_true), theta_mu[0.8,0.9]mrad with clones",120,-0.0004,0.0004);
h_res_mu_pre_clones.at(9)=new TH1D("res9_mu_pre_clones", "(thmu_rec-thmu_true), theta_mu[0.9,1]mrad with clones",120,-0.0004,0.0004);
h_res_mu_pre_clones.at(10)=new TH1D("res10_mu_pre_clones", "(thmu_rec-thmu_true), theta_mu[1,2]mrad with clones",120,-0.0004,0.0004);
h_res_mu_pre_clones.at(11)=new TH1D("res11_mu_pre_clones", "(thmu_rec-thmu_true), theta_mu[2,3]mrad with clones",120,-0.0004,0.0004);
h_res_mu_pre_clones.at(12)=new TH1D("res12_mu_pre_clones", "(thmu_rec-thmu_true), theta_mu[3,4]mrad with clones",120,-0.0004,0.0004);
h_res_mu_pre_clones.at(13)=new TH1D("res13_mu_pre_clones", "(thmu_rec-thmu_true), theta_mu[4,5]mrad with clones",120,-0.0004,0.0004);



TH1D *h_res_vrtx_muin=new TH1D("h_res_vrtx_muin", "(thmu_in_rec-thmu_in_true) ",120,-0.0004,0.0004);

	double r_wnorm[6]={20.786765103274643,33.313091221576336,42.396733790329876,50.584815206143652,61.828110400735824,106.88513370134392};




for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
                cbmsim->GetEntry(i);
                if(i%10000 == 0) cout<<"Entry "<<i<<endl;

int code_mu_in=-99;
int code_e=-99;
int code_mu=-99;

        TVector3 p_muin_MC;
        TVector3 p_mu_MC;
        TVector3 p_e_MC;
        Double_t the_gen, thmu_gen,theX_gen,theY_gen,thmuX_gen,thmuY_gen;

           int hit_modXmu=0; int hit_modXe=0;
           int hit_modYmu=0; int hit_modYe=0;
           int stereo_mu=0; int stereo_e=0;
     	   int hit_modXmuin=0;
	   int hit_modYmuin=0;
	   int stereo_muin=0;

// Checking if in the MCTracks container there are elastic particles and if they are in acceptance (4 hits in X modules (== 2 stubs, as TrackerPoints give the hit per sensor) + 4 hits in Y modules + at least 2 hits in UV)

        for(int n = 0; n < MCTrack->GetEntries(); n++) {
         const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
         if(MCTr->interactionID()==0 and MCTr->pdgCode()==-13) {code_mu_in=n; p_muin_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); p_muin_MC.Unit();
        for(int s=0; s<TrackerPoints->GetEntries(); s++)
                         {const MUonETrackerPoint *TrackerPt = static_cast<const MUonETrackerPoint*>(TrackerPoints->At(s));
         if(TrackerPt->trackPDGCode()==-13 and TrackerPt->trackID()==n and TrackerPt->stationID()==0){
                                                                                                 if(TrackerPt->moduleID()==0 or TrackerPt->moduleID()==4) hit_modXmuin++;
                                                                                                 else if(TrackerPt->moduleID()==1 or TrackerPt->moduleID()==5) hit_modYmuin++;
                                                                                                 else if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_muin++;}
                }
	}
         if(MCTr->interactionID()==45 and MCTr->pdgCode()==11) {code_e=n; p_e_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); p_e_MC.Unit(); the_gen=p_muin_MC.Angle(p_e_MC);
                                                                                 theX_gen=MCTr->ax();
                                                                                 theY_gen=MCTr->ay();

        for(int s=0; s<TrackerPoints->GetEntries(); s++)
                         {const MUonETrackerPoint *TrackerPt = static_cast<const MUonETrackerPoint*>(TrackerPoints->At(s));
         if(TrackerPt->trackPDGCode()==11 and TrackerPt->trackID()==code_e and TrackerPt->stationID()==1){
                                                                                                 if(TrackerPt->moduleID()==0 or TrackerPt->moduleID()==4) hit_modXe++;
                                                                                                 else if(TrackerPt->moduleID()==1 or TrackerPt->moduleID()==5) hit_modYe++;
                                                                                                 else if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_e++;}
        	}
	 }
         if(MCTr->interactionID()==45 and MCTr->pdgCode()==-13) {code_mu=n; p_mu_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); p_mu_MC.Unit(); thmu_gen=p_muin_MC.Angle(p_mu_MC);
                                                                                 thmuX_gen=MCTr->ax();
                                                                                 thmuY_gen=MCTr->ay();
        for(int s=0; s<TrackerPoints->GetEntries(); s++)
                         {const MUonETrackerPoint *TrackerPt = static_cast<const MUonETrackerPoint*>(TrackerPoints->At(s));
	 if(TrackerPt->trackPDGCode()==-13 and TrackerPt->trackID()==code_mu and TrackerPt->stationID()==1){
                                                                                                 if(TrackerPt->moduleID()==0 or TrackerPt->moduleID()==4) hit_modXmu++;
                                                                                                 else if(TrackerPt->moduleID()==1 or TrackerPt->moduleID()==5) hit_modYmu++;
                                                                                                 else if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_mu++;}
		}
	 }
	}

// Look at reconstruction if events are reconstructible (all three particles with necessary hits to be potentially reconstructed)

	if(code_mu_in!=-99 and code_e!=-99 and code_mu!=-99 and hit_modXmuin==4 and hit_modYmuin==4 and stereo_muin>1 and hit_modXe==4 and hit_modYe==4 and stereo_e>1 and hit_modXmu==4 and hit_modYmu==4 and stereo_mu>1){


double wnorm=0.;

// wnorm is a number needed for normalization when we use different mesmer sample together (like in this case, 6 sample in differen kinematic region of the electron). We also register an idx for plots dipendently on the mu/el angle for the single event

int idx=99; int idx_mu=99;

if(the_gen>=0. and the_gen<0.001){idx=0;}
if(the_gen>=0.001 and the_gen<0.002){idx=1;}
if(the_gen>=0.002 and the_gen<0.003){idx=2;}
if(the_gen>=0.003 and the_gen<0.004){idx=3;}
if(the_gen>=0.004 and the_gen<0.005){idx=4;}
if(the_gen>=0.005 and the_gen<0.006){idx=5;}
if(the_gen>=0.006 and the_gen<0.007){idx=6;}
if(the_gen>=0.007 and the_gen<0.008){idx=7;}
if(the_gen>=0.008 and the_gen<0.009){idx=8;}
if(the_gen>=0.009 and the_gen<0.010){idx=9;}
if(the_gen>=0.010 and the_gen<0.015){idx=10;}
if(the_gen>=0.015 and the_gen<0.020){idx=11;}
if(the_gen>=0.020 and the_gen<0.025){idx=12;}
if(the_gen>=0.025 and the_gen<=0.033){idx=13;}



if(thmu_gen>=0. and thmu_gen<0.0001)idx_mu=0;
if(thmu_gen>=0.0001 and thmu_gen<0.0002)idx_mu=1;
if(thmu_gen>=0.0002 and thmu_gen<0.0003)idx_mu=2;
if(thmu_gen>=0.0003 and thmu_gen<0.0004)idx_mu=3;
if(thmu_gen>=0.0004 and thmu_gen<0.0005)idx_mu=4;
if(thmu_gen>=0.0005 and thmu_gen<0.0006)idx_mu=5;
if(thmu_gen>=0.0006 and thmu_gen<0.0007)idx_mu=6;
if(thmu_gen>=0.0007 and thmu_gen<0.0008)idx_mu=7;
if(thmu_gen>=0.0008 and thmu_gen<0.0009)idx_mu=8;
if(thmu_gen>=0.0009 and thmu_gen<0.001)idx_mu=9;
if(thmu_gen>=0.001 and thmu_gen<0.002)idx_mu=10;
if(thmu_gen>=0.002 and thmu_gen<0.003)idx_mu=11;
if(thmu_gen>=0.003 and thmu_gen<0.004)idx_mu=12;
if(thmu_gen>=0.004 and thmu_gen<=0.005)idx_mu=13;


if(index==0){wnorm=r_wnorm[0];}
else if(index==1){wnorm=r_wnorm[1];}
else if(index==2){wnorm=r_wnorm[2];}
else if(index==3){wnorm=r_wnorm[3];}
else if(index==4){wnorm=r_wnorm[4];}
else if(index==5){wnorm=r_wnorm[5];}


Int_t yes2=0;
Double_t th_inx,th_iny,x0_in,y0_in;
Double_t chi2_muin=999.;
Int_t stubs_muin=0;
Int_t e=0;
Int_t mu=0;
TVector3 p_muin,p_e,p_mu;
Double_t the_rec,theX_rec,theY_rec,thmu_rec,thmuX_rec,thmuY_rec,the_rec_vrtx,thmu_rec_vrtx;

// Define best vertex and its chi2

MUonERecoOutputVertex vrtx = ReconstructionOutput->bestVertex();
Double_t chi=vrtx.chi2perDegreeOfFreedom();


// Look at reconstructed track container to see if outgoing mu and e are reconstructed
vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
int other=0;
int mu_in=0;
vector<double> quality_e; quality_e.reserve(5);
vector<double> quality_mu; quality_mu.reserve(5);

         for (auto&& track : tracks) {

        if(code_mu_in==track.linkedTrackID() and track.sector()==0){
	mu_in++;
         th_inx=track.xSlope();
         th_iny=track.ySlope();
         x0_in=track.x0();
         y0_in=track.y0();
         chi2_muin=track.chi2perDegreeOfFreedom();
         stubs_muin=track.hits().size();
        p_muin.SetXYZ(th_inx,th_iny,1.0);
        p_muin=p_muin.Unit();
                        }
        else if(track.processIDofLinkedTrack()==45 and track.sector()==1)
                {
                 if(code_e==track.linkedTrackID()) {yes2++; e++; theX_rec=track.xSlope(); theY_rec=track.ySlope();p_e.SetXYZ(theX_rec,theY_rec,1.0); p_e=p_e.Unit();
							if(track.fractionOfHitsSharedWithLinkedTrack()>=0.65)the_rec=p_e.Angle(p_muin);
                                                    quality_e.push_back(track.fractionOfHitsSharedWithLinkedTrack());}

                 if(code_mu==track.linkedTrackID()) {yes2++; mu++;thmuX_rec=track.xSlope(); thmuY_rec=track.ySlope();p_mu.SetXYZ(thmuX_rec,thmuY_rec,1.0); p_mu=p_mu.Unit();
							if(track.fractionOfHitsSharedWithLinkedTrack()>=0.65)thmu_rec=p_mu.Angle(p_muin);
                                                    quality_mu.push_back(track.fractionOfHitsSharedWithLinkedTrack());}
                }
        else if(track.processIDofLinkedTrack()!=45 and track.sector()==1){other++;}

         }//end of for cycle


// Check that the incoming muon is well reconstructed
if(stubs_muin>=5 and chi2_muin<5){

			// Evaluate resolution incoming muon
                        h_res_vrtx_muin->Fill(p_muin.Theta()-p_muin_MC.Theta(),MesmerEvent->wgt_LO*wnorm);

// Check that within reconstructed e/mu there are good quality tracks PRE-VRTX
if(find_if(quality_e.begin(),quality_e.end(),[](double i){return i>=0.65;})!=end(quality_e) and find_if(quality_mu.begin(),quality_mu.end(),[](double i){return i>=0.65;})!=end(quality_mu)){


// If 1 e and 1 mu are reconstructed with also a vrtx (chi!=0 -> if vrtx is not reconstructed chi2==0) are present, then plot resolution PRE-VRTX (without clones)
	if(chi!=0 and e==1 and mu==1){
                        h_res_el_pre.at(idx)->Fill(the_rec-the_gen,MesmerEvent->wgt_LO*wnorm);
                        h_res_mu_pre.at(idx_mu)->Fill(thmu_rec-thmu_gen,MesmerEvent->wgt_LO*wnorm);
			}
// If at least 1 mu and at least 1 el with also a vrtx (chi!=0 -> if vrtx is not reconstructed chi2==0) are present, then plot resolution PRE-VRTX (including clones)
	if(chi!=0 and e>=1 and mu>=1){
                        h_res_el_pre_clones.at(idx)->Fill(the_rec-the_gen,MesmerEvent->wgt_LO*wnorm);
                        h_res_mu_pre_clones.at(idx_mu)->Fill(thmu_rec-thmu_gen,MesmerEvent->wgt_LO*wnorm);}

}

// Check that the e/mu tracks associated to the reconstructed vrtx are good quality tracks POST-VRTX
if(vrtx.outgoingMuon().fractionOfHitsSharedWithLinkedTrack()>=0.65 and vrtx.outgoingElectron().fractionOfHitsSharedWithLinkedTrack()>=0.65){

//Evaluate resolution of post-vrtx track when just 1 mu and 1 el are reconstructed (without clones). Use MC truth to associate correctly the track (do correct PID and correcting wrong fairmu PID -happening for small angles-)
	if(chi!=0 and e==1 and mu==1){
			if(vrtx.outgoingMuon().linkedTrackID()==code_mu)h_res_vrtx_mu.at(idx_mu)->Fill(vrtx.muonTheta()-thmu_gen,MesmerEvent->wgt_LO*wnorm);
			if(vrtx.outgoingMuon().linkedTrackID()==code_e)h_res_vrtx_el.at(idx)->Fill(vrtx.muonTheta()-the_gen,MesmerEvent->wgt_LO*wnorm);
                        if(vrtx.outgoingElectron().linkedTrackID()==code_mu)h_res_vrtx_mu.at(idx_mu)->Fill(vrtx.electronTheta()-thmu_gen,MesmerEvent->wgt_LO*wnorm);
                        if(vrtx.outgoingElectron().linkedTrackID()==code_e)h_res_vrtx_el.at(idx)->Fill(vrtx.electronTheta()-the_gen,MesmerEvent->wgt_LO*wnorm);}

//Evaluate resolution of post-vrtx track when at least 1 mu and at least 1 el are reconstructed (including clones) Use MC truth to associate correctly the track (do correct PID and correcting wrong fairmu PID -happening for small angles-)
	if(chi!=0 and e>=1 and mu>=1){
                        if(vrtx.outgoingMuon().linkedTrackID()==code_mu)h_res_vrtx_mu_clones.at(idx_mu)->Fill(vrtx.muonTheta()-thmu_gen,MesmerEvent->wgt_LO*wnorm);
                        if(vrtx.outgoingMuon().linkedTrackID()==code_e)h_res_vrtx_el_clones.at(idx)->Fill(vrtx.muonTheta()-the_gen,MesmerEvent->wgt_LO*wnorm);
                        if(vrtx.outgoingElectron().linkedTrackID()==code_mu)h_res_vrtx_mu_clones.at(idx_mu)->Fill(vrtx.electronTheta()-thmu_gen,MesmerEvent->wgt_LO*wnorm);
                        if(vrtx.outgoingElectron().linkedTrackID()==code_e)h_res_vrtx_el_clones.at(idx)->Fill(vrtx.electronTheta()-the_gen,MesmerEvent->wgt_LO*wnorm);}
			}//if quality vrtx


                }//if mu_in
        }//if generated
}//for


TCanvas n0("n0","n0",2100,3500);
n0.Divide(3,5);
for(int m=0; m<NBINS; m++){
n0.cd(m+1);
h_res_vrtx_el.at(m)->SetMinimum(1.);
h_res_vrtx_el.at(m)->Draw("hist");
h_res_el_pre.at(m)->SetLineColor(kRed+1);
h_res_el_pre.at(m)->Draw("hist same");
h_res_el_pre.at(m)->SaveAs(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/res_prevrtx_el_%d_%dhit_LO_index%d_allNewFeatures.root",static_cast<char>(m), static_cast<char>(nhits), static_cast<char>(index)));
h_res_vrtx_el.at(m)->SaveAs(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/res_vrtx_el_%d_%dhit_LO_index%d_allNewFeatures.root",static_cast<char>(m), static_cast<char>(nhits), static_cast<char>(index)));
//gPad->SetLogy();
}
n0.SaveAs(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/h_res_pre_el_%dhit_LO_index%d_allNewFeatures.pdf",static_cast<char>(nhits),static_cast<char>(index)));


TCanvas n1("n1","n1",2100,3500);
n1.Divide(3,5);
for(int m=0; m<NBINS; m++){
n1.cd(m+1);
h_res_vrtx_el_clones.at(m)->Draw("hist");
h_res_el_pre_clones.at(m)->SetLineColor(kRed+1);
h_res_el_pre_clones.at(m)->Draw("hist same");
h_res_el_pre_clones.at(m)->SaveAs(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/res_prevrtx_clones_el_%d_%dhit_LO_index%d_allNewFeatures.root",static_cast<char>(m), static_cast<char>(nhits), static_cast<char>(index)));
h_res_vrtx_el_clones.at(m)->SaveAs(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/res_vrtx_clones_el_%d_%dhit_LO_index%d_allNewFeatures.root",static_cast<char>(m), static_cast<char>(nhits), static_cast<char>(index)));
}
n1.SaveAs(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/h_res_vrtx_el_clones_%dhitLO_index%d_allNewFeatures.pdf",static_cast<char>(nhits),static_cast<char>(index)));


TCanvas n2("n2","n2",2100,3500);
n2.Divide(3,5);
for(int m=0; m<NBINS_mu; m++){
n2.cd(m+1);
h_res_vrtx_mu.at(m)->SetMinimum(1.);
h_res_vrtx_mu.at(m)->Draw("hist");
h_res_mu_pre.at(m)->SetLineColor(kRed+1);
h_res_mu_pre.at(m)->Draw("hist same");
h_res_mu_pre.at(m)->SaveAs(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/res_prevrtx_AngleMu_%d_%dhit_LO_index%d_allNewFeatures.root",static_cast<char>(m), static_cast<char>(nhits), static_cast<char>(index)));
h_res_vrtx_mu.at(m)->SaveAs(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/res_vrtx_AngleMu_%d_%dhit_LO_index%d_allNewFeatures.root",static_cast<char>(m), static_cast<char>(nhits), static_cast<char>(index)));
//gPad->SetLogy();
}
n2.SaveAs(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/h_res_pre_mu_%dhit_LO_index%d_allNewFeatures.pdf",static_cast<char>(nhits),static_cast<char>(index)));

TCanvas n3("n3","n3",2100,3500);
n3.Divide(3,5);
for(int m=0; m<NBINS_mu; m++){
n3.cd(m+1);
h_res_vrtx_mu_clones.at(m)->Draw("hist");
h_res_mu_pre_clones.at(m)->SetLineColor(kRed+1);
h_res_mu_pre_clones.at(m)->Draw("hist same");
h_res_mu_pre_clones.at(m)->SaveAs(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/res_prevrtx_clones_AngleMu_%d_%dhit_LO_index%d_allNewFeatures.root",static_cast<char>(m),static_cast<char>(nhits),static_cast<char>(index)));
h_res_vrtx_mu_clones.at(m)->SaveAs(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/res_vrtx_clones_AngleMu_%d_%dhit_LO_index%d_allNewFeatures.root",static_cast<char>(m),static_cast<char>(nhits),static_cast<char>(index)));
}
n3.SaveAs(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/h_res_vrtx_mu_clones_%dhit_LO_index%d_allNewFeatures.pdf",static_cast<char>(nhits),static_cast<char>(index)));

TCanvas in("in","in",700,700);
h_res_vrtx_muin->Draw();
in.SaveAs(Form("/home/espedica/macros_fairmu/clean_codes/separate/validation/h_res_pmuin_%dhit_LO_index%d_allNewFeatures.pdf",static_cast<char>(nhits),static_cast<char>(index)));




}



