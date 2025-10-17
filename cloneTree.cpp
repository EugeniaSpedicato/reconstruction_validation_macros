#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TStyle.h>



void cloneTree()
{

TFile *f=new TFile("/mnt/raid10/DATA/espedica/fairmu/gen_digi/minbias_SIM-DIGI_8Mevents.root","RECREATE");


TFile *inputfile = new TFile("/eos/user/e/ehess/tr_1.7.24/tr2023_promptanalysis/GiovanniA/job_24-09-27_20:26:44_master_8Mevents_gnn1/outputs/run_100/skim_interactions_singlemu_MC_run_100_inputs.root");
   TTree *cbmsim = (TTree*)inputfile->Get("cbmsim");

/*TChain * cbmsim = new TChain("cbmsim");
cbmsim->Add("/eos/user/e/ehess/tr_1.7.24/tr2023_promptanalysis/GiovanniA/job_24-08-01_14:57:47_master_14.7.24_10Mevents/outputs/run_100/skim_interactions_singlemu_MC_run_100_inputs_master.root");
cbmsim->Add("/eos/user/e/ehess/tr_1.7.24/tr2023_promptanalysis/GiovanniA/job_24-09-25_18:37:27_master_11Mevents_gnn1/outputs/run_100/skim_interactions_singlemu_MC_run_100_inputs.root");
cbmsim->Add("/eos/user/e/ehess/tr_1.7.24/tr2023_promptanalysis/GiovanniA/job_24-09-27_20:26:44_master_8Mevents_gnn1/outputs/run_100/skim_interactions_singlemu_MC_run_100_inputs.root");
*/
/*cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_5-10mrad_1M_2hitFirstModules_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_10-15mrad_1M_2hitFirstModules_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_15-20mrad_1M_2hitFirstModules_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_20-25mrad_1M_2hitFirstModules_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_25-32mrad_1M_2hitFirstModules_NOoutchi2_1M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/theta_32-inf_mrad_1M_2hitFirstModules_NOoutchi2_reassign.root");*/

//cbmsim->SetBranchStatus("MCTracks", 0);
//cbmsim->SetBranchStatus("TrackerStubs", 0);
cbmsim->SetBranchStatus("Bend*", 0);
cbmsim->SetBranchStatus("Bx*", 0);
cbmsim->SetBranchStatus("Link*", 0);
cbmsim->SetBranchStatus("SuperID*", 0);
cbmsim->SetBranchStatus("LocalX*", 0);
cbmsim->SetBranchStatus("LocalY*", 0);

//cbmsim->SetBranchStatus("BestVertex*", 0);
//cbmsim->SetBranchStatus("AdaptiveFitter*", 0);
//cbmsim->SetBranchStatus("IsReconstructed", 0);
//cbmsim->SetBranchStatus("SourceEventNumber", 0);
//cbmsim->SetBranchStatus("IsMC", 0);
//cbmsim->SetBranchStatus("TotalEventEnergy", 0);
//cbmsim->SetBranchStatus("BestVertex", 0);
//cbmsim->SetBranchStatus("ReconstructionOutput", 0);

        f->cd();
        TList* tl = (TList*)inputfile->Get("BranchList");
        tl->Write("BranchList", 1);


   auto newtree = cbmsim->CloneTree();
   newtree->Print();
   newtree->Write();
   f->Close();

}
