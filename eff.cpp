
#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TStyle.h>


void eff(){
TFile *f_thmu_NLO = TFile::Open("GM_NLOvsLO/50outchi2_mu_eff_NLO.root");
TFile *f_thmu_LO = TFile::Open("GM_NLOvsLO/50outchi2_mu_eff_LO.root");

TFile *f_thmu_NLOc = TFile::Open("GM_NLOvsLO/50outchi2_mu_eff_NLO_02cut.root");
TFile *f_thmu_LOc = TFile::Open("GM_NLOvsLO/50outchi2_mu_eff_LO_02cut.root");


TH1::SetDefaultSumw2(kTRUE);

TH1D *h_theta_mu_NLO=(TH1D*)f_thmu_NLO->Get("theta_mu");
TH1D *h_theta_mu_LO=(TH1D*)f_thmu_LO->Get("theta_mu");
TH1D *h_theta_mu_NLOc=(TH1D*)f_thmu_NLOc->Get("theta_mu");
TH1D *h_theta_mu_LOc=(TH1D*)f_thmu_LOc->Get("theta_mu");

cout << "50outchi2_mu_eff_NLO.root " << endl;
for(int b=1; b<11; b++){
cout << b << " ) " << h_theta_mu_NLO->GetBinContent(b) << " +- " << h_theta_mu_NLO->GetBinError(b)<< " ----> " << h_theta_mu_NLO->GetBinError(b)/h_theta_mu_NLO->GetBinContent(b)*100 << "%" << endl;
}
cout << "50outchi2_mu_eff_LO.root " << endl;
for(int b=1; b<11; b++){
cout << b << " ) " << h_theta_mu_LO->GetBinContent(b) << " +- " << h_theta_mu_LO->GetBinError(b)<< " ---> " << h_theta_mu_LO->GetBinError(b)/h_theta_mu_LO->GetBinContent(b)*100 << "%" <<endl ;
}

cout << "50outchi2_mu_eff_NLO_02cut.root " << endl;
for(int b=1; b<11; b++){
cout << b << " ) " << h_theta_mu_NLOc->GetBinContent(b) << " +- " << h_theta_mu_NLOc->GetBinError(b)<< " ----> " << h_theta_mu_NLOc->GetBinError(b)/h_theta_mu_NLOc->GetBinContent(b)*100 << "%" << endl;
}
cout << "50outchi2_mu_eff_LO_02cut.root " << endl;
for(int b=1; b<11; b++){
cout << b << " ) " << h_theta_mu_LOc->GetBinContent(b) << " +- " << h_theta_mu_LOc->GetBinError(b)<< " ---> " << h_theta_mu_LOc->GetBinError(b)/h_theta_mu_LOc->GetBinContent(b)*100 << "%" <<endl ;
}


}
