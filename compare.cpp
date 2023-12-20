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
#include "TMatrix.h"
#include "TSystemDirectory.h"
#include <TStyle.h>
#include <yaml-cpp/yaml.h>
#include "geometry.h"

using namespace std;


void RealDataAnalyzer(){

 	TFile *fin_6= new TFile("chi2_muin_6.root");
        TFile *fout_6= new TFile("chi2_muout_6.root");
        TFile *fin_12= new TFile("chi2_muin_12.root");
        TFile *fout_12= new TFile("chi2_muout_12.root");


TH1F * h_in_6 = new TH1F("h_in_6","chi2 mu_in 6+6" , 200, 0., 100.);
h_in_6 = (TH1F*)fin_6->Get("h_chi20");
TH1F * h_out_6 = new TH1F("h_out_6","chi2 mu_out 6+6" , 200, 0., 100.);
h_out_6 = (TH1F*)fout_6->Get("h_chi21");



TH1F * h_in_12 = new TH1F("h_in_12","chi2 mu_in 12" , 200, 0., 100.);
h_in_12 = (TH1F*)fin_12->Get("h_chi20");
TH1F * h_out_12 = new TH1F("h_out_12","chi2 mu_out 12" , 200, 0., 100.);
h_out_12 = (TH1F*)fout_12->Get("h_chi21");

TCanvas n1("n1","n1",700,700);
n1.Divide(1,2);
n1.cd(1);
h_in_6->Draw();
h_in_12->SetLineColor(kRed);
h_in_12->Draw("same");

n1.cd(2);
h_out_6->Draw();
h_out_12->SetLineColor(kRed);
h_out_12->Draw("same");

n1.SaveAs("compare.pdf");


}
