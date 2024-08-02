
#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TStyle.h>


void single_divide(){

TFile *f_RD = TFile::Open("d_eff_RD_prev.root");
TFile *f_MC = TFile::Open("d_eff_MC_prev.root");

TFile *f_the_MC = TFile::Open("theta_e_MC_prev.root");
TFile *f_thmu_MC = TFile::Open("theta_mu_MC_prev.root");

TFile *f_thmu_RD = TFile::Open("theta_mu_RD_prev.root");
TFile *f_the_RD = TFile::Open("theta_e_RD_prev.root");

/*TFile *f_2D_RD = TFile::Open("comparison_RDMC/2D_RD_prev.root");
TFile *f_2D_MC = TFile::Open("comparison_RDMC/2D_MC_prev.root");*/

TH1::SetDefaultSumw2(kTRUE);


TH1D *h_RD=(TH1D*)f_RD->Get("d_eff_real");
//TH1D * h1 = (TH1D*) h_RD->Clone();
//h1->Divide(h_MC);



TH1D *h_theta_e_MC=(TH1D*)f_the_MC->Get("theta_e");
TH1D *h_theta_mu_MC=(TH1D*)f_thmu_MC->Get("theta_mu");


TH1D *h_theta_e_RD=(TH1D*)f_the_RD->Get("theta_e");
TH1D *h_theta_mu_RD=(TH1D*)f_thmu_RD->Get("theta_mu");

/*TH2D *h_2D_RD=(TH2D*)f_2D_RD->Get("h2D");
TH2D *h_2D_MC=(TH2D*)f_2D_MC->Get("h2D");*/



/*h_theta_e_RD->SetBins(8,0.,0.032);
h_theta_e_MC->SetBins(8,0.,0.032);

h_theta_mu_RD->SetBins(9,0.0005,0.0029);
h_theta_mu_MC->SetBins(9,0.0005,0.0029);
*/

TH1::SetDefaultSumw2(kTRUE);

   Double_t cross_section = 1338.66391322781;//316.39947555367;
//   Double_t entries = 679760;//670327;//261456;
   Double_t entries = 597766;//prev

TCanvas a1("a1","a1",1000,700);
a1.Divide(2,2);
a1.cd(1);
cout << "theta el. reco:" << endl;
h_theta_e_MC->Scale(9.3792854e+08* 5.5*1E+23*3.*1E-30 * cross_section/entries);
h_theta_e_MC->SetMinimum(1.);
h_theta_e_MC->Draw("E");
h_theta_e_RD->SetLineColor(kPink+10);
h_theta_e_RD->Draw("E same");
gPad->SetLogy();
gStyle->SetOptStat(0);

a1.cd(2);

cout << "theta mu. reco:" << endl;
h_theta_mu_MC->Scale(9.3792854e+08* 5.5*1E+23*3.*1E-30 * cross_section/entries);
h_theta_mu_MC->SetMinimum(1.);
h_theta_mu_MC->GetXaxis()->SetRangeUser(0.,0.0025);
h_theta_mu_RD->GetXaxis()->SetRangeUser(0.,0.0025);
h_theta_mu_MC->Draw("E");
h_theta_mu_RD->SetLineColor(kPink+10);
h_theta_mu_RD->Draw("E same");
//gPad->SetLogy();
gStyle->SetOptStat(0);


TH1D * h2 = (TH1D*) h_theta_mu_RD->Clone();
h2->Divide(h_theta_mu_MC);
TH1D * h3 = (TH1D*) h_theta_e_RD->Clone();
h3->Divide(h_theta_e_MC);
a1.cd(3);
h3->SetMinimum(0.);
h3->Draw("E");
a1.cd(4);
h2->SetMinimum(0.);
h2->Draw("E");
//gPad->SetLogy();
a1.SaveAs("comparison_RDMC/theta_MC_reco_RD_prev.pdf");
/*
TCanvas a2("a2","a2",2800,4800);
a2.Divide(1,3);
a2.cd(1);
h_2D_MC->SetTitle("#theta el -vs- #theta mu MC");
h_2D_MC->Draw("COLZ");
a2.cd(2);
h_2D_RD->SetTitle("#theta el -vs- #theta mu RD");
h_2D_RD->Draw("COLZ");
a2.cd(3);


TH2D * h5 = (TH2D*) h_2D_RD->Clone();
h5->Divide(h_2D_MC);


Int_t nx = h5->GetNbinsX();
Int_t ny = h5->GetNbinsY();
for (Int_t i=1; i<nx+1; i++) {
for (Int_t j=1; j<ny+1; j++) {
if (h5->GetBinError(i,j)/h5->GetBinContent(i,j)*100>=5) h5->SetBinContent(i,j,0);}}
  Int_t n=13;
  Int_t *colors = new Int_t[n];
//  for (Int_t i = 0; i < n; i++){if(i<5) colors[i] = kRed+i; else if(i>=5 and i<10) colors[i] = kOrange+i-5; else colors[i] = kYellow+i;}
  colors[0] = kRed-10;
  colors[1] = kRed-10;
  colors[2] = kRed-10;
  colors[3] = kOrange+6;
  colors[4] = kOrange+6;
  colors[5] = kOrange+6;
  colors[6] = kOrange+10;
  colors[7] = kOrange+10;
  colors[8] = kOrange+10;
  colors[9] = kViolet-9;
  colors[10] = kViolet-9;
  colors[11] = kAzure+6;
  colors[12] = kAzure;

  gStyle->SetPalette(n, colors);

  // #levels <= #colors + 1 (notes: +-3.4e38 = +-FLT_MAX; +1.17e-38 = +FLT_MIN)
  Double_t levels[] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.00, 10.0, 10., 30., 60};

h5->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
// gStyle->SetPalette(kIsland);
h5->SetTitle("#theta el -vs- #theta mu ratio RD/MC");
gStyle->SetPaintTextFormat(".2f");
h5->Draw("COLZ TEXT");
a2.SaveAs("comparison_RDMC/comparison_2D_prev.pdf");

TCanvas a3("a3","a3",1400,700);
a3.Divide(2,1);
a3.cd(1);
auto p1 =h5->ProfileX();
p1->SetTitle("Efficiency as a function of #theta el");
p1->Draw();
a3.cd(2);
auto p2 =h5->ProfileY();
p2->SetTitle("Efficiency as a function of #theta mu");
p2->Draw();
a3.SaveAs("comparison_RDMC/profile_prev.pdf");
*/



}

