#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TStyle.h>


void divide(){

TH1::SetDefaultSumw2(kTRUE);
int const NBINS = 7;

TFile *f_RD = TFile::Open("comparison_RDMC/d_eff_RD_z.root");
TFile *f_MC = TFile::Open("comparison_RDMC/d_eff_MC_z.root");

array<TFile*,NBINS> f_thmu_MC,f_the_MC, f_2D_MC,f_op_MC;

for(int m=0; m<NBINS; m++)
{f_the_MC.at(m) = TFile::Open(Form("comparison_RDMC/theta_e_MC_%d_z.root",static_cast<char>(m)));
 f_thmu_MC.at(m) = TFile::Open(Form("comparison_RDMC/theta_mu_MC_%d_z.root",static_cast<char>(m)));
 f_2D_MC.at(m) = TFile::Open(Form("comparison_RDMC/2D_MC_%d_z.root",static_cast<char>(m))); 
 f_op_MC.at(m) = TFile::Open(Form("comparison_RDMC/opening_MC_%d_z.root",static_cast<char>(m)));}

TFile *f_thmu_RD = TFile::Open("comparison_RDMC/theta_mu_RD_z.root");
TFile *f_the_RD = TFile::Open("comparison_RDMC/theta_e_RD_z.root");
TFile *f_2D_RD = TFile::Open("comparison_RDMC/2D_RD_z.root");
TFile *f_op_RD = TFile::Open("comparison_RDMC/opening_RD_z.root");


TH1::SetDefaultSumw2(kTRUE);


TH1D *h_RD=(TH1D*)f_RD->Get("d_eff_real");
//TH1D * h1 = (TH1D*) h_RD->Clone();
//h1->Divide(h_MC);


array<TH1D*,NBINS> h_theta_mu_MC,h_theta_e_MC,h_opening_MC;
array<TH2D*,NBINS> h_2D_MC;

for(int m=0; m<NBINS; m++)
{
h_theta_e_MC.at(m)=(TH1D*)f_the_MC.at(m)->Get(Form("theta_e%d",static_cast<char>(m)));
h_theta_mu_MC.at(m)=(TH1D*)f_thmu_MC.at(m)->Get(Form("theta_mu%d",static_cast<char>(m)));
h_2D_MC.at(m)=(TH2D*)f_2D_MC.at(m)->Get(Form("h2D%d",static_cast<char>(m)));
h_opening_MC.at(m)=(TH1D*)f_op_MC.at(m)->Get(Form("h_opening%d",static_cast<char>(m)));
}

TH1D *h_theta_e_RD=(TH1D*)f_the_RD->Get("theta_e");
TH1D *h_theta_mu_RD=(TH1D*)f_thmu_RD->Get("theta_mu");
TH2D *h_2D_RD=(TH2D*)f_2D_RD->Get("h2D");
TH1D *h_opening_RD=(TH1D*)f_op_RD->Get("h_opening");

/*h_theta_e_RD->SetBins(8,0.,0.032);
h_theta_e_MC->SetBins(8,0.,0.032);

h_theta_mu_RD->SetBins(6,0.,0.003);
for(int m=0; m<NBINS; m++){h_theta_mu_MC.at(m)->SetBins(6,0.,0.003);}
*/
   Double_t cross_section[NBINS] = {19.93244038235,32.68221635672,40.34387346048,51.24449139646,62.78779293344,108.45556192659,1022.04019743297};

//   Double_t entries[NBINS] = {3.38196e+06,196868,132707,117330,100338,124103,498447}; //new sel
   Double_t entries[NBINS] = {3.02526e+06,164901,111828,97681.8,85351.9,105260,421551};
//   Double_t entries[NBINS] = {3.13392e+06,169007,115901,101213,90838.3,108235,435498}; //old sel

for(int m=0; m<6; m++){
	 h_theta_e_MC.at(m)->Scale(9.3792854e+08* 5.5*1E+23*3.*1E-30 * cross_section[m]/entries[m]);
         h_theta_mu_MC.at(m)->Scale(9.3792854e+08* 5.5*1E+23*3.*1E-30 * cross_section[m]/entries[m]);
	 h_2D_MC.at(m)->Scale(9.3792854e+08* 5.5*1E+23*3.*1E-30 * cross_section[m]/entries[m]);
	 h_opening_MC.at(m)->Scale(9.3792854e+08* 5.5*1E+23*3.*1E-30 * cross_section[m]/entries[m]);
}

TCanvas a("a","a",700,700);
for(int m=1; m<NBINS; m++){h_opening_MC.at(0)->Add(h_opening_MC.at(m));}
a.Divide(1,2);
a.cd(1);
h_opening_MC.at(0)->SetMinimum(1.);
//h_opening_MC.at(0)->Rebin(2);
//h_opening_RD->Rebin(2);
h_opening_RD->SetMinimum(1.);
h_opening_MC.at(0)->Draw("hist");
h_opening_RD->SetLineColor(kPink+10);
h_opening_RD->Draw("hist same");
a.cd(2);

TH1D * h = (TH1D*) h_opening_RD->Clone();
h->Divide(h_opening_MC.at(0));
//h->Rebin(2);
h->SetMaximum(1.);
h->SetMinimum(0.);
h->Draw("E");
a.SaveAs("comparison_RDMC/opening_reb_z.pdf");

TCanvas a1("a1","a1",1000,700);
a1.Divide(2,2);
a1.cd(1);
cout << "theta el. reco:" << endl;
//h_theta_e_MC->Scale(9.3792854e+08*5.5*1E+23*3.*1E-30*315.44638/6.79333e+07);//26563993,605675    137720.00    1.3242886e+09   6.79333e+07 9.3792854e+08
h_theta_e_MC.at(0)->SetMinimum(1.);
for(int m=1; m<NBINS; m++){h_theta_e_MC.at(0)->Add(h_theta_e_MC.at(m));}
h_theta_e_MC.at(0)->Draw("E");
h_theta_e_RD->SetLineColor(kPink+10);
h_theta_e_RD->Draw("E same");
gPad->SetLogy();
gStyle->SetOptStat(0);
a1.cd(2);
TH1D * h3 = (TH1D*) h_theta_e_RD->Clone();
h3->Divide(h_theta_e_MC.at(0));
h3->SetMinimum(0.);
h3->Draw("E");
a1.cd(3);
cout << "theta mu. reco:" << endl;
//h_theta_e_MC->Scale(9.3792854e+08*5.5*1E+23*3.*1E-30*315.44638/6.79333e+07);//26563993,605675    137720.00    1.3242886e+09   6.79333e+07 9.3792854e+08
h_theta_mu_MC.at(0)->SetMinimum(1.);
h_theta_mu_MC.at(0)->GetXaxis()->SetRangeUser(0.,0.0025);
for(int m=1; m<NBINS; m++){h_theta_mu_MC.at(0)->Add(h_theta_mu_MC.at(m));}
h_theta_mu_MC.at(0)->Draw("E");
h_theta_mu_RD->GetXaxis()->SetRangeUser(0.,0.0025);
h_theta_mu_RD->SetLineColor(kPink+10);
h_theta_mu_RD->Draw("E same");
gStyle->SetOptStat(0);
gPad->SetLogy();
a1.cd(4);
TH1D * h4 = (TH1D*) h_theta_mu_RD->Clone();
h4->Divide(h_theta_mu_MC.at(0));
h4->SetMinimum(0.);
h4->Draw("E");
a1.SaveAs("comparison_RDMC/theta_MC_reco_RD_add_z.pdf");


TCanvas a2("a2","a2",2800,4800);
a2.Divide(1,3);
for(int m=1; m<NBINS; m++){h_2D_MC.at(0)->Add(h_2D_MC.at(m));}

a2.cd(1);
h_2D_MC.at(0)->SetTitle("#theta el -vs- #theta mu MC");
//h_2D_MC.at(0)->RebinY(2);
//h_2D_MC.at(0)->RebinX(2);
h_2D_MC.at(0)->Draw("COLZ");
a2.cd(2);
h_2D_RD->SetTitle("#theta el -vs- #theta mu RD");
//h_2D_RD->RebinX(2);
//h_2D_RD->RebinY(2);
h_2D_RD->Draw("COLZ");
a2.cd(3);


TH2D * h5 = (TH2D*) h_2D_RD->Clone();
h5->Divide(h_2D_MC.at(0));


Int_t nx = h5->GetNbinsX();
Int_t ny = h5->GetNbinsY();
for (Int_t i=1; i<nx+1; i++) {
for (Int_t j=1; j<ny+1; j++) {
if (h5->GetBinError(i,j)/h5->GetBinContent(i,j)*100>=5) h5->SetBinContent(i,j,0);}}

for(int a=1; a<h5->GetNbinsX(); a++){
 for(int b=1; b<h5->GetNbinsY(); b++){
	cout << a<<"," <<b << " ) " << h5->GetBinContent(a,b) << " +- " << h5->GetBinError(a,b)<< " ----> " << h5->GetBinError(a,b)/h5->GetBinContent(a,b)*100 << "%" << endl;
 }
}

  //Int_t colors[] = {0, 1, 2, 3, 4, 5, 6}; // #colors >= #levels - 1
  //gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
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
a2.SaveAs("comparison_RDMC/comparison_2D_z.pdf");

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
a3.SaveAs("comparison_RDMC/profile_z.pdf");


}

