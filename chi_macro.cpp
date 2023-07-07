#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TSystemDirectory.h"
#include <TStyle.h>

using namespace std;

void chi_macro()
{

TFile s("vrtx_chi_sig02gev_norm_ovrf.root");
TFile s2("vrtx_chi_sig2gev_norm_ovrf.root");
TFile p("vrtx_chi_PP_norm_ovrflw.root");
TFile sm("vrtx_chi_sigminb_norm_ovrflw.root");

TH1F *sig_mes_02gev=(TH1F*)s.Get("htmp_s");
TH1F *sig_mes_2gev=(TH1F*)s2.Get("htmp_s");

sig_mes_02gev->SetTitle("Vertex quality signal MESMER  Ee>0.2GeV (blue) and bkg minbias (red) with acceptance cut th_mu<35mrad");
sig_mes_2gev->SetTitle("Vertex quality signal MESMER  Ee>2GeV (blue) and bkg minbias (red) with acceptance cut th_mu<35mrad");

TH1F *bkg_min=(TH1F*)p.Get("h1");
TH1F *sig_min=(TH1F*)sm.Get("h1_sm_o");
sig_min->SetTitle("Vertex quality signal minbias (blue) and bkg minbias (red) with acceptance cut th_mu<35mrad");

TCanvas c1("c1","c1",700,700);
c1.Divide(1,3);
c1.cd(1);
sig_mes_02gev->SetMaximum(1);
sig_mes_02gev->Draw("hist");
bkg_min->SetLineColor(kRed);
bkg_min->Draw("hist same");
gPad->SetLogy();
gStyle->SetOptStat(0);
c1.cd(2);
sig_mes_2gev->SetMaximum(1);
sig_mes_2gev->Draw("hist");
bkg_min->SetLineColor(kRed);
bkg_min->Draw("hist same");
gPad->SetLogy();
gStyle->SetOptStat(0);
c1.cd(3);
sig_min->SetMaximum(1);
sig_min->Draw("hist");
bkg_min->SetLineColor(kRed);
bkg_min->Draw("hist same");
gPad->SetLogy();
gStyle->SetOptStat(0);
c1.SaveAs("vrtxquality_sigbkg.pdf");


TFile s_a("vrtx_chi_sig02gev_norm_after02_ovrf.root");
TFile s2_a("vrtx_chi_sig2gev_norm_after02_ovrf.root");
TFile p_a("vrtx_chi_PP_norm_after02_ovrflw.root");
TFile sm_a("vrtx_chi_sigminb_norm_after02_ovrflw.root");

TH1F *sig_mes_02gev_a=(TH1F*)s_a.Get("h2s");
TH1F *sig_mes_2gev_a=(TH1F*)s2_a.Get("h2s");

sig_mes_02gev_a->SetTitle("Vertex quality signal MESMER  Ee>0.2GeV (blue) and bkg minbias (red) with acceptance cut th_mu<35mrad and th_mu>0.2mrad");
sig_mes_2gev_a->SetTitle("Vertex quality signal MESMER  Ee>2GeV (blue) and bkg minbias (red) with acceptance cut th_mu<35mrad and th_mu>0.2mrad");

TH1F *bkg_min_a=(TH1F*)p_a.Get("h2");
TH1F *sig_min_a=(TH1F*)sm_a.Get("htmp_sm2_o");
sig_min_a->SetTitle("Vertex quality signal minbias (blue) and bkg minbias (red) with acceptance cut th_mu<35mrad and th_mu>0.2mrad");

TCanvas c2("c2","c2",700,700);
c2.Divide(1,3);
c2.cd(1);
sig_mes_02gev_a->SetMaximum(1);
sig_mes_02gev_a->Draw("hist");
bkg_min_a->SetLineColor(kRed);
bkg_min_a->Draw("hist same");
gPad->SetLogy();
gStyle->SetOptStat(0);
c2.cd(2);
sig_mes_2gev_a->SetMaximum(1);
sig_mes_2gev_a->Draw("hist");
bkg_min_a->SetLineColor(kRed);
bkg_min_a->Draw("hist same");
gPad->SetLogy();
gStyle->SetOptStat(0);
c2.cd(3);
sig_min_a->SetLineColor(kBlue);
sig_min_a->SetMaximum(1);
sig_min_a->Draw("hist");
bkg_min_a->SetLineColor(kRed);
bkg_min_a->Draw("hist same");
gPad->SetLogy();
gStyle->SetOptStat(0);
c2.SaveAs("vrtxquality_sigbkg_after02.pdf");


TFile ac_s("h_aco_sig_02gev.root");
TFile ac_s2("h_aco_sig_2gev.root");
TFile ac_p("h_aco_PP.root");
TFile ac_sm("h_aco_sig_min.root");

TH1F *h_ac_s=(TH1F*)ac_s.Get("h_acoN");
TH1F *h_ac_s2=(TH1F*)ac_s2.Get("h_acoN");
TH1F *h_ac_p=(TH1F*)ac_p.Get("h_acoN");
TH1F *h_ac_sm=(TH1F*)ac_sm.Get("h_acoN");

h_ac_s->Rebin(4);
h_ac_s2->Rebin(4);
h_ac_p->Rebin(4);
h_ac_sm->Rebin(4);

h_ac_s->SetTitle("Acoplanarity sig MESEMER Ee>0.2GeV (blu) and bkg minbias (red) with acceptance cut th_mu<35mrad");
h_ac_s2->SetTitle("Acoplanarity sig MESMER Ee>2GeV (blu) and bkg minbias (red) with acceptance cut th_mu<35mrad");
h_ac_sm->SetTitle("Acoplanarity sig minbias (green) and bkg minbias (red) with acceptance cut th_mu<35mrad");

TCanvas a1("a1","a1",700,700);
a1.Divide(2,4);
a1.cd(1);
//h_ac_s->SetMinimum(0);
h_ac_s->Draw("hist");
h_ac_p->SetLineColor(kRed);
h_ac_p->Draw("hist same");
a1.cd(2);
//h_ac_s->SetMinimum(0);
h_ac_s->Draw("hist");
h_ac_p->SetLineColor(kRed);
h_ac_p->Draw("hist same");
gPad->SetLogy();
a1.cd(3);
//h_ac_s2->SetMinimum(0);
h_ac_s2->Draw("hist");
h_ac_p->SetLineColor(kRed);
h_ac_p->Draw("hist same");
a1.cd(4);
//h_ac_s2->SetMinimum(0);
h_ac_s2->Draw("hist");
h_ac_p->SetLineColor(kRed);
h_ac_p->Draw("hist same");
gPad->SetLogy();
a1.cd(5);
h_ac_sm->Draw("hist");
h_ac_p->SetLineColor(kRed);
h_ac_p->Draw("hist same");
a1.cd(6);
h_ac_sm->Draw("hist");
h_ac_p->SetLineColor(kRed);
h_ac_p->Draw("hist same");
gPad->SetLogy();
a1.cd(7);
//h_ac_s->SetMinimum(0);
h_ac_sm->SetLineColor(kRed);
h_ac_sm->Draw("hist");
h_ac_s->Draw("hist same");
a1.cd(8);
//h_ac_s->SetMinimum(0);
h_ac_sm->SetTitle("Acoplanarity sig minbias (green) and sig MESMER (blu) with acceptance cut th_mu<35mrad");
h_ac_sm->SetLineColor(8);
h_ac_sm->Draw("hist");
h_ac_s->Draw("hist same");
gPad->SetLogy();

a1.SaveAs("h_acoN.pdf");



TFile ac_sa("h_aco_sig_02gev_after02.root");
TFile ac_s2a("h_aco_sig_2gev_after02.root");
TFile ac_pa("h_aco_PP_after02.root");
TFile ac_sma("h_aco_sig_min_after02.root");

TH1F *h_ac_sa=(TH1F*)ac_sa.Get("h_acoN_02");
TH1F *h_ac_s2a=(TH1F*)ac_s2a.Get("h_acoN_02");
TH1F *h_ac_pa=(TH1F*)ac_pa.Get("h_acoN_02");
TH1F *h_ac_sma=(TH1F*)ac_sma.Get("h_acoN_02");

h_ac_sa->Rebin(4);
h_ac_s2a->Rebin(4);
h_ac_pa->Rebin(4);
h_ac_sma->Rebin(4);

h_ac_sa->SetTitle("Acoplanarity sig MESEMER Ee>0.2GeV (blu) and bkg minbias (red) with acceptance cut th_mu<35mrad and th_mu>0.2mrad");
h_ac_s2a->SetTitle("Acoplanarity sig MESMER Ee>2GeV (blu) and bkg minbias (red) with acceptance cut th_mu<35mrad and th_mu>0.2mrad");
h_ac_sma->SetTitle("Acoplanarity sig minbias (green) and bkg minbias (red) with acceptance cut th_mu<35mrad and th_mu>0.2mrad");


TCanvas a2("a2","a2",700,700);
a2.Divide(2,4);
a2.cd(1);
//h_ac_s->SetMinimum(0);
h_ac_sa->Draw("hist");
h_ac_pa->SetLineColor(kRed);
h_ac_pa->Draw("hist same");
a2.cd(2);
//h_ac_s->SetMinimum(0);
h_ac_sa->Draw("hist");
h_ac_pa->SetLineColor(kRed);
h_ac_pa->Draw("hist same");
gPad->SetLogy();
a2.cd(3);
//h_ac_s2->SetMinimum(0);
h_ac_s2a->Draw("hist");
h_ac_pa->SetLineColor(kRed);
h_ac_pa->Draw("hist same");
a2.cd(4);
//h_ac_s2->SetMinimum(0);
h_ac_s2a->Draw("hist");
h_ac_pa->SetLineColor(kRed);
h_ac_pa->Draw("hist same");
gPad->SetLogy();
a2.cd(5);
h_ac_sma->Draw("hist");
h_ac_pa->SetLineColor(kRed);
h_ac_pa->Draw("hist same");
a2.cd(6);
h_ac_sma->Draw("hist");
h_ac_pa->SetLineColor(kRed);
h_ac_pa->Draw("hist same");
gPad->SetLogy();
a2.cd(7);
//h_ac_s->SetMinimum(0);
h_ac_sma->SetLineColor(kRed);
h_ac_sma->Draw("hist");
h_ac_sa->Draw("hist same");
a2.cd(8);
//h_ac_s->SetMinimum(0);
h_ac_sma->SetTitle("Acoplanarity sig minbias (green) and sig MESMER (blu) with acceptance cut th_mu<35mrad");
h_ac_sma->SetLineColor(8);
h_ac_sma->Draw("hist");
h_ac_sa->Draw("hist same");
gPad->SetLogy();

a2.SaveAs("h_acoN_after02.pdf");



}
