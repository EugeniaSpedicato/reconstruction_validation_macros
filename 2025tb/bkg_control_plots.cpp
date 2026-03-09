#include <iostream>
#include <fstream>
#include <string>
#include <regex>
#include <filesystem>
#include <vector>
#include <algorithm>
#include "TROOT.h"
#include "TLegend.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TStyle.h>

namespace fs = std::filesystem;

//void control_plots(int run,int hits){

int main(int argc, char* argv[]) {

//TFile* f_0 = TFile::Open(Form("txt/run%i/root_single_muon_interaction_0_%s_nhits%i.root",run,file.c_str(),hits));
//TFile* f_1 = TFile::Open(Form("txt/run%i/root_single_muon_interaction_1_%s_nhits%i.root",run,file.c_str(),hits));


    fs::path folder = argv[1];
    std::vector<fs::path> files_0;
    std::vector<fs::path> files_1;

int run = std::stoi(argv[2]);
int hits = std::stoi(argv[3]);

    for (auto& entry : fs::directory_iterator(folder)) {
        if (fs::is_regular_file(entry.path()) && entry.path().extension() == ".root" and entry.path().string().find("possible") != std::string::npos) {
            if(entry.path().string().find("single_muon_interaction_0") != std::string::npos){files_0.push_back(entry.path());}
	    else if(entry.path().string().find("single_muon_interaction_1") != std::string::npos){files_1.push_back(entry.path());}
        }
    }


std::vector<TH1D*> h_nhits_pre_0(files_0.size(), nullptr);
std::vector<TH1D*> h_nhits_post_0(files_0.size(), nullptr);
std::vector<TH1D*> h_ntracks_post_0(files_0.size(), nullptr);
std::vector<TH1D*> h_chi2_vrtx_0(files_0.size(), nullptr);
std::vector<TH1D*> h_Dz_vrtx_0(files_0.size(), nullptr);
std::vector<TH1D*> h_aco_0(files_0.size(), nullptr);
std::vector<TH1D*> h_th_min_0(files_0.size(), nullptr);
std::vector<TH1D*> h_th_max_0(files_0.size(), nullptr);
std::vector<TH1D*> h_IP_0(files_0.size(), nullptr);
std::vector<TH1D*> h_mf_0(files_0.size(), nullptr);
std::vector<TH1D*> h_mf_mod_0(files_0.size(), nullptr);
std::vector<TH1D*> h_nhits_pre_1(files_1.size(), nullptr);
std::vector<TH1D*> h_nhits_post_1(files_1.size(), nullptr);
std::vector<TH1D*> h_ntracks_post_1(files_1.size(), nullptr);
std::vector<TH1D*> h_chi2_vrtx_1(files_1.size(), nullptr);
std::vector<TH1D*> h_Dz_vrtx_1(files_1.size(), nullptr);
std::vector<TH1D*> h_aco_1(files_1.size(), nullptr);
std::vector<TH1D*> h_th_min_1(files_1.size(), nullptr);
std::vector<TH1D*> h_th_max_1(files_1.size(), nullptr);
std::vector<TH1D*> h_IP_1(files_1.size(), nullptr);
std::vector<TH1D*> h_mf_1(files_1.size(), nullptr);
std::vector<TH1D*> h_mf_mod_1(files_1.size(), nullptr);



std::vector<TH1D*> h_nhits_pre_0_sig(files_0.size(), nullptr);
std::vector<TH1D*> h_nhits_post_0_sig(files_0.size(), nullptr);
std::vector<TH1D*> h_ntracks_post_0_sig(files_0.size(), nullptr);
std::vector<TH1D*> h_chi2_vrtx_0_sig(files_0.size(), nullptr);
std::vector<TH1D*> h_Dz_vrtx_0_sig(files_0.size(), nullptr);
std::vector<TH1D*> h_aco_0_sig(files_0.size(), nullptr);
std::vector<TH1D*> h_th_min_0_sig(files_0.size(), nullptr);
std::vector<TH1D*> h_th_max_0_sig(files_0.size(), nullptr);
std::vector<TH1D*> h_IP_0_sig(files_0.size(), nullptr);
std::vector<TH1D*> h_mf_0_sig(files_0.size(), nullptr);
std::vector<TH1D*> h_mf_mod_0_sig(files_0.size(), nullptr);
std::vector<TH1D*> h_nhits_pre_1_sig(files_1.size(), nullptr);
std::vector<TH1D*> h_nhits_post_1_sig(files_1.size(), nullptr);
std::vector<TH1D*> h_ntracks_post_1_sig(files_1.size(), nullptr);
std::vector<TH1D*> h_chi2_vrtx_1_sig(files_1.size(), nullptr);
std::vector<TH1D*> h_Dz_vrtx_1_sig(files_1.size(), nullptr);
std::vector<TH1D*> h_aco_1_sig(files_1.size(), nullptr);
std::vector<TH1D*> h_th_min_1_sig(files_1.size(), nullptr);
std::vector<TH1D*> h_th_max_1_sig(files_1.size(), nullptr);
std::vector<TH1D*> h_IP_1_sig(files_1.size(), nullptr);
std::vector<TH1D*> h_mf_1_sig(files_1.size(), nullptr);
std::vector<TH1D*> h_mf_mod_1_sig(files_1.size(), nullptr);



for(int f=0; f<files_0.size(); f++){


TFile* f_0 = TFile::Open(files_0.at(f).string().c_str(),"READ");

    if (!f_0 ||f_0->IsZombie()) {
        std::cerr << "Errore nell'aprire il file ROOT!\n";
        return 1;
    }

h_nhits_pre_0.at(f)=(TH1D*)f_0->Get("h_nhits_pre_0");
h_nhits_post_0.at(f)=(TH1D*)f_0->Get("h_nhits_post_0");
h_ntracks_post_0.at(f)=(TH1D*)f_0->Get("h_ntracks_post_0");
h_chi2_vrtx_0.at(f)=(TH1D*)f_0->Get("h_chi2_vrtx_0");
h_Dz_vrtx_0.at(f)=(TH1D*)f_0->Get("h_Dz_vrtx_0");
h_aco_0.at(f)=(TH1D*)f_0->Get("h_aco_0");
h_th_min_0.at(f)=(TH1D*)f_0->Get("h_th_min_0");
h_th_max_0.at(f)=(TH1D*)f_0->Get("h_th_max_0");
h_IP_0.at(f)=(TH1D*)f_0->Get("h_IP_0");
h_mf_0.at(f)=(TH1D*)f_0->Get("h_mf_0");
h_mf_mod_0.at(f)=(TH1D*)f_0->Get("h_mf_mod_0");

h_nhits_pre_0_sig.at(f)=(TH1D*)f_0->Get("h_nhits_pre_0_sig");
h_nhits_post_0_sig.at(f)=(TH1D*)f_0->Get("h_nhits_post_0_sig");
h_ntracks_post_0_sig.at(f)=(TH1D*)f_0->Get("h_ntracks_post_0_sig");
h_chi2_vrtx_0_sig.at(f)=(TH1D*)f_0->Get("h_chi2_vrtx_0_sig");
h_Dz_vrtx_0_sig.at(f)=(TH1D*)f_0->Get("h_Dz_vrtx_0_sig");
h_aco_0_sig.at(f)=(TH1D*)f_0->Get("h_aco_0_sig");
h_th_min_0_sig.at(f)=(TH1D*)f_0->Get("h_th_min_0_sig");
h_th_max_0_sig.at(f)=(TH1D*)f_0->Get("h_th_max_0_sig");
h_IP_0_sig.at(f)=(TH1D*)f_0->Get("h_IP_0_sig");
h_mf_0_sig.at(f)=(TH1D*)f_0->Get("h_mf_0_sig");
h_mf_mod_0_sig.at(f)=(TH1D*)f_0->Get("h_mf_mod_0_sig");

        if(f>0){
h_nhits_pre_0.at(0)->Add(h_nhits_pre_0.at(f));
h_nhits_post_0.at(0)->Add(h_nhits_post_0.at(f));
h_ntracks_post_0.at(0)->Add(h_ntracks_post_0.at(f));
h_chi2_vrtx_0.at(0)->Add(h_chi2_vrtx_0.at(f));
h_Dz_vrtx_0.at(0)->Add(h_Dz_vrtx_0.at(f));
h_aco_0.at(0)->Add(h_aco_0.at(f));
h_th_min_0.at(0)->Add(h_th_min_0.at(f));
h_th_max_0.at(0)->Add(h_th_max_0.at(f));
h_IP_0.at(0)->Add(h_IP_0.at(f));
h_mf_0.at(0)->Add(h_mf_0.at(f));
h_mf_mod_0.at(0)->Add(h_mf_mod_0.at(f));


h_nhits_pre_0_sig.at(0)->Add(h_nhits_pre_0_sig.at(f));
h_nhits_post_0_sig.at(0)->Add(h_nhits_post_0_sig.at(f));
h_ntracks_post_0_sig.at(0)->Add(h_ntracks_post_0_sig.at(f));
h_chi2_vrtx_0_sig.at(0)->Add(h_chi2_vrtx_0_sig.at(f));
h_Dz_vrtx_0_sig.at(0)->Add(h_Dz_vrtx_0_sig.at(f));
h_aco_0_sig.at(0)->Add(h_aco_0_sig.at(f));
h_th_min_0_sig.at(0)->Add(h_th_min_0_sig.at(f));
h_th_max_0_sig.at(0)->Add(h_th_max_0_sig.at(f));
//h_IP_0_sig.at(0)->Add(h_IP_0_sig.at(f));
h_mf_0_sig.at(0)->Add(h_mf_0_sig.at(f));
h_mf_mod_0_sig.at(0)->Add(h_mf_mod_0_sig.at(f));
	}

}

for(int f=0; f<files_1.size(); f++){


TFile* f_0 = TFile::Open(files_1.at(f).string().c_str(),"READ");


    if (!f_0 ||f_0->IsZombie()) {
        std::cerr << "Errore nell'aprire il file ROOT!\n";
        return 1;
    }



h_nhits_pre_1.at(f)=(TH1D*)f_0->Get("h_nhits_pre_1");
h_nhits_post_1.at(f)=(TH1D*)f_0->Get("h_nhits_post_1");
h_ntracks_post_1.at(f)=(TH1D*)f_0->Get("h_ntracks_post_1");
h_chi2_vrtx_1.at(f)=(TH1D*)f_0->Get("h_chi2_vrtx_1");
h_Dz_vrtx_1.at(f)=(TH1D*)f_0->Get("h_Dz_vrtx_1");
h_aco_1.at(f)=(TH1D*)f_0->Get("h_aco_1");
h_th_min_1.at(f)=(TH1D*)f_0->Get("h_th_min_1");
h_th_max_1.at(f)=(TH1D*)f_0->Get("h_th_max_1");
h_IP_1.at(f)=(TH1D*)f_0->Get("h_IP_1");
h_mf_1.at(f)=(TH1D*)f_0->Get("h_mf_1");
h_mf_mod_1.at(f)=(TH1D*)f_0->Get("h_mf_mod_1");

h_nhits_pre_1_sig.at(f)=(TH1D*)f_0->Get("h_nhits_pre_1_sig");
h_nhits_post_1_sig.at(f)=(TH1D*)f_0->Get("h_nhits_post_1_sig");
h_ntracks_post_1_sig.at(f)=(TH1D*)f_0->Get("h_ntracks_post_1_sig");
h_chi2_vrtx_1_sig.at(f)=(TH1D*)f_0->Get("h_chi2_vrtx_1_sig");
h_Dz_vrtx_1_sig.at(f)=(TH1D*)f_0->Get("h_Dz_vrtx_1_sig");
h_aco_1_sig.at(f)=(TH1D*)f_0->Get("h_aco_1_sig");
h_th_min_1_sig.at(f)=(TH1D*)f_0->Get("h_th_min_1_sig");
h_th_max_1_sig.at(f)=(TH1D*)f_0->Get("h_th_max_1_sig");
//h_IP_1_sig.at(f)=(TH1D*)f_0->Get("h_IP_1_sig");
h_mf_1_sig.at(f)=(TH1D*)f_0->Get("h_mf_1_sig");
h_mf_mod_1_sig.at(f)=(TH1D*)f_0->Get("h_mf_mod_1_sig");


	if(f>0){
h_nhits_pre_1.at(0)->Add(h_nhits_pre_1.at(f));
h_nhits_post_1.at(0)->Add(h_nhits_post_1.at(f));
h_ntracks_post_1.at(0)->Add(h_ntracks_post_1.at(f));
h_chi2_vrtx_1.at(0)->Add(h_chi2_vrtx_1.at(f));
h_Dz_vrtx_1.at(0)->Add(h_Dz_vrtx_1.at(f));
h_aco_1.at(0)->Add(h_aco_1.at(f));
h_th_min_1.at(0)->Add(h_th_min_1.at(f));
h_th_max_1.at(0)->Add(h_th_max_1.at(f));
h_IP_1.at(0)->Add(h_IP_1.at(f));
h_mf_1.at(0)->Add(h_mf_1.at(f));
h_mf_mod_1.at(0)->Add(h_mf_mod_1.at(f));


h_nhits_pre_1_sig.at(0)->Add(h_nhits_pre_1_sig.at(f));
h_nhits_post_1_sig.at(0)->Add(h_nhits_post_1_sig.at(f));
h_ntracks_post_1_sig.at(0)->Add(h_ntracks_post_1_sig.at(f));
h_chi2_vrtx_1_sig.at(0)->Add(h_chi2_vrtx_1_sig.at(f));
h_Dz_vrtx_1_sig.at(0)->Add(h_Dz_vrtx_1_sig.at(f));
h_aco_1_sig.at(0)->Add(h_aco_1_sig.at(f));
h_th_min_1_sig.at(0)->Add(h_th_min_1_sig.at(f));
h_th_max_1_sig.at(0)->Add(h_th_max_1_sig.at(f));
//h_IP_1_sig.at(0)->Add(h_IP_1_sig.at(f));
h_mf_1_sig.at(0)->Add(h_mf_1_sig.at(f));
h_mf_mod_1_sig.at(0)->Add(h_mf_mod_1_sig.at(f));
		}
}


h_nhits_pre_0.at(0)->Scale(1./h_nhits_pre_0.at(0)->Integral());
h_nhits_post_0.at(0)->Scale(1./h_nhits_post_0.at(0)->Integral());
h_ntracks_post_0.at(0)->Scale(1./h_ntracks_post_0.at(0)->Integral());
h_chi2_vrtx_0.at(0)->Scale(1./h_chi2_vrtx_0.at(0)->Integral());
h_Dz_vrtx_0.at(0)->Scale(1./h_Dz_vrtx_0.at(0)->Integral());
h_aco_0.at(0)->Scale(1./h_aco_0.at(0)->Integral());
h_th_min_0.at(0)->Scale(1./h_th_min_0.at(0)->Integral());
h_th_max_0.at(0)->Scale(1./h_th_max_0.at(0)->Integral());
h_IP_0.at(0)->Scale(1./h_IP_0.at(0)->Integral());
h_mf_0.at(0)->Scale(1./h_mf_0.at(0)->Integral());
h_mf_mod_0.at(0)->Scale(1./h_mf_mod_0.at(0)->Integral());


h_nhits_pre_1.at(0)->Scale(1./h_nhits_pre_1.at(0)->Integral());
h_nhits_post_1.at(0)->Scale(1./h_nhits_post_1.at(0)->Integral());
h_ntracks_post_1.at(0)->Scale(1./h_ntracks_post_1.at(0)->Integral());
h_chi2_vrtx_1.at(0)->Scale(1./h_chi2_vrtx_1.at(0)->Integral());
h_Dz_vrtx_1.at(0)->Scale(1./h_Dz_vrtx_1.at(0)->Integral());
h_aco_1.at(0)->Scale(1./h_aco_1.at(0)->Integral());
h_th_min_1.at(0)->Scale(1./h_th_min_1.at(0)->Integral());
h_th_max_1.at(0)->Scale(1./h_th_max_1.at(0)->Integral());
h_IP_1.at(0)->Scale(1./h_IP_1.at(0)->Integral());
h_mf_1.at(0)->Scale(1./h_mf_1.at(0)->Integral());
h_mf_mod_1.at(0)->Scale(1./h_mf_mod_1.at(0)->Integral());



h_nhits_pre_0_sig.at(0)->Scale(1./h_nhits_pre_0_sig.at(0)->Integral());
h_nhits_post_0_sig.at(0)->Scale(1./h_nhits_post_0_sig.at(0)->Integral());
h_ntracks_post_0_sig.at(0)->Scale(1./h_ntracks_post_0_sig.at(0)->Integral());
h_chi2_vrtx_0_sig.at(0)->Scale(1./h_chi2_vrtx_0_sig.at(0)->Integral());
h_Dz_vrtx_0_sig.at(0)->Scale(1./h_Dz_vrtx_0_sig.at(0)->Integral());
h_aco_0_sig.at(0)->Scale(1./h_aco_0_sig.at(0)->Integral());
h_th_min_0_sig.at(0)->Scale(1./h_th_min_0_sig.at(0)->Integral());
h_th_max_0_sig.at(0)->Scale(1./h_th_max_0_sig.at(0)->Integral());
//h_IP_0_sig.at(0)->Scale(1./h_IP_0_sig.at(0)->Integral());
h_mf_0_sig.at(0)->Scale(1./h_mf_0_sig.at(0)->Integral());
h_mf_mod_0_sig.at(0)->Scale(1./h_mf_mod_0_sig.at(0)->Integral());


h_nhits_pre_1_sig.at(0)->Scale(1./h_nhits_pre_1_sig.at(0)->Integral());
h_nhits_post_1_sig.at(0)->Scale(1./h_nhits_post_1_sig.at(0)->Integral());
h_ntracks_post_1_sig.at(0)->Scale(1./h_ntracks_post_1_sig.at(0)->Integral());
h_chi2_vrtx_1_sig.at(0)->Scale(1./h_chi2_vrtx_1_sig.at(0)->Integral());
h_Dz_vrtx_1_sig.at(0)->Scale(1./h_Dz_vrtx_1_sig.at(0)->Integral());
h_aco_1_sig.at(0)->Scale(1./h_aco_1_sig.at(0)->Integral());
h_th_min_1_sig.at(0)->Scale(1./h_th_min_1_sig.at(0)->Integral());
h_th_max_1_sig.at(0)->Scale(1./h_th_max_1_sig.at(0)->Integral());
//h_IP_1_sig.at(0)->Scale(1./h_IP_1_sig.at(0)->Integral());
h_mf_1_sig.at(0)->Scale(1./h_mf_1_sig.at(0)->Integral());
h_mf_mod_1_sig.at(0)->Scale(1./h_mf_mod_1_sig.at(0)->Integral());

auto leg1=new TLegend(0.65, 0.30, 0.90, 0.45);
leg1->AddEntry(h_nhits_pre_0.at(0),"Target 0");
leg1->AddEntry(h_nhits_pre_1.at(0),"Target 1");

TCanvas gen("gen","gen",1500,1500);
gen.Divide(3,3);
gen.cd(1);
h_nhits_pre_0.at(0)->Draw("hist");
h_nhits_pre_1.at(0)->SetLineColor(kRed);
h_nhits_pre_1.at(0)->Draw("hist same");
gPad->SetLogy();
leg1->Draw();
gen.cd(2);
h_nhits_post_0.at(0)->Draw("hist");
h_nhits_post_1.at(0)->SetLineColor(kRed);
h_nhits_post_1.at(0)->Draw("hist same");
gPad->SetLogy();
gen.cd(3);
h_ntracks_post_0.at(0)->Draw("hist");
h_ntracks_post_1.at(0)->SetLineColor(kRed);
h_ntracks_post_1.at(0)->Draw("hist same");
gen.cd(4);
h_chi2_vrtx_0.at(0)->Draw("hist");
h_chi2_vrtx_1.at(0)->SetLineColor(kRed);
h_chi2_vrtx_1.at(0)->Draw("hist same");
gen.cd(5);
h_Dz_vrtx_0.at(0)->Draw("hist");
h_Dz_vrtx_1.at(0)->SetLineColor(kRed);
h_Dz_vrtx_1.at(0)->Draw("hist same");
gen.cd(6);
h_aco_0.at(0)->Rebin(2);
h_aco_1.at(0)->Rebin(2);
h_aco_0.at(0)->Draw("hist");
h_aco_1.at(0)->SetLineColor(kRed);
h_aco_1.at(0)->Draw("hist same");
gen.cd(7);
h_th_min_0.at(0)->Draw("hist");
h_th_min_1.at(0)->SetLineColor(kRed);
h_th_min_1.at(0)->Draw("hist same");
gen.cd(8);
h_th_max_0.at(0)->Draw("hist");
h_th_max_1.at(0)->SetLineColor(kRed);
h_th_max_1.at(0)->Draw("hist same");
gen.cd(9);
h_IP_0.at(0)->Draw("hist");
h_IP_1.at(0)->SetLineColor(kRed);
h_IP_1.at(0)->Draw("hist same");
gen.SaveAs(Form("plots/run%i/bkg_postFid_nhits%i.pdf",run,hits));

TCanvas pe("pe","pe",1500,1000);
pe.Divide(2,1);
pe.cd(1);
h_mf_0.at(0)->Draw("hist");
h_mf_1.at(0)->SetLineColor(kRed);
h_mf_1.at(0)->Draw("hist same");
leg1->Draw();
gPad->SetLogy();
pe.cd(2);
h_mf_mod_0.at(0)->Draw("hist");
h_mf_mod_1.at(0)->SetLineColor(kRed);
h_mf_mod_1.at(0)->Draw("hist same");
pe.SaveAs(Form("plots/run%i/bkg_preEl_nhits%i.pdf",run,hits));






TCanvas gen_sig("gen_sig","gen_sig",1500,1500);
gen_sig.Divide(3,3);
gen_sig.cd(1);
h_nhits_pre_0_sig.at(0)->Draw("hist");
h_nhits_pre_1_sig.at(0)->SetLineColor(kRed);
h_nhits_pre_1_sig.at(0)->Draw("hist same");
leg1->Draw();
gPad->SetLogy();
gen_sig.cd(2);
h_nhits_post_0_sig.at(0)->Draw("hist");
h_nhits_post_1_sig.at(0)->SetLineColor(kRed);
h_nhits_post_1_sig.at(0)->Draw("hist same");
gPad->SetLogy();
gen_sig.cd(3);
h_ntracks_post_0_sig.at(0)->Draw("hist");
h_ntracks_post_1_sig.at(0)->SetLineColor(kRed);
h_ntracks_post_1_sig.at(0)->Draw("hist same");
gen_sig.cd(4);
h_chi2_vrtx_0_sig.at(0)->Draw("hist");
h_chi2_vrtx_1_sig.at(0)->SetLineColor(kRed);
h_chi2_vrtx_1_sig.at(0)->Draw("hist same");
gen_sig.cd(5);
h_Dz_vrtx_0_sig.at(0)->Draw("hist");
h_Dz_vrtx_1_sig.at(0)->SetLineColor(kRed);
h_Dz_vrtx_1_sig.at(0)->Draw("hist same");
gen_sig.cd(6);
h_aco_0_sig.at(0)->Draw("hist");
h_aco_1_sig.at(0)->SetLineColor(kRed);
h_aco_1_sig.at(0)->Draw("hist same");
gen_sig.cd(7);
h_th_min_0_sig.at(0)->Draw("hist");
h_th_min_1_sig.at(0)->SetLineColor(kRed);
h_th_min_1_sig.at(0)->Draw("hist same");
gen_sig.cd(8);
h_th_max_0_sig.at(0)->Draw("hist");
h_th_max_1_sig.at(0)->SetLineColor(kRed);
h_th_max_1_sig.at(0)->Draw("hist same");
/*gen_sig.cd(9);
h_IP_0_sig.at(0)->Draw("hist");
h_IP_1_sig.at(0)->SetLineColor(kRed);
h_IP_1_sig.at(0)->Draw("hist same");*/
gen_sig.SaveAs(Form("plots/run%i/sig_postFid_nhits%i.pdf",run,hits));

TCanvas pe_sig("pe_sig","pe_sig",1500,1000);
pe_sig.Divide(2,1);
pe_sig.cd(1);
h_mf_0_sig.at(0)->Draw("hist");
h_mf_1_sig.at(0)->SetLineColor(kRed);
h_mf_1_sig.at(0)->Draw("hist same");
gPad->SetLogy();
leg1->Draw();
pe_sig.cd(2);
h_mf_mod_0_sig.at(0)->Draw("hist");
h_mf_mod_1_sig.at(0)->SetLineColor(kRed);
h_mf_mod_1_sig.at(0)->Draw("hist same");
pe_sig.SaveAs(Form("plots/run%i/sig_preEl_nhits%i.pdf",run,hits));


auto leg2=new TLegend(0.65, 0.30, 0.90, 0.45);
leg2->AddEntry(h_nhits_pre_0_sig.at(0),"Signal candidate");
leg2->AddEntry(h_nhits_pre_0.at(0),"Bkg candidate");

TCanvas gen_sig0("gen_sig0","gen_sig0",1500,1500);
gen_sig0.Divide(3,3);
gen_sig0.cd(1);
h_nhits_pre_0_sig.at(0)->SetLineColor(kRed);
h_nhits_pre_0_sig.at(0)->Draw("hist");
h_nhits_pre_0.at(0)->Draw("hist same");
leg2->Draw();
gPad->SetLogy();
gen_sig0.cd(2);
h_nhits_post_0_sig.at(0)->SetLineColor(kRed);
h_nhits_post_0_sig.at(0)->Draw("hist");
h_nhits_post_0.at(0)->Draw("hist same");
gPad->SetLogy();
gen_sig0.cd(3);
h_ntracks_post_0_sig.at(0)->SetLineColor(kRed);
h_ntracks_post_0_sig.at(0)->Draw("hist");
h_ntracks_post_0.at(0)->Draw("hist same");
gen_sig0.cd(4);
h_chi2_vrtx_0_sig.at(0)->SetLineColor(kRed);
h_chi2_vrtx_0_sig.at(0)->Draw("hist");
h_chi2_vrtx_0.at(0)->Draw("hist same");
gen_sig0.cd(5);
h_Dz_vrtx_0_sig.at(0)->SetLineColor(kRed);
h_Dz_vrtx_0_sig.at(0)->Draw("hist");
h_Dz_vrtx_0.at(0)->Draw("hist same");
gen_sig0.cd(6);
h_aco_0_sig.at(0)->SetLineColor(kRed);
h_aco_0_sig.at(0)->Draw("hist");
h_aco_0.at(0)->Draw("hist same");
gen_sig0.cd(7);
h_th_min_0.at(0)->GetXaxis()->SetRangeUser(0.,0.001);
h_th_min_0.at(0)->Draw("hist");
h_th_min_0_sig.at(0)->SetLineColor(kRed);
h_th_min_0_sig.at(0)->Draw("hist same");
gen_sig0.cd(8);
h_th_max_0.at(0)->Draw("hist");
h_th_max_0_sig.at(0)->SetLineColor(kRed);
h_th_max_0_sig.at(0)->Draw("hist same");
/*gen_sig0.cd(9);
h_IP_0.at(0)->Draw("hist");
h_IP_0_sig.at(0)->SetLineColor(kRed);
h_IP_0_sig.at(0)->Draw("hist same");*/
gen_sig0.SaveAs(Form("plots/run%i/tar0_postFid_nhits%i.pdf",run,hits));

TCanvas pe_sig0("pe_sig0","pe_sig0",1500,1000);
pe_sig0.Divide(2,1);
pe_sig0.cd(1);
h_mf_0_sig.at(0)->SetLineColor(kRed);
h_mf_0_sig.at(0)->Draw("hist");
h_mf_0.at(0)->Draw("hist same");
leg2->Draw();
gPad->SetLogy();
pe_sig0.cd(2);
h_mf_mod_0.at(0)->Draw("hist");
h_mf_mod_0_sig.at(0)->SetLineColor(kRed);
h_mf_mod_0_sig.at(0)->Draw("hist same");
pe_sig0.SaveAs(Form("plots/run%i/tar0_preEl_nhits%i.pdf",run,hits));



TCanvas gen_sig1("gen_sig1","gen_sig1",1500,1500);
gen_sig1.Divide(3,3);
gen_sig1.cd(1);
h_nhits_pre_1.at(0)->SetLineColor(kBlue);
h_nhits_pre_1_sig.at(0)->Draw("hist");
h_nhits_pre_1.at(0)->Draw("hist same");
leg2->Draw();
gPad->SetLogy();
gen_sig1.cd(2);
h_nhits_post_1.at(0)->Draw("hist");
h_nhits_post_1.at(0)->SetLineColor(kBlue);
h_nhits_post_1_sig.at(0)->Draw("hist same");
gPad->SetLogy();
gen_sig1.cd(3);
h_ntracks_post_1.at(0)->Draw("hist");
h_ntracks_post_1.at(0)->SetLineColor(kBlue);
h_ntracks_post_1_sig.at(0)->Draw("hist same");
gen_sig1.cd(4);
h_chi2_vrtx_1.at(0)->SetLineColor(kBlue);
h_chi2_vrtx_1_sig.at(0)->Draw("hist");
h_chi2_vrtx_1.at(0)->Draw("hist same");
gen_sig1.cd(5);
h_Dz_vrtx_1.at(0)->Draw("hist");
h_Dz_vrtx_1_sig.at(0)->SetLineColor(kBlue);
h_Dz_vrtx_1_sig.at(0)->Draw("hist same");
gen_sig1.cd(6);
h_aco_1_sig.at(0)->SetLineColor(kBlue);
h_aco_1_sig.at(0)->Draw("hist");
h_aco_1.at(0)->Draw("hist same");
gen_sig1.cd(7);
h_th_min_1.at(0)->GetXaxis()->SetRangeUser(0.,0.001);
h_th_min_1.at(0)->Draw("hist");
h_th_min_1.at(0)->SetLineColor(kBlue);
h_th_min_1_sig.at(0)->Draw("hist same");
gen_sig1.cd(8);
h_th_max_1.at(0)->Draw("hist");
h_th_max_1.at(0)->SetLineColor(kBlue);
h_th_max_1_sig.at(0)->Draw("hist same");
/*gen_sig1.cd(9);
h_IP_1.at(0)->Draw("hist");
h_IP_1_sig.at(0)->SetLineColor(kBlue);
h_IP_1_sig.at(0)->Draw("hist same");*/
gen_sig1.SaveAs(Form("plots/run%i/tar1_postFid_nhits%i.pdf",run,hits));

TCanvas pe_sig1("pe_sig1","pe_sig1",1500,1000);
pe_sig1.Divide(2,1);
pe_sig1.cd(1);
h_mf_1.at(0)->SetLineColor(kBlue);
h_mf_1_sig.at(0)->Draw("hist");
h_mf_1.at(0)->Draw("hist same");
gPad->SetLogy();
leg2->Draw();
pe_sig1.cd(2);
h_mf_mod_1.at(0)->Draw("hist");
h_mf_mod_1.at(0)->SetLineColor(kBlue);
h_mf_mod_1_sig.at(0)->Draw("hist same");
pe_sig1.SaveAs(Form("plots/run%i/tar1_preEl_nhits%i.pdf",run,hits));


return 0;

}
