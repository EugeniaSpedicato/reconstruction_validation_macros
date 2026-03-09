#include <iostream>
#include <fstream>
#include <string>
#include <regex>
#include <filesystem>
#include <vector>
#include <algorithm>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TStyle.h>

namespace fs = std::filesystem;

//void control_plots(int run,int hits){

int main(int argc, char* argv[]) {

//TFile* f_0 = TFile::Open(Form("txt/run%s/root_single_int_0_%s.root",run,file.c_str(),hits));
//TFile* f_1 = TFile::Open(Form("txt/run%s/root_single_int_1_%s.root",run,file.c_str(),hits));


    fs::path folder = argv[1];
    std::vector<fs::path> files_0;
    std::vector<fs::path> files_1;

char* run = argv[2];

    for (auto& entry : fs::directory_iterator(folder)) {
        if (fs::is_regular_file(entry.path()) && entry.path().extension() == ".root" and entry.path().string().find("zin3") != std::string::npos and entry.path().string().find("possible") == std::string::npos and entry.path().string().find("preVrtx") == std::string::npos and entry.path().string().find("smearing") == std::string::npos) {
            if(entry.path().string().find("single_int_0") != std::string::npos){files_0.push_back(entry.path());}
	    else if(entry.path().string().find("single_int_1") != std::string::npos){files_1.push_back(entry.path());}
        }
    }


std::vector<TH1D*> h_x_0(files_0.size(), nullptr);
std::vector<TH1D*> h_y_0(files_0.size(), nullptr);
std::vector<TH1D*> h_thx_0(files_0.size(), nullptr);
std::vector<TH1D*> h_thy_0(files_0.size(), nullptr);
std::vector<TH1D*> h_nhits_pre_0(files_0.size(), nullptr);
std::vector<TH1D*> h_nhits_post_0(files_0.size(), nullptr);
std::vector<TH1D*> h_ntracks_post_0(files_0.size(), nullptr);
std::vector<TH1D*> h_chi2_vrtx_0(files_0.size(), nullptr);
std::vector<TH1D*> h_Dz_vrtx_0(files_0.size(), nullptr);
std::vector<TH1D*> h_aco_0(files_0.size(), nullptr);
std::vector<TH1D*> h_th_min_0(files_0.size(), nullptr);
std::vector<TH1D*> h_th_max_0(files_0.size(), nullptr);
std::vector<TH1D*> h_preE_th_min_0(files_0.size(), nullptr);
std::vector<TH1D*> h_preE_th_max_0(files_0.size(), nullptr);
std::vector<TH1D*> h_elastic_0(files_0.size(), nullptr);
std::vector<TH1D*> h_el_1(files_0.size(), nullptr);
std::vector<TH1D*> h_E_th_min_0(files_0.size(), nullptr);
std::vector<TH1D*> h_E_th_max_0(files_0.size(), nullptr);
std::vector<TH1D*> h_final_nhits_pre_0(files_0.size(), nullptr);
std::vector<TH1D*> h_final_nhits_post_0(files_0.size(), nullptr);
std::vector<TH1D*> h_final_ntracks_post_0(files_0.size(), nullptr);
std::vector<TH1D*> h_final_chi2_vrtx_0(files_0.size(), nullptr);
std::vector<TH1D*> h_final_Dz_vrtx_0(files_0.size(), nullptr);
std::vector<TH1D*> h_final_aco_0(files_0.size(), nullptr);
std::vector<TH1D*> h_final_chi2_tr0_0(files_0.size(), nullptr);
std::vector<TH1D*> h_final_chi2_trMax_0(files_0.size(), nullptr);
std::vector<TH1D*> h_final_chi2_trMin_0(files_0.size(), nullptr);

std::vector<TH1D*> h_x_1(files_1.size(), nullptr);
std::vector<TH1D*> h_y_1(files_1.size(), nullptr);
std::vector<TH1D*> h_thx_1(files_1.size(), nullptr);
std::vector<TH1D*> h_thy_1(files_1.size(), nullptr);
std::vector<TH1D*> h_nhits_pre_1(files_1.size(), nullptr);
std::vector<TH1D*> h_nhits_post_1(files_1.size(), nullptr);
std::vector<TH1D*> h_ntracks_post_1(files_1.size(), nullptr);
std::vector<TH1D*> h_chi2_vrtx_1(files_1.size(), nullptr);
std::vector<TH1D*> h_Dz_vrtx_1(files_1.size(), nullptr);
std::vector<TH1D*> h_aco_1(files_1.size(), nullptr);
std::vector<TH1D*> h_th_min_1(files_1.size(), nullptr);
std::vector<TH1D*> h_th_max_1(files_1.size(), nullptr);
std::vector<TH1D*> h_preE_th_min_1(files_1.size(), nullptr);
std::vector<TH1D*> h_preE_th_max_1(files_1.size(), nullptr);
std::vector<TH1D*> h_elastic_1(files_1.size(), nullptr);
std::vector<TH1D*> h_el_2(files_1.size(), nullptr);
std::vector<TH1D*> h_E_th_min_1(files_1.size(), nullptr);
std::vector<TH1D*> h_E_th_max_1(files_1.size(), nullptr);
std::vector<TH1D*> h_final_nhits_pre_1(files_1.size(), nullptr);
std::vector<TH1D*> h_final_nhits_post_1(files_1.size(), nullptr);
std::vector<TH1D*> h_final_ntracks_post_1(files_1.size(), nullptr);
std::vector<TH1D*> h_final_chi2_vrtx_1(files_1.size(), nullptr);
std::vector<TH1D*> h_final_Dz_vrtx_1(files_1.size(), nullptr);
std::vector<TH1D*> h_final_aco_1(files_1.size(), nullptr);
std::vector<TH1D*> h_final_chi2_tr0_1(files_1.size(), nullptr);
std::vector<TH1D*> h_final_chi2_trMax_1(files_1.size(), nullptr);
std::vector<TH1D*> h_final_chi2_trMin_1(files_1.size(), nullptr);


for(int f=0; f<files_0.size(); f++){


TFile* f_0 = TFile::Open(files_0.at(f).string().c_str(),"READ");

    if (!f_0 ||f_0->IsZombie()) {
        std::cerr << "Errore nell'aprire il file ROOT!\n";
        return 1;
    }

h_x_0.at(f)=(TH1D*)f_0->Get("h_x_0");
h_y_0.at(f)=(TH1D*)f_0->Get("h_y_0");
h_thx_0.at(f)=(TH1D*)f_0->Get("h_thx_0");
h_thy_0.at(f)=(TH1D*)f_0->Get("h_thy_0");
h_nhits_pre_0.at(f)=(TH1D*)f_0->Get("h_nhits_pre_0");
h_nhits_post_0.at(f)=(TH1D*)f_0->Get("h_nhits_post_0");
h_ntracks_post_0.at(f)=(TH1D*)f_0->Get("h_ntracks_post_0");
h_chi2_vrtx_0.at(f)=(TH1D*)f_0->Get("h_chi2_vrtx_0");
h_Dz_vrtx_0.at(f)=(TH1D*)f_0->Get("h_Dz_vrtx_0");
h_aco_0.at(f)=(TH1D*)f_0->Get("h_aco_0");
h_th_min_0.at(f)=(TH1D*)f_0->Get("h_th_min_0");
h_th_max_0.at(f)=(TH1D*)f_0->Get("h_th_max_0");
h_preE_th_min_0.at(f)=(TH1D*)f_0->Get("h_preE_th_min_0");
h_preE_th_max_0.at(f)=(TH1D*)f_0->Get("h_preE_th_max_0");
h_elastic_0.at(f)=(TH1D*)f_0->Get("h_elastic_0");
h_el_1.at(f)=(TH1D*)f_0->Get("h_el_1");
h_E_th_min_0.at(f)=(TH1D*)f_0->Get("h_E_th_min_0");
h_E_th_max_0.at(f)=(TH1D*)f_0->Get("h_E_th_max_0");
h_final_nhits_pre_0.at(f)=(TH1D*)f_0->Get("h_final_nhits_pre_0");
h_final_nhits_post_0.at(f)=(TH1D*)f_0->Get("h_final_nhits_post_0");
h_final_ntracks_post_0.at(f)=(TH1D*)f_0->Get("h_final_ntracks_post_0");
h_final_chi2_vrtx_0.at(f)=(TH1D*)f_0->Get("h_final_chi2_vrtx_0");
h_final_Dz_vrtx_0.at(f)=(TH1D*)f_0->Get("h_final_Dz_vrtx_0");
h_final_aco_0.at(f)=(TH1D*)f_0->Get("h_final_aco_0");
h_final_chi2_tr0_0.at(f)=(TH1D*)f_0->Get("h_final_chi2_tr0_0");
h_final_chi2_trMax_0.at(f)=(TH1D*)f_0->Get("h_final_chi2_trMax_0");
h_final_chi2_trMin_0.at(f)=(TH1D*)f_0->Get("h_final_chi2_trMin_0");

        if(f>0){
h_x_0.at(0)->Add(h_x_0.at(f));
h_y_0.at(0)->Add(h_y_0.at(f));
h_thx_0.at(0)->Add(h_thx_0.at(f));
h_thy_0.at(0)->Add(h_thy_0.at(f));
h_nhits_pre_0.at(0)->Add(h_nhits_pre_0.at(f));
h_nhits_post_0.at(0)->Add(h_nhits_post_0.at(f));
h_ntracks_post_0.at(0)->Add(h_ntracks_post_0.at(f));
h_chi2_vrtx_0.at(0)->Add(h_chi2_vrtx_0.at(f));
h_Dz_vrtx_0.at(0)->Add(h_Dz_vrtx_0.at(f));
h_aco_0.at(0)->Add(h_aco_0.at(f));
h_th_min_0.at(0)->Add(h_th_min_0.at(f));
h_th_max_0.at(0)->Add(h_th_max_0.at(f));
h_preE_th_min_0.at(0)->Add(h_preE_th_min_0.at(f));
h_preE_th_max_0.at(0)->Add(h_preE_th_max_0.at(f));
h_elastic_0.at(0)->Add(h_elastic_0.at(f));
h_el_1.at(0)->Add(h_el_1.at(f));
h_E_th_min_0.at(0)->Add(h_E_th_min_0.at(f));
h_E_th_max_0.at(0)->Add(h_E_th_max_0.at(f));
h_final_nhits_pre_0.at(0)->Add(h_final_nhits_pre_0.at(f));
h_final_nhits_post_0.at(0)->Add(h_final_nhits_post_0.at(f));
h_final_ntracks_post_0.at(0)->Add(h_final_ntracks_post_0.at(f));
h_final_chi2_vrtx_0.at(0)->Add(h_final_chi2_vrtx_0.at(f));
h_final_Dz_vrtx_0.at(0)->Add(h_final_Dz_vrtx_0.at(f));
h_final_aco_0.at(0)->Add(h_final_aco_0.at(f));
h_final_chi2_tr0_0.at(0)->Add(h_final_chi2_tr0_0.at(f));
h_final_chi2_trMax_0.at(0)->Add(h_final_chi2_trMax_0.at(f));
h_final_chi2_trMin_0.at(0)->Add(h_final_chi2_trMin_0.at(f));
	}

}

for(int f=0; f<files_1.size(); f++){


TFile* f_0 = TFile::Open(files_1.at(f).string().c_str(),"READ");


    if (!f_0 ||f_0->IsZombie()) {
        std::cerr << "Errore nell'aprire il file ROOT!\n";
        return 1;
    }


h_x_1.at(f)=(TH1D*)f_0->Get("h_x_1");
h_y_1.at(f)=(TH1D*)f_0->Get("h_y_1");
h_thx_1.at(f)=(TH1D*)f_0->Get("h_thx_1");
h_thy_1.at(f)=(TH1D*)f_0->Get("h_thy_1");
h_nhits_pre_1.at(f)=(TH1D*)f_0->Get("h_nhits_pre_1");
h_nhits_post_1.at(f)=(TH1D*)f_0->Get("h_nhits_post_1");
h_ntracks_post_1.at(f)=(TH1D*)f_0->Get("h_ntracks_post_1");
h_chi2_vrtx_1.at(f)=(TH1D*)f_0->Get("h_chi2_vrtx_1");
h_Dz_vrtx_1.at(f)=(TH1D*)f_0->Get("h_Dz_vrtx_1");
h_aco_1.at(f)=(TH1D*)f_0->Get("h_aco_1");
h_th_min_1.at(f)=(TH1D*)f_0->Get("h_th_min_1");
h_th_max_1.at(f)=(TH1D*)f_0->Get("h_th_max_1");
h_preE_th_min_1.at(f)=(TH1D*)f_0->Get("h_preE_th_min_1");
h_preE_th_max_1.at(f)=(TH1D*)f_0->Get("h_preE_th_max_1");
h_elastic_1.at(f)=(TH1D*)f_0->Get("h_elastic_1");
h_el_2.at(f)=(TH1D*)f_0->Get("h_el_2");
h_E_th_min_1.at(f)=(TH1D*)f_0->Get("h_E_th_min_1");
h_E_th_max_1.at(f)=(TH1D*)f_0->Get("h_E_th_max_1");
h_final_nhits_pre_1.at(f)=(TH1D*)f_0->Get("h_final_nhits_pre_1");
h_final_nhits_post_1.at(f)=(TH1D*)f_0->Get("h_final_nhits_post_1");
h_final_ntracks_post_1.at(f)=(TH1D*)f_0->Get("h_final_ntracks_post_1");
h_final_chi2_vrtx_1.at(f)=(TH1D*)f_0->Get("h_final_chi2_vrtx_1");
h_final_Dz_vrtx_1.at(f)=(TH1D*)f_0->Get("h_final_Dz_vrtx_1");
h_final_aco_1.at(f)=(TH1D*)f_0->Get("h_final_aco_1");
h_final_chi2_tr0_1.at(f)=(TH1D*)f_0->Get("h_final_chi2_tr0_1");
h_final_chi2_trMax_1.at(f)=(TH1D*)f_0->Get("h_final_chi2_trMax_1");
h_final_chi2_trMin_1.at(f)=(TH1D*)f_0->Get("h_final_chi2_trMin_1");

	if(f>0){
h_x_1.at(0)->Add(h_x_1.at(f));
h_y_1.at(0)->Add(h_y_1.at(f));
h_thx_1.at(0)->Add(h_thx_1.at(f));
h_thy_1.at(0)->Add(h_thy_1.at(f));
h_nhits_pre_1.at(0)->Add(h_nhits_pre_1.at(f));
h_nhits_post_1.at(0)->Add(h_nhits_post_1.at(f));
h_ntracks_post_1.at(0)->Add(h_ntracks_post_1.at(f));
h_chi2_vrtx_1.at(0)->Add(h_chi2_vrtx_1.at(f));
h_Dz_vrtx_1.at(0)->Add(h_Dz_vrtx_1.at(f));
h_aco_1.at(0)->Add(h_aco_1.at(f));
h_th_min_1.at(0)->Add(h_th_min_1.at(f));
h_th_max_1.at(0)->Add(h_th_max_1.at(f));
h_preE_th_min_1.at(0)->Add(h_preE_th_min_1.at(f));
h_preE_th_max_1.at(0)->Add(h_preE_th_max_1.at(f));
h_elastic_1.at(0)->Add(h_elastic_1.at(f));
h_el_2.at(0)->Add(h_el_2.at(f));
h_E_th_min_1.at(0)->Add(h_E_th_min_1.at(f));
h_E_th_max_1.at(0)->Add(h_E_th_max_1.at(f));
h_final_nhits_pre_1.at(0)->Add(h_final_nhits_pre_1.at(f));
h_final_nhits_post_1.at(0)->Add(h_final_nhits_post_1.at(f));
h_final_ntracks_post_1.at(0)->Add(h_final_ntracks_post_1.at(f));
h_final_chi2_vrtx_1.at(0)->Add(h_final_chi2_vrtx_1.at(f));
h_final_Dz_vrtx_1.at(0)->Add(h_final_Dz_vrtx_1.at(f));
h_final_aco_1.at(0)->Add(h_final_aco_1.at(f));
h_final_chi2_tr0_1.at(0)->Add(h_final_chi2_tr0_1.at(f));
h_final_chi2_trMax_1.at(0)->Add(h_final_chi2_trMax_1.at(f));
h_final_chi2_trMin_1.at(0)->Add(h_final_chi2_trMin_1.at(f));

		}
}



h_x_0.at(0)->Scale(1./h_x_0.at(0)->Integral());
h_y_0.at(0)->Scale(1./h_y_0.at(0)->Integral());
h_thx_0.at(0)->Scale(1./h_thx_0.at(0)->Integral());
h_thy_0.at(0)->Scale(1./h_thy_0.at(0)->Integral());
h_nhits_pre_0.at(0)->Scale(1./h_nhits_pre_0.at(0)->Integral());
h_nhits_post_0.at(0)->Scale(1./h_nhits_post_0.at(0)->Integral());
h_ntracks_post_0.at(0)->Scale(1./h_ntracks_post_0.at(0)->Integral());
h_chi2_vrtx_0.at(0)->Scale(1./h_chi2_vrtx_0.at(0)->Integral());
h_Dz_vrtx_0.at(0)->Scale(1./h_Dz_vrtx_0.at(0)->Integral());
h_aco_0.at(0)->Scale(1./h_aco_0.at(0)->Integral());
h_th_min_0.at(0)->Scale(1./h_th_min_0.at(0)->Integral());
h_th_max_0.at(0)->Scale(1./h_th_max_0.at(0)->Integral());
h_preE_th_min_0.at(0)->Scale(1./h_preE_th_min_0.at(0)->Integral());
h_preE_th_max_0.at(0)->Scale(1./h_preE_th_max_0.at(0)->Integral());
h_elastic_0.at(0)->Scale(1./h_elastic_0.at(0)->Integral());
h_el_1.at(0)->Scale(1./h_el_1.at(0)->Integral());
h_E_th_min_0.at(0)->Scale(1./h_E_th_min_0.at(0)->Integral());
h_E_th_max_0.at(0)->Scale(1./h_E_th_max_0.at(0)->Integral());
h_final_nhits_pre_0.at(0)->Scale(1./h_final_nhits_pre_0.at(0)->Integral());
h_final_nhits_post_0.at(0)->Scale(1./h_final_nhits_post_0.at(0)->Integral());
h_final_ntracks_post_0.at(0)->Scale(1./h_final_ntracks_post_0.at(0)->Integral());
h_final_chi2_vrtx_0.at(0)->Scale(1./h_final_chi2_vrtx_0.at(0)->Integral());
h_final_Dz_vrtx_0.at(0)->Scale(1./h_final_Dz_vrtx_0.at(0)->Integral());
h_final_aco_0.at(0)->Scale(1./h_final_aco_0.at(0)->Integral());
h_final_chi2_tr0_0.at(0)->Scale(1./h_final_chi2_tr0_0.at(0)->Integral());
h_final_chi2_trMax_0.at(0)->Scale(1./h_final_chi2_trMax_0.at(0)->Integral());
h_final_chi2_trMin_0.at(0)->Scale(1./h_final_chi2_trMin_0.at(0)->Integral());

h_x_1.at(0)->Scale(1./h_x_1.at(0)->Integral());
h_y_1.at(0)->Scale(1./h_y_1.at(0)->Integral());
h_thx_1.at(0)->Scale(1./h_thx_1.at(0)->Integral());
h_thy_1.at(0)->Scale(1./h_thy_1.at(0)->Integral());
h_nhits_pre_1.at(0)->Scale(1./h_nhits_pre_1.at(0)->Integral());
h_nhits_post_1.at(0)->Scale(1./h_nhits_post_1.at(0)->Integral());
h_ntracks_post_1.at(0)->Scale(1./h_ntracks_post_1.at(0)->Integral());
h_chi2_vrtx_1.at(0)->Scale(1./h_chi2_vrtx_1.at(0)->Integral());
h_Dz_vrtx_1.at(0)->Scale(1./h_Dz_vrtx_1.at(0)->Integral());
h_aco_1.at(0)->Scale(1./h_aco_1.at(0)->Integral());
h_th_min_1.at(0)->Scale(1./h_th_min_1.at(0)->Integral());
h_th_max_1.at(0)->Scale(1./h_th_max_1.at(0)->Integral());
h_preE_th_min_1.at(0)->Scale(1./h_preE_th_min_1.at(0)->Integral());
h_preE_th_max_1.at(0)->Scale(1./h_preE_th_max_1.at(0)->Integral());
h_elastic_1.at(0)->Scale(1./h_elastic_1.at(0)->Integral());
h_el_2.at(0)->Scale(1./h_el_2.at(0)->Integral());
h_E_th_min_1.at(0)->Scale(1./h_E_th_min_1.at(0)->Integral());
h_E_th_max_1.at(0)->Scale(1./h_E_th_max_1.at(0)->Integral());
h_final_nhits_pre_1.at(0)->Scale(1./h_final_nhits_pre_1.at(0)->Integral());
h_final_nhits_post_1.at(0)->Scale(1./h_final_ntracks_post_1.at(0)->Integral());
h_final_ntracks_post_1.at(0)->Scale(1./h_final_ntracks_post_1.at(0)->Integral());
h_final_chi2_vrtx_1.at(0)->Scale(1./h_final_chi2_vrtx_1.at(0)->Integral());
h_final_Dz_vrtx_1.at(0)->Scale(1./h_final_Dz_vrtx_1.at(0)->Integral());
h_final_aco_1.at(0)->Scale(1./h_final_aco_1.at(0)->Integral());
h_final_chi2_tr0_1.at(0)->Scale(1./h_final_chi2_tr0_1.at(0)->Integral());
h_final_chi2_trMax_1.at(0)->Scale(1./h_final_chi2_trMax_1.at(0)->Integral());
h_final_chi2_trMin_1.at(0)->Scale(1./h_final_chi2_trMin_1.at(0)->Integral());


TCanvas pos("pos","pos",1500,1500);
pos.Divide(2,2);
pos.cd(1);
h_x_0.at(0)->Draw("hist");
h_x_1.at(0)->SetLineColor(kRed);
h_x_1.at(0)->Draw("hist sames");
pos.cd(2);
h_y_0.at(0)->Draw("hist");
h_y_1.at(0)->SetLineColor(kRed);
h_y_1.at(0)->Draw("hist sames");
pos.cd(3);
h_thx_0.at(0)->Draw("hist");
h_thx_1.at(0)->SetLineColor(kRed);
h_thx_1.at(0)->Draw("hist sames");
pos.cd(4);
h_thy_0.at(0)->Draw("hist");
h_thy_1.at(0)->SetLineColor(kRed);
h_thy_1.at(0)->Draw("hist sames");
pos.SaveAs(Form("plots/run%s/zin3_preFid_postion.pdf",run));

TCanvas gen("gen","gen",1500,1500);
gen.Divide(3,3);
gen.cd(1);
h_nhits_pre_0.at(0)->Draw("hist");
h_nhits_pre_1.at(0)->SetLineColor(kRed);
h_nhits_pre_1.at(0)->Draw("hist sames");
gen.cd(2);
h_nhits_post_0.at(0)->Draw("hist");
h_nhits_post_1.at(0)->SetLineColor(kRed);
h_nhits_post_1.at(0)->Draw("hist sames");
gen.cd(3);
h_ntracks_post_0.at(0)->Draw("hist");
h_ntracks_post_1.at(0)->SetLineColor(kRed);
h_ntracks_post_1.at(0)->Draw("hist sames");
gen.cd(4);
h_chi2_vrtx_0.at(0)->Draw("hist");
h_chi2_vrtx_1.at(0)->SetLineColor(kRed);
h_chi2_vrtx_1.at(0)->Draw("hist sames");
gen.cd(5);
h_Dz_vrtx_0.at(0)->Draw("hist");
h_Dz_vrtx_1.at(0)->SetLineColor(kRed);
h_Dz_vrtx_1.at(0)->Draw("hist sames");
gen.cd(6);
h_aco_0.at(0)->Draw("hist");
h_aco_1.at(0)->SetLineColor(kRed);
h_aco_1.at(0)->Draw("hist sames");
gen.cd(7);
h_th_min_0.at(0)->Draw("hist");
h_th_min_1.at(0)->SetLineColor(kRed);
h_th_min_1.at(0)->Draw("hist sames");
gen.cd(8);
h_th_max_0.at(0)->Draw("hist");
h_th_max_1.at(0)->SetLineColor(kRed);
h_th_max_1.at(0)->Draw("hist sames");
gen.SaveAs(Form("plots/run%s/zin3_postFid.pdf",run));

TCanvas pe("pe","pe",1500,1500);
pe.Divide(3,3);
pe.cd(1);
h_preE_th_min_0.at(0)->GetXaxis()->SetRangeUser(0.,0.002);
h_preE_th_min_1.at(0)->GetXaxis()->SetRangeUser(0.,0.002);
h_preE_th_min_0.at(0)->Draw("hist");
h_preE_th_min_1.at(0)->SetLineColor(kRed);
h_preE_th_min_1.at(0)->Draw("hist sames");
pe.cd(2);
h_preE_th_max_0.at(0)->GetXaxis()->SetRangeUser(0.005,0.02);
h_preE_th_max_1.at(0)->GetXaxis()->SetRangeUser(0.005,0.02);
h_preE_th_max_0.at(0)->Draw("hist");
h_preE_th_max_1.at(0)->SetLineColor(kRed);
h_preE_th_max_1.at(0)->Draw("hist sames");
pe.cd(3);
h_elastic_0.at(0)->Draw("hist");
h_elastic_1.at(0)->SetLineColor(kRed);
h_elastic_1.at(0)->Draw("hist sames");


h_preE_th_min_0.at(0)->Rebin(5);
h_preE_th_min_1.at(0)->Rebin(5);
h_preE_th_max_0.at(0)->Rebin(4);
h_preE_th_max_1.at(0)->Rebin(4);
pe.cd(4);
TH1D * h1 = (TH1D*) h_preE_th_min_0.at(0)->Clone();
h1->Divide(h_preE_th_min_1.at(0));
h1->GetXaxis()->SetRangeUser(0.,0.002);
h1->SetMinimum(0.2);
h1->SetMaximum(2);
h1->SetTitle("#theta_min Tar_0/Tar_1");
h1->Draw();

pe.cd(5);
TH1D * h2 = (TH1D*) h_preE_th_max_0.at(0)->Clone();
h2->Divide(h_preE_th_max_1.at(0));
h2->GetXaxis()->SetRangeUser(0.005,0.02);
h2->SetMinimum(0.2);
h2->SetMaximum(2);
h2->SetTitle("#theta_max Tar_0/Tar_1");
h2->Draw();

pe.cd(7);
h_el_1.at(0)->Draw("COLZ");
pe.cd(8);
h_el_2.at(0)->Draw("COLZ");

pe.SaveAs(Form("plots/run%s/zin3_preEl.pdf",run));
pe.SaveAs(Form("plots/run%s/zin3_preEl.root",run));



TCanvas fe("fe","fe",1000,1000);
fe.Divide(2,2);
fe.cd(1);
h_E_th_min_0.at(0)->GetXaxis()->SetRangeUser(0.,0.002);
h_E_th_min_1.at(0)->GetXaxis()->SetRangeUser(0.,0.002);
h_E_th_min_0.at(0)->Draw("hist");
h_E_th_min_1.at(0)->SetLineColor(kRed);
h_E_th_min_1.at(0)->Draw("hist sames");
fe.cd(2);
h_E_th_max_0.at(0)->GetXaxis()->SetRangeUser(0.005,0.02);
h_E_th_max_1.at(0)->GetXaxis()->SetRangeUser(0.005,0.02);
h_E_th_max_0.at(0)->Draw("hist");
h_E_th_max_1.at(0)->SetLineColor(kRed);
h_E_th_max_1.at(0)->Draw("hist sames");

h_E_th_min_0.at(0)->Rebin(5);
h_E_th_min_1.at(0)->Rebin(5);
h_E_th_max_0.at(0)->Rebin(4);
h_E_th_max_1.at(0)->Rebin(4);

fe.cd(3);
TH1D * h1_el = (TH1D*) h_E_th_min_0.at(0)->Clone();
h1_el->Divide(h_E_th_min_1.at(0));
h1_el->GetXaxis()->SetRangeUser(0.,0.002);
h1_el->SetMinimum(0.2);
h1_el->SetMaximum(2);
h1_el->SetTitle("#theta_min Tar_0/Tar_1");
h1_el->Draw();

fe.cd(4);
TH1D * h2_el = (TH1D*) h_E_th_max_0.at(0)->Clone();
h2_el->Divide(h_E_th_max_1.at(0));
h2_el->GetXaxis()->SetRangeUser(0.005,0.02);
h2_el->SetMinimum(0.2);
h2_el->SetMaximum(2);
h2_el->SetTitle("#theta_max Tar_0/Tar_1");
h2_el->Draw();

fe.SaveAs(Form("plots/run%s/zin3_finalEl.pdf",run));



TCanvas gen_fin("gen_fin","gen_fin",2000,1500);
gen_fin.Divide(3,3);
gen_fin.cd(1);
h_final_nhits_pre_0.at(0)->Draw("hist");
h_final_nhits_pre_1.at(0)->SetLineColor(kRed);
h_final_nhits_pre_1.at(0)->Draw("hist sames");
gen_fin.cd(2);
h_final_nhits_post_0.at(0)->Draw("hist");
h_final_nhits_post_1.at(0)->SetLineColor(kRed);
h_final_nhits_post_1.at(0)->Draw("hist sames");
gen_fin.cd(3);
h_final_ntracks_post_0.at(0)->Draw("hist");
h_final_ntracks_post_1.at(0)->SetLineColor(kRed);
h_final_ntracks_post_1.at(0)->Draw("hist sames");
gen_fin.cd(4);
h_final_chi2_vrtx_0.at(0)->Draw("hist");
h_final_chi2_vrtx_1.at(0)->SetLineColor(kRed);
h_final_chi2_vrtx_1.at(0)->Draw("hist sames");
gen_fin.cd(5);
h_final_Dz_vrtx_0.at(0)->Draw("hist");
h_final_Dz_vrtx_1.at(0)->SetLineColor(kRed);
h_final_Dz_vrtx_1.at(0)->Draw("hist sames");
gen_fin.cd(6);
h_final_aco_0.at(0)->Draw("hist");
h_final_aco_1.at(0)->SetLineColor(kRed);
h_final_aco_1.at(0)->Draw("hist sames");
gen_fin.cd(7);
h_final_chi2_tr0_0.at(0)->Draw("hist");
h_final_chi2_tr0_1.at(0)->SetLineColor(kRed);
h_final_chi2_tr0_1.at(0)->Draw("hist sames");
gen_fin.cd(8);
h_final_chi2_trMax_0.at(0)->Draw("hist");
h_final_chi2_trMax_1.at(0)->SetLineColor(kRed);
h_final_chi2_trMax_1.at(0)->Draw("hist sames");
gen_fin.cd(9);
h_final_chi2_trMin_0.at(0)->Draw("hist");
h_final_chi2_trMin_1.at(0)->SetLineColor(kRed);
h_final_chi2_trMin_1.at(0)->Draw("hist sames");


gen_fin.SaveAs(Form("plots/run%s/zin3_final.pdf",run));
gen_fin.SaveAs(Form("plots/run%s/zin3_final.root",run));



return 0;

}
