#include <iostream>
#include <fstream>
#include <string>
#include <regex>
#include <filesystem>
#include <vector>
#include <algorithm>

//#include "TApplication.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TPaveStats.h"

namespace fs = std::filesystem;
/*
long extract_number(const std::string& filename) {
	std::regex re("-(\\d+)"); std::smatch match;
	if (std::regex_search(filename, match, re))
	 { return std::stoi(match[1].str()); }
	return -1;
}
*/

long extract_number(const std::string& filename) {
//    std::regex re(".*-\\d+-(\\d+)_nhits\\d+\\.txt$");
    std::regex re(".*-(\\d+)_nhits\\d+\\.txt$");
    std::smatch match;
    if (std::regex_match(filename, match, re)) {
        return std::stol(match[1].str());
    }
    return -1;
}

int main(int argc, char* argv[]) {

    if (argc < 2) {
        std::cerr << "Uso: " << argv[0] << " <cartella_con_txt>" << std::endl;
        return 1;
    }

    fs::path folder = argv[1];
    if (!fs::exists(folder) || !fs::is_directory(folder)) {
        std::cerr << "Cartella non valida: " << folder << std::endl;
        return 1;
    }

    std::vector<fs::path> files;

    for (auto& entry : fs::directory_iterator(folder)) {
        if (fs::is_regular_file(entry.path()) && entry.path().extension() == ".txt" and entry.path().string().find("wip11x") != std::string::npos) {
            files.push_back(entry.path());
        }
    }

    // Ordina per numero dopo "muedaq01-"
    std::sort(files.begin(), files.end(), [](const fs::path& a, const fs::path& b) {
        return extract_number(a.filename().string()) < extract_number(b.filename().string());
    });

    // Regex per i valori
    std::regex pattern_fid_tr0(R"(Ratio fiducial tr0 over entries:\s+([0-9]*\.?[0-9]+))");
    std::regex pattern_fid_tr1(R"(Ratio fiducial tr1 over entries:\s+([0-9]*\.?[0-9]+))");

    // vettori per TGraph
    std::vector<double> x0_fid; // numero estratto dal nome
    std::vector<double> y0_fid; // valore Ratio good fid
    std::vector<double> x1_fid; // numero estratto dal nome
    std::vector<double> y1_fid; // valore Ratio good fid

    // Regex per i valori
    std::regex pattern_vrtx_tr0(R"(Ratio any vrtx tr0 over fiducial:\s+([0-9]*\.?[0-9]+))");
    std::regex pattern_vrtx_tr1(R"(Ratio any vrtx tr1 over fiducial:\s+([0-9]*\.?[0-9]+))");

    // vettori per TGraph
    std::vector<double> x0_vrtx; // numero estratto dal nome
    std::vector<double> y0_vrtx; // valore Ratio any vrtx
    std::vector<double> x1_vrtx; // numero estratto dal nome
    std::vector<double> y1_vrtx; // valore Ratio any vrtx

    // Regex per i valori
//    std::regex pattern_preEl_tr0(R"(Ratio pre_elastic tr0 over best vrtx:\s+([0-9]*\.?[0-9]+))");
//    std::regex pattern_preEl_tr1(R"(Ratio pre_elastic tr1 over best vrtx:\s+([0-9]*\.?[0-9]+))");
    std::regex pattern_preEl_tr0(R"(Ratio pre_elastic tr0 over fiducial:\s+([0-9]*\.?[0-9]+))");
    std::regex pattern_preEl_tr1(R"(Ratio pre_elastic tr1 over fiducial:\s+([0-9]*\.?[0-9]+))");

    // vettori per TGraph
    std::vector<double> x0_preEl; // numero estratto dal nome
    std::vector<double> y0_preEl; // valore Ratio any vrtx
    std::vector<double> x1_preEl; // numero estratto dal nome
    std::vector<double> y1_preEl; // valore Ratio any vrtx

    // Regex per i valori
//    std::regex pattern_el_tr0(R"(Ratio elastic tr0 over pre_elastic:\s+([0-9]*\.?[0-9]+))");
//    std::regex pattern_el_tr1(R"(Ratio elastic tr1 over pre_elastic:\s+([0-9]*\.?[0-9]+))");
    std::regex pattern_el_tr0(R"(Ratio elastic tr0 over fiducial:\s+([0-9]*\.?[0-9]+))");
    std::regex pattern_el_tr1(R"(Ratio elastic tr1 over fiducial:\s+([0-9]*\.?[0-9]+))");

    // vettori per TGraph
    std::vector<double> x0_el; // numero estratto dal nome
    std::vector<double> y0_el; // valore Ratio any vrtx
    std::vector<double> x1_el; // numero estratto dal nome
    std::vector<double> y1_el; // valore Ratio any vrtx



    std::regex pattern_num_fid_tr0(R"(Events with fiducial 0\s+([0-9]*\.?[0-9]+))");
    std::regex pattern_num_fid_tr1(R"(Events with fiducial 1\s+([0-9]*\.?[0-9]+))");
//    std::regex pattern_num_fid_tr0(R"(Events with reco vrtx in st01\s+([0-9]*\.?[0-9]+))");
//    std::regex pattern_num_fid_tr1(R"(Events with reco vrtx in st12\s+([0-9]*\.?[0-9]+))");

    // vettori per TGraph
    std::vector<double> x0_num_fid; // numero estratto dal nome
    std::vector<double> y0_num_fid; // valore Ratio any vrtx
    std::vector<double> x1_num_fid; // numero estratto dal nome
    std::vector<double> y1_num_fid; // valore Ratio any vrtx

    std::regex pattern_num_el_tr0(R"(Events elastic 0\s+([0-9]*\.?[0-9]+))");
    std::regex pattern_num_el_tr1(R"(Events elastic 1\s+([0-9]*\.?[0-9]+))");

    // vettori per TGraph
    std::vector<double> x0_num_el; // numero estratto dal nome
    std::vector<double> y0_num_el; // valore Ratio any vrtx
    std::vector<double> x1_num_el; // numero estratto dal nome
    std::vector<double> y1_num_el; // valore Ratio any vrtx


    std::regex pattern_num_bkg_tr0(R"(Events with IP>0.2cm and op_angle<2mrad tar0\s+([0-9]*\.?[0-9]+))");
    std::regex pattern_num_bkg_tr1(R"(Events with IP>0.2cm and op_angle<2mrad tar1\s+([0-9]*\.?[0-9]+))");

    // vettori per TGraph
    std::vector<double> x0_num_bkg; // numero estratto dal nome
    std::vector<double> y0_num_bkg; // valore Ratio any vrtx
    std::vector<double> x1_num_bkg; // numero estratto dal nome
    std::vector<double> y1_num_bkg; // valore Ratio any vrtx

    std::regex pattern_num_sig_tr0(R"(Events with IP<0.2cm and op_angle>2mrad tar0\s+([0-9]*\.?[0-9]+))");
    std::regex pattern_num_sig_tr1(R"(Events with IP<0.2cm and op_angle>2mrad tar1\s+([0-9]*\.?[0-9]+))");

    // vettori per TGraph
    std::vector<double> x0_num_sig; // numero estratto dal nome
    std::vector<double> y0_num_sig; // valore Ratio any vrtx
    std::vector<double> x1_num_sig; // numero estratto dal nome
    std::vector<double> y1_num_sig; // valore Ratio any vrtx




    for (const auto& file : files) {
        std::ifstream infile(file);
        if (!infile) {
            std::cerr << "Errore nell'apertura del file: " << file << std::endl;
            continue;
        }

        std::string line;
        std::smatch match;
        bool trovato = false;

        std::string filename = file.filename().string();
	std::time_t num = extract_number(filename);

//std::cout <<"num " << num << std::endl;
        if (filename.find("single_muon_interaction_0") != std::string::npos) {
            while (std::getline(infile, line)) {
                if (std::regex_search(line, match, pattern_num_fid_tr0) && match.size() > 1) {
                    //std::cout << filename << " -> Valore num_fid tr0: " << match[1] << std::endl;
                    x0_num_fid.push_back(num);
                    y0_num_fid.push_back(std::stod(match[1].str()));
                    trovato = true;
                    continue;
                }
                else if (std::regex_search(line, match, pattern_num_el_tr0) && match.size() > 1) {
                    //std::cout << filename << " -> Valore num_el tr0: " << match[1] << std::endl;
                    x0_num_el.push_back(num);
                    y0_num_el.push_back(std::stod(match[1].str()));
                    trovato = true;
                    continue;
                }
                else if (std::regex_search(line, match, pattern_num_bkg_tr0) && match.size() > 1) {
                    //std::cout << filename << " -> Valore num_bkg tr0: " << match[1] << std::endl;
                    x0_num_bkg.push_back(num);
                    y0_num_bkg.push_back(std::stod(match[1].str()));
                    trovato = true;
                    continue;
                }
                else if (std::regex_search(line, match, pattern_num_sig_tr0) && match.size() > 1) {
                    //std::cout << filename << " -> Valore num_sig tr0: " << match[1] << std::endl;
                    x0_num_sig.push_back(num);
                    y0_num_sig.push_back(std::stod(match[1].str()));
                    trovato = true;
                    continue;
                }
                else if (std::regex_search(line, match, pattern_fid_tr0) && match.size() > 1) {
//                    std::cout << filename << " -> Valore fid tr0: " << match[1] << std::endl;
                    x0_fid.push_back(num);
                    y0_fid.push_back(std::stod(match[1].str()));
                    trovato = true;
                    continue;
                }
		else if (std::regex_search(line, match, pattern_vrtx_tr0) && match.size() > 1) {
//                    std::cout << filename << " -> Valore vrtx tr0: " << match[1] << std::endl;
                    x0_vrtx.push_back(num);
                    y0_vrtx.push_back(std::stod(match[1].str()));
                    trovato = true;
                    continue;
                }
                else if (std::regex_search(line, match, pattern_preEl_tr0) && match.size() > 1) {
//                    std::cout << filename << " -> Valore preEl tr0: " << match[1] << std::endl;
                    x0_preEl.push_back(num);
                    y0_preEl.push_back(std::stod(match[1].str()));
                    trovato = true;
                    continue;
                }
                else if (std::regex_search(line, match, pattern_el_tr0) && match.size() > 1) {
//                    std::cout << filename << " -> Valore el tr0: " << match[1] << std::endl;
                    x0_el.push_back(num);
                    y0_el.push_back(std::stod(match[1].str()));
                    trovato = true;
                    break;
                }
            }
        } else if (filename.find("single_muon_interaction_1") != std::string::npos) {
            while (std::getline(infile, line)) {
                if (std::regex_search(line, match, pattern_num_fid_tr1) && match.size() > 1) {
                    //std::cout << filename << " -> Valore num_fid tr1: " << match[1] << std::endl;
                    x1_num_fid.push_back(num);
                    y1_num_fid.push_back(std::stod(match[1].str()));
                    trovato = true;
                    continue;
                }
                else if (std::regex_search(line, match, pattern_num_el_tr1) && match.size() > 1) {
                    //std::cout << filename << " -> Valore num_el tr1: " << match[1] << std::endl;
                    x1_num_el.push_back(num);
                    y1_num_el.push_back(std::stod(match[1].str()));
                    trovato = true;
                    continue;
                }
                else if (std::regex_search(line, match, pattern_num_bkg_tr1) && match.size() > 1) {
                    //std::cout << filename << " -> Valore num_bkg tr1: " << match[1] << std::endl;
                    x1_num_bkg.push_back(num);
                    y1_num_bkg.push_back(std::stod(match[1].str()));
                    trovato = true;
                    continue;
                }
                else if (std::regex_search(line, match, pattern_num_sig_tr1) && match.size() > 1) {
                    //std::cout << filename << " -> Valore num_sig tr1: " << match[1] << std::endl;
                    x1_num_sig.push_back(num);
                    y1_num_sig.push_back(std::stod(match[1].str()));
                    trovato = true;
                    continue;
                }
                else if (std::regex_search(line, match, pattern_fid_tr1) && match.size() > 1) {
//                    std::cout << filename << " -> Valore fid t10: " << match[1] << std::endl;
                    x1_fid.push_back(num);
                    y1_fid.push_back(std::stod(match[1].str()));
                    trovato = true;
                    continue;
                }
                else if (std::regex_search(line, match, pattern_vrtx_tr1) && match.size() > 1) {
//                    std::cout << filename << " -> Valore vrtx tr1: " << match[1] << std::endl;
                    x1_vrtx.push_back(num);
                    y1_vrtx.push_back(std::stod(match[1].str()));
                    trovato = true;
                    continue;
                }
                else if (std::regex_search(line, match, pattern_preEl_tr1) && match.size() > 1) {
//                    std::cout << filename << " -> Valore preEl tr1: " << match[1] << std::endl;
                    x1_preEl.push_back(num);
                    y1_preEl.push_back(std::stod(match[1].str()));
                    trovato = true;
                    continue;
                }
                else if (std::regex_search(line, match, pattern_el_tr1) && match.size() > 1) {
//                    std::cout << filename << " -> Valore el tr1: " << match[1] << std::endl;
                    x1_el.push_back(num);
                    y1_el.push_back(std::stod(match[1].str()));
                    trovato = true;
                    break;
                }
            }
        }
 }


double tot_fid_0=0.;
double tot_fid_1=0.;
double tot_el_0=0.;
double tot_el_1=0.;
double tot_bkg_0=0.;
double tot_bkg_1=0.;
double tot_sig_0=0.;
double tot_sig_1=0.;

for(int h=0; h<y0_num_fid.size(); h++)tot_fid_0+=y0_num_fid.at(h);
for(int h=0; h<y1_num_fid.size(); h++)tot_fid_1+=y1_num_fid.at(h);
for(int h=0; h<y0_num_el.size(); h++)tot_el_0+=y0_num_el.at(h);
for(int h=0; h<y1_num_el.size(); h++)tot_el_1+=y1_num_el.at(h);
for(int h=0; h<y0_num_bkg.size(); h++)tot_bkg_0+=y0_num_bkg.at(h);
for(int h=0; h<y1_num_bkg.size(); h++)tot_bkg_1+=y1_num_bkg.at(h);
for(int h=0; h<y0_num_sig.size(); h++)tot_sig_0+=y0_num_sig.at(h);
for(int h=0; h<y1_num_sig.size(); h++)tot_sig_1+=y1_num_sig.at(h);

std::cout<<"Elastic vrtx0 " << tot_el_0 << " +- " << sqrt(tot_el_0) << std::endl;
std::cout<<"Elastic vrtx1 " << tot_el_1 << " +- " << sqrt(tot_el_1) << std::endl;

std::cout<<"FIducial vrtx0 " << tot_fid_0 << " +- " << sqrt(tot_fid_0) << std::endl;
std::cout<<"FIducial vrtx1 " << tot_fid_1 << " +- " << sqrt(tot_fid_1) << std::endl;


std::cout<<"Total ratio el/fid vrtx0 " << tot_el_0/tot_fid_0 << " +- " <<
sqrt( ( 1./tot_fid_0 * sqrt(tot_el_0) )*(1./tot_fid_0 * sqrt(tot_el_0)) + ( tot_el_0/(tot_fid_0*tot_fid_0) * sqrt(tot_fid_0) )*( tot_el_0/(tot_fid_0*tot_fid_0) * sqrt(tot_fid_0) ) )
<< std::endl;
std::cout<<"Total ratio el/fid vrtx1 " << tot_el_1/tot_fid_1 << " +- " <<
sqrt( ( 1./tot_fid_1 * sqrt(tot_el_1) )*(1./tot_fid_1 * sqrt(tot_el_1)) + ( tot_el_1/(tot_fid_1*tot_fid_1) * sqrt(tot_fid_1) )*( tot_el_1/(tot_fid_1*tot_fid_1) * sqrt(tot_fid_1) ) )
<< std::endl;

std::cout<<"Total ratio bkg/fid vrtx0 " << tot_bkg_0/tot_fid_0 << " +- " <<
sqrt( ( 1./tot_fid_0 * sqrt(tot_bkg_0) )*(1./tot_fid_0 * sqrt(tot_bkg_0)) + ( tot_bkg_0/(tot_fid_0*tot_fid_0) * sqrt(tot_fid_0) )*( tot_bkg_0/(tot_fid_0*tot_fid_0) * sqrt(tot_fid_0) ) )
<< std::endl;
std::cout<<"Total ratio bkg/fid vrtx1 " << tot_bkg_1/tot_fid_1 << " +- " <<
sqrt( ( 1./tot_fid_1 * sqrt(tot_bkg_1) )*(1./tot_fid_1 * sqrt(tot_bkg_1)) + ( tot_bkg_1/(tot_fid_1*tot_fid_1) * sqrt(tot_fid_1) )*( tot_bkg_1/(tot_fid_1*tot_fid_1) * sqrt(tot_fid_1) ) )
<< std::endl;

std::cout<<"Total ratio sig/fid vrtx0 " << tot_sig_0/tot_fid_0 << " +- " <<
sqrt( ( 1./tot_fid_0 * sqrt(tot_sig_0) )*(1./tot_fid_0 * sqrt(tot_sig_0)) + ( tot_sig_0/(tot_fid_0*tot_fid_0) * sqrt(tot_fid_0) )*( tot_sig_0/(tot_fid_0*tot_fid_0) * sqrt(tot_fid_0) ) )
<< std::endl;
std::cout<<"Total ratio sig/fid vrtx1 " << tot_sig_1/tot_fid_1 << " +- " <<
sqrt( ( 1./tot_fid_1 * sqrt(tot_sig_1) )*(1./tot_fid_1 * sqrt(tot_sig_1)) + ( tot_sig_1/(tot_fid_1*tot_fid_1) * sqrt(tot_fid_1) )*( tot_sig_1/(tot_fid_1*tot_fid_1) * sqrt(tot_fid_1) ) )
<< std::endl;




TH1D* h_fid0 = new TH1D("h_fid0","fiducial/online_trigger vrtx0",200,0.,1.);
TH1D* h_fid1 = new TH1D("h_fid1","fiducial/online_trigger vrtx1",200,0.,1.);

TH1D* h_vrtx0 = new TH1D("h_vrtx0","reconstructed/fiducial vrtx0",200,0.,1.);
TH1D* h_vrtx1 = new TH1D("h_vrtx1","reconstructed/fiducial vrtx1",200,0.,1.);

/*

TH1D* h_preEl0 = new TH1D("h_preEl0","pre_elastic/any_vrtx1",200,0.,1.);
TH1D* h_preEl1 = new TH1D("h_preEl1","pre_elastic/any_vrtx0",200,0.,1.);

TH1D* h_el0 = new TH1D("h_el0","elastic/pre_elastic vrtx1",200,0.,1.);
TH1D* h_el1 = new TH1D("h_el1","elastic/pre_elastic vrtx0",200,0.,1.);
*/


TH1D* h_preEl0 = new TH1D("h_preEl0","pre_elastic/fiducial vrtx1",200,0.,1.);
TH1D* h_preEl1 = new TH1D("h_preEl1","pre_elastic/fiducial vrtx0",200,0.,1.);

TH1D* h_el0 = new TH1D("h_el0","elastic/fiducial vrtx1",200,0.,1.);
TH1D* h_el1 = new TH1D("h_el1","elastic/fiducial vrtx0",200,0.,1.);


for(int h=0; h<y0_fid.size(); h++)h_fid0->Fill(y0_fid.at(h));
for(int h=0; h<y1_fid.size(); h++)h_fid1->Fill(y1_fid.at(h));

for(int h=0; h<y0_vrtx.size(); h++)h_vrtx0->Fill(y0_vrtx.at(h));
for(int h=0; h<y1_vrtx.size(); h++)h_vrtx1->Fill(y1_vrtx.at(h));

for(int h=0; h<y0_preEl.size(); h++)h_preEl0->Fill(y0_preEl.at(h));
for(int h=0; h<y1_preEl.size(); h++)h_preEl1->Fill(y1_preEl.at(h));

for(int h=0; h<y0_el.size(); h++)h_el0->Fill(y0_el.at(h));
for(int h=0; h<y1_el.size(); h++)h_el1->Fill(y1_el.at(h));

/*
    TCanvas* b = new TCanvas("b", "Hist", 4000, 3000);
    b->Divide(2,2);
    b->cd(1);
h_fid0->SetLineColor(kRed);
h_fid0->Draw("hist");
h_fid1->SetLineColor(kOrange);
h_fid1->Draw("hist sames");
    b->cd(2);
h_vrtx0->SetLineColor(kBlue);
h_vrtx0->Draw("hist");
h_vrtx1->SetLineColor(kCyan);
h_vrtx1->Draw("hist same");
    b->cd(3);
h_preEl0->SetLineColor(kGreen+2);
h_preEl0->Draw("hist");
h_preEl1->SetLineColor(kGreen);
h_preEl1->Draw("hist same");
    b->cd(4);
h_el0->SetLineColor(kPink);
h_el0->Draw("hist");
h_el1->SetLineColor(kPink+2);
h_el1->Draw("hist same");

   b->SaveAs("meanRate.pdf");
*/

    TGraph* g = new TGraph(x0_fid.size(), x0_fid.data(), y0_fid.data());
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kRed);
    g->SetTitle("fiducial/online_trigger vrtx0");
//    g->Draw("A*");

    TGraph* g1 = new TGraph(x1_fid.size(), x1_fid.data(), y1_fid.data());
    g1->SetMarkerStyle(20);
    g1->SetMarkerColor(kOrange);
    g1->SetTitle("fiducial/online_trigger vrtx1");
//    g1->Draw("A*");

    TGraph* g2 = new TGraph(x0_vrtx.size(), x0_vrtx.data(), y0_vrtx.data());
    g2->SetMarkerStyle(20);
    g2->SetMarkerColor(kBlue);
    g2->SetTitle("any vrtx0/fiducial");
//    g2->Draw("A*");

    TGraph* g12 = new TGraph(x1_vrtx.size(), x1_vrtx.data(), y1_vrtx.data());
    g12->SetMarkerStyle(20);
    g12->SetMarkerColor(kCyan);
    g12->SetTitle("any vrtx1/fiducial");
//    g12->Draw("A*");

    TGraph* g3 = new TGraph(x0_preEl.size(), x0_preEl.data(), y0_preEl.data());
    g3->SetMarkerStyle(20);
    g3->SetMarkerColor(kGreen);
//    g3->SetTitle("pre_elastic/any_vrtx0");
    g3->SetTitle("pre_elastic/fiducial vrtx0");
//    g3->Draw("A*");

    TGraph* g13 = new TGraph(x1_preEl.size(), x1_preEl.data(), y1_preEl.data());
    g13->SetMarkerStyle(20);
    g13->SetMarkerColor(kGreen+2);
//    g13->SetTitle("pre_elastic/any_vrtx1");
    g13->SetTitle("pre_elastic/fiducial vrtx1");
//    g13->Draw("A*");

    TGraph* g4 = new TGraph(x0_el.size(), x0_el.data(), y0_el.data());
    g4->SetMarkerStyle(20);
    g4->SetMarkerColor(kPink);
//    g4->SetTitle("elastic/pre_elastic vrtx0");
    g4->SetTitle("elastic/fiducial vrtx0");
//    g4->Draw("A*");

    TGraph* g14 = new TGraph(x1_el.size(), x1_el.data(), y1_el.data());
    g14->SetMarkerStyle(20);
    g14->SetMarkerColor(kPink+2);
//    g14->SetTitle("elastic/pre_elastic vrtx1");
    g14->SetTitle("elastic/fiducial vrtx1");
//    g14->Draw("A*");

    TCanvas* c = new TCanvas("c", "Graph", 2200, 1000);
      c->Divide(2,1);
        c->cd(1);
TMultiGraph *mgx = new TMultiGraph();
mgx->Add(g,"A*");
mgx->Add(g1,"A*");
mgx->Add(g2,"A*");
mgx->Add(g12,"A*");
mgx->Add(g3,"A*");
mgx->Add(g13,"A*");
mgx->Add(g4,"A*");
mgx->Add(g14,"A*");
mgx->SetMaximum(0.5);
mgx->SetMinimum(0.);
mgx->Draw("A* ");
mgx->SetTitle("Rates vs time");
    mgx->GetXaxis()->SetLabelOffset(0.03);
    mgx->GetXaxis()->SetTimeDisplay(1);
    mgx->GetXaxis()->SetTimeFormat("#splitline{%H:%M}{%d/%m}");
    mgx->GetXaxis()->SetTimeOffset(0,"CEST"); 
    mgx->GetXaxis()->LabelsOption("v");
gPad->BuildLegend(0.55,0.6,0.85,0.8);
//gPad->SetLogy();


	c->cd(2);

TMultiGraph *mgx2 = new TMultiGraph();
mgx2->Add(g3,"A*");
mgx2->Add(g13,"A*");
mgx2->Add(g4,"A*");
mgx2->Add(g14,"A*");
mgx2->SetMaximum(0.01);
mgx2->SetMinimum(0.);
mgx2->Draw("A* ");
mgx2->SetTitle("Rates vs time (zoom)");
    mgx2->GetXaxis()->SetLabelOffset(0.03);
    mgx2->GetXaxis()->SetTimeDisplay(1);
    mgx2->GetXaxis()->SetTimeFormat("#splitline{%H:%M}{%d/%m}");
    mgx2->GetXaxis()->SetTimeOffset(0,"CEST"); 
    mgx2->GetXaxis()->LabelsOption("v");
gPad->BuildLegend(0.55,0.6,0.85,0.8);

        c->SaveAs(Form("wip11x_ratevstime_fiducial_run%s.pdf",argv[2]));



return 0;
}
