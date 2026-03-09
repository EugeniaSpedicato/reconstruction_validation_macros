#include <iostream>
#include <fstream>
#include <string>
#include <regex>
#include <filesystem>
#include <vector>
#include <algorithm>
#include "TGraph.h"
#include "TH1D.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TPaveStats.h"



void graph_rate_bkg(){

std::array<double,20> x{1.26,1.31,0.60,0.63,0.63,0.49,0.47,0.45,0.21,0.57,0.63,0.62,0.43,0.36,1.11,0.53,0.51,0.43,0.66,0.57};
std::array<double,20> y_0{0.40,0.40,0.53,1.39,1.80,1.86,1.86,1.83,3.11,4.11,4.06,4.17,2.56,4.90,2.61,4.04,4.21};
std::array<double,20> y_1{0.70,0.71,0.049,0.13,0.17,0.10,0.10,0.10,1.29,0.65,0.65,0.66,0.71,1.95,8.22,1.72,1.77};

 TGraph* tar0 = new TGraph(x.size(), x.data(), y_0.data());
    tar0->SetMarkerStyle(22);
    tar0->SetMarkerSize(1.5);
    tar0->SetMarkerColor(9);
    tar0->SetTitle("Target 0");
 TGraph* tar1 = new TGraph(x.size(), x.data(), y_1.data());
    tar1->SetMarkerStyle(20);
    tar1->SetMarkerSize(1.5);
    tar1->SetMarkerColor(46);
    tar1->SetTitle("Target 1");

TCanvas g("g","g",1300,700);
TMultiGraph *mgx = new TMultiGraph();
mgx->Add(tar0,"AP");
mgx->Add(tar1,"AP");
mgx->Draw("AP ");
mgx->SetTitle("Population rate (bkg/sig) VS tar0/tar1 elastic rate");
mgx->GetYaxis()->SetTitle("bkg/sig");
mgx->GetXaxis()->SetTitle("elastic rate tar0/tar1");
gPad->BuildLegend(0.35,0.6,0.55,0.8);
g.SaveAs("rateVSpopulations.pdf");





}
