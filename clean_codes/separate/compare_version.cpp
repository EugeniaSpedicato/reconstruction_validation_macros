#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TStyle.h>


void compare_version(){


/*TFile * f_el_pre = TFile::Open("/home/espedica/macros_fairmu/clean_codes/separate/validation/results/sigma_el_pre_2hit.root");
TFile * f_el_post = TFile::Open("/home/espedica/macros_fairmu/clean_codes/separate/validation/results/sigma_el_post_2hit.root");
TFile * f_mu_pre = TFile::Open("/home/espedica/macros_fairmu/clean_codes/separate/validation/results/sigma_mu_pre_2hit.root");
TFile * f_mu_post = TFile::Open("/home/espedica/macros_fairmu/clean_codes/separate/validation/results/sigma_mu_post_2hit.root");


//TFile *f_el_1 = TFile::Open("/home/espedica/macros_fairmu/clean_codes/separate/validation/results/el_1.root");
//TFile *f_mu_1 = TFile::Open("/home/espedica/macros_fairmu/clean_codes/separate/validation/results/mu_1.root");
TFile *f_el_2 = TFile::Open("/home/espedica/macros_fairmu/clean_codes/separate/validation/results/el_2.root");
TFile *f_mu_2 = TFile::Open("/home/espedica/macros_fairmu/clean_codes/separate/validation/results/mu_2.root");
*/

TFile * f_el_pre = TFile::Open("/home/espedica/macros_fairmu/snakemake/plots/results/sigma_el_pre_2hit_first2.root");
TFile * f_el_post = TFile::Open("/home/espedica/macros_fairmu/snakemake/plots/results/sigma_el_post_2hit_first2.root");
TFile * f_mu_pre = TFile::Open("/home/espedica/macros_fairmu/snakemake/plots/results/sigma_mu_pre_2hit_first2.root");
TFile * f_mu_post = TFile::Open("/home/espedica/macros_fairmu/snakemake/plots/results/sigma_mu_post_2hit_first2.root");


//TFile *f_el_1 = TFile::Open("/home/espedica/macros_fairmu/snakemake/plots/results/el_1.root");
//TFile *f_mu_1 = TFile::Open("/home/espedica/macros_fairmu/snakemake/plots/results/mu_1.root");
TFile *f_el_2 = TFile::Open("/home/espedica/macros_fairmu/snakemake/plots/results/el_2_first2.root");
TFile *f_mu_2 = TFile::Open("/home/espedica/macros_fairmu/snakemake/plots/results/mu_2_first2.root");

/*
TFile * f_el_pre_MCcorr = TFile::Open("/home/espedica/macros_fairmu/clean_codes/separate/validation/results/sigma_el_pre_2hit_MCcorr.root");
TFile * f_el_post_MCcorr = TFile::Open("/home/espedica/macros_fairmu/clean_codes/separate/validation/results/sigma_el_post_2hit_MCcorr.root");
TFile * f_mu_pre_MCcorr = TFile::Open("/home/espedica/macros_fairmu/clean_codes/separate/validation/results/sigma_mu_pre_2hit_MCcorr.root");
TFile * f_mu_post_MCcorr = TFile::Open("/home/espedica/macros_fairmu/clean_codes/separate/validation/results/sigma_mu_post_2hit_MCcorr.root");


TFile *f_el_1_MCcorr = TFile::Open("/home/espedica/fair_install/instFairRoot/share/MUonE/macros/proposal/el_1_MCcorr.root");
TFile *f_mu_1_MCcorr = TFile::Open("/home/espedica/fair_install/instFairRoot/share/MUonE/macros/proposal/mu_1_MCcorr.root");
TFile *f_el_2_MCcorr = TFile::Open("/home/espedica/fair_install/instFairRoot/share/MUonE/macros/proposal/el_2_MCcorr.root");
TFile *f_mu_2_MCcorr = TFile::Open("/home/espedica/fair_install/instFairRoot/share/MUonE/macros/proposal/mu_2_MCcorr.root");*/

// /home/espedica/macros_fairmu/snakemake/plots/results/

TFile * f_el_pre_MCcorr =TFile::Open("/home/espedica/macros_fairmu/snakemake/plots/results/sigma_el_pre_2hit_first2vertexSplitting.root");
TFile * f_el_post_MCcorr = TFile::Open("/home/espedica/macros_fairmu/snakemake/plots/results/sigma_el_post_2hit_first2vertexSplitting.root");
TFile * f_mu_pre_MCcorr = TFile::Open("/home/espedica/macros_fairmu/snakemake/plots/results/sigma_mu_pre_2hit_first2vertexSplitting.root");
TFile * f_mu_post_MCcorr = TFile::Open("/home/espedica/macros_fairmu/snakemake/plots/results/sigma_mu_post_2hit_first2vertexSplitting.root");


TFile *f_el_2_MCcorr = TFile::Open("/home/espedica/macros_fairmu/snakemake/plots/results/el_2_first2vertexSplitting.root");
TFile *f_mu_2_MCcorr = TFile::Open("/home/espedica/macros_fairmu/snakemake/plots/results/mu_2_first2vertexSplitting.root");


TH1::SetDefaultSumw2(kTRUE);

TGraph *el_pre_MCcorr=(TGraph*)f_el_pre_MCcorr->Get("el_pre");
TGraph *mu_pre_MCcorr=(TGraph*)f_mu_pre_MCcorr->Get("mu_pre");
TGraph *el_post_MCcorr=(TGraph*)f_el_post_MCcorr->Get("el_post");
TGraph *mu_post_MCcorr=(TGraph*)f_mu_post_MCcorr->Get("mu_post");


TGraph *el_pre=(TGraph*)f_el_pre->Get("el_pre");
TGraph *mu_pre=(TGraph*)f_mu_pre->Get("mu_pre");
TGraph *el_post=(TGraph*)f_el_post->Get("el_post");
TGraph *mu_post=(TGraph*)f_mu_post->Get("mu_post");



//TGraph *el_1_MCcorr=(TGraph*)f_el_1_MCcorr->Get("el_1");
//TGraph *mu_1_MCcorr=(TGraph*)f_mu_1_MCcorr->Get("mu_1");
TGraph *el_2_MCcorr=(TGraph*)f_el_2_MCcorr->Get("el_2");
TGraph *mu_2_MCcorr=(TGraph*)f_mu_2_MCcorr->Get("mu_2");


//TGraph *el_1=(TGraph*)f_el_1->Get("el_1");
//TGraph *mu_1=(TGraph*)f_mu_1->Get("mu_1");
TGraph *el_2=(TGraph*)f_el_2->Get("el_2");
TGraph *mu_2=(TGraph*)f_mu_2->Get("mu_2");

 /*auto legend_mu1 = new TLegend(0.75,0.15,0.9,0.3);
legend_mu1->AddEntry(mu_1,"WP_14 default","LEP");
legend_mu1->AddEntry(mu_1_MCcorr,"WP_14 vertexSplitting","LEP");//Reco-fixes GA","LEP");
*/

 auto legend_mu2 = new TLegend(0.75,0.15,0.9,0.3);
legend_mu2->AddEntry(mu_2,"WP_14 default","LEP");
legend_mu2->AddEntry(mu_2_MCcorr,"WP_14 vertexSplitting","LEP");//Reco-fixes GA","LEP");

/*
 auto legend_e1 = new TLegend(0.75,0.75,0.9,0.9);
legend_e1->AddEntry(el_1,"WP_14 default","LEP");
legend_e1->AddEntry(el_1_MCcorr,"WP_14 vertexSplitting","LEP");//Reco-fixes GA","LEP");
*/
 auto legend_e2 = new TLegend(0.75,0.75,0.9,0.9);
legend_e2->AddEntry(el_2,"WP_14 default","LEP");
legend_e2->AddEntry(el_2_MCcorr,"WP_14 vertexSplitting","LEP");//Reco-fixes GA","LEP");


TCanvas d("d","d",2100,1400);
d.Divide(1,2);
/*d.cd(1);
TMultiGraph *mg1 = new TMultiGraph();
el_1->SetMarkerColor(kRed);
el_1_MCcorr->SetMarkerColor(kBlue);
mg1->Add(el_1,"A*");
mg1->Add(el_1_MCcorr,"A*");
mg1->Draw("A*");
legend_e1->Draw();
mg1->SetTitle("Sigma difference el 1 hit shared");
mg1->GetYaxis()->SetTitle("(#sigma_post - #sigma_pre)/#sigma_pre");
mg1->GetXaxis()->SetTitle("#theta_#el (mrad)");
*/
d.cd(1);
TMultiGraph *mg2 = new TMultiGraph();
el_2->SetMarkerColor(kRed);
el_2_MCcorr->SetMarkerColor(kBlue);
mg2->Add(el_2,"A*");
mg2->Add(el_2_MCcorr,"A*");
mg2->Draw("A*");
legend_e2->Draw();
mg2->SetTitle("Sigma difference el 2 hit shared");
mg2->GetYaxis()->SetTitle("(#sigma_post - #sigma_pre)/#sigma_pre");
mg2->GetXaxis()->SetTitle("#theta_#el (mrad)");

/*d.cd(3);
TMultiGraph *mg3 = new TMultiGraph();
mu_1->SetMarkerColor(kRed);
mu_1_MCcorr->SetMarkerColor(kBlue);
mg3->Add(mu_1,"A*");
mg3->Add(mu_1_MCcorr,"A*");
mg3->Draw("A*");
legend_mu1->Draw();

mg3->SetTitle("Sigma difference #mu 1 hit shared");
mg3->GetYaxis()->SetTitle("(#sigma_post - #sigma_pre)/#sigma_pre");
mg3->GetXaxis()->SetTitle("#theta_#mu (mrad)");
*/
d.cd(2);
TMultiGraph *mg4 = new TMultiGraph();
mu_2->SetMarkerColor(kRed);
mu_2_MCcorr->SetMarkerColor(kBlue);
mg4->Add(mu_2,"A*");
mg4->Add(mu_2_MCcorr,"A*");
mg4->Draw("A*");
legend_mu2->Draw();
mg4->SetTitle("Sigma difference #mu 2 hit shared");
mg4->GetYaxis()->SetTitle("(#sigma_post - #sigma_pre)/#sigma_pre");
mg4->GetXaxis()->SetTitle("#theta_#mu (mrad)");

d.SaveAs("/home/espedica/macros_fairmu/clean_codes/separate/validation/results/diff_first2vertexSplitting.pdf");



 auto legend_mu1pre = new TLegend(0.75,0.15,0.9,0.3);
legend_mu1pre->AddEntry(mu_pre,"pre WP_14","LEP");
legend_mu1pre->AddEntry(mu_pre_MCcorr,"pre WP_14 vertexSplitting","LEP");//Reco-fixes GA","LEP");

 auto legend_mu2post = new TLegend(0.75,0.15,0.9,0.3);
legend_mu2post->AddEntry(mu_post,"post WP_14","LEP");
legend_mu2post->AddEntry(mu_post_MCcorr,"post WP_14 vertexSplitting","LEP");//Reco-fixes GA","LEP");


 auto legend_e1pre = new TLegend(0.75,0.15,0.9,0.3);
legend_e1pre->AddEntry(el_pre,"pre WP_14","LEP");
legend_e1pre->AddEntry(el_pre_MCcorr,"pre WP_14 vertexSplitting","LEP");//Reco-fixes GA","LEP");

 auto legend_e2post = new TLegend(0.75,0.15,0.9,0.3);
legend_e2post->AddEntry(el_post,"post WP_14","LEP");
legend_e2post->AddEntry(el_post_MCcorr,"post WP_14 vertexSplitting","LEP");//Reco-fixes GA","LEP");


auto el_pre_diff = new TGraphErrors();
auto el_post_diff = new TGraphErrors();
auto mu_pre_diff = new TGraphErrors();
auto mu_post_diff = new TGraphErrors();

   const Int_t NBINS_mu = 14;
   Double_t edges_mu[NBINS_mu + 1] = {0.00005,0.00015,0.00025,0.00035,0.00045,0.00055,0.00065,0.00075,0.00085,0.00095,0.0015,0.0025,0.0035,0.0045};

   const Int_t NBINS = 14;
   Double_t edges_el[NBINS + 1] = {0.0005,0.0015,0.0025,0.0035,0.0045,0.0055,0.0065,0.0075,0.0085,0.0095, 0.0125, 0.0175, 0.0225, 0.0285};


for(int m=0; m< el_pre->GetN(); m++){
el_pre_diff->SetPoint(m,edges_el[m]*1000,el_pre->GetPointY(m)-el_pre_MCcorr->GetPointY(m));
el_post_diff->SetPoint(m,edges_el[m]*1000,el_post->GetPointY(m)-el_post_MCcorr->GetPointY(m));
}

for(int m=0; m< mu_pre->GetN(); m++){
mu_pre_diff->SetPoint(m,edges_mu[m]*1000,mu_pre->GetPointY(m)-mu_pre_MCcorr->GetPointY(m));
mu_post_diff->SetPoint(m,edges_mu[m]*1000,mu_post->GetPointY(m)-mu_post_MCcorr->GetPointY(m));
}



TCanvas d2("d2","d2",2100,1400);
d2.Divide(2,2);
d2.cd(1);
TMultiGraph *mg1pre = new TMultiGraph();
el_pre->SetMarkerColor(kRed);
el_pre_MCcorr->SetMarkerColor(kBlue);
mg1pre->Add(el_pre,"A*");
mg1pre->Add(el_pre_MCcorr,"A*");
mg1pre->Draw("A*");
legend_e1pre->Draw();
mg1pre->SetTitle("Pre sigma el 2 hit shared");
mg1pre->GetYaxis()->SetTitle("#sigma_pre");
mg1pre->GetXaxis()->SetTitle("#theta_#el (mrad)");
mg1pre->SetMinimum(0.0001);

/*el_pre_diff->SetTitle("Pre-Pre_MCcorr sigma el 2 hit shared");
el_pre_diff->SetMinimum(-0.008);
el_pre_diff->SetMaximum(0.08);
el_pre_diff->Draw("A*");
*///gPad->SetLogy();


d2.cd(2);
TMultiGraph *mg2post = new TMultiGraph();
el_post->SetMarkerColor(kRed);
el_post_MCcorr->SetMarkerColor(kBlue);
mg2post->Add(el_post,"A*");
mg2post->Add(el_post_MCcorr,"A*");
mg2post->Draw("A*");
legend_e2post->Draw();
mg2post->SetTitle("Post sigma el 2 hit shared");
mg2post->GetYaxis()->SetTitle("#sigma_post");
mg2post->GetXaxis()->SetTitle("#theta_#el (mrad)");
mg2post->SetMinimum(0.0001);

/*el_post_diff->SetTitle("Post-Post_MCcorr sigma el 2 hit shared");
el_post_diff->SetMinimum(-0.008);
el_post_diff->SetMaximum(0.08);
el_post_diff->Draw("A*");
*///gPad->SetLogy();


d2.cd(3);
TMultiGraph *mg3pre = new TMultiGraph();
mu_pre->SetMarkerColor(kRed);
mu_pre_MCcorr->SetMarkerColor(kBlue);
mg3pre->Add(mu_pre,"A*");
mg3pre->Add(mu_pre_MCcorr,"A*");
mg3pre->Draw("A*");
legend_mu1pre->Draw();
mg3pre->SetMinimum(0.0001);
mg3pre->SetTitle("Pre sigma #mu 2 hit shared");
mg3pre->GetYaxis()->SetTitle("#sigma_pre");
mg3pre->GetXaxis()->SetTitle("#theta_#mu (mrad)");
/*mu_pre_diff->SetTitle("Pre-Pre_MCcorr sigma mu 2 hit shared");
mu_pre_diff->SetMinimum(-0.001);
mu_pre_diff->SetMaximum(0.02);
mu_pre_diff->Draw("A*");
//gPad->SetLogy();
*/

d2.cd(4);
TMultiGraph *mg4post = new TMultiGraph();
mu_post->SetMarkerColor(kRed);
mu_post_MCcorr->SetMarkerColor(kBlue);
mg4post->Add(mu_post,"A*");
mg4post->Add(mu_post_MCcorr,"A*");
mg4post->Draw("A*");
legend_mu2post->Draw();
mg4post->SetTitle("Post sigma #mu 2 hit shared");
mg4post->GetYaxis()->SetTitle("#sigma_post");
mg4post->GetXaxis()->SetTitle("#theta_#mu (mrad)");
mg4post->SetMinimum(0.0001);
/*mu_post_diff->SetTitle("Post-Post_MCcorr sigma mu 2 hit shared");
mu_post_diff->SetMinimum(-0.001);
mu_post_diff->SetMaximum(0.02);
mu_post_diff->Draw("A*");
*///gPad->SetLogy();

d2.SaveAs("/home/espedica/macros_fairmu/clean_codes/separate/validation/results/sigma_overlap_first2vertexSplitting.pdf");

/*
for(int i=0; i< el_2->GetN(); i++){
cout << "diff el_2 " << el_2->GetPointY(i) << endl;;
cout << "el pre " << el_pre->GetPointY(i) << endl;
cout << "el post " << el_post->GetPointY(i) << endl;
cout << "<><><><><><>"<<endl;
cout << "diff el_2_MCcorr " << el_2_MCcorr->GetPointY(i) << endl;;
cout << "el_MCcorr pre " << el_pre_MCcorr->GetPointY(i) << endl;
cout << "el_MCcorr post " << el_post_MCcorr->GetPointY(i) << endl;
cout << "<.>.<.>.<.>.<.>.<.>.<.>"<<endl;
cout <<endl;}
cout<<"_________________________"<<endl;
for(int i=0; i< mu_2->GetN(); i++){
cout << "diff mu_2 " << mu_2->GetPointY(i) << endl;;
cout << "mu pre " << mu_pre->GetPointY(i) << endl;
cout << "mu post " << mu_post->GetPointY(i) << endl;
cout << "<><><><><><>"<<endl;
cout << "diff mu_2_MCcorr " << mu_2_MCcorr->GetPointY(i) << endl;;
cout << "mu_MCcorr pre " << mu_pre_MCcorr->GetPointY(i) << endl;
cout << "mu_MCcorr post " << mu_post_MCcorr->GetPointY(i) << endl;
cout << "<.>.<.>.<.>.<.>.<.>.<.>"<<endl;
cout <<endl;
}*/


}


