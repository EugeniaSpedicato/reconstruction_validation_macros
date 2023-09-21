#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TGraph.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TSystemDirectory.h"
#include <TStyle.h>

using namespace std;

void RealDataAnalyzer(int numb){


        TFile *inputfile = new TFile("/mnt/raid10/DATA/espedica/fairmu/dataReconstruction_run3234_3235_1.root");
//TRMesmer_3cm.root");//dataReconstruction_try_sigma.root");//nrrowBP_20GeV_1shared.root");//trPoints,$
        TTree* cbmsim = (TTree*) inputfile->Get("cbmsim");


        MUonERecoOutput *ReconstructionOutput = 0;
        cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);


double thmu=0; double the=0;
double signal=0; double reco=0; double reco1=0; double more_reco=0; double reco0=0;

auto gx = new TGraph();
auto gy = new TGraph();
auto gx_v1 = new TGraph();
auto gy_v1 = new TGraph();

std::vector<double> qx_in;
std::vector<double> mx_in;

std::vector<double> qy_in;
std::vector<double> my_in;

std::vector<double> qx;
std::vector<double> mx;

std::vector<double> qy;
std::vector<double> my;

std::vector<double> qx_in_v;
std::vector<double> mx_in_v;

std::vector<double> qy_in_v;
std::vector<double> my_in_v;

std::vector<double> qx_v;
std::vector<double> mx_v;

std::vector<double> qy_v;
std::vector<double> my_v;

          auto trackXatZ = [](double q, double m,double z) {

                return q + (z ) * m;
            };

            auto trackYatZ = [](double q, double m,double z) {

                return q + (z ) * m;
            };


// U and V from strip number to cm in local reference frame
                auto rotationU=[&](double seedClusterCenterStrip){ double posStripU = seedClusterCenterStrip;///2. -1;
                                      return ( (posStripU - 1016/2.)*0.009 +0.009/2.);};
                auto rotationV=[&](double seedClusterCenterStrip){ double posStripV = seedClusterCenterStrip;///2. -1;
                                      return ( + (posStripV - 1016/2.)*0.009 +0.009/2.);};
// U and V rotation in global XY frame
                auto newX=[](double angle, double U, double V){return cos(angle)*U + sin(angle)*V;};
                auto newY=[](double angle, double U, double V){return -sin(angle)*U + cos(angle)*V;};

// X and Y rotation in local UV frame
                auto newU=[](double angle, auto X, auto Y){return cos(angle)*X - sin(angle)*Y;};
                auto newV=[](double angle, auto X, auto Y){return sin(angle)*X + cos(angle)*Y;};

double x=0.;double y=0.;double z=0.;

for(Long64_t i = numb; i < numb+1; i++) {
//for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {

		cbmsim->GetEntry(i);
		if(i%1000 == 0) cout<<"Entry "<<i<<endl;

cout<<"Entry "<< cbmsim->GetEntry(i)<<endl;

vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();

cout << "tracks.size() " << tracks.size() << endl;

 if(tracks.size()>2){

    int sec0=0; int sec1=0;

    for(int j=0; j<tracks.size();j++)
    {
     	if(tracks.at(j).sector()==0) sec0++;
        if(tracks.at(j).sector()==1) sec1++;
        }

cout << "Track size 3 and " << sec0 << ", " << sec1 << endl;



MUonERecoOutputVertex vrtx = ReconstructionOutput->bestVertex();

x=vrtx.x();
y=vrtx.y();
z=vrtx.z();

 MUonERecoOutputTrack mu_in = vrtx.incomingMuon();
 MUonERecoOutputTrack mu_out = vrtx.outgoingMuon();
 MUonERecoOutputTrack e_out = vrtx.outgoingElectron();

cout << "vrtx chi2 " << vrtx.chi2perDegreeOfFreedom() << endl;

        gx_v1->SetPoint(gx_v1->GetN(),5,vrtx.x());
        gy_v1->SetPoint(gy_v1->GetN(),5,vrtx.y());

qx_in_v.push_back(mu_in.x0());
mx_in_v.push_back(mu_in.xSlope());
qy_in_v.push_back(mu_in.y0());
my_in_v.push_back(mu_in.ySlope());


qx_v.push_back(mu_out.x0());
mx_v.push_back(mu_out.xSlope());
qy_v.push_back(mu_out.y0());
my_v.push_back(mu_out.ySlope());

qx_v.push_back(e_out.x0());
mx_v.push_back(e_out.xSlope());
qy_v.push_back(e_out.y0());
my_v.push_back(e_out.ySlope());



std::array<std::vector<double>,6> position;//={{{0.,0.},{0.,0.},{0.,0.},{0.,0.},{0.,0.},{0.,0.}}};
std::array<std::vector<double>,6> position_in;//={{{0.,0.},{0.,0.},{0.,0.},{0.,0.},{0.,0.},{0.,0.}}};

double px0;
double py0;

for(int j=0; j<tracks.size();j++)
{

if(tracks.at(j).sector()==0 and sec0==1)
{
std::vector<MUonERecoOutputHit> hits_=tracks.at(j).hits();
for(int h=0;h<hits_.size();h++){
position_in.at(hits_.at(h).moduleID()).push_back(hits_.at(h).positionPerpendicular()); 
                }

qx_in.push_back(tracks.at(j).x0());
mx_in.push_back(tracks.at(j).xSlope());
qy_in.push_back(tracks.at(j).y0());
my_in.push_back(tracks.at(j).ySlope());

        cout << "post vrtx: trackY in Z " << trackYatZ(mu_in.y0(),mu_in.ySlope(),0.) << " VS " << y << endl;
        cout << "pre vrtx: trackY in Z " << trackYatZ(tracks.at(j).y0(),tracks.at(j).ySlope(),z) << " VS " << y << endl;

}
if(tracks.at(j).sector()==1 and sec1>0) //and tracks.at(0).processIDofLinkedTrack()==45 and tracks.at(0).linkedTrackID()!=tracks.at(1).linkedTrackID()){
{
std::vector<MUonERecoOutputHit> hits_=tracks.at(j).hits();

qx.push_back(tracks.at(j).x0());
mx.push_back(tracks.at(j).xSlope());
qy.push_back(tracks.at(j).y0());
my.push_back(tracks.at(j).ySlope());

for(int h=0;h<hits_.size();h++){
position.at(hits_.at(h).moduleID()).push_back(hits_.at(h).positionPerpendicular());
cout << hits_.at(h).moduleID() << ") hits_.at(h).positionPerpendicular() " << hits_.at(h).positionPerpendicular() << endl;
		}
	}
}

if(position_in.at(0).size()!=0){
for(int h=0;h<position_in.at(0).size(); h++){ gx->SetPoint(gx->GetN(),0,position_in.at(0).at(h));} 
}
if(position_in.at(1).size()!=0){
for(int h=0;h<position_in.at(1).size(); h++){ gy->SetPoint(gy->GetN(),0,position_in.at(1).at(h));} 
}


if(position_in.at(2).size()==position_in.at(3).size()!=0){
for(int h=0;h<position_in.at(2).size(); h++){   gx->SetPoint(gx->GetN(),2,newX(45,-position_in.at(2).at(h),position_in.at(3).at(h)));
						gy->SetPoint(gy->GetN(),2,newY(45,-position_in.at(2).at(h),position_in.at(3).at(h)));} 
}

if(position_in.at(4).size()!=0){
for(int h=0;h<position_in.at(4).size(); h++){ gx->SetPoint(gx->GetN(),4,position_in.at(4).at(h));} 
}
if(position_in.at(5).size()!=0){
for(int h=0;h<position_in.at(5).size(); h++){ gy->SetPoint(gy->GetN(),4,position_in.at(5).at(h));} 
}


/*if(position.at(2).size()!=0 and position.at(2).size()==position.at(3).size())
{
for(int h=0;h<position.at(2).size(); h++){
	double X = newX(45,-position.at(2).at(h),position.at(3).at(h));
 	double Y = newY(45,-position.at(2).at(h),position.at(3).at(h));
gx->SetPoint(gx->GetN(),8,X);
gy->SetPoint(gy->GetN(),8,Y);
	}
}
}
*/

if(position.at(0).size()!=0){
for(int h=0;h<position.at(0).size(); h++){ gx->SetPoint(gx->GetN(),6,position.at(0).at(h));} 
}

/*if(position.at(2).size()!=0){
for(int h=0;h<position.at(2).size(); h++){ gx->SetPoint(gx->GetN(),8,position.at(2).at(h));} 
}*/

if(position.at(4).size()!=0){
for(int h=0;h<position.at(4).size(); h++){ gx->SetPoint(gx->GetN(),10,position.at(4).at(h));} 
}

if(position.at(1).size()!=0){
for(int h=0;h<position.at(1).size(); h++){ gy->SetPoint(gy->GetN(),6,position.at(1).at(h));} 
}

/*if(position.at(3).size()!=0){
for(int h=0;h<position.at(3).size(); h++){ gy->SetPoint(gy->GetN(),8,position.at(3).at(h));} 
}*/

if(position.at(5).size()!=0){
for(int h=0;h<position.at(5).size(); h++){ gy->SetPoint(gy->GetN(),10,position.at(5).at(h));} 
}

if(position.at(2).size()==position.at(3).size()!=0){
for(int h=0;h<position.at(2).size(); h++){   gx->SetPoint(gx->GetN(),8,newX(45,-position.at(2).at(h),position.at(3).at(h)));
                                                gy->SetPoint(gy->GetN(),8,newY(45,-position.at(2).at(h),position.at(3).at(h)));} 
}

}
cout << "---------------------"<<endl;
} //end of general for


TCanvas a2("a2","a2",1400,1400);
a2.Divide(2,2);
a2.cd(1);

TLine *t = new TLine(5,-5,5,5);
t->SetLineWidth(4);
t->SetLineColor(47);
TLine *l0 = new TLine(0,-5,0,5);
TLine *l1 = new TLine(2,-5,2,5);
TLine *l2 = new TLine(4,-5,4,5);
TLine *l3 = new TLine(6,-5,6,5);
TLine *l4 = new TLine(8,-5,8,5);
TLine *l5 = new TLine(10,-5,10,5);

gx->SetMinimum(-6);
gx->SetMaximum(6);
gx->SetMarkerColor(kRed);
gx->SetTitle("sig track stubs");


TMultiGraph *mgx = new TMultiGraph();
mgx->SetMinimum(-6);
mgx->SetMaximum(6);
mgx->Add(gx,"A*");
mgx->Draw("A* ");
mgx->SetTitle("X projection");


if(qx_in.size()!=0 and mx_in.size()!=0)
{
        for(int c=0; c<qx_in.size(); c++)
               {if(c==0)
		{TLine* lx = new TLine(0, trackXatZ(qx_in.at(c),mx_in.at(c),810.), 5, trackXatZ(qx_in.at(c),mx_in.at(c), z));
                 lx->SetLineColor(kBlue+c);
                 lx->Draw("same");
		}
		else if(c==1){
		 TLine* lx = new TLine(5, trackXatZ(qx_in.at(c),mx_in.at(c), 912.7), 10, trackXatZ(qx_in.at(c),mx_in.at(c), 912.7+89.9218));
                 lx->SetLineColor(kBlue+c);
                 lx->Draw("same");
			} 
       }
}

if(qx.size()!=0 and mx.size()!=0)
{ cout << qx.size() << " qx"<<endl;
        for(int c=0; c<qx.size(); c++)
               {TLine* lx = new TLine(5, trackXatZ(qx.at(c),mx.at(c), 912.7), 10, trackXatZ(qx.at(c),mx.at(c), 912.7+89.9218));
                lx->SetLineColor(kRed+c);
                //lx->SetTitle("sig track");
                lx->Draw("same");
        }
}

gPad->BuildLegend(0.25,0.15,0.25,0.15);
l0->Draw("same");
l1->Draw("same");
l2->Draw("same");
t->Draw("same");
l3->Draw("same");
l4->Draw("same");
l5->Draw("same");

a2.cd(2);

gx_v1->SetMinimum(-6);
gx_v1->SetMaximum(6);
gx_v1->SetMarkerColor(kOrange);
gx_v1->SetTitle("vertex");
gx_v1->SetMarkerStyle(29);
gx_v1->SetMarkerSize(5);

gx->SetMinimum(-6);
gx->SetMaximum(6);
gx->SetMarkerColor(kRed);
gx->SetTitle("sig track stubs");

TMultiGraph *mgx_v = new TMultiGraph();
mgx_v->SetMinimum(-6);
mgx_v->SetMaximum(6);
mgx_v->Add(gx_v1,"A*");
mgx_v->Add(gx,"A*");
mgx_v->Draw("A* ");
mgx_v->SetTitle("X projection vrtx fit");



if(qx_in_v.size()!=0 and mx_in_v.size()!=0)
{
        for(int c=0; c<qx_in_v.size(); c++)
               {if(c==0)
                {//TLine* lx = new TLine(0, trackXatZ(qx_in_v.at(c),mx_in_v.at(c), 0.), 5, trackXatZ(qx_in_v.at(c),mx_in_v.at(c),z-810.));
		 TLine* lx = new TLine(0, trackXatZ(qx_in_v.at(c),mx_in_v.at(c), -93.7693), 5, trackXatZ(qx_in_v.at(c),mx_in_v.at(c),0.));
                 lx->SetLineColor(kBlue+c);
                 lx->Draw("same");
        cout << "trackX in Z " << trackXatZ(qx_in_v.at(c),mx_in_v.at(c),0.) << " VS " << x << endl;

                }
                else if(c==1){
                 TLine* lx = new TLine(5, trackXatZ(qx_in_v.at(c),mx_in_v.at(c), 0.), 10, trackXatZ(qx_in_v.at(c),mx_in_v.at(c),93.7693));
                 lx->SetLineColor(kBlue+c);
                 lx->Draw("same");
                        } 
       }
}

if(qx_v.size()!=0 and mx_v.size()!=0)
{ cout << qx_v.size() << " qx_v"<<endl;
        for(int c=0; c<qx_v.size(); c++)
               {TLine* lx = new TLine(5, trackXatZ(qx_v.at(c),mx_v.at(c), 0.), 10, trackXatZ(qx_v.at(c),mx_v.at(c), 93.7693));
                lx->SetLineColor(kOrange+c);
                //lx->SetTitle("bkg track");
                lx->Draw("same");
        }
}

gPad->BuildLegend(0.25,0.15,0.25,0.15);

l0->Draw("same");
l1->Draw("same");
l2->Draw("same");
t->Draw("same");
l3->Draw("same");
l4->Draw("same");
l5->Draw("same");

a2.cd(3);

gy->SetMinimum(-6);
gy->SetMaximum(6);
gy->SetMarkerColor(kRed);

gy->SetTitle("sig track stubs");

TMultiGraph *mg = new TMultiGraph();
mg->SetMinimum(-6);
mg->SetMaximum(6);
mg->Add(gy,"A*");
mg->Draw("A* ");
mg->SetTitle("Y projection");

if(qy_in.size()!=0 and my_in.size()!=0)
{
        for(int c=0; c<qy_in.size(); c++)
               {if(c==0)
                {TLine* ly = new TLine(0, trackYatZ(qy_in.at(c),my_in.at(c), 811), 5, trackYatZ(qy_in.at(c),my_in.at(c), z));
                 ly->SetLineColor(kBlue+c);
                 ly->Draw("same");
	        cout << "pre vrtx: trackY in Z " << trackYatZ(qy_in.at(c),my_in.at(c),z) << " VS " << y << endl;

                }
                else if(c==1){ 
                 TLine* ly = new TLine(5, trackYatZ(qy_in.at(c),my_in.at(c), 912.7), 10, trackYatZ(qy_in.at(c),my_in.at(c), 912.7+89.9218));
                 ly->SetLineColor(kBlue+c);
                 ly->Draw("same");
                        } 
        }
}

if(qy.size()!=0 and my.size()!=0)
{ cout << qy.size() << " qy"<<endl;
        for(int c=0; c<qy.size(); c++)
               {TLine* ly = new TLine(5, trackYatZ(qy.at(c),my.at(c), 912.7), 10, trackYatZ(qy.at(c),my.at(c), 912.7+89.9218));
                ly->SetLineColor(kRed+c);
                //ly->SetTitle("sig track");
                ly->Draw("same");
	}
}

gPad->BuildLegend(0.25,0.15,0.25,0.15);
l0->Draw("same");
l1->Draw("same");
l2->Draw("same");
t->Draw("same");
l3->Draw("same");
l4->Draw("same");
l5->Draw("same");

a2.cd(4);

gy_v1->SetMinimum(-6);
gy_v1->SetMaximum(6);
gy_v1->SetMarkerColor(kOrange);
gy_v1->SetTitle("vertex");
gy_v1->SetMarkerStyle(29);
gy_v1->SetMarkerSize(5);

gy->SetMinimum(-6);
gy->SetMaximum(6);
gy->SetMarkerColor(kRed);
gy->SetTitle("sig track stubs");

TMultiGraph *mgy_v = new TMultiGraph();
mgy_v->SetMinimum(-6);
mgy_v->SetMaximum(6);
mgy_v->Add(gy_v1,"A*");
mgy_v->Add(gy,"A*");
mgy_v->Draw("A* ");
mgy_v->SetTitle("Y projection vrtx fit");



if(qy_in_v.size()!=0 and my_in_v.size()!=0)
{
        for(int c=0; c<qy_in_v.size(); c++)
               {if(c==0)
                {TLine* ly = new TLine(0, trackYatZ(qy_in_v.at(c),my_in_v.at(c), -93.7693), 5, trackYatZ(qy_in_v.at(c),my_in_v.at(c),0.));
                 ly->SetLineColor(kBlue+c);
                 ly->Draw("same");
	cout << "post vrtx: trackY in Z " << trackYatZ(qy_in_v.at(c),my_in_v.at(c),0.) << " VS " << y << endl;
                }
                else if(c==1){
                 TLine* ly = new TLine(5, trackYatZ(qy_in_v.at(c),my_in_v.at(c), 0.), 10, trackYatZ(qy_in_v.at(c),my_in_v.at(c), 93.7693));
                 ly->SetLineColor(kBlue+c);
                 ly->Draw("same");
                        }
	}
}

if(qy_v.size()!=0 and my_v.size()!=0)
{ cout << qy_v.size() << " qy_v"<<endl;
        for(int c=0; c<qy_v.size(); c++)
               {TLine* ly = new TLine(5, trackYatZ(qy_v.at(c),my_v.at(c), 0.), 10, trackYatZ(qy_v.at(c),my_v.at(c), 93.7693));
                ly->SetLineColor(kOrange+c);
                //lx->SetTitle("bkg track");
                ly->Draw("same");
        }
}

gPad->BuildLegend(0.25,0.15,0.25,0.15);

l0->Draw("same");
l1->Draw("same");
l2->Draw("same");
t->Draw("same");
l3->Draw("same");
l4->Draw("same");
l5->Draw("same");

a2.SaveAs("pdf_exercise.pdf");

}

