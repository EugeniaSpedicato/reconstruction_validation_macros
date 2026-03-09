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

void MC_tb2025_ev_draw(int numb, int nhits , int version, int interaction){

double z1=0.;
double z2=0.;
double z3=0.;
double tar_z1=0.;
double tar_z2=0.;
double other_z1=0;
double other_z2=0;


	const char *name;

        if(nhits==0 and version==2025 and interaction==0){name="/mnt/raid10/DATA/espedica/fairmu/TB2025/reco/snakemake/range_5_10_0hit_defaultWiP0176IdealGeom.root";}
//name="/mnt/raid10/DATA/espedica/fairmu/TB2025/run32/single_muon_interaction_1/muedaq02-1754197878_0hit_WiP_17_6.root";}
//name="/mnt/raid10/DATA/espedica/fairmu/TB2025/run32/single_muon_interaction_0/muedaq02-1754238276_0hit_WiP_17_6.root";}


        TFile *inputfile = new TFile(name);
        TTree* cbmsim = (TTree*) inputfile->Get("cbmsim");

//IF MC
        TFile *f2 = new TFile("/eos/user/r/rpilato/shared_with_Eugenia/MCSamples_validation_WiP_1.1.x/idealGeometry/generated_5_10.root");
        TTree *tree2 = (TTree*)f2->Get("cbmsim");
        tree2->SetEntries(cbmsim->GetEntries());
           cbmsim->AddFriend(tree2);


        TClonesArray *TrackerStubs = 0;
	cbmsim->SetBranchAddress("TrackerStubs", &TrackerStubs);
              MUonERecoOutputFull *ReconstructionOutput = 0;
        cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);

    Double_t clusterEnergy = 0.0;
    cbmsim->SetBranchAddress("ReconstructedCalorimeterCluster.ClusterEnergy", &clusterEnergy);
    Double_t xCentroid = 0.0;
    cbmsim->SetBranchAddress("ReconstructedCalorimeterCluster.XCentroid", &xCentroid);
    Double_t yCentroid = 0.0;
    cbmsim->SetBranchAddress("ReconstructedCalorimeterCluster.YCentroid", &yCentroid);

double thmu=0; double the=0;
double signal=0; double reco=0; double reco1=0; double more_reco=0; double reco0=0;

auto gx = new TGraph();
auto gy = new TGraph();
auto gx_s = new TGraph();
auto gy_s = new TGraph();
auto gx_v1 = new TGraph();
auto gy_v1 = new TGraph();
auto guv = new TGraph();
auto guv_t = new TGraph();
auto g_mod_u = new TGraph();
auto g_mod_u_t = new TGraph();
auto g_mod_v = new TGraph();
auto g_mod_v_t = new TGraph();

std::vector<double> qx_in;
std::vector<double> mx_in;

std::vector<double> qy_in;
std::vector<double> my_in;

std::vector<double> qx;
std::vector<double> mx;

std::vector<double> qy;
std::vector<double> my;

std::vector<double> qx_other;
std::vector<double> mx_other;

std::vector<double> qy_other;
std::vector<double> my_other;

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
//std::array<double,6>={912.7+18.0218,912.7+21.8693,912.7+55.3635,912.7+56.6205,912.7+89.9218,912.7+92.3722};

int sec0=0;
int sec1=0;
int sec2=0;

for(Long64_t i = numb; i < numb+1; i++) {
//for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
cout << i << " vs " << numb << endl;
		cbmsim->GetEntry(i);
//		if(i%1000 == 0) cout<<"Entry "<<i<<endl;

cout<<"Entry "<< cbmsim->GetEntry(i)<<endl;

           int hit_modXmu=0; int hit_modXe=0;
           int hit_modYmu=0; int hit_modYe=0;
           int stereo_mu=0; int stereo_e=0;


vector<MUonERecoOutputTrackFull> tracks = ReconstructionOutput->reconstructedTracks();

cout << "tracks.size() " << tracks.size() << endl;
cout << "clusterEnergy " << clusterEnergy << endl;
cout << "xCentroid " << xCentroid << endl;
cout << "yCentroid " << yCentroid << endl;
 if(tracks.size()>0){


    for(int j=0; j<tracks.size();j++)
    {
     	if(tracks.at(j).sector()==0) sec0++;
        if(tracks.at(j).sector()==1) sec1++;
        if(tracks.at(j).sector()==2) sec2++;
        }

cout << "Track size " << tracks.size() << " and st0: " << sec0 << ", st1: " << sec1 << ", st2: " << sec2 << endl;



MUonERecoOutputVertexFull vrtx = ReconstructionOutput->bestVertex();

x=vrtx.xKinematicFit();
y=vrtx.yKinematicFit();
z=vrtx.zKinematicFit();

 MUonERecoOutputTrackFull mu_in = vrtx.incomingMuon();
 MUonERecoOutputTrackFull mu_out = vrtx.outgoingMuon();
 MUonERecoOutputTrackFull e_out = vrtx.outgoingElectron();

double chi=vrtx.chi2perDegreeOfFreedom();
cout << "vrtx chi2 " << vrtx.chi2perDegreeOfFreedom() << endl;
cout << "vrtx Z " << vrtx.zPositionFit() << endl;
        gx_v1->SetPoint(gx_v1->GetN(),5,vrtx.xKinematicFit());
        gy_v1->SetPoint(gy_v1->GetN(),5,vrtx.yKinematicFit());

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




std::array<std::array<std::vector<double>,6>,5> position;//={{{0.,0.},{0.,0.},{0.,0.},{0.,0.},{0.,0.},{0.,0.}}};
std::array<std::array<std::vector<double>,6>,5> position_other;//={{{0.,0.},{0.,0.},{0.,0.},{0.,0.},{0.,0.},{0.,0.}}};
std::array<std::array<std::vector<double>,6>,5> position_in;//={{{0.,0.},{0.,0.},{0.,0.},{0.,0.},{0.,0.},{0.,0.}}};

std::array<std::vector<double>,6> position_s;
std::array<std::vector<double>,6> position_s1;
std::array<std::vector<double>,6> position_s2;
double px0;
double py0;

int code_e,code_mu;

int stub0=0;
int stub1=0;
int stub2=0;

vector<MUonERecoOutputHitFull> stubs_=ReconstructionOutput->reconstructedHits();
for(int s=0; s<stubs_.size(); s++){

if(stubs_.at(s).stationID()==0)stub0++;
else if(stubs_.at(s).stationID()==1)stub1++;
else if(stubs_.at(s).stationID()==2)stub2++;

if(chi!=0){
	if( (vrtx.stationIndex()==1 and stubs_.at(s).stationID()==0) or (vrtx.stationIndex()==2 and stubs_.at(s).stationID()==1))position_s.at(stubs_.at(s).moduleID()).push_back(stubs_.at(s).positionPerpendicular());
	if((vrtx.stationIndex()==1 and stubs_.at(s).stationID()==1) or (vrtx.stationIndex()==2 and stubs_.at(s).stationID()==2)){position_s1.at(stubs_.at(s).moduleID()).push_back(stubs_.at(s).positionPerpendicular());
				cout << stubs_.at(s).moduleID()<< ") stub posPerp " << stubs_.at(s).positionPerpendicular() << endl;}
	if((vrtx.stationIndex()==1 and stubs_.at(s).stationID()==2) or (vrtx.stationIndex()==2 and stubs_.at(s).stationID()==0)){position_s2.at(stubs_.at(s).moduleID()).push_back(stubs_.at(s).positionPerpendicular());
				cout << stubs_.at(s).moduleID()<< ") other stub posPerp " << stubs_.at(s).positionPerpendicular() << endl;}
 }
else if(chi==0){
	if(stubs_.at(s).stationID()==0)position_s.at(stubs_.at(s).moduleID()).push_back(stubs_.at(s).positionPerpendicular());
	else if(stubs_.at(s).stationID()==1)position_s1.at(stubs_.at(s).moduleID()).push_back(stubs_.at(s).positionPerpendicular());
	else if(stubs_.at(s).stationID()==2)position_s2.at(stubs_.at(s).moduleID()).push_back(stubs_.at(s).positionPerpendicular());
	}
}

cout << "stub0 " << stub0 << endl;
cout << "stub1 " << stub1 << endl;
cout << "stub2 " << stub2 << endl;



std::array<int,12> trigger;
std::array<int,6> other;
int nt=0;
int nt_other=0;
int nt_0=0;

cout << "vrtx.stationIndex() " << vrtx.stationIndex() << endl;
if(vrtx.stationIndex()==1){ trigger={0,1,2,3,4,5,6,7,8,9,10,11}; other={12,13,14,15,16,17}; nt=sec1; nt_other=sec2; nt_0=sec0;
			    z1=550; z2=667.3-2.7; z3=784.6; tar_z1=5;tar_z2=10; other_z1=11; other_z2=16;}
if(vrtx.stationIndex()==2){ trigger={6,7,8,9,10,11,12,13,14,15,16,17}; other={0,1,2,3,4,5,}; nt=sec2; nt_other=sec0; nt_0=sec1;
			    z1=667.3; z2=784.6-3.9; z3=550; tar_z1=11;tar_z2=16; other_z1=0; other_z2=5;}
else{trigger={0,1,2,3,4,5,6,7,8,9,10,11}; other={12,13,14,15,16,17}; nt=sec1; nt_other=sec2; nt_0=sec0;
                            z1=550; z2=667.3-2.7; z3=784.6; tar_z1=5;tar_z2=10; other_z1=11; other_z2=16;}

int n0=0;
int n=0;
int n1=0;


for(int j=0; j<tracks.size();j++)
{
if(chi==0){
	if(tracks.at(j).sector()==0 and sec0>0){

std::vector<MUonERecoOutputHitFull> hits_=tracks.at(j).hits();
for(int h=0;h<hits_.size();h++){
position_in.at(n0).at(hits_.at(h).moduleID()).push_back(hits_.at(h).positionPerpendicular()); 
                }

	qx_in.push_back(tracks.at(j).x0());
	mx_in.push_back(tracks.at(j).xSlope());
	qy_in.push_back(tracks.at(j).y0());
	my_in.push_back(tracks.at(j).ySlope());
	n0++;
        cout << "post vrtx: trackY in Z " << trackYatZ(mu_in.y0(),mu_in.ySlope(),0.) << " VS " << y << endl;
        cout << "pre vrtx: trackY in Z " << trackYatZ(tracks.at(j).y0(),tracks.at(j).ySlope(),912.7) << " VS " << y << endl;
	}
	else if(tracks.at(j).sector()==1 and sec1>0){

	std::vector<MUonERecoOutputHitFull> hits_=tracks.at(j).hits();

	qx.push_back(tracks.at(j).x0());
	mx.push_back(tracks.at(j).xSlope());
	qy.push_back(tracks.at(j).y0());
	my.push_back(tracks.at(j).ySlope());

for(int h=0;h<hits_.size();h++){
position.at(n).at(hits_.at(h).moduleID()).push_back(hits_.at(h).positionPerpendicular());
cout << hits_.at(h).moduleID() << ") stub track " << j << " posPerp " << hits_.at(h).positionPerpendicular() << endl;
                }
	n++;
	}
	else if(tracks.at(j).sector()==2 and sec2>0){

        std::vector<MUonERecoOutputHitFull> hits_=tracks.at(j).hits();

	for(int h=0;h<hits_.size();h++){
	position_other.at(n).at(hits_.at(h).moduleID()).push_back(hits_.at(h).positionPerpendicular());
                }

	qx_other.push_back(tracks.at(j).x0());
	mx_other.push_back(tracks.at(j).xSlope());
	qy_other.push_back(tracks.at(j).y0());
	my_other.push_back(tracks.at(j).ySlope());

	n1++;
	}
}
else if(vrtx.stationIndex()==1)
 {
	if(tracks.at(j).sector()==0 and sec0==1){
std::vector<MUonERecoOutputHitFull> hits_=tracks.at(j).hits();
for(int h=0;h<hits_.size();h++){
position_in.at(n0).at(hits_.at(h).moduleID()).push_back(hits_.at(h).positionPerpendicular());
}
	qx_in.push_back(tracks.at(j).x0());
	mx_in.push_back(tracks.at(j).xSlope());
	qy_in.push_back(tracks.at(j).y0());
	my_in.push_back(tracks.at(j).ySlope());
        n0++;
        cout << "post vrtx: trackY in Z " << trackYatZ(mu_in.y0(),mu_in.ySlope(),0.) << " VS " << y << endl;
        cout << "pre vrtx: trackY in Z " << trackYatZ(tracks.at(j).y0(),tracks.at(j).ySlope(),912.7) << " VS " << y << endl;
	}
	else if(tracks.at(j).sector()==1 and sec1>0){

        std::vector<MUonERecoOutputHitFull> hits_=tracks.at(j).hits();
for(int h=0;h<hits_.size();h++){
position.at(n).at(hits_.at(h).moduleID()).push_back(hits_.at(h).positionPerpendicular());
cout << hits_.at(h).moduleID() << ") stub track " << j << " posPerp " << hits_.at(h).positionPerpendicular() << endl;
                }

	cout << j << ") int ID " <<tracks.at(j).processIDofLinkedTrack()<< endl;
	cout << j << ") hits_size " <<hits_.size()<< endl;
	cout << j << ") chi_2 " <<tracks.at(j).chi2perDegreeOfFreedom()<< endl;
	qx.push_back(tracks.at(j).x0());
	mx.push_back(tracks.at(j).xSlope());
	qy.push_back(tracks.at(j).y0());
	my.push_back(tracks.at(j).ySlope());

	n++;
	}
	else if(tracks.at(j).sector()==2 and sec2>0){

        std::vector<MUonERecoOutputHitFull> hits_=tracks.at(j).hits();

        for(int h=0;h<hits_.size();h++){
        position_other.at(n).at(hits_.at(h).moduleID()).push_back(hits_.at(h).positionPerpendicular());
                }

	qx_other.push_back(tracks.at(j).x0());
	mx_other.push_back(tracks.at(j).xSlope());
	qy_other.push_back(tracks.at(j).y0());
	my_other.push_back(tracks.at(j).ySlope());

	n1++;
	}
 }
 else if(vrtx.stationIndex()==2)
 {
	if(tracks.at(j).sector()==1 and sec0==1){

std::vector<MUonERecoOutputHitFull> hits_=tracks.at(j).hits();
for(int h=0;h<hits_.size();h++){
position_in.at(n0).at(hits_.at(h).moduleID()).push_back(hits_.at(h).positionPerpendicular());
}
	qx_in.push_back(tracks.at(j).x0());
	mx_in.push_back(tracks.at(j).xSlope());
	qy_in.push_back(tracks.at(j).y0());
	my_in.push_back(tracks.at(j).ySlope());
        n0++;
        cout << "post vrtx: trackY in Z " << trackYatZ(mu_in.y0(),mu_in.ySlope(),0.) << " VS " << y << endl;
        cout << "pre vrtx: trackY in Z " << trackYatZ(tracks.at(j).y0(),tracks.at(j).ySlope(),z2) << " VS " << y << endl;
	}
	else if(tracks.at(j).sector()==2 and sec2>0){

        std::vector<MUonERecoOutputHitFull> hits_=tracks.at(j).hits();
for(int h=0;h<hits_.size();h++){
position.at(n).at(hits_.at(h).moduleID()).push_back(hits_.at(h).positionPerpendicular());
cout << hits_.at(h).moduleID() << ") stub track " << j << " posPerp " << hits_.at(h).positionPerpendicular() << endl;
                }

	cout << j << ") int ID " <<tracks.at(j).processIDofLinkedTrack()<< endl;
	cout << j << ") hits_size " <<hits_.size()<< endl;
	cout << j << ") chi_2 " <<tracks.at(j).chi2perDegreeOfFreedom()<< endl;
	qx.push_back(tracks.at(j).x0());
	mx.push_back(tracks.at(j).xSlope());
	qy.push_back(tracks.at(j).y0());
	my.push_back(tracks.at(j).ySlope());

	n++;
	}
        else if(tracks.at(j).sector()==0 and sec0>0){

        std::vector<MUonERecoOutputHitFull> hits_=tracks.at(j).hits();

        for(int h=0;h<hits_.size();h++){
        position_other.at(n).at(hits_.at(h).moduleID()).push_back(hits_.at(h).positionPerpendicular());
                }

        qx_other.push_back(tracks.at(j).x0());
        mx_other.push_back(tracks.at(j).xSlope());
        qy_other.push_back(tracks.at(j).y0());
        my_other.push_back(tracks.at(j).ySlope());

	n1++;
        }
 }
}


	if(position_s.at(0).size()!=0){
	for(int h=0;h<position_s.at(0).size(); h++){ gx_s->SetPoint(gx_s->GetN(),trigger.at(0),position_s.at(0).at(h));}
	}
	if(position_s.at(1).size()!=0){
	for(int h=0;h<position_s.at(1).size(); h++){ gy_s->SetPoint(gy_s->GetN(),trigger.at(0),position_s.at(1).at(h));} 
	}

if(position_s.at(2).size()!=0 and position_s.at(3).size()!=0){
for(int h=0;h<position_s.at(2).size(); h++){
        for(int h3=0;h3<position_s.at(3).size(); h3++){
						gx_s->SetPoint(gx_s->GetN(),trigger.at(2),newX(45,-position_s.at(2).at(h),position_s.at(3).at(h3)));
                                                gy_s->SetPoint(gy_s->GetN(),trigger.at(2),newY(45,-position_s.at(2).at(h),position_s.at(3).at(h3)));
		}
	}
}


if(position_s.at(4).size()!=0){
for(int h=0;h<position_s.at(4).size(); h++){ gx_s->SetPoint(gx_s->GetN(),trigger.at(4),position_s.at(4).at(h));} 
}
if(position_s.at(5).size()!=0){
for(int h=0;h<position_s.at(5).size(); h++){ gy_s->SetPoint(gy_s->GetN(),trigger.at(4),position_s.at(5).at(h));} 
}



if(position_s1.at(0).size()!=0){
for(int h=0;h<position_s1.at(0).size(); h++){ gx_s->SetPoint(gx_s->GetN(),trigger.at(6),position_s1.at(0).at(h));} 
}
if(position_s1.at(1).size()!=0){
for(int h=0;h<position_s1.at(1).size(); h++){ gy_s->SetPoint(gy_s->GetN(),trigger.at(6),position_s1.at(1).at(h));} 
}

if(position_s1.at(2).size()!=0 and position_s1.at(3).size()!=0){//cout<<"Position all hits in the station second"<<endl;
for(int h=0;h<position_s1.at(2).size(); h++){
	for(int h3=0;h3<position_s1.at(3).size(); h3++){
						gx_s->SetPoint(gx_s->GetN(),trigger.at(8),newX(45,-position_s1.at(2).at(h),position_s1.at(3).at(h3)));
                                                gy_s->SetPoint(gy_s->GetN(),trigger.at(8),newY(45,-position_s1.at(2).at(h),position_s1.at(3).at(h3)));
//						cout << "new X(uv) " << newX(45,-position_s1.at(2).at(h),position_s1.at(3).at(h3)) << endl; cout << "new Y(uv) " <<newY(45,-position_s1.at(2).at(h),position_s1.at(3).at(h3)) << endl;} //cout <<"----------"<<endl;
		}
	}
}

if(position_s1.at(4).size()!=0){
for(int h=0;h<position_s1.at(4).size(); h++){ gx_s->SetPoint(gx_s->GetN(),trigger.at(10),position_s1.at(4).at(h));} 
}
if(position_s1.at(5).size()!=0){
for(int h=0;h<position_s1.at(5).size(); h++){ gy_s->SetPoint(gy_s->GetN(),trigger.at(10),position_s1.at(5).at(h));} 
}


if(position_s2.at(0).size()!=0){
for(int h=0;h<position_s2.at(0).size(); h++){ gx_s->SetPoint(gx_s->GetN(),other.at(0),position_s2.at(0).at(h));} 
}

if(position_s2.at(1).size()!=0){
for(int h=0;h<position_s2.at(1).size(); h++){ gy_s->SetPoint(gy_s->GetN(),other.at(0),position_s2.at(1).at(h));} 
}

if(position_s2.at(2).size()!=0 and position_s2.at(3).size()!=0){
for(int h=0;h<position_s2.at(2).size(); h++){
	for(int h3=0;h3<position_s2.at(3).size(); h3++){
						gx_s->SetPoint(gx_s->GetN(),other.at(2),newX(45,-position_s2.at(2).at(h),position_s2.at(3).at(h3)));
                                                gy_s->SetPoint(gy_s->GetN(),other.at(2),newY(45,-position_s2.at(2).at(h),position_s2.at(3).at(h3)));
						}
	}
}

if(position_s2.at(4).size()!=0){
for(int h=0;h<position_s2.at(4).size(); h++){ gx_s->SetPoint(gx_s->GetN(),other.at(4),position_s2.at(4).at(h));} 
}
if(position_s2.at(5).size()!=0){
for(int h=0;h<position_s2.at(5).size(); h++){ gy_s->SetPoint(gy_s->GetN(),other.at(4),position_s2.at(5).at(h));} 
}

for(int j=0; j<nt_0; j++){
	if(position_in.at(j).at(0).size()!=0){
	for(int h=0;h<position_in.at(j).at(0).size(); h++){ gx->SetPoint(gx->GetN(),trigger.at(0),position_in.at(j).at(0).at(h));} 
	}
	if(position_in.at(j).at(1).size()!=0){
	for(int h=0;h<position_in.at(j).at(1).size(); h++){ gy->SetPoint(gy->GetN(),trigger.at(0),position_in.at(j).at(1).at(h));} 
	}


	if(position_in.at(j).at(2).size()==position_in.at(j).at(3).size()!=0){
	for(int h=0;h<position_in.at(j).at(2).size(); h++){   gx->SetPoint(gx->GetN(),trigger.at(2),newX(45,-position_in.at(j).at(2).at(h),position_in.at(j).at(3).at(h)));
						gy->SetPoint(gy->GetN(),trigger.at(2),newY(45,-position_in.at(j).at(2).at(h),position_in.at(j).at(3).at(h)));} 
	}

	if(position_in.at(j).at(4).size()!=0){
	for(int h=0;h<position_in.at(j).at(4).size(); h++){ gx->SetPoint(gx->GetN(),trigger.at(4),position_in.at(j).at(4).at(h));} 
	}
	if(position_in.at(j).at(5).size()!=0){
	for(int h=0;h<position_in.at(j).at(5).size(); h++){ gy->SetPoint(gy->GetN(),trigger.at(4),position_in.at(j).at(5).at(h));} 
	}
}


for(int j=0; j<nt; j++){

if(position.at(j).at(0).size()!=0){
for(int h=0;h<position.at(j).at(0).size(); h++){ gx->SetPoint(gx->GetN(),trigger.at(6),position.at(j).at(0).at(h));} 
}

if(position.at(j).at(4).size()!=0){
for(int h=0;h<position.at(j).at(4).size(); h++){ gx->SetPoint(gx->GetN(),trigger.at(10),position.at(j).at(4).at(h));} 
}

if(position.at(j).at(1).size()!=0){
for(int h=0;h<position.at(j).at(1).size(); h++){ gy->SetPoint(gy->GetN(),trigger.at(6),position.at(j).at(1).at(h));} 
}

if(position.at(j).at(5).size()!=0){
for(int h=0;h<position.at(j).at(5).size(); h++){ gy->SetPoint(gy->GetN(),trigger.at(10),position.at(j).at(5).at(h));} 
}


if(position.at(j).at(2).size()!=0 and position.at(j).at(3).size()!=0){
for(int h=0;h<position.at(j).at(2).size(); h++){
        for(int h3=0;h3<position.at(j).at(3).size(); h3++){
                                                gx->SetPoint(gx->GetN(),trigger.at(8),newX(45,-position.at(j).at(2).at(h),position.at(j).at(3).at(h3)));
                                                gy->SetPoint(gy->GetN(),trigger.at(8),newY(45,-position.at(j).at(2).at(h),position.at(j).at(3).at(h3)));//cout << "track X(uv) " << newX(45,-position.at(j).at(2).at(h),position.at(j).at(3).at(h3)) << endl;//cout << "track Y(uv) " <<newY(45,-position.at(j).at(2).at(h),position.at(j).at(3).at(h3)) << endl; //cout << "------------------"<<endl;
		}
        }
 }

}//end sec0 for




for(int j=0; j<nt_other; j++){


if(position_other.at(j).at(0).size()!=0){
for(int h=0;h<position_other.at(j).at(0).size(); h++){ gx->SetPoint(gx->GetN(),other.at(0),position_other.at(j).at(0).at(h));}
}

if(position_other.at(j).at(4).size()!=0){
for(int h=0;h<position_other.at(j).at(4).size(); h++){ gx->SetPoint(gx->GetN(),other.at(4),position_other.at(j).at(4).at(h));}
}

if(position_other.at(j).at(1).size()!=0){
for(int h=0;h<position_other.at(j).at(1).size(); h++){ gy->SetPoint(gy->GetN(),other.at(0),position_other.at(j).at(1).at(h));}
}

if(position_other.at(j).at(5).size()!=0){
for(int h=0;h<position_other.at(j).at(5).size(); h++){ gy->SetPoint(gy->GetN(),other.at(4),position_other.at(j).at(5).at(h));}
}

if(position_other.at(j).at(2).size()!=0 and position_other.at(j).at(3).size()!=0){
for(int h=0;h<position_other.at(j).at(2).size(); h++){
        for(int h3=0;h3<position_other.at(j).at(3).size(); h3++){
                                                gx->SetPoint(gx->GetN(),other.at(2),newX(45,-position_other.at(j).at(2).at(h),position_other.at(j).at(3).at(h3)));
                                                gy->SetPoint(gy->GetN(),other.at(2),newY(45,-position_other.at(j).at(2).at(h),position_other.at(j).at(3).at(h3)));//cout << "track X(uv) " << newX(45,-position_other.at(j).at(2).at(h),position_other.at(j).at(3).at(h3)) << endl;//cout << "track Y(uv) " <<newY(45,-position_other.at(j).at(2).at(h),position_other.at(j).at(3).at(h3)) << endl; //cout << "------------------"<<endl;
		}
        }
 }

}//end sec0 for






}
cout << "---------------------"<<endl;
} //end of general for






std::vector<TGraph*> guv_loc_t;
std::vector<TGraph*> g_mod_u_loc_t;
std::vector<TGraph*> g_mod_v_loc_t;

TCanvas a2("a2","a2",1400,2100);
a2.Divide(2,3);
a2.cd(1);

TLine *t = new TLine(5,-5,5,5);
t->SetLineWidth(4);
t->SetLineColor(47);
TLine *t1 = new TLine(11,-5,11,5);
t1->SetLineWidth(4);
t1->SetLineColor(47);
TLine *l0 = new TLine(0,-5,0,5);
TLine *l1 = new TLine(2,-5,2,5);
TLine *l2 = new TLine(4,-5,4,5);
TLine *l3 = new TLine(6,-5,6,5);
TLine *l4 = new TLine(8,-5,8,5);
TLine *l5 = new TLine(10,-5,10,5);
TLine *l6 = new TLine(12,-5,12,5);
TLine *l7 = new TLine(14,-5,14,5);
TLine *l8 = new TLine(16,-5,16,5);
//TLine *l9 = new TLine(18,-5,18,5);
//TLine *l10 = new TLine(20,-5,20,5);


gx_s->SetMinimum(-6);
gx_s->SetMaximum(6);
gx_s->SetMarkerColor(kBlue);
gx_s->SetTitle("stubs");

gx->SetMinimum(-6);
gx->SetMaximum(6);
gx->SetMarkerColor(kRed);
gx->SetTitle("sig track stubs");


TMultiGraph *mgx = new TMultiGraph();
mgx->SetMinimum(-6);
mgx->SetMaximum(6);
mgx->Add(gx_s,"A*");
mgx->Add(gx,"A*");
mgx->Draw("A* ");
mgx->SetTitle("X projection");


if(qx_in.size()!=0 and mx_in.size()!=0)
{
        for(int c=0; c<qx_in.size(); c++)//               {if(c==0)
		{TLine* lx = new TLine(tar_z1-5, trackXatZ(qx_in.at(c),mx_in.at(c),z1), tar_z1, trackXatZ(qx_in.at(c),mx_in.at(c), z2));
                 lx->SetLineColor(kBlue+c);
                 lx->Draw("same");
//		}
/*		else if(c==1){
		 TLine* lx = new TLine(tar_z2-5, trackXatZ(qx_in.at(c),mx_in.at(c), z2), tar_z2, trackXatZ(qx_in.at(c),mx_in.at(c), z2+88.5227));
                 lx->SetLineColor(kBlue+c);
                 lx->Draw("same");
			}*/
       }
}

//RECO TRACKS
if(qx.size()!=0 and mx.size()!=0)
{ cout << qx.size() << " qx"<<endl;


        for(int c=0; c<qx.size(); c++)
               {
		cout << "linea " << c << endl;
		TLine* lx = new TLine(tar_z1, trackXatZ(qx.at(c),mx.at(c), z2), tar_z2, trackXatZ(qx.at(c),mx.at(c), z2+88.5227));
                lx->SetLineColor(kRed+c);
                //lx->SetTitle(Form("sig track %d",c));
                lx->Draw("same");
        }
}

if(qx_other.size()!=0 and mx_other.size()!=0)
{ cout << qx_other.size() << " qx"<<endl;


        for(int c=0; c<qx_other.size(); c++)
               {
                cout << "linea " << c << endl;
                TLine* lx = new TLine(other_z1, trackXatZ(qx_other.at(c),mx_other.at(c), z3), other_z2, trackXatZ(qx_other.at(c),mx_other.at(c), z3+88.5227));
                lx->SetLineColor(kRed+c);
                //lx->SetTitle(Form("sig track %d",c));
                lx->Draw("same");
        }
}

gPad->BuildLegend(0.25,0.15,0.25,0.15);
l0->Draw("same");
l1->Draw("same");
l2->Draw("same");
t->Draw("same");
t1->Draw("same");
l3->Draw("same");
l4->Draw("same");
l5->Draw("same");
l6->Draw("same");
l7->Draw("same");
l8->Draw("same");
//l9->Draw("same");
//l10->Draw("same");

a2.cd(2);

gx_v1->SetMinimum(-6);
gx_v1->SetMaximum(6);
gx_v1->SetMarkerColor(kOrange);
gx_v1->SetTitle("vertex");
gx_v1->SetMarkerStyle(29);
gx_v1->SetMarkerSize(5);

gx_s->SetMinimum(-6);
gx_s->SetMaximum(6);
gx_s->SetMarkerColor(kBlue);
gx_s->SetTitle("stubs");

gx->SetMinimum(-6);
gx->SetMaximum(6);
gx->SetMarkerColor(kRed);
gx->SetTitle("sig track stubs");

TMultiGraph *mgx_v = new TMultiGraph();
mgx_v->SetMinimum(-6);
mgx_v->SetMaximum(6);
mgx_v->Add(gx_v1,"A*");
mgx_v->Add(gx_s,"A*");
mgx_v->Add(gx,"A*");
mgx_v->Draw("A* ");
mgx_v->SetTitle("X projection vrtx fit");



if(qx_in_v.size()!=0 and mx_in_v.size()!=0)
{
        for(int c=0; c<qx_in_v.size(); c++)
               {if(c==0)
                {//TLine* lx = new TLine(0, trackXatZ(qx_in_v.at(c),mx_in_v.at(c), 0.), 5, trackXatZ(qx_in_v.at(c),mx_in_v.at(c),z-810.));
		 TLine* lx = new TLine(tar_z1-5, trackXatZ(qx_in_v.at(c),mx_in_v.at(c), -88.5227), tar_z1, trackXatZ(qx_in_v.at(c),mx_in_v.at(c),0.));
                 lx->SetLineColor(kBlue+c);
                 lx->Draw("same");
        cout << "trackX in Z " << trackXatZ(qx_in_v.at(c),mx_in_v.at(c),0.) << " VS " << x << endl;

                }
                else if(c==1){
                 TLine* lx = new TLine(tar_z2-5, trackXatZ(qx_in_v.at(c),mx_in_v.at(c), 0.), tar_z2, trackXatZ(qx_in_v.at(c),mx_in_v.at(c),88.5227));
                 lx->SetLineColor(kBlue+c);
                 lx->Draw("same");
                        } 
       }
}

if(qx_v.size()!=0 and mx_v.size()!=0)
{ cout << qx_v.size() << " qx_v"<<endl;
        for(int c=0; c<qx_v.size(); c++)
               {TLine* lx = new TLine(tar_z1, trackXatZ(qx_v.at(c),mx_v.at(c), 0.), tar_z2, trackXatZ(qx_v.at(c),mx_v.at(c), 88.5227));
                lx->SetLineColor(kOrange+c);
                //lx->SetTitle("bkg track");
                lx->Draw("same");
        }
}

if(qx_other.size()!=0 and mx_other.size()!=0)
{ cout << qx_other.size() << " qx"<<endl;


        for(int c=0; c<qx_other.size(); c++)
               {
                cout << "linea " << c << endl;
                TLine* lx = new TLine(other_z1, trackXatZ(qx_other.at(c),mx_other.at(c), z3), other_z2, trackXatZ(qx_other.at(c),mx_other.at(c), z3+88.5227));
                lx->SetLineColor(kRed+c);
                //lx->SetTitle(Form("sig track %d",c));
                lx->Draw("same");
        }
}

gPad->BuildLegend(0.25,0.15,0.25,0.15);

l0->Draw("same");
l1->Draw("same");
l2->Draw("same");
t->Draw("same");
t1->Draw("same");
l3->Draw("same");
l4->Draw("same");
l5->Draw("same");
l6->Draw("same");
l7->Draw("same");
l8->Draw("same");
//l9->Draw("same");
//l10->Draw("same");

a2.cd(3);

gy_s->SetMinimum(-6);
gy_s->SetMaximum(6);
gy_s->SetMarkerColor(kBlue);
gy_s->SetTitle("stubs");

gy->SetMinimum(-6);
gy->SetMaximum(6);
gy->SetMarkerColor(kRed);
gy->SetTitle("sig track stubs");


TMultiGraph *mg = new TMultiGraph();
mg->SetMinimum(-6);
mg->SetMaximum(6);
mg->Add(gy_s,"A*");
mg->Add(gy,"A*");
mg->Draw("A* ");
mg->SetTitle("Y projection");

if(qy_in.size()!=0 and my_in.size()!=0)
{
        for(int c=0; c<qy_in.size(); c++)////               {if(c==0)
                {TLine* ly = new TLine(tar_z1-5, trackYatZ(qy_in.at(c),my_in.at(c), z1), tar_z1, trackYatZ(qy_in.at(c),my_in.at(c), z2));
                 ly->SetLineColor(kBlue+c);
                 ly->Draw("same");
	        cout << "pre vrtx: trackY in Z " << trackYatZ(qy_in.at(c),my_in.at(c),912.7) << " VS " << y << endl;

/*                }
                else if(c==1){ 
                 TLine* ly = new TLine(tar_z2-5, trackYatZ(qy_in.at(c),my_in.at(c), z2), tar_z2, trackYatZ(qy_in.at(c),my_in.at(c), z2+92.3722));
                 ly->SetLineColor(kBlue+c);
                 ly->Draw("same");
                        }*/
        }
}

if(qy.size()!=0 and my.size()!=0)
{ cout << qy.size() << " qy"<<endl;
        for(int c=0; c<qy.size(); c++)
               {TLine* ly = new TLine(tar_z1, trackYatZ(qy.at(c),my.at(c), z2), tar_z2, trackYatZ(qy.at(c),my.at(c), z2+92.3722));
                ly->SetLineColor(kRed+c);
                //ly->SetTitle(Form("sig track %d",c));
                ly->Draw("same");
	}
}

if(qy_other.size()!=0 and my_other.size()!=0)
{ cout << qy_other.size() << " qy"<<endl;


        for(int c=0; c<qy_other.size(); c++)
               {
                cout << "linea " << c << endl;
                TLine* ly = new TLine(other_z1, trackYatZ(qy_other.at(c),my_other.at(c), z3), other_z2, trackYatZ(qy_other.at(c),my_other.at(c), z3+92.3722));
                ly->SetLineColor(kRed+c);
                //ly->SetTitle(Form("sig track %d",c));
                ly->Draw("same");
        }
}


gPad->BuildLegend(0.25,0.15,0.25,0.15);
l0->Draw("same");
l1->Draw("same");
l2->Draw("same");
t->Draw("same");
t1->Draw("same");
l3->Draw("same");
l4->Draw("same");
l5->Draw("same");
l6->Draw("same");
l7->Draw("same");
l8->Draw("same");
//l9->Draw("same");
//l10->Draw("same");

a2.cd(4);

gy_v1->SetMinimum(-6);
gy_v1->SetMaximum(6);
gy_v1->SetMarkerColor(kOrange);
gy_v1->SetTitle("vertex");
gy_v1->SetMarkerStyle(29);
gy_v1->SetMarkerSize(5);

gy_s->SetMinimum(-6);
gy_s->SetMaximum(6);
gy_s->SetMarkerColor(kBlue);
gy_s->SetTitle("stubs");

gy->SetMinimum(-6);
gy->SetMaximum(6);
gy->SetMarkerColor(kRed);
gy->SetTitle("sig track stubs");


TMultiGraph *mgy_v = new TMultiGraph();
mgy_v->SetMinimum(-6);
mgy_v->SetMaximum(6);
mgy_v->Add(gy_v1,"A*");
mgy_v->Add(gy_s,"A*");
mgy_v->Add(gy,"A*");
mgy_v->Draw("A* ");
mgy_v->SetTitle("Y projection vrtx fit");



if(qy_in_v.size()!=0 and my_in_v.size()!=0)
{
        for(int c=0; c<qy_in_v.size(); c++)
               {if(c==0)
                {TLine* ly = new TLine(tar_z1-5, trackYatZ(qy_in_v.at(c),my_in_v.at(c), -92.3722), tar_z1, trackYatZ(qy_in_v.at(c),my_in_v.at(c),0.));
                 ly->SetLineColor(kBlue+c);
                 ly->Draw("same");
	cout << "post vrtx: trackY in Z " << trackYatZ(qy_in_v.at(c),my_in_v.at(c),0.) << " VS " << y << endl;
                }
                else if(c==1){
                 TLine* ly = new TLine(tar_z2-5, trackYatZ(qy_in_v.at(c),my_in_v.at(c), 0.), tar_z2, trackYatZ(qy_in_v.at(c),my_in_v.at(c), 92.3722));
                 ly->SetLineColor(kBlue+c);
                 ly->Draw("same");
                        }
	}
}

if(qy_v.size()!=0 and my_v.size()!=0)
{ cout << qy_v.size() << " qy_v"<<endl;
        for(int c=0; c<qy_v.size(); c++)
               {TLine* ly = new TLine(tar_z1, trackYatZ(qy_v.at(c),my_v.at(c), 0.), tar_z2, trackYatZ(qy_v.at(c),my_v.at(c), 92.3722));
                ly->SetLineColor(kOrange+c);
                //lx->SetTitle("bkg track");
                ly->Draw("same");
        }
}
if(qy_other.size()!=0 and my_other.size()!=0)
{ cout << qy_other.size() << " qy"<<endl;


        for(int c=0; c<qy_other.size(); c++)
               {
                cout << "linea " << c << endl;
                TLine* ly = new TLine(other_z1, trackYatZ(qy_other.at(c),my_other.at(c), z3), other_z2, trackYatZ(qy_other.at(c),my_other.at(c), z3+92.3722));
                ly->SetLineColor(kRed+c);
                //ly->SetTitle(Form("sig track %d",c));
                ly->Draw("same");
        }
}
gPad->BuildLegend(0.25,0.15,0.25,0.15);

t->Draw("same");
t1->Draw("same");
l0->Draw("same");
l1->Draw("same");
l2->Draw("same");
l3->Draw("same");
l4->Draw("same");
l5->Draw("same");
l6->Draw("same");
l7->Draw("same");
l8->Draw("same");
//l9->Draw("same");
//l10->Draw("same");

a2.SaveAs(Form("pdf_exercise_%d_version%d_int%i.pdf",nhits,version,interaction));

}

