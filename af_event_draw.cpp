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

	TFile *inputfile = new TFile("TRMesmer_2_2stations.root");//nrrowBP_20GeV_1shared.root");//trPoints,20GeV.root");//exampleProductionJob.root");narrowBP.root
        TTree* cbmsim = (TTree*) inputfile->Get("cbmsim");

        TClonesArray *MCTrack = 0;
        TClonesArray *SignalTracks = 0;
        TClonesArray *TrackerStripDigis = 0;
        TClonesArray *TrackerPoints = 0;
        TClonesArray *TrackerStubs = 0;
        MuE::Event *MesmerEvent = 0;
        MUonERecoOutput *ReconstructionOutput = 0;

        cbmsim->SetBranchAddress("MCTrack", &MCTrack);
        cbmsim->SetBranchAddress("SignalTracks", &SignalTracks);
        cbmsim->SetBranchAddress("TrackerPoints", &TrackerPoints);
//        cbmsim->SetBranchAddress("TrackerStripDigis", &TrackerStripDigis);
        cbmsim->SetBranchAddress("TrackerStubs", &TrackerStubs);
//        cbmsim->SetBranchAddress("MesmerEvent", &MesmerEvent);
        cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);


double the_gen=0; double thmu_gen=0; double thmu=0; double the=0;
double signal=0; double reco=0; double reco1=0; double more_reco=0; double reco0=0;

auto gx = new TGraph();
auto gx_no = new TGraph();
auto gx_else = new TGraph();
auto gx_v = new TGraph();
auto gy_v = new TGraph();
auto gx_v1 = new TGraph();
auto gy_v1 = new TGraph();
auto gy = new TGraph();
auto gy_no = new TGraph();
auto gy_else = new TGraph();

int yes_e=0;int yes_mu=0; int yes2=0;
int point_mu=0; int point_el=0;
int code_mu=-99; int code_e=-99;
int TrackIdreco=-99;

std::vector<double> qx_in;
std::vector<double> mx_in;

std::vector<double> qy_in;
std::vector<double> my_in;

std::vector<double> qx;
std::vector<double> mx;

std::vector<double> qy;
std::vector<double> my;

std::vector<double> qx_else;
std::vector<double> mx_else;

std::vector<double> qy_else;
std::vector<double> my_else;

std::vector<double> qx_in_v;
std::vector<double> mx_in_v;

std::vector<double> qy_in_v;
std::vector<double> my_in_v;

std::vector<double> qx_v;
std::vector<double> mx_v;

std::vector<double> qy_v;
std::vector<double> my_v;

double z;
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



for(Long64_t i = numb; i < numb+1; i++) {//for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
		cbmsim->GetEntry(i);
		if(i%1000 == 0) cout<<"Entry "<<i<<endl;
	TVector3 pmuin,pe,pmu;

	for(int n = 0; n < MCTrack->GetEntries(); n++) {
	 const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
         if(MCTr->interactionID()==45 and MCTr->pdgCode()==11) code_e=n;
         if(MCTr->interactionID()==45 and MCTr->pdgCode()==13) code_mu=n;
	}

	for(int n = 0; n < SignalTracks->GetEntries(); n++) {
	const MUonETrack *SigTracks = static_cast<const MUonETrack*>(SignalTracks->At(n));
	 if(SignalTracks->GetEntries()>1 and SigTracks->interactionID()==45 )
	 { 
           int last_modXmu=0; int last_modXe=0; //= qx+mx*(189.9218);
           int last_modYmu=0; int last_modYe=0; // = qy+my*(193.7693);
	   int stereo_mu=0; int stereo_e=0;

                         for(int s=0; s<TrackerPoints->GetEntries(); s++)
                         {const MUonETrackerPoint *TrackerPt = static_cast<const MUonETrackerPoint*>(TrackerPoints->At(s));
                          if(TrackerPt->trackPDGCode()==11 and SigTracks->pdgCode()==11 and TrackerPt->trackID()==code_e){ point_el++; 
                                                        cout << s << ") point_el "<< point_el << " mod" << TrackerPt->moduleID() << endl;
                                                                                                 if(TrackerPt->moduleID()==4) last_modXe++; 
                                                                                                 if(TrackerPt->moduleID()==5) last_modYe++;
                                                                                                 if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_e++;}
                          if(TrackerPt->trackPDGCode()==13 and SigTracks->pdgCode()==13 and TrackerPt->trackID()==code_mu){ point_mu++;
                                                        cout << s << ") point_mu "<< point_mu << " mod "<<  TrackerPt->moduleID() << endl;
                                                                                                 if(TrackerPt->moduleID()==4) last_modXmu++;
                                                                                                 if(TrackerPt->moduleID()==5) last_modYmu++;
                                                                                                 if(TrackerPt->moduleID()==2 or TrackerPt->moduleID()==3) stereo_mu++;}
                         }

           if(SigTracks->pdgCode()==11 and last_modXe==2 and last_modYe==2 and stereo_e>1){ //and abs(SigTracks->startX())<4 and abs(SigTracks->startY())<4 and abs(last_modX)<4 and abs(last_modY)<4){ 
		yes_e=1;}//if(the_sdr<0.035)yes_e=1;}
 	   if(SigTracks->pdgCode()==13 and last_modXmu==2 and last_modYmu==2 and stereo_mu>1){//and abs(SigTracks->startX())<4 and abs(SigTracks->startY())<4 and abs(last_modX)<4 and abs(last_modY)<4){ 
		yes_mu=1;}//if(thmu_sdr>0.0001)yes_mu=1;}
	 }
	}
 if(yes_e!=1 or yes_mu!=1) cout << "NOT RECONSTRUCTIBLE" << endl;

double the_rec;

 if(yes_e==1 and yes_mu==1){

cout << "RECONSTRUCTIBLE" << endl;

//vector<MUonERecoOutputVertex> vrtx = ReconstructionOutput->reconstructedVertices();
vector<MUonERecoOutputAdaptiveFitterVertex> vrtx = ReconstructionOutput->adaptiveFitterVertices();
z=vrtx.at(0).z();

for(int j=0; j<vrtx.size();j++)
{

if(vrtx.at(j).stationIndex()==1) //and tracks.at(0).processIDofLinkedTrack()==45 and tracks.at(0).linkedTrackID()!=tracks.at(1).linkedTrackID()){
{
///////
vector<MUonERecoOutputTrack> tr_ad = vrtx.at(j).tracks();
/////

cout << "vrtx.size() " << vrtx.size() << endl;
        gx_v->SetPoint(gx_v->GetN(),6,vrtx.at(j).x());
	gy_v->SetPoint(gy_v->GetN(),6,vrtx.at(j).y());
/*
 MUonERecoOutputTrack mu_in = vrtx.at(j).incomingMuon();
 MUonERecoOutputTrack mu_out = vrtx.at(j).outgoingMuon();
 MUonERecoOutputTrack e_out = vrtx.at(j).outgoingElectron();
*/
 MUonERecoOutputTrack mu_in;
 MUonERecoOutputTrack mu_out;
 MUonERecoOutputTrack e_out;

   if(vrtx.at(j).stationIndex()==1) //and tracks.at(0).processIDofLinkedTrack()==45 and tracks.at(0).linkedTrackID()!=tracks.at(1).linkedTrackID()){
        {


	for(int v=0; v<vrtx.at(j).tracks().size(); v++)
        {
         if(tr_ad.at(v).processIDofLinkedTrack()==45 and vrtx.at(j).tracks().size()>=2 and tr_ad.at(v).sector()==1 and code_mu==tr_ad.at(v).linkedTrackID()) {mu_out = tr_ad.at(v); cout << mu_out.xSlope() <<endl;}
         if(tr_ad.at(v).processIDofLinkedTrack()==45 and vrtx.at(j).tracks().size()>=2 and tr_ad.at(v).sector()==1 and code_e==tr_ad.at(v).linkedTrackID()) e_out =  tr_ad.at(v);
        }
}

if(j==0){

        gx_v1->SetPoint(gx_v1->GetN(),6,vrtx.at(j).x());
        gy_v1->SetPoint(gy_v1->GetN(),6,vrtx.at(j).y());

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
		}
       }
}



std::array<std::vector<double>,6> position;//={{{0.,0.},{0.,0.},{0.,0.},{0.,0.},{0.,0.},{0.,0.}}};
std::array<std::vector<double>,6> position_else;//={{{0.,0.},{0.,0.},{0.,0.},{0.,0.},{0.,0.},{0.,0.}}};

vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
int vtx = ReconstructionOutput->adaptiveFitterVerticesMultiplicity();

int sig=0;
vector<double> the_rec_v;
the_rec_v.reserve(5);
vector<double> chi_min;
chi_min.reserve(5);

cout << "Track size " << tracks.size() << endl;
double px0;
double py0;
for(int j=0; j<tracks.size();j++)
{

if(tracks.at(j).processIDofLinkedTrack()==0 and tracks.at(j).sector()==0)
{
std::vector<MUonERecoOutputTrackHit> hits_=tracks.at(j).hits();
for(int h=0;h<hits_.size();h++){
if(hits_.at(h).moduleID()==0) px0=hits_.at(h).position();
if(hits_.at(h).moduleID()==1) py0=hits_.at(h).position();
}

qx_in.push_back(tracks.at(j).x0());
mx_in.push_back(tracks.at(j).xSlope());
qy_in.push_back(tracks.at(j).y0());
my_in.push_back(tracks.at(j).ySlope());

}
if(tracks.at(j).processIDofLinkedTrack()==45) { TrackIdreco=tracks.at(j).linkedTrackID();}
if(tracks.at(j).processIDofLinkedTrack()==45 and tracks.size()>=2 and tracks.at(j).sector()==1) //and tracks.at(0).processIDofLinkedTrack()==45 and tracks.at(0).linkedTrackID()!=tracks.at(1).linkedTrackID()){
{yes2++;
 int sum = tracks.at(j).numberOfXProjectionHits() + tracks.at(j).numberOfYProjectionHits() + tracks.at(j).numberOfStereoHits();
 cout << "TrackerStubs " << TrackerStubs->GetEntries() << endl;
 cout << "TrackID " << tracks.at(j).linkedTrackID() << endl;
 cout << "#%hitsshared " <<tracks.at(j).percentageOfHitsSharedWithLinkedTrack() << " and sum of hits " << sum << endl;
 cout << "chi2perDegreeOfFreedom " << tracks.at(j).chi2perDegreeOfFreedom() << endl;

	if(tracks.at(j).linkedTrackID()==code_e) {the_rec_v.push_back(sqrt(tracks.at(j).xSlope() * tracks.at(j).xSlope() + tracks.at(j).ySlope() * tracks.at(j).ySlope())); chi_min.push_back(tracks.at(j).chi2perDegreeOfFreedom());}

std::vector<MUonERecoOutputTrackHit> hits_=tracks.at(j).hits();

qx.push_back(tracks.at(j).x0());
mx.push_back(tracks.at(j).xSlope());
qy.push_back(tracks.at(j).y0());
my.push_back(tracks.at(j).ySlope());

for(int h=0;h<hits_.size();h++){
position.at(hits_.at(h).moduleID()).push_back(hits_.at(h).position());
		}
	}

if(tracks.at(j).processIDofLinkedTrack()!=45 or tracks.size()<2) {

if(tracks.at(j).processIDofLinkedTrack()!=0)
{cout << " Id!= 45 ---> " << tracks.at(j).processIDofLinkedTrack() << endl;
qx_else.push_back(tracks.at(j).x0());
mx_else.push_back(tracks.at(j).xSlope());
qy_else.push_back(tracks.at(j).y0());
my_else.push_back(tracks.at(j).ySlope());
std::vector<MUonERecoOutputTrackHit> hits_=tracks.at(j).hits();
for(int h=0;h<hits_.size();h++){
position_else.at(hits_.at(h).moduleID()).push_back(hits_.at(h).position());
	                }
		}
	}
 else if(tracks.at(j).processIDofLinkedTrack()==0 and tracks.at(j).sector()==1)
	{
	 qx_in.push_back(tracks.at(j).x0());
	 mx_in.push_back(tracks.at(j).xSlope());
	 qy_in.push_back(tracks.at(j).y0());
	 my_in.push_back(tracks.at(j).ySlope());
	}

}
if(the_rec_v.size()!=0){ auto it = min_element(chi_min.begin(),chi_min.end()); the_rec = the_rec_v.at(std::distance(chi_min.begin(), it));}


if(yes2==1){
cout <<"NOT RECONSTRUCTED"<<endl;
}

std::array<std::vector<double>,6> position_no;

for(int s=0; s<TrackerStubs->GetEntries(); s++){
double pos; 
         const MUonETrackerStub *TrS = static_cast<const MUonETrackerStub*>(TrackerStubs->At(s));
	 if(TrS->moduleID()==0 | TrS->moduleID()==4){pos = rotationU(TrS->seedClusterCenterStrip()); pos = pos*cos(0.233); position_no.at(TrS->moduleID()).push_back(pos);}
         if(TrS->moduleID()==1 | TrS->moduleID()==5){pos = rotationV(TrS->seedClusterCenterStrip()); pos = pos*cos(0.233); position_no.at(TrS->moduleID()).push_back(pos);}
         if(TrS->moduleID()==2){ pos = rotationU(TrS->seedClusterCenterStrip()); position_no.at(TrS->moduleID()).push_back(pos);}
         if(TrS->moduleID()==3){ pos = rotationV(TrS->seedClusterCenterStrip()); position_no.at(TrS->moduleID()).push_back(pos);}
	}
gx_no->SetPoint(gx_no->GetN(),0,px0);
gy_no->SetPoint(gy_no->GetN(),0,py0);

for(int h=0;h<position_no.at(0).size(); h++){ gx_no->SetPoint(gx_no->GetN(),8,position_no.at(0).at(h));}
for(int h=0;h<position_no.at(4).size(); h++){ gx_no->SetPoint(gx_no->GetN(),12,position_no.at(4).at(h));}

/*
if(position_no.at(2).size()!=0 and position_no.at(2).size()==position_no.at(3).size())
{
for(int h=0;h<position_no.at(2).size(); h++){
        double X = newX(-45,position_no.at(2).at(h),position_no.at(3).at(h));
        double Y = newY(45,position_no.at(2).at(h),position_no.at(3).at(h));
gx_no->SetPoint(gx_no->GetN(),10,X);
gy_no->SetPoint(gy_no->GetN(),10,Y);
        }
}
else if(position_no.at(2).size()!=position_no.at(3).size()){

if(position_no.at(2).size()!=0){for(int h=0;h<position_no.at(2).size(); h++){ gx_no->SetPoint(gx_no->GetN(),10,position_no.at(2).at(h));} }
if(position_no.at(3).size()!=0){for(int h=0;h<position_no.at(3).size(); h++){ gy_no->SetPoint(gy_no->GetN(),10,position_no.at(3).at(h));} }
}*/

for(int h=0;h<position_no.at(1).size(); h++){ gy_no->SetPoint(gy_no->GetN(),8,position_no.at(1).at(h));}
for(int h=0;h<position_no.at(5).size(); h++){ gy_no->SetPoint(gy_no->GetN(),12,position_no.at(5).at(h));}

if(position.at(0).size()!=0){
for(int h=0;h<position.at(0).size(); h++){ gx->SetPoint(gx->GetN(),8,position.at(0).at(h));} 
}

if(position.at(2).size()!=0 and position.at(2).size()==position.at(3).size())
{
for(int h=0;h<position.at(2).size(); h++){
	double X = newX(135,position.at(2).at(h),position.at(3).at(h));
 	double Y = newY(45,position.at(2).at(h),position.at(3).at(h));
gx->SetPoint(gx->GetN(),10,X);
gy->SetPoint(gy->GetN(),10,Y);
	}
}
else if(position.at(2).size()!=position.at(3).size()){

if(position.at(2).size()!=0){for(int h=0;h<position.at(2).size(); h++){ gx->SetPoint(gx->GetN(),10,position.at(2).at(h));} }
if(position.at(3).size()!=0){for(int h=0;h<position.at(3).size(); h++){ gy->SetPoint(gy->GetN(),10,position.at(3).at(h));} }
}

if(position.at(4).size()!=0){
for(int h=0;h<position.at(4).size(); h++){ gx->SetPoint(gx->GetN(),12,position.at(4).at(h));} 
}
if(position.at(1).size()!=0){
for(int h=0;h<position.at(1).size(); h++){ gy->SetPoint(gy->GetN(),8,position.at(1).at(h));} 
}
if(position.at(5).size()!=0){
for(int h=0;h<position.at(5).size(); h++){ gy->SetPoint(gy->GetN(),12,position.at(5).at(h));} 
}

if(position_else.at(0).size()!=0){
for(int h=0;h<position_else.at(0).size(); h++){ gx_else->SetPoint(gx_else->GetN(),8,position_else.at(0).at(h));} 
}

if(position_else.at(2).size()!=0 and position_else.at(2).size()==position_else.at(3).size()) 
{
for(int h=0;h<position_else.at(2).size(); h++){ 
        double X = newX(-45,position_else.at(2).at(h),position_else.at(3).at(h));
        double Y = newY(45,position_else.at(2).at(h),position_else.at(3).at(h));
gx_else->SetPoint(gx_else->GetN(),10,X);
gy_else->SetPoint(gy_else->GetN(),10,Y);
        }
}
else if(position_else.at(2).size()!=position_else.at(3).size()){

if(position_else.at(2).size()!=0){for(int h=0;h<position_else.at(2).size(); h++){ gx_else->SetPoint(gx_else->GetN(),10,position_else.at(2).at(h));} }
if(position_else.at(3).size()!=0){for(int h=0;h<position_else.at(3).size(); h++){ gy_else->SetPoint(gy_else->GetN(),10,position_else.at(3).at(h));} }
}

if(position_else.at(4).size()!=0){
for(int h=0;h<position_else.at(4).size(); h++){ gx_else->SetPoint(gx_else->GetN(),12,position_else.at(4).at(h));} 
}
if(position_else.at(1).size()!=0){
for(int h=0;h<position_else.at(1).size(); h++){ gy_else->SetPoint(gy_else->GetN(),8,position_else.at(1).at(h));} 
}
if(position_else.at(5).size()!=0){
for(int h=0;h<position_else.at(5).size(); h++){ gy_else->SetPoint(gy_else->GetN(),12,position_else.at(5).at(h));} 
}

}
cout << "---------------------"<<endl;
yes_e=0;yes_mu=0;
code_e=-99;code_mu=-99;
point_el=0;point_mu=0;
TrackIdreco=-99;
yes2=0;
} //end of general for

double ratio =reco/signal;
double ratio1 =reco1/signal;
double ratio0 =reco0/signal;
double ratioM =more_reco/signal;

cout << "Su " << signal << " eventi di segnale, " << reco << " sono ricostruiti, con un rapporto del " << ratio*100 << "%"<< endl;
cout << "Su " << signal << " eventi di segnale con piu tracce (stesso id), " << more_reco << " sono ricostruiti, con un rapporto del " << ratioM*100 << "%"<< endl;
cout << "Su " << signal << " eventi di segnale con 0 tracce di segnale reco, " << reco0 << ", con un rapporto del " << ratio0*100 << "%"<< endl;
cout << "Su " << signal << " eventi di con 1 sola traccia di segnale reco, " << reco1 << ", con un rapporto del " << ratio1*100 << "%"<< endl;


TCanvas a2("a2","a2",1400,1400);
a2.Divide(2,2);
a2.cd(1);

TLine *t = new TLine(6,-5,6,5);
t->SetLineWidth(4);
t->SetLineColor(47);
TLine *l0 = new TLine(7,-5,7,5);
TLine *l1 = new TLine(8,-5,8,5);
TLine *l2 = new TLine(9,-5,9,5);
TLine *l3 = new TLine(10,-5,10,5);
TLine *l4 = new TLine(11,-5,11,5);
TLine *l5 = new TLine(12,-5,12,5);

gx_v->SetMinimum(-6);
gx_v->SetMaximum(6);
gx_v->SetMarkerColor(kOrange);
gx_v->SetTitle("vertex");
gx_v->SetMarkerStyle(29);
gx_v->SetMarkerSize(5);

gx->SetMinimum(-6);
gx->SetMaximum(6);
gx->SetMarkerColor(kRed);
gx->SetTitle("sig track stubs");

gx_no->SetMinimum(-6);
gx_no->SetMaximum(6);
gx_no->SetMarkerColor(kBlue);
gx_no->SetTitle("digitalizer stubs");

gx_else->SetMinimum(-6);
gx_else->SetMaximum(6);
gx_else->SetMarkerColor(kGreen);
gx_else->SetTitle("bkg track stubs");


TMultiGraph *mgx = new TMultiGraph();
mgx->SetMinimum(-6);
mgx->SetMaximum(6);
mgx->Add(gx_v,"A*");
mgx->Add(gx_no,"A*");
mgx->Add(gx,"A*");
mgx->Add(gx_else,"A*");
mgx->Draw("A* ");
mgx->SetTitle("X projection");


if(qx_in.size()!=0 and mx_in.size()!=0)
{
        for(int c=0; c<qx_in.size(); c++)
               {if(c==0)
		{TLine* lx = new TLine(0, trackXatZ(qx_in.at(c),mx_in.at(c), 931.7), 6, trackXatZ(qx_in.at(c),mx_in.at(c), 931.7+89.9218));
                 lx->SetLineColor(kBlue+c);
                 lx->Draw("same");
		}
		else if(c==1){
		 TLine* lx = new TLine(6, trackXatZ(qx_in.at(c),mx_in.at(c), 1031.7), 12, trackXatZ(qx_in.at(c),mx_in.at(c), 1031.7+89.9218));
                 lx->SetLineColor(kBlue+c);
                 lx->Draw("same");
			} 
       }
}

if(qx.size()!=0 and mx.size()!=0)
{ cout << qx.size() << " qx"<<endl;
        for(int c=0; c<qx.size(); c++)
               {TLine* lx = new TLine(6, trackXatZ(qx.at(c),mx.at(c), 1031.7), 12, trackXatZ(qx.at(c),mx.at(c), 1031.7+89.9218));
                lx->SetLineColor(kRed+c);
                //lx->SetTitle("sig track");
                lx->Draw("same");
        }
}

if(qx_else.size()!=0 and mx_else.size()!=0)
{ cout << qx_else.size() << " qx_else"<<endl;
        for(int c=0; c<qx_else.size(); c++)
               {TLine* lx = new TLine(6, trackXatZ(qx_else.at(c),mx_else.at(c), 1031.7), 12, trackXatZ(qx_else.at(c),mx_else.at(c), 1031.7+89.9218));
                lx->SetLineColor(kGreen+c);
                //lx->SetTitle("bkg track");
                lx->Draw("same");
        }
}
gPad->BuildLegend(0.25,0.15,0.25,0.15);

t->Draw("same");
l1->Draw("same");
//l1->Draw("same");
l3->Draw("same");
//l3->Draw("same");
l5->Draw("same");
//l5->Draw("same");

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

gx_no->SetMinimum(-6);
gx_no->SetMaximum(6);
gx_no->SetMarkerColor(kBlue);
gx_no->SetTitle("digitalizer stubs");

gx_else->SetMinimum(-6);
gx_else->SetMaximum(6);
gx_else->SetMarkerColor(kGreen);
gx_else->SetTitle("bkg track stubs");

TMultiGraph *mgx_v = new TMultiGraph();
mgx_v->SetMinimum(-6);
mgx_v->SetMaximum(6);
mgx_v->Add(gx_v,"A*");
mgx_v->Add(gx_no,"A*");
mgx_v->Add(gx,"A*");
mgx_v->Add(gx_else,"A*");
mgx_v->Draw("A* ");
mgx_v->SetTitle("X projection vrtx fit");



if(qx_in_v.size()!=0 and mx_in_v.size()!=0)
{
        for(int c=0; c<qx_in_v.size(); c++)
               {if(c==0)
                {TLine* lx = new TLine(0, trackXatZ(qx_in_v.at(c),mx_in_v.at(c), 0.), 6, trackXatZ(qx_in_v.at(c),mx_in_v.at(c),89.9218));
                 lx->SetLineColor(kBlue+c);
                 lx->Draw("same");
                }
                else if(c==1){
                 TLine* lx = new TLine(6, trackXatZ(qx_in_v.at(c),mx_in_v.at(c), 0.), 12, trackXatZ(qx_in_v.at(c),mx_in_v.at(c),89.9218));
                 lx->SetLineColor(kBlue+c);
                 lx->Draw("same");
                        } 
       }
}

if(qx_v.size()!=0 and mx_v.size()!=0)
{ cout << qx_v.size() << " qx_v"<<endl;
        for(int c=0; c<qx_v.size(); c++)
               {TLine* lx = new TLine(6, trackXatZ(qx_v.at(c),mx_v.at(c), z), 12, trackXatZ(qx_v.at(c),mx_v.at(c),z+89.9218));
                lx->SetLineColor(kOrange+c);
                //lx->SetTitle("bkg track");
                lx->Draw("same");
        }
}

gPad->BuildLegend(0.25,0.15,0.25,0.15);

t->Draw("same");
l1->Draw("same");
//l1->Draw("same");
l3->Draw("same");
//l3->Draw("same");
l5->Draw("same");
//l5->Draw("same");

a2.cd(3);

gy_v->SetMinimum(-6);
gy_v->SetMaximum(6);
gy_v->SetMarkerStyle(29);
gy_v->SetMarkerColor(kOrange);
gy_v->SetTitle("vertex");
gy_v->SetMarkerSize(5);

gy->SetMinimum(-6);
gy->SetMaximum(6);
gy->SetMarkerColor(kRed);

gy_no->SetMinimum(-6);
gy_no->SetMaximum(6);
gy_no->SetMarkerColor(kBlue);

gy_else->SetMinimum(-6);
gy_else->SetMaximum(6);
gy_else->SetMarkerColor(kGreen);

gy->SetTitle("sig track stubs");
gy_no->SetTitle("digitalizer stubs");
gy_else->SetTitle("bkg track stubs");

TMultiGraph *mg = new TMultiGraph();
mg->SetMinimum(-6);
mg->SetMaximum(6);
mg->Add(gy_v,"A*");
mg->Add(gy_no,"A*");
mg->Add(gy,"A*");
mg->Add(gy_else,"A*");
mg->Draw("A* ");
mg->SetTitle("Y projection");

if(qy_in.size()!=0 and my_in.size()!=0)
{
        for(int c=0; c<qy_in.size(); c++)
               {if(c==0)
                {TLine* ly = new TLine(0, trackYatZ(qy_in.at(c),my_in.at(c), 931.7), 6, trackYatZ(qy_in.at(c),my_in.at(c), 931.7+89.9218));
                 ly->SetLineColor(kBlue+c);
                 ly->Draw("same");
                } 
                else if(c==1){ 
                 TLine* ly = new TLine(6, trackYatZ(qy_in.at(c),my_in.at(c), 1031.7), 12, trackYatZ(qy_in.at(c),my_in.at(c), 1031.7+89.9218));
                 ly->SetLineColor(kBlue+c);
                 ly->Draw("same");
                        } 
        }
}

if(qy.size()!=0 and my.size()!=0)
{ cout << qy.size() << " qy"<<endl;
        for(int c=0; c<qy.size(); c++)
               {TLine* ly = new TLine(6, trackYatZ(qy.at(c),my.at(c), 1031.7), 12, trackYatZ(qy.at(c),my.at(c), 1031.7+89.9218));
                ly->SetLineColor(kRed+c);
                //ly->SetTitle("sig track");
                ly->Draw("same");
	}
}

if(qy_else.size()!=0 and my_else.size()!=0)
{ cout << qy_else.size() << " qy_else"<<endl;
	for(int c=0; c<qy_else.size(); c++)
               {TLine* ly = new TLine(6, trackYatZ(qy_else.at(c),my_else.at(c), 1031.7), 12, trackYatZ(qy_else.at(c),my_else.at(c), 1031.7+89.9218));
                ly->SetLineColor(kGreen+c);
		//ly->SetTitle("bkg track");
                ly->Draw("same");
	}
}
gPad->BuildLegend(0.25,0.15,0.25,0.15);
t->Draw("same");
l1->Draw("same");
//l1->Draw("same");
l3->Draw("same");
//l3->Draw("same");
l5->Draw("same");
//l5->Draw("same");

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

gy_no->SetMinimum(-6);
gy_no->SetMaximum(6);
gy_no->SetMarkerColor(kBlue);
gy_no->SetTitle("digitalizer stubs");

gy_else->SetMinimum(-6);
gy_else->SetMaximum(6);
gy_else->SetMarkerColor(kGreen);
gy_else->SetTitle("bkg track stubs");

TMultiGraph *mgy_v = new TMultiGraph();
mgy_v->SetMinimum(-6);
mgy_v->SetMaximum(6);
mgy_v->Add(gy_v,"A*");
mgy_v->Add(gy_no,"A*");
mgy_v->Add(gy,"A*");
mgy_v->Add(gy_else,"A*");
mgy_v->Draw("A* ");
mgy_v->SetTitle("Y projection vrtx fit");



if(qy_in_v.size()!=0 and my_in_v.size()!=0)
{
        for(int c=0; c<qy_in_v.size(); c++)
		{if(c==0)
                {TLine* ly = new TLine(0, trackYatZ(qy_in_v.at(c),my_in_v.at(c), 0.), 6, trackYatZ(qy_in_v.at(c),my_in_v.at(c), 89.9218));
                 ly->SetLineColor(kBlue+c);
                 ly->Draw("same");
                }
                else if(c==1){
                 TLine* lx = new TLine(6, trackXatZ(qx_in_v.at(c),mx_in_v.at(c), 0.), 12, trackXatZ(qx_in_v.at(c),mx_in_v.at(c), 89.9218));
                 lx->SetLineColor(kBlue+c);
                 lx->Draw("same");
                        }
	}
}

if(qy_v.size()!=0 and my_v.size()!=0)
{ cout << qy_v.size() << " qy_v"<<endl;
        for(int c=0; c<qy_v.size(); c++)
               {TLine* ly = new TLine(6, trackYatZ(qy_v.at(c),my_v.at(c), z), 12, trackYatZ(qy_v.at(c),my_v.at(c), z+89.9218));
                ly->SetLineColor(kOrange+c);
                //lx->SetTitle("bkg track");
                ly->Draw("same");
        }
}

gPad->BuildLegend(0.25,0.15,0.25,0.15);

t->Draw("same");
l1->Draw("same");
//l1->Draw("same");
l3->Draw("same");
//l3->Draw("same");
l5->Draw("same");

a2.SaveAs("pdf_exercise.pdf");

}


