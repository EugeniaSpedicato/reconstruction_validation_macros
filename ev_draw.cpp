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

void ev_draw(int numb){

double z_mod[6]={911.2+18.0218,911.2+21.8693,911.2+55.3635,911.2+56.6205,911.2+89.9218,911.2+93.7693};
double alpha[6]={0  *TMath::DegToRad(), 90 *TMath::DegToRad(),135*TMath::DegToRad(), 45 *TMath::DegToRad(),0  *TMath::DegToRad(), 90 *TMath::DegToRad()};

        TFile *inputfile = new TFile("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_0hit.root");
        TTree* cbmsim = (TTree*) inputfile->Get("cbmsim");

/*   TFile *inputfile = new TFile("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_commit_2f4e96f4_MCsignal_bestConfig_0hit.root");
   TFile *f2 = new TFile("/mnt/raid10/DATA/espedica/fairmu/gen_digi/commit_2f4e96f4_MCsignal_SIM-DIGI.root");

   TTree *cbmsim = (TTree*)inputfile->Get("cbmsim");
   TTree *t2 = (TTree*)f2->Get("cbmsim");
        t2->SetEntries(cbmsim->GetEntries());
   cbmsim->AddFriend(t2);
*/


//        TClonesArray *TrackerPoints = 0;
//        cbmsim->SetBranchAddress("TrackerPoints", &TrackerPoints);
        TClonesArray *TrackerStubs = 0;
	cbmsim->SetBranchAddress("TrackerStubs", &TrackerStubs);
        MUonERecoOutput *ReconstructionOutput = 0;
        cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);
       TClonesArray *MCTrack = 0;
        cbmsim->SetBranchAddress("MCTrack", &MCTrack);


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
//std::array<double,6>={912.7+18.0218,912.7+21.8693,912.7+55.3635,912.7+56.6205,912.7+89.9218,912.7+93.7693};
    int sec0=0; int sec1=0;

for(Long64_t i = numb; i < numb+1; i++) {
//for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {

		cbmsim->GetEntry(i);
//		if(i%1000 == 0) cout<<"Entry "<<i<<endl;

cout<<"Entry "<< cbmsim->GetEntry(i)<<endl;

           int hit_modXmu=0; int hit_modXe=0;
           int hit_modYmu=0; int hit_modYe=0;
           int stereo_mu=0; int stereo_e=0;


vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();

cout << "tracks.size() " << tracks.size() << endl;

int sec1=0;
 if(tracks.size()>0){


    for(int j=0; j<tracks.size();j++)
    {
     	if(tracks.at(j).sector()==0) sec0++;
        if(tracks.at(j).sector()==1) sec1++;
        }

cout << "Track size " << tracks.size() << " and " << sec0 << ", " << sec1 << endl;



MUonERecoOutputVertex vrtx = ReconstructionOutput->bestVertex();

x=vrtx.xKinematicFit();
y=vrtx.yKinematicFit();
z=vrtx.zKinematicFit();

 MUonERecoOutputTrack mu_in = vrtx.incomingMuon();
 MUonERecoOutputTrack mu_out = vrtx.outgoingMuon();
 MUonERecoOutputTrack e_out = vrtx.outgoingElectron();

cout << "vrtx chi2 " << vrtx.chi2perDegreeOfFreedom() << endl;

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



std::array<std::vector<double>,6> position;//={{{0.,0.},{0.,0.},{0.,0.},{0.,0.},{0.,0.},{0.,0.}}};
std::array<std::vector<double>,6> position_in;//={{{0.,0.},{0.,0.},{0.,0.},{0.,0.},{0.,0.},{0.,0.}}};
std::array<std::vector<double>,6> position_s;
std::array<std::vector<double>,6> position_s1;
double px0;
double py0;

int code_e,code_mu;

/*	for(int n = 0; n < MCTrack->GetEntries(); n++) {
	 const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
         if(MCTr->interactionID()==45 and MCTr->pdgCode()==11) {code_e=n;}
         if(MCTr->interactionID()==45 and MCTr->pdgCode()==-13) {code_mu=n;}
}*/

/*for(int s=0; s<TrackerPoints->GetEntries(); s++)
                         {const MUonETrackerPoint *TrackerPt = static_cast<const MUonETrackerPoint*>(TrackerPoints->At(s));

if(TrackerPt->trackPDGCode()==-13 and TrackerPt->trackID()==code_mu and TrackerPt->stationID()==1)
				{
				  TVector3 exiting=TrackerPt->exitingPositionGlobalCoordinates();

				  if(TrackerPt->moduleID()==0 or TrackerPt->moduleID()==4) {position_s1.at(TrackerPt->moduleID()).push_back(exiting.X());
											cout <<TrackerPt->moduleID()<< ") mu " << exiting.X() <<endl;}
				  if(TrackerPt->moduleID()==1 or TrackerPt->moduleID()==5) {position_s1.at(TrackerPt->moduleID()).push_back(exiting.Y());
                                                                                        cout <<TrackerPt->moduleID()<< ") mu " << exiting.Y() <<endl;}
				}
if(TrackerPt->trackPDGCode()==11 and TrackerPt->trackID()==code_e and TrackerPt->stationID()==1)
{
                                  TVector3 exiting=TrackerPt->exitingPositionGlobalCoordinates();
                                  if(TrackerPt->moduleID()==0 or TrackerPt->moduleID()==4){ position_s1.at(TrackerPt->moduleID()).push_back(exiting.X());
                                                                                        cout <<TrackerPt->moduleID()<< ") el " << exiting.X() <<endl;}
                                  if(TrackerPt->moduleID()==1 or TrackerPt->moduleID()==5){ position_s1.at(TrackerPt->moduleID()).push_back(exiting.Y());
                                                                                        cout <<TrackerPt->moduleID()<< ") el " << exiting.Y() <<endl;}
}



			 }

*/


/*for(int t=0; t<TrackerStubs->GetEntries(); t++)
                         {const MUonETrackerStub *stubs = static_cast<const MUonETrackerStub*>(TrackerStubs->At(t));
	if(stubs->stationID()==1){ double stub_pos=(stubs->seedClusterCenterStrip() + 0.5 + 0.5 * stubs->bend()) * 9.144 / 1016 - 0.5 * 9.144;
 			  if(stubs->moduleID()==0) {position_s1.at(stubs->moduleID()).push_back(stub_pos);}
                          if( stubs->moduleID()==4) {position_s1.at(stubs->moduleID()).push_back(stub_pos);}
                          if(stubs->moduleID()==1) {position_s1.at(stubs->moduleID()).push_back(stub_pos*sin(90));}
			  if(stubs->moduleID()==5) {position_s1.at(stubs->moduleID()).push_back(stub_pos*sin(90));}
			}
	else if(stubs->stationID()==0){ double stub_pos=(stubs->seedClusterCenterStrip() + 0.5 + 0.5 * stubs->bend()) * 9.144 / 1016 - 0.5 * 9.144;
                          if(stubs->moduleID()==0) {position_s.at(stubs->moduleID()).push_back(stub_pos);}
                          if( stubs->moduleID()==4) {position_s.at(stubs->moduleID()).push_back(stub_pos);}
                          if(stubs->moduleID()==1) {position_s.at(stubs->moduleID()).push_back(stub_pos*sin(90));}
                          if(stubs->moduleID()==5) {position_s.at(stubs->moduleID()).push_back(stub_pos*sin(90));}
                        }
}
*/

vector<MUonERecoOutputHit> stubs_=ReconstructionOutput->reconstructedHits();
for(int s=0; s<stubs_.size(); s++){
if(stubs_.at(s).stationID()==0)position_s.at(stubs_.at(s).moduleID()).push_back(stubs_.at(s).positionPerpendicular());
if(stubs_.at(s).stationID()==1){position_s1.at(stubs_.at(s).moduleID()).push_back(stubs_.at(s).positionPerpendicular());
				cout << stubs_.at(s).moduleID()<< ") stub posPerp " << stubs_.at(s).positionPerpendicular() << endl;}
}

// position stub all event UV local frame
if(position_s1.at(2).size()!=0){
for(int h=0;h<position_s1.at(2).size(); h++){ guv->SetPoint(guv->GetN(),2,position_s1.at(2).at(h)); g_mod_u->SetPoint(g_mod_u->GetN(),0.,position_s1.at(2).at(h));} 
}
if(position_s1.at(3).size()!=0){
for(int h=0;h<position_s1.at(3).size(); h++){ guv->SetPoint(guv->GetN(),3,position_s1.at(3).at(h)); g_mod_v->SetPoint(g_mod_v->GetN(),position_s1.at(3).at(h),0.);} 
}



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
        cout << "pre vrtx: trackY in Z " << trackYatZ(tracks.at(j).y0(),tracks.at(j).ySlope(),912.7) << " VS " << y << endl;

}
if(tracks.at(j).sector()==1 and sec1>=0) //and tracks.at(0).processIDofLinkedTrack()==45 and tracks.at(0).linkedTrackID()!=tracks.at(1).linkedTrackID()){
{
std::vector<MUonERecoOutputHit> hits_=tracks.at(j).hits();
cout << j << ") int ID " <<tracks.at(j).processIDofLinkedTrack()<< endl;
cout << j << ") hits_size " <<hits_.size()<< endl;
qx.push_back(tracks.at(j).x0());
mx.push_back(tracks.at(j).xSlope());
qy.push_back(tracks.at(j).y0());
my.push_back(tracks.at(j).ySlope());

for(int h=0;h<hits_.size();h++){
position.at(hits_.at(h).moduleID()).push_back(hits_.at(h).positionPerpendicular());
cout << hits_.at(h).moduleID() << ") stub track " << j << " posPerp " << hits_.at(h).positionPerpendicular() << endl;
		}
	}
}

// position stub used by tracks UV local frame
if(position.at(2).size()!=0){
for(int h=0;h<position.at(2).size(); h++){ guv_t->SetPoint(guv_t->GetN(),2,position.at(2).at(h)); g_mod_u_t->SetPoint(g_mod_u_t->GetN(),0.,position.at(2).at(h));} 
}
if(position.at(3).size()!=0){
for(int h=0;h<position.at(3).size(); h++){ guv_t->SetPoint(guv_t->GetN(),3,position.at(3).at(h)); g_mod_v_t->SetPoint(g_mod_v_t->GetN(),position.at(3).at(h),0.);} 
}





if(position_s.at(0).size()!=0){
for(int h=0;h<position_s.at(0).size(); h++){ gx_s->SetPoint(gx_s->GetN(),0,position_s.at(0).at(h));} 
}
if(position_s.at(1).size()!=0){
for(int h=0;h<position_s.at(1).size(); h++){ gy_s->SetPoint(gy_s->GetN(),0,position_s.at(1).at(h));} 
}


if(position_s.at(2).size()==position_s.at(3).size()!=0){
for(int h=0;h<position_s.at(2).size(); h++){    gx_s->SetPoint(gx_s->GetN(),2,newX(45,-position_s.at(2).at(h),position_s.at(3).at(h)));
                                                gy_s->SetPoint(gy_s->GetN(),2,newY(45,-position_s.at(2).at(h),position_s.at(3).at(h)));} 
}


if(position_s.at(4).size()!=0){
for(int h=0;h<position_s.at(4).size(); h++){ gx_s->SetPoint(gx_s->GetN(),4,position_s.at(4).at(h));} 
}
if(position_s.at(5).size()!=0){
for(int h=0;h<position_s.at(5).size(); h++){ gy_s->SetPoint(gy_s->GetN(),4,position_s.at(5).at(h));} 
}



if(position_s1.at(0).size()!=0){
for(int h=0;h<position_s1.at(0).size(); h++){ gx_s->SetPoint(gx_s->GetN(),6,position_s1.at(0).at(h));} 
}
if(position_s1.at(1).size()!=0){
for(int h=0;h<position_s1.at(1).size(); h++){ gy_s->SetPoint(gy_s->GetN(),6,position_s1.at(1).at(h));} 
}

/*
cout <<"position_s1.at(2).size() " << position_s1.at(2).size() << " and position_s.at(3).size() " << position_s.at(3).size() << endl;
if(position_s1.at(2).size()==position_s.at(3).size()!=0){
cout<<"quiii"<<endl;
for(int h=0;h<position_s1.at(2).size(); h++){   gx_s->SetPoint(gx_s->GetN(),8,newX(45,-position_s1.at(2).at(h),position_s1.at(3).at(h)));
						cout << "new X(uv) " << newX(45,-position_s1.at(2).at(h),position_s1.at(3).at(h)) << endl;
                                                gy_s->SetPoint(gy_s->GetN(),8,newY(45,-position_s1.at(2).at(h),position_s1.at(3).at(h))); 
						cout << "new Y(uv) " <<newY(45,-position_s1.at(2).at(h),position_s1.at(3).at(h)) << endl;} 
}*/



if(position_s1.at(4).size()!=0){
for(int h=0;h<position_s1.at(4).size(); h++){ gx_s->SetPoint(gx_s->GetN(),10,position_s1.at(4).at(h));} 
}
if(position_s1.at(5).size()!=0){
for(int h=0;h<position_s1.at(5).size(); h++){ gy_s->SetPoint(gy_s->GetN(),10,position_s1.at(5).at(h));} 
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

/*
if(position.at(2).size()==position.at(3).size()!=0){
for(int h=0;h<position.at(2).size(); h++){   gx->SetPoint(gx->GetN(),8,newX(45,-position.at(2).at(h),position.at(3).at(h)));
                                                gy->SetPoint(gy->GetN(),8,newY(45,-position.at(2).at(h),position.at(3).at(h)));} 
}*/

}
cout << "---------------------"<<endl;
} //end of general for


std::vector <TGraph*> guv_loc_t(sec1);
std::vector <TGraph*> g_mod_u_loc_t(sec1);
std::vector <TGraph*> g_mod_v_loc_t(sec1);

TCanvas a2("a2","a2",1400,2100);
a2.Divide(2,3);
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

gx_s->SetMinimum(-6);
gx_s->SetMaximum(6);
gx_s->SetMarkerColor(kBlue);
gx_s->SetTitle("stubs");

TMultiGraph *mgx = new TMultiGraph();
mgx->SetMinimum(-6);
mgx->SetMaximum(6);
mgx->Add(gx,"A*");
mgx->Add(gx_s,"A*");
mgx->Draw("A* ");
mgx->SetTitle("X projection");


if(qx_in.size()!=0 and mx_in.size()!=0)
{
        for(int c=0; c<qx_in.size(); c++)
               {if(c==0)
		{TLine* lx = new TLine(0, trackXatZ(qx_in.at(c),mx_in.at(c),810.), 5, trackXatZ(qx_in.at(c),mx_in.at(c), 912.7));
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

//RECO TRACKS
if(qx.size()!=0 and mx.size()!=0)
{ cout << qx.size() << " qx"<<endl;


        for(int c=0; c<qx.size(); c++)
               {
          double tracklocxy2 = (qx.at(c) + mx.at(c)*z_mod[2])*cos(alpha[2]) +
                               (qy.at(c) + my.at(c)*z_mod[2])*sin(alpha[2]);
          double tracklocxy3 = (qx.at(c) + mx.at(c)*z_mod[3])*cos(alpha[3]) +
                               (qy.at(c) + my.at(c)*z_mod[3])*sin(alpha[3]);

guv_loc_t.push_back(new TGraph());
g_mod_u_loc_t.push_back(new TGraph());
g_mod_v_loc_t.push_back(new TGraph());

guv_loc_t.at(c)->SetPoint(guv_loc_t.at(c)->GetN(),2,tracklocxy2);
guv_loc_t.at(c)->SetPoint(guv_loc_t.at(c)->GetN(),3,tracklocxy3);

g_mod_u_loc_t.at(c)->SetPoint(g_mod_u_loc_t.at(c)->GetN(),0.,tracklocxy2);
g_mod_v_loc_t.at(c)->SetPoint(g_mod_v_loc_t.at(c)->GetN(),tracklocxy3,0.);

cout << "Local pos (posPerp) of the track " << c <<" in U " << tracklocxy2 << endl;
cout << "Local pos (posPerp) of the track " << c <<" in V " << tracklocxy3 << endl;

		cout << "linea " << c << endl;
		TLine* lx = new TLine(5, trackXatZ(qx.at(c),mx.at(c), 912.7), 10, trackXatZ(qx.at(c),mx.at(c), 912.7+89.9218));
                lx->SetLineColor(kRed+c);
                //lx->SetTitle(Form("sig track %d",c));
                lx->Draw("same");
        }
}

//MC TRACKS
/*        for(int n = 0; n < MCTrack->GetEntries(); n++) {
         const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));

         if(MCTr->interactionID()==45 and MCTr->pdgCode()==11) {
                                                                TLine* lxe = new TLine(5, trackXatZ(MCTr->bx(),MCTr->ax(), 912.7), 10, trackXatZ(MCTr->bx(),MCTr->ax(), 912.7+89.9218));
                                                                lxe->SetLineColor(kGreen);
                                                                lxe->Draw("same");}
         if(MCTr->interactionID()==45 and MCTr->pdgCode()==-13) {
                                                                TLine* lxmu = new TLine(5, trackXatZ(MCTr->bx(),MCTr->ax(), 912.7), 10, trackXatZ(MCTr->bx(),MCTr->ax(), 912.7+89.9218));
                                                                lxmu->SetLineColor(kGreen+1);
                                                                lxmu->Draw("same");}
        }
*/
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

gx_s->SetMinimum(-6);
gx_s->SetMaximum(6);
gx_s->SetMarkerColor(kBlue);
gx_s->SetTitle("stubs");

TMultiGraph *mgx_v = new TMultiGraph();
mgx_v->SetMinimum(-6);
mgx_v->SetMaximum(6);
mgx_v->Add(gx_v1,"A*");
mgx_v->Add(gx,"A*");
mgx_v->Add(gx_s,"A*");
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

gy_s->SetMinimum(-6);
gy_s->SetMaximum(6);
gy_s->SetMarkerColor(kBlue);
gy_s->SetTitle("stubs");

TMultiGraph *mg = new TMultiGraph();
mg->SetMinimum(-6);
mg->SetMaximum(6);
mg->Add(gy,"A*");
mg->Add(gy_s,"A*");
mg->Draw("A* ");
mg->SetTitle("Y projection");

if(qy_in.size()!=0 and my_in.size()!=0)
{
        for(int c=0; c<qy_in.size(); c++)
               {if(c==0)
                {TLine* ly = new TLine(0, trackYatZ(qy_in.at(c),my_in.at(c), 811), 5, trackYatZ(qy_in.at(c),my_in.at(c), 912.7));
                 ly->SetLineColor(kBlue+c);
                 ly->Draw("same");
	        cout << "pre vrtx: trackY in Z " << trackYatZ(qy_in.at(c),my_in.at(c),912.7) << " VS " << y << endl;

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
                //ly->SetTitle(Form("sig track %d",c));
                ly->Draw("same");
	}
}


/*        for(int n = 0; n < MCTrack->GetEntries(); n++) {
         const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));

         if(MCTr->interactionID()==45 and MCTr->pdgCode()==11) {
                                                                TLine* lye = new TLine(5, trackYatZ(MCTr->by(),MCTr->ay(), 912.7), 10, trackYatZ(MCTr->by(),MCTr->ay(), 912.7+89.9218));
								lye->SetLineColor(kGreen);
								lye->Draw("same");
								cout <<"Ee " << MCTr->energy() << endl;

 for(int s=0; s<TrackerPoints->GetEntries(); s++){
  const MUonETrackerPoint *TrackerPt = static_cast<const MUonETrackerPoint*>(TrackerPoints->At(s));
  if(TrackerPt->trackID()==n and TrackerPt->stationID()==1){if(TrackerPt->moduleID()==5)cout << "mod 5 " << TrackerPt->exitingPositionGlobalCoordinates().Y();}
 }

								}
         if(MCTr->interactionID()==45 and MCTr->pdgCode()==-13) {
                                                                TLine* lymu = new TLine(5, trackYatZ(MCTr->by(),MCTr->ay(), 912.7), 10, trackYatZ(MCTr->by(),MCTr->ay(), 912.7+89.9218));
                                                                lymu->SetLineColor(kGreen+1);
                                                                lymu->Draw("same");
								}
	}
*/

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

gy_s->SetMinimum(-6);
gy_s->SetMaximum(6);
gy_s->SetMarkerColor(kBlue);
gy_s->SetTitle("stubs");

TMultiGraph *mgy_v = new TMultiGraph();
mgy_v->SetMinimum(-6);
mgy_v->SetMaximum(6);
mgy_v->Add(gy_v1,"A*");
mgy_v->Add(gy,"A*");
mgy_v->Add(gy_s,"A*");
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

a2.cd(5);

guv->SetMinimum(-6);
guv->SetMaximum(6);
guv->SetMarkerColor(kRed);
guv->SetTitle("stubs");

guv_t->SetMinimum(-6);
guv_t->SetMaximum(6);
guv_t->SetMarkerColor(kBlue);
guv_t->SetTitle("track stubs");

for(int c=0; c<guv_loc_t.size();c++){
guv_loc_t.at(c)->SetMinimum(-6);
guv_loc_t.at(c)->SetMaximum(6);
guv_loc_t.at(c)->SetMarkerColor(kGreen+c);
guv_loc_t.at(c)->SetMarkerStyle(29);
guv_loc_t.at(c)->SetTitle("loc track pos");
}

TMultiGraph *mguv = new TMultiGraph();
mguv->SetMinimum(-6);
mguv->SetMaximum(6);
mguv->Add(guv,"A*");
mguv->Add(guv_t,"A*");
for(int c=0; c<guv_loc_t.size();c++) mguv->Add(guv_loc_t.at(c),"A*");
mguv->Draw("A* ");
mguv->SetTitle("UV local position");

gPad->BuildLegend(0.25,0.15,0.25,0.15);

TLine *l0_uv = new TLine(0,-5,0,5);
TLine *l1_uv = new TLine(1,-5,1,5);
TLine *l2_uv = new TLine(2,-5,2,5);
TLine *l3_uv = new TLine(3,-5,3,5);
TLine *l4_uv = new TLine(4,-5,4,5);
TLine *l5_uv = new TLine(5,-5,5,5);
TLine *t2 = new TLine(-1,-5,-1,5);
t2->SetLineWidth(4);
t2->SetLineColor(47);

l0_uv->Draw("same");
l1_uv->Draw("same");
l2_uv->Draw("same");
l3_uv->Draw("same");
l4_uv->Draw("same");
l5_uv->Draw("same");
t2->Draw("same");

mguv->GetXaxis()->SetLimits(-2.,6.);


a2.cd(6);

g_mod_u->SetMinimum(-6);
g_mod_u->SetMaximum(6);
g_mod_u->SetMarkerColor(kRed);
g_mod_u->SetTitle("stubs");

g_mod_u_t->SetMinimum(-6);
g_mod_u_t->SetMaximum(6);
g_mod_u_t->SetMarkerColor(kBlue);
g_mod_u_t->SetTitle("track stubs");

for(int c=0; c<g_mod_u_loc_t.size();c++){
g_mod_u_loc_t.at(c)->SetMinimum(-6);
g_mod_u_loc_t.at(c)->SetMaximum(6);
g_mod_u_loc_t.at(c)->SetMarkerColor(kPink+c);
g_mod_u_loc_t.at(c)->SetMarkerStyle(29);
g_mod_u_loc_t.at(c)->SetTitle("loc track pos");
}


g_mod_v->SetMinimum(-6);
g_mod_v->SetMaximum(6);
g_mod_v->SetMarkerColor(kRed);
g_mod_v->SetTitle("stubs");

g_mod_v_t->SetMinimum(-6);
g_mod_v_t->SetMaximum(6);
g_mod_v_t->SetMarkerColor(kBlue);
g_mod_v_t->SetTitle("track stubs");


for(int c=0; c<g_mod_v_loc_t.size();c++){
g_mod_v_loc_t.at(c)->SetMinimum(-6);
g_mod_v_loc_t.at(c)->SetMaximum(6);
g_mod_v_loc_t.at(c)->SetMarkerColor(kPink+c);
g_mod_v_loc_t.at(c)->SetMarkerStyle(29);
g_mod_v_loc_t.at(c)->SetTitle("loc track pos");
}

TMultiGraph *mguv_mod = new TMultiGraph();
mguv_mod->SetMinimum(-5);
mguv_mod->SetMaximum(5);
mguv_mod->Add(g_mod_u,"A*");
mguv_mod->Add(g_mod_u_t,"A*");
for(int c=0; c<g_mod_u_loc_t.size();c++) mguv_mod->Add(g_mod_u_loc_t.at(c),"A*");
mguv_mod->Add(g_mod_v,"A*");
mguv_mod->Add(g_mod_v_t,"A*");
for(int c=0; c<g_mod_v_loc_t.size();c++) mguv_mod->Add(g_mod_v_loc_t.at(c),"A*");
mguv_mod->Draw("A* ");
mguv_mod->SetTitle("UV local position");

mguv_mod->GetXaxis()->SetLimits(-5.,5);

a2.SaveAs("pdf_exercise.pdf");

}

