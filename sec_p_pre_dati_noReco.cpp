#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TSystemDirectory.h"
#include <TStyle.h>


using namespace std;

                auto pos_on_track = [](double q, double m, double z){return (q + m*z);};

int sec_p_pre_dati_noReco(string version, int nhits,string bend,string mc){

  int nthreads = 6;
  ROOT::EnableImplicitMT(nthreads);


TChain * cbmsim = new TChain("cbmsim");
TChain * cbmsim_g = new TChain("cbmsim");

if(version=="default" and nhits==2){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_1.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_2.root");
}
else if(version=="default" and nhits==0 and bend=="Bend" and mc=="data"){
cout <<"si"<<endl;
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_0hit.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_0hit_1.root");

/*cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_uvWindow045_0hit.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_uvWindow045_0hit_1.root");*/

//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_UVwindow015_0hit.root");
}
else if(version=="default" and nhits==0 and bend=="noBend" and mc=="data"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_noBend_0hit.root");
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_noBend_0hit_1.root");
}
else if(version=="ricMis" and nhits==2 and bend=="Bend" and mc=="data"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_commit_b2ed7c3b.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_commit_b2ed7c3b_1.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_commit_b2ed7c3b_2.root");
}
else if(version=="ricMis" and nhits==0 and bend=="Bend" and mc=="data"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_commit_b2ed7c3b_0hit.root");}
else if(version=="default" and nhits==0 and bend=="Bend" and mc=="mc"){
          cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_commit_2f4e96f4_MCsignal_bestConfig_0hit.root");
          cbmsim_g->Add("/mnt/raid10/DATA/espedica/fairmu/gen_digi/commit_2f4e96f4_MCsignal_SIM-DIGI.root");
          cbmsim->AddFriend(cbmsim_g);}
else if(version=="ricMis" and nhits==0 and bend=="Bend" and mc=="mc"){
          cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/commit_b2ed7c3b_MCsignal_bestConfig.root");
          cbmsim_g->Add("/mnt/raid10/DATA/espedica/fairmu/gen_digi/commit_b2ed7c3b_MCsignal_SIM-DIGI.root");
          cbmsim->AddFriend(cbmsim_g);}


ROOT::TTreeProcessorMT tp1(*cbmsim,nthreads);

      MUonERecoOutput *ReconstructionOutput = 0;



   ROOT::TThreadedObject<TH1D> residualU_critic("residualU_critic","residual U when track distance <0.2cm and track_e no Ustub",200,-1.,1.);
   ROOT::TThreadedObject<TH1D> residualU_no_critic("residualU_no_critic","residual U when track distance <0.2cm",400,-2.,2.);

   ROOT::TThreadedObject<TH1D> residualU("residualU","residual U when track distance >0.2cm",200,-1.,1.);
   ROOT::TThreadedObject<TH1D> residualU_no("residualU_no","residual U when track distance >0.2cm and track_e no Ustub",400,-2.,2.);

   ROOT::TThreadedObject<TH1D> residual0("residual0","residual 0",200,-1.,1.);
   ROOT::TThreadedObject<TH1D> residual1("residual1","residual 1",200,-1.,1.);
   ROOT::TThreadedObject<TH1D> residual2("residual2","residual 2",200,-1.,1.);
   ROOT::TThreadedObject<TH1D> residual3("residual3","residual 3",200,-1.,1.);
   ROOT::TThreadedObject<TH1D> residual4("residual4","residual 4",200,-1.,1.);
   ROOT::TThreadedObject<TH1D> residual5("residual5","residual 5",200,-1.,1.);


   ROOT::TThreadedObject<TH1D> h_dist0("h_dist0", "distance stub electron - stub muon in cm mod0",200,-1.,1.);
   ROOT::TThreadedObject<TH1D> h_dist1("h_dist1", "distance stub electron - stub muon in cm mod1",200,-1.,1.);
   ROOT::TThreadedObject<TH1D> h_dist2("h_dist2", "distance stub electron - stub muon in cm mod2",200,-1.,1.);
   ROOT::TThreadedObject<TH1D> h_dist3("h_dist3", "distance stub electron - stub muon in cm mod3",200,-1.,1.);
   ROOT::TThreadedObject<TH1D> h_dist4("h_dist4", "distance stub electron - stub muon in cm mod4",200,-1.,1.);
   ROOT::TThreadedObject<TH1D> h_dist5("h_dist5", "distance stub electron - stub muon in cm mod5",200,-1.,1.);

   ROOT::TThreadedObject<TH1D> h_hist_distance0("h_hist_distance0","Distance between stub1 and stub0 in module0",200,-1.,1.);
   ROOT::TThreadedObject<TH1D> h_hist_distance_zoom0("h_hist_distance_zoom0","Distance between stub0 and stub1 in module0 zoomed",80,-0.4,0.4);
   ROOT::TThreadedObject<TH1D> h_hist_distance2("h_hist_distance2","Distance between stub1 and stub0 in module",200,-1.,1.);
   ROOT::TThreadedObject<TH1D> h_hist_distance_zoom2("h_hist_distance_zoom2","Distance between stub0 and stub1 in module2 zoomed",80,-0.4,0.4);

   ROOT::TThreadedObject<TH1D> h_nstubs_track_e("h_nstubs_track_e","nstubs track e",8,0,8);
   ROOT::TThreadedObject<TH1D> h_nstubs_track_mu("h_nstubs_track_mu","nstubs track mu",8,0,8);
   ROOT::TThreadedObject<TH1D> reco_trackel("reco_trackel","# recco track el", 6,0,6);
   ROOT::TThreadedObject<TH1D> reco_trackmu("reco_trackmu","# recco track mu", 6,0,6);

auto trackatZ = [](double q, double m,double z) {return q + (z ) * m;};
double z_mod[6]={911.2+18.0218,911.2+21.8693,911.2+55.3635,911.2+56.6205,911.2+89.9218,911.2+93.7693};
double alpha[6]={0  *TMath::DegToRad(), 90 *TMath::DegToRad(), 135*TMath::DegToRad(), 45 *TMath::DegToRad(),0  *TMath::DegToRad(), 90 *TMath::DegToRad()};

 auto myFunction = [&](TTreeReader &myReader) {

     TTreeReaderValue<std::vector<MUonERecoOutputTrack>> RVtracks(myReader, "ReconstructedTracks");
     TTreeReaderValue<MUonERecoOutputVertex> vrtx(myReader, "BestVertex");
     TTreeReaderValue<std::vector<MUonERecoOutputHit>> RVstubs(myReader, "ReconstructedHits");

     while (myReader.Next()) {

double reco_v=0.; double more_reco_v=0.; double reco0_v=0.;
int yes2=0; int yes_v=0;
int code_mu=-99; int code_e=-99; int code_mu_in=-99;
double z_fix=912.7;


 double chi=vrtx->chi2perDegreeOfFreedom();

int yes_mu=0;
int yes_e=0;
double th_inx,th_iny,x0_in,y0_in;
double chi2_muin;
int sec0=0;
int sec1=0;
int stubs_muin=0.;
double th_muin=0.;


double the_rec=-99.;
double thmu_rec=-99.;
TVector3 p_e,p_mu,p_muin;

        std::vector<MUonERecoOutputHit> hits_mu=vrtx->outgoingMuon().hits();
        std::vector<double> pos_mu; pos_mu.resize(6);
        for(int p=0;p<hits_mu.size();p++)pos_mu.push_back(hits_mu.at(p).position());
        std::vector<MUonERecoOutputHit> hits_e=vrtx->outgoingElectron().hits();
        std::vector<double> pos_e; pos_e.resize(6);
        for(int p=0;p<hits_e.size();p++)pos_e.push_back(hits_e.at(p).position());


         auto tracks = *RVtracks;
         for (auto&& track : tracks) {
        if(track.sector()==0) sec0++;
        if(track.sector()==1) sec1++;
        }

 MUonERecoOutputTrack mu_in_v = vrtx->incomingMuon();
 MUonERecoOutputTrack mu_out_v = vrtx->outgoingMuon();
 MUonERecoOutputTrack e_out_v = vrtx->outgoingElectron();
        TVector3 p_muin_v(mu_in_v.xSlope(),mu_in_v.ySlope(),1.0);
        TVector3 p_mu_v(mu_out_v.xSlope(),mu_out_v.ySlope(),1.0);
        TVector3 p_e_v(e_out_v.xSlope(),e_out_v.ySlope(),1.0);
        p_e_v.Unit();p_mu_v.Unit();p_muin_v.Unit();

MUonERecoOutputTrack t_mu;
MUonERecoOutputTrack t_e;
int  size1=0;
        std::array<std::vector<float>,6> posT{{std::vector<float>(0.0)}};

        std::array<std::vector<float>,6> posTZ{{std::vector<float>(0.0)}};


         for (auto&& track : tracks) {
        std::vector<MUonERecoOutputHit> hits_=track.hits();

        if(track.sector()==0 and sec0==1){
	stubs_muin=hits_.size();
        th_inx=track.xSlope();
        th_iny=track.ySlope();
        x0_in=track.x0();
        y0_in=track.y0();
        chi2_muin=track.chi2perDegreeOfFreedom();
        p_muin.SetXYZ(th_inx,th_iny,1.0);
        p_muin=p_muin.Unit();
        th_muin=p_muin.Theta();
                        }
	if(track.sector()==1)
	{

        std::array<std::vector<float>,6> pos{{std::vector<float>(0.0)}};

        for(int p=0;p<hits_.size();p++){pos.at(hits_.at(p).moduleID()).push_back(hits_.at(p).position());posTZ.at(hits_.at(p).moduleID()).push_back(hits_.at(p).z());}
	size1=hits_.size();

         if(mu_out_v.index()==track.index()){p_mu.SetXYZ(track.xSlope(),track.ySlope(),1.0); p_mu.Unit(); thmu_rec=p_mu.Angle(p_muin);t_mu=track;}
        else if(e_out_v.index()==track.index()){p_e.SetXYZ(track.xSlope(),track.ySlope(),1.0); p_e.Unit(); the_rec=p_e.Angle(p_muin);t_e=track;}

	posT=pos;

	yes2++;}
}

double posxIN=99.;
double posyIN=99.;

int stub0 = 0;
int stub1 = 0;

         auto stubs = *RVstubs;
         std::array<int,6> module_st1;module_st1.fill({0});
        std::array<std::vector<float>,6> localX{{std::vector<float>(0.0)}};
        std::array<std::vector<float>,6> localZ{{std::vector<float>(0.0)}};
         for (auto&& stub : stubs) {
                if(stub.stationID()==0){stub0++; if(stub.moduleID()==4){posxIN=stub.position(); } else if(stub.moduleID()==5){posyIN=stub.position();}   }
                if(stub.stationID()==1){stub1++; module_st1.at(stub.moduleID())+=1;int link= stub.moduleID(); localX.at(link).push_back(stub.position());
					localZ.at(link).push_back(stub.z());}
        }


 if(sec0==1 and stubs_muin==6 and abs(posxIN)<=1.5 and abs(posyIN)<=1.5 and chi2_muin<=2 and stub0==6 and th_muin<0.004){// and stub1<=15){

bool allmod=std::all_of(std::begin(module_st1), std::end(module_st1), [](int i){return i==2;});
double u=99.;


if(allmod and stub1==12 and chi==0 and size1==6){
double tilt[6]={0.233,0.233,0.,0.,0.233,0.233};
        std::array<std::vector<float>,6> not_used{{std::vector<float>(0.0)}};
        std::array<std::vector<float>,6> not_usedZ{{std::vector<float>(0.0)}};

for(int m=0;m<6;m++){
        double real_pos1= cos(tilt[m])*localX.at(m).at(1);
        double real_pos0= cos(tilt[m])*localX.at(m).at(0);
                            if(m==0){ h_hist_distance0->Fill(real_pos1-real_pos0);
                             	      h_hist_distance_zoom0->Fill(real_pos1-real_pos0);
                             }

                            if(m==2){ h_hist_distance2->Fill(real_pos1-real_pos0);
                                      //h_hist_distance_zoom2->Fill(real_pos1-real_pos0);
					if(abs(real_pos1-real_pos0)<0.2 and the_rec!=99)reco_trackel->Fill(sec1);
					if(abs(real_pos1-real_pos0)<0.2 and thmu_rec!=99)reco_trackmu->Fill(sec1);
					u=real_pos1-real_pos0;
	                             }
                } //end for

int aiut=0;
for(int p=0;p<6;p++){
/*cout << "0 " << localX.at(0).size() << endl;
cout << "1 " << localX.at(1).size() << endl;
cout << "4 " << localX.at(4).size() << endl;
cout << "5 " << localX.at(5).size() << endl;

cout << "localX.at(p).at(0) " << localX.at(p).at(0) <<endl;
cout<< "posT.at(p).at(0) " << posT.at(p).at(0) << endl;
cout << "localX.at(p).at(1) " << localX.at(p).at(1) <<endl;
*/


if(localX.at(p).at(0)!=posT.at(p).at(0)){not_used.at(p).push_back(localX.at(p).at(0));not_usedZ.at(p).push_back(localZ.at(p).at(0));}
else if(localX.at(p).at(1)!=posT.at(p).at(0)){not_used.at(p).push_back(localX.at(p).at(1));not_usedZ.at(p).push_back(localZ.at(p).at(1));}
else{aiut=1;}// cout << "aiuto"<<endl;}
}

if(aiut==0){
/*cout << "0 " << not_used.at(0).size() << endl;
cout << "1 " << not_used.at(1).size() << endl;
cout << "4 " << not_used.at(4).size() << endl;
cout << "5 " << not_used.at(5).size() << endl;
*/

  double mx=( posT.at(4).at(0)-posT.at(0).at(0) )/( posTZ.at(4).at(0)-posTZ.at(0).at(0) );
  double my=( posT.at(5).at(0)-posT.at(1).at(0) )/( posTZ.at(5).at(0)-posTZ.at(1).at(0) );
  double qx=( posTZ.at(4).at(0)*posT.at(0).at(0) - posTZ.at(0).at(0)*posT.at(4).at(0) )/( posTZ.at(4).at(0)-posTZ.at(0).at(0) );
  double qy=( posTZ.at(5).at(0)*posT.at(1).at(0) - posTZ.at(1).at(0)*posT.at(5).at(0) )/( posTZ.at(5).at(0)-posTZ.at(1).at(0) );

          double tracklocxy2 = (qx + mx*posTZ.at(2).at(0))*cos(alpha[2]) +
                               (qy + my*posTZ.at(2).at(0))*sin(alpha[2]);

          double tracklocxy3 = (qx + mx*posTZ.at(3).at(0))*cos(alpha[3]) +
                               (qy + my*posTZ.at(3).at(0))*sin(alpha[3]);


 residual0->Fill(posT.at(0).at(0)- pos_on_track(qx,mx,posTZ.at(0).at(0)) );
 residual1->Fill(posT.at(1).at(0)- pos_on_track(qy,my,posTZ.at(1).at(0)) );
 residual2->Fill(posT.at(2).at(0)- tracklocxy2 );
 residual3->Fill(posT.at(3).at(0)- tracklocxy3 );
 residual4->Fill(posT.at(4).at(0)- pos_on_track(qx,mx,posTZ.at(4).at(0)) );
 residual5->Fill(posT.at(5).at(0)- pos_on_track(qy,my,posTZ.at(5).at(0)) );

/*
  double mx=( not_used.at(4).at(0)-not_used.at(0).at(0) )/( not_usedZ.at(4).at(0)-not_usedZ.at(0).at(0) );
  double my=( not_used.at(5).at(0)-not_used.at(1).at(0) )/( not_usedZ.at(5).at(0)-not_usedZ.at(1).at(0) );
  double qx=( not_usedZ.at(4).at(0)*not_used.at(0).at(0) - not_usedZ.at(0).at(0)*not_used.at(4).at(0) )/( not_usedZ.at(4).at(0)-not_usedZ.at(0).at(0) );
  double qy=( not_usedZ.at(5).at(0)*not_used.at(1).at(0) - not_usedZ.at(1).at(0)*not_used.at(5).at(0) )/( not_usedZ.at(5).at(0)-not_usedZ.at(1).at(0) );

          double tracklocxy2 = (qx + mx*not_usedZ.at(2).at(0))*cos(alpha[2]) +
                               (qy + my*not_usedZ.at(2).at(0))*sin(alpha[2]);

          double tracklocxy3 = (qx + mx*not_usedZ.at(3).at(0))*cos(alpha[3]) +
                               (qy + my*not_usedZ.at(3).at(0))*sin(alpha[3]);


 residual0->Fill(not_used.at(0).at(0)- pos_on_track(qx,mx,not_usedZ.at(0).at(0)) );
 residual1->Fill(not_used.at(1).at(0)- pos_on_track(qy,my,not_usedZ.at(1).at(0)) );
 residual2->Fill(not_used.at(2).at(0)- tracklocxy2 );
 residual3->Fill(not_used.at(3).at(0)- tracklocxy3 );
 residual4->Fill(not_used.at(4).at(0)- pos_on_track(qx,mx,not_usedZ.at(4).at(0)) );
 residual5->Fill(not_used.at(5).at(0)- pos_on_track(qy,my,not_usedZ.at(5).at(0)) );
*/
h_hist_distance_zoom2->Fill(u);
 }


}//if chi==0
//else if(size1!=6){cout << "size1 " << size1 <<endl;}
if(chi!=0 and the_rec!=-99 and thmu_rec!=-99){


                                                double dotProduct_v = p_mu.Dot(p_e);
                                                TVector3 crossProduct_v = p_mu.Cross(p_e);
                                                double T_v = p_muin.Dot(crossProduct_v);
                                                TVector3 im_v= p_muin.Cross(p_mu);
                                                TVector3 ie_v= p_muin.Cross(p_e);
                                                T_v = T_v>0? 1:-1;
                                                double acoplanarity_v= T_v*(TMath::Pi()- acos( ((im_v).Dot(ie_v))/(im_v.Mag()*ie_v.Mag()) ));

int mod2=0;
double mod2_mu=99.;
for(int h=0;h<t_mu.hits().size(); h++){if(t_mu.hits().at(h).moduleID()==2)mod2_mu=t_mu.hits().at(h).position();}

if(allmod and stub1==12){//and thmu_rec>0.0005 and thmu_rec<=0.002 and the_rec<0.010 and the_rec>=0.003 and thmu_rec<=Elastic2+0.0002 and thmu_rec>=Elastic2-0.0002 and p_mu.Angle(p_e)>0.005){

 if(yes2>=2){


          double tracklocxy0 = (t_e.x0() + t_e.xSlope()*z_mod[2])*cos(alpha[2]) +
                               (t_e.y0() + t_e.ySlope()*z_mod[2])*sin(alpha[2]);



std::array<double,6> t_e_hits; t_e_hits.fill({-99.});
 std::array<double,6> t_mu_hits; t_mu_hits.fill({-99.});

 for(int h=0; h<t_mu.hits().size(); h++){
  if(t_mu.hits().at(h).stationID()==1) t_mu_hits.at(t_mu.hits().at(h).moduleID())=t_mu.hits().at(h).position();
 }

 for(int h=0; h<t_e.hits().size(); h++){
  if(t_e.hits().at(h).stationID()==1) t_e_hits.at(t_e.hits().at(h).moduleID())=t_e.hits().at(h).position();
 }


//if(abs(u)<0.2 ){
        h_nstubs_track_e->Fill(t_e.hits().size());
        h_nstubs_track_mu->Fill(t_mu.hits().size());
for(int m=0;m<hits_e.size();m++){
	if(t_e_hits.at(2)!=-99. and hits_e.at(m).moduleID()==2) residualU_critic->Fill(hits_e.at(m).perpendicularResiduum());}
//	if(t_e_hits.at(2)!=-99.) residualU_critic->Fill(tracklocxy0-t_e_hits.at(2));
        if(t_e_hits.at(2)==-99. and abs(tracklocxy0-localX.at(2).at(0)) < abs(tracklocxy0-localX.at(2).at(1))) residualU_no_critic->Fill(tracklocxy0-localX.at(2).at(0));
        else if(t_e_hits.at(2)==-99. and abs(tracklocxy0-localX.at(2).at(1)) < abs(tracklocxy0-localX.at(2).at(0))) residualU_no_critic->Fill(tracklocxy0-localX.at(2).at(1));
/* if(t_e.hits().size()==6 and t_e_hits.at(2)!=-99.) residualU_critic->Fill(tracklocxy0-t_e_hits.at(2));
 else if( t_e.hits().size()==5 and t_e_hits.at(2)==-99. and t_mu_hits.at(2)!=localX.at(2).at(0) ) residualU_no_critic->Fill(tracklocxy0-localX.at(2).at(0));
 else if( t_e.hits().size()==5 and t_e_hits.at(2)==-99 and t_mu_hits.at(2)!=localX.at(2).at(1) ) residualU_no_critic->Fill(tracklocxy0-localX.at(2).at(1));*/

//}

//else if(abs(u)>0.2 ){
        if(t_e_hits.at(2)!=-99.) residualU->Fill(tracklocxy0-t_e_hits.at(2));
	if(t_e_hits.at(2)==-99. and abs(tracklocxy0-localX.at(2).at(0)) < abs(tracklocxy0-localX.at(2).at(1))) residualU_no->Fill(tracklocxy0-localX.at(2).at(0));
	else if(t_e_hits.at(2)==-99. and abs(tracklocxy0-localX.at(2).at(1)) < abs(tracklocxy0-localX.at(2).at(0))) residualU_no->Fill(tracklocxy0-localX.at(2).at(1));

 /*if(t_e.hits().size()==6 and t_e_hits.at(2)!=-99.){residualU->Fill(tracklocxy0-t_e_hits.at(2));}
 else if(t_e.hits().size()==5 and t_e_hits.at(2)==-99. and t_mu_hits.at(2)!=localX.at(2).at(0))residualU_no->Fill(tracklocxy0-localX.at(2).at(0));
 else if(t_e.hits().size()==5  and t_e_hits.at(2)==-99. and t_mu_hits.at(2)!=localX.at(2).at(1))residualU_no->Fill(tracklocxy0-localX.at(2).at(1));*/
//}




 if( t_e_hits.at(0)!=-99. and t_mu_hits.at(0)!=-99.) h_dist0->Fill(t_e_hits.at(0)-t_mu_hits.at(0));
 if( t_e_hits.at(1)!=-99. and t_mu_hits.at(1)!=-99.) h_dist1->Fill(t_e_hits.at(1)-t_mu_hits.at(1));
 if( t_e_hits.at(2)!=-99. and t_mu_hits.at(2)!=-99.) h_dist2->Fill(t_e_hits.at(2)-t_mu_hits.at(2));
 if( t_e_hits.at(3)!=-99. and t_mu_hits.at(3)!=-99.) h_dist3->Fill(t_e_hits.at(3)-t_mu_hits.at(3));
 if( t_e_hits.at(4)!=-99. and t_mu_hits.at(4)!=-99.) h_dist4->Fill(t_e_hits.at(4)-t_mu_hits.at(4));
 if( t_e_hits.at(5)!=-99. and t_mu_hits.at(5)!=-99.) h_dist5->Fill(t_e_hits.at(5)-t_mu_hits.at(5));

				}//yes2>=2
			}//aco e chi
		}//chi!=0
	}//mu_in
yes2=0;
 } //end of general while
}; //end of ySlopefunction

  tp1.Process(myFunction);

auto residual0M = residual0.Merge();
auto residual1M = residual1.Merge();
auto residual2M = residual2.Merge();
auto residual3M = residual3.Merge();
auto residual4M = residual4.Merge();
auto residual5M = residual5.Merge();

auto h_hist_distance0M = h_hist_distance0.Merge();
auto h_hist_distance_zoom0M = h_hist_distance_zoom0.Merge();
auto h_hist_distance2M = h_hist_distance2.Merge();
auto h_hist_distance_zoom2M = h_hist_distance_zoom2.Merge();

auto h_nstubs_track_eM=h_nstubs_track_e.Merge();
auto h_nstubs_track_muM=h_nstubs_track_mu.Merge();

  auto h_dist0M=h_dist0.Merge();
  auto h_dist1M=h_dist1.Merge();
  auto h_dist2M=h_dist2.Merge();
  auto h_dist3M=h_dist3.Merge();
  auto h_dist4M=h_dist4.Merge();
  auto h_dist5M=h_dist5.Merge();

  auto residualU_criticM = residualU_critic.Merge();
  auto residualU_no_criticM =residualU_no_critic.Merge();
  auto residualUM_post = residualU.Merge();
  auto residualU_noM_post =residualU_no.Merge();

  auto reco_trackelM=reco_trackel.Merge();
  auto reco_trackmuM=reco_trackmu.Merge();

TCanvas r("r","r",1400,2100);
r.Divide(2,3);
r.cd(1);
residual0M->Draw("hist");
r.cd(2);
residual1M->Draw("hist");
r.cd(3);
residual2M->Draw("hist");
r.cd(4);
residual3M->Draw("hist");
r.cd(5);
residual4M->Draw("hist");
r.cd(6);
residual5M->Draw("hist");
r.SaveAs(Form("comparison_RDMC/%s_residuals_track_%s_%dhit_%s.pdf",version.c_str(),bend.c_str(),nhits,mc.c_str()));


TCanvas n1("n1","n1",1000,1000);
n1.Divide(2,2);
n1.cd(1);
h_hist_distance0M->Draw("hist");
n1.cd(2);
h_hist_distance2M->Draw("hist");
n1.cd(3);
//h_hist_distance_zoom0M->Draw("hist");
reco_trackmuM->Draw("hist");
//reco_trackelM->SetLineColor(kRed);
//reco_trackelM->Draw("hist same");
n1.cd(4);
h_hist_distance_zoom2M->Draw("hist");
n1.SaveAs(Form("comparison_RDMC/%s_distance_stubs_zoom_absU_chi0_%s_%dhit_%s.pdf",version.c_str(),bend.c_str(),nhits,mc.c_str()));



TCanvas t("t","t",700,700);
t.Divide(2,2);
t.cd(1);
residualU_criticM->Draw("hist");
t.cd(2);
//residualU_no_criticM->Draw("hist");
t.cd(3);
residualUM_post->Draw("hist");
t.cd(4);
//residualU_noM_post->Draw("hist");
//h_2d_preMerged->Draw();
t.SaveAs(Form("comparison_RDMC/%s_RD_muin_absU_chi0_%s_%dhit_%s.pdf",version.c_str(),bend.c_str(),nhits,mc.c_str()));
/*
TCanvas d("d","d",1400,2100);
d.Divide(2,3);
d.cd(1);
h_dist0M->Draw("hist");
h_dist0M->SaveAs(Form("comparison_RDMC/%s_distance_MOD0_absU_chi0_%s_%dhit.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));
d.cd(2);
h_dist1M->Draw("hist");
h_dist1M->SaveAs(Form("comparison_RDMC/%s_distance_MOD1_absU_chi0_%s_%dhit.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));
d.cd(3);
h_dist2M->Draw("hist");
h_dist2M->SaveAs(Form("comparison_RDMC/%s_distance_MOD2_absU_chi0_%s_%dhit.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));
d.cd(4);
h_dist3M->Draw("hist");
h_dist3M->SaveAs(Form("comparison_RDMC/%s_distance_MOD3_absU_chi0_%s_%dhit.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));
d.cd(5);
h_dist4M->Draw("hist");
h_dist4M->SaveAs(Form("comparison_RDMC/%s_distance_MOD4_absU_chi0_%s_%dhit.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));
d.cd(6);
h_dist5M->Draw("hist");
h_dist5M->SaveAs(Form("comparison_RDMC/%s_distance_MOD5_absU_chi0_%s_%dhit.root",version.c_str(),bend.c_str(),nhits,mc.c_str()));

d.SaveAs(Form("comparison_RDMC/%s_distance_absU_chi0_%s_%dhit_%s.pdf",version.c_str(),bend.c_str(),nhits,mc.c_str()));
*/

TCanvas ns("ns","ns",700,700);
ns.Divide(1,2);
ns.cd(1);
h_nstubs_track_eM->Draw("hist");
ns.cd(2);
h_nstubs_track_muM->Draw("hist");
ns.SaveAs(Form("comparison_RDMC/%s_nstubs_track_absU_chi0_%s_%dhit_%s.pdf",version.c_str(),bend.c_str(),nhits,mc.c_str()));

return 0;
}

