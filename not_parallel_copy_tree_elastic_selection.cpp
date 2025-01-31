
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

int not_parallel_copy_tree_elastic_selection(string version, int nhits,string bend,string mc){


TChain * cbmsim = new TChain("cbmsim");
TChain * cbmsim_g = new TChain("cbmsim");

if(version=="default" and nhits==2 and mc=="data"){
cout <<"si 2 hit default"<<endl;
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_1.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_2.root");
}
else if(version=="default" and nhits==0 and bend=="Bend" and mc=="data"){
cout <<"si 0 hit default"<<endl;
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_0hit.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_0hit_1.root");

/*cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_uvWindow045_0hit.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_uvWindow045_0hit_1.root");*/

//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_UVwindow015_0hit.root");
}
else if(version=="default" and nhits==1 and bend=="Bend" and mc=="data"){
cout <<"si 1 hit default"<<endl;
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_1hit.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_1hit_1.root");
}
else if(version=="default" and nhits==0 and bend=="noBend" and mc=="data"){
cout <<"si 0 hit default no bend"<<endl;
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_noBend_0hit.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_noBend_0hit_1.root");
}
else if(version=="ricMis" and nhits==2 and bend=="Bend" and mc=="data"){
cout <<"si 2 hit ricMis"<<endl;
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_commit_b2ed7c3b.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_commit_b2ed7c3b_1.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_commit_b2ed7c3b_2.root");
}
else if(version=="wip12" and nhits==0 and bend=="Bend" and mc=="data"){
cout <<"si 0 hit wip12"<<endl;
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_commit_9e0035f0_0hit.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_commit_9e0035f0_0hit_1.root");
}
else if(version=="wip11" and nhits==0 and bend=="Bend" and mc=="data"){
cout <<"si 0 hit wip11 old skim e old reco"<<endl;

cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_3234-3235_new_1hit_skimGA_inputs_0.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_3234-3235_new_1hit_skimGA_inputs_1.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_3234-3235_new_1hit_skimGA_inputs_2.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_3234-3235_new_1hit_skimGA_inputs_3.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_3234-3235_new_1hit_skimGA_inputs_4.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_3234-3235_new_1hit_skimGA_inputs_5.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_3234-3235_new_1hit_skimGA_inputs_6.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_3234-3235_new_1hit_skimGA_inputs_7.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_3234-3235_new_1hit_skimGA_inputs_8.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_3234-3235_new_1hit_skimGA_inputs_9.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_3234-3235_new_1hit_skimGA_inputs_10.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_3234-3235_new_1hit_skimGA_inputs_11.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_3234-3235_new_1hit_skimGA_inputs_12.root");

}
else if(version=="noskim" and nhits==0 and bend=="Bend" and mc=="data"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_noskim_run6_1.root");
}
else if(version=="oldskim" and nhits==0 and bend=="Bend" and mc=="data"){
cout <<"si 0 hit master pero su old skim"<<endl;
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_run6_decemberSkim_bestConfig_0hit.root");
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_OLDskim_run6.root");
}
else if(version=="oldskim_wip11" and nhits==0 and bend=="Bend" and mc=="data"){
cout <<"si 0 hit wip11 pero su old skim"<<endl;
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_run6_decemberSkim_bestConfig_wip11_0hit_5M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_run6_decemberSkim_bestConfig_wip11_0hit.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_run6_decemberSkim_bestConfig_wip11_0hit_1.root");
}
else if(version=="oldskim_wip11" and nhits==1 and bend=="Bend" and mc=="data"){
cout <<"si 1 hit wip11 pero su old skim"<<endl;
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_run6_decemberSkim_bestConfig_wip11_1hit_5M.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_run6_decemberSkim_bestConfig_wip11_1hit.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_run6_decemberSkim_bestConfig_wip11_1hit_1.root");
}
else if(version=="ricMis_noAlign" and nhits==0 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_commit_ricMis_0hit_noalignment_5M.root");
}
else if(version=="uOpposit" and nhits==0 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_0hit_uOpposit.root");
}
else if(version=="default" and nhits==-1 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_minus1hit.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_minus1hit_1.root");
}
else if(version=="alin" and nhits==0 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_0hit_alin.root");
}
else if(version=="cloneremovaltest" and nhits==2 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_2hit_cloneremovaltest.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_2hit_cloneremovaltest_1.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_2hit_cloneremovaltest_2.root");
}
else if(version=="run8" and nhits==-1 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/commit_6c846c96_MCsignal_RECO_m1hit_run8_partial.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/commit_6c846c96_MCsignal_RECO_m1hit_run8_partial_1.root");
}
else if(version=="default_zTarFix" and nhits==-1 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/commit_6c846c96_MCsignal_RECO_m1hit.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/commit_6c846c96_MCsignal_RECO_m1hit_1.root");
}
else if(version=="run17" and nhits==-1 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/commit_6c846c96_MCsignal_RECO_m1hit_run17_partial.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/commit_6c846c96_MCsignal_RECO_m1hit_run17_partial_1.root");
}
else if(version=="run17" and nhits==2 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/commit_6c846c96_MCsignal_RECO_2hit_run17_partial.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/commit_6c846c96_MCsignal_RECO_2hit_run17_partial_1.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/commit_6c846c96_MCsignal_RECO_2hit_run17_partial_2.root");
}
else if(version=="default_chi2out50" and nhits==2 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/commit_6c846c96_MCsignal_RECO_2hit_run6_chi2out50.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/commit_6c846c96_MCsignal_RECO_2hit_run6_chi2out50_1.root");
}
else if(version=="default_bb35b5de" and nhits==2 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/commit_bb35b5de_RECO_2hit_run6.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/commit_bb35b5de_RECO_2hit_run6_1.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/commit_bb35b5de_RECO_2hit_run6_2.root");
}
else if(version=="default_bb35b5de" and nhits==-1 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/commit_bb35b5de_RECO_m1hit_run6.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/commit_bb35b5de_RECO_m1hit_run6_1.root");
}
else if(version=="default_chi2out50" and nhits==0 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/commit_bb35b5de_RECO_0hit_run6_chi2out50_partial.root");
}
else if(version=="default_bb35b5de" and nhits==0 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/commit_bb35b5de_RECO_0hit_run6.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/commit_bb35b5de_RECO_0hit_run6_1.root");
}
else if(version=="my_modifica" and nhits==0 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/all_run6_nuovo_algo_0hit.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/all_run6_nuovo_algo_0hit_1.root");
}
else if(version=="my_modifica_reassign" and nhits==2 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/all_run6_nuovo_algo_2hit_reassign.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/all_run6_nuovo_algo_2hit_reassign_1.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/all_run6_nuovo_algo_2hit_reassign_2.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/all_run6_nuovo_algo_2hit_reassign_3.root");
}


      MUonERecoOutput *ReconstructionOutput = 0;
        cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);

std::vector<uint64_t> nID_v_0;
std::vector<uint64_t> nID_v_1;
std::vector<uint64_t> nID_v_2;



auto trackatZ = [](double q, double m,double z) {return q + (z ) * m;};
double z_mod[6]={911.2+18.0218,911.2+21.8693,911.2+55.3635,911.2+56.6205,911.2+89.9218,911.2+93.7693};
double alpha[6]={0  *TMath::DegToRad(), 90 *TMath::DegToRad(),135*TMath::DegToRad(), 45 *TMath::DegToRad(),0  *TMath::DegToRad(), 90 *TMath::DegToRad()};




for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
                cbmsim->GetEntry(i);
                if(i%1000000 == 0) cout<<"Entry "<<i<<endl;


double reco_v=0.; double more_reco_v=0.; double reco0_v=0.;
int yes2=0; int yes_v=0;
int code_mu=-99; int code_e=-99; int code_mu_in=-99;
double z_fix=912.7;

MUonERecoOutputVertex vrtx = ReconstructionOutput->bestVertex();

 double chi=vrtx.chi2perDegreeOfFreedom();

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

        std::vector<MUonERecoOutputHit> hits_mu=vrtx.outgoingMuon().hits();
        std::vector<double> pos_mu; pos_mu.resize(6);
        for(int p=0;p<hits_mu.size();p++)pos_mu.push_back(hits_mu.at(p).position());
        std::vector<MUonERecoOutputHit> hits_e=vrtx.outgoingElectron().hits();
        std::vector<double> pos_e; pos_e.resize(6);
        for(int p=0;p<hits_e.size();p++)pos_e.push_back(hits_e.at(p).position());


vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();

         for (auto&& track : tracks) {
        if(track.sector()==0) sec0++;
        if(track.sector()==1) sec1++;
        }

 MUonERecoOutputTrack mu_in_v = vrtx.incomingMuon();
 MUonERecoOutputTrack mu_out_v = vrtx.outgoingMuon();
 MUonERecoOutputTrack e_out_v = vrtx.outgoingElectron();
        TVector3 p_muin_v(mu_in_v.xSlope(),mu_in_v.ySlope(),1.0);
        TVector3 p_mu_v(mu_out_v.xSlope(),mu_out_v.ySlope(),1.0);
        TVector3 p_e_v(e_out_v.xSlope(),e_out_v.ySlope(),1.0);
        p_e_v.Unit();p_mu_v.Unit();p_muin_v.Unit();

MUonERecoOutputTrack t_mu;
MUonERecoOutputTrack t_e;

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
	std::vector<double> pos; pos.resize(6);

        for(int p=0;p<hits_.size();p++){pos.push_back(hits_.at(p).position());
      //if(hits_.at(p).moduleID()==0 or hits_.at(p).moduleID()==4) residualX_pre->Fill(hits_.at(p).position()-trackatZ(track.x0(),track.xSlope(),hits_.at(p).z()));
      //if(hits_.at(p).moduleID()==1 or hits_.at(p).moduleID()==5) residualY_pre->Fill(hits_.at(p).position()-trackatZ(track.y0(),track.ySlope(),hits_.at(p).z()));
	}

         if(std::equal(pos.begin(),pos.end(),pos_mu.begin())){p_mu.SetXYZ(track.xSlope(),track.ySlope(),1.0); p_mu.Unit(); thmu_rec=p_mu.Angle(p_muin);t_mu=track;}
        else if(std::equal(pos.begin(),pos.end(),pos_e.begin())){p_e.SetXYZ(track.xSlope(),track.ySlope(),1.0); p_e.Unit(); the_rec=p_e.Angle(p_muin);t_e=track;}



	yes2++;}
}

double posxIN=99.;//pos_on_track(x0_in,th_inx,z_fix);
double posyIN=99.;//pos_on_track(y0_in,th_iny,z_fix);

int stub0 = 0;
int stub1 = 0;

std::vector<MUonERecoOutputHit> stubs=ReconstructionOutput->reconstructedHits();

         std::array<int,6> module_st1_2;module_st1_2.fill({0});
         std::array<int,6> module_st1;module_st1.fill({0});
        std::array<std::vector<float>,6> localX{{std::vector<float>(0.0)}};
         for (auto&& stub : stubs) {
                if(stub.stationID()==0){stub0++; if(stub.moduleID()==4){posxIN=stub.position(); } else if(stub.moduleID()==5){posyIN=stub.position();}   }
                if(stub.stationID()==1){stub1++; module_st1_2.at(stub.moduleID())+=1; module_st1.at(stub.moduleID())=1;int link= stub.moduleID(); localX.at(link).push_back(stub.position());}
                //if(stub.stationID()==1){stub1++; module_st1_2.at(stub.moduleID())+=1;int link= stub.moduleID(); localX.at(link).push_back(stub.position());}
        }



 if(sec0==1 and stubs_muin==6 and abs(posxIN)<=1.5 and abs(posyIN)<=1.5 and chi2_muin<=2 and stub0==6 and th_muin<0.004){// and stub1<=15){


//bool allmod=std::all_of(std::begin(module_st1), std::end(module_st1), [](int i){return i==1;});
bool allmod2=std::all_of(std::begin(module_st1_2), std::end(module_st1_2), [](int i){return i==2;});
double u=99.;
double v=99.;


  if(allmod2 and stub1==12){
double tilt[6]={0.233,0.233,0.,0.,0.233,0.233};
for(int m=0;m<6;m++){
        double real_pos1= cos(tilt[m])*localX.at(m).at(1);
        double real_pos0= cos(tilt[m])*localX.at(m).at(0);
                            if(m==2){
                             u=real_pos1-real_pos0;
                             }

                            if(m==3){
				v=real_pos1-real_pos0;
                             }
                } //end for

}

if(chi!=0 and the_rec!=-99 and thmu_rec!=-99){


                                                double dotProduct_v = p_mu.Dot(p_e);
                                                TVector3 crossProduct_v = p_mu.Cross(p_e);
                                                double T_v = p_muin.Dot(crossProduct_v);
                                                TVector3 im_v= p_muin.Cross(p_mu);
                                                TVector3 ie_v= p_muin.Cross(p_e);
                                                T_v = T_v>0? 1:-1;
                                                double acoplanarity_v= T_v*(TMath::Pi()- acos( ((im_v).Dot(ie_v))/(im_v.Mag()*ie_v.Mag()) ));


double Elastic=0.5109989461*0.001*((1+(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*cos(the_rec)*cos(the_rec))/(1-(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*cos(the_rec)*cos(the_rec)));
double Elastic2=asin( (sin(the_rec)*sqrt(Elastic*Elastic-0.5109989461*0.001*0.5109989461*0.001))/sqrt( (160+0.5109989461*0.001-Elastic)*(160+0.5109989461*0.001-Elastic)-105.6583745 *0.001*105.6583745 *0.001 ) );
//double Elastic3=asin( (sin(Elastic2)*sqrt( (160+0.5109989461*0.001-Elastic)*(160+0.5109989461*0.001-Elastic)-105.6583745 *0.001*105.6583745 *0.001 ))/sqrt(Elastic*Elastic-0.5109989461*0.001*0.5109989461*0.001) );

bool allmod=std::all_of(std::begin(module_st1), std::end(module_st1), [](int i){return i==1;});

int mod2=0;
double mod2_mu=99.;

//  if(allmod and abs(acoplanarity_v)<=0.4 and chi<20 and thmu_rec>0.0002 and stub1<=20 and thmu_rec<=Elastic2+0.0002 and thmu_rec>=Elastic2-0.0002 and the_rec<0.025 and the_rec>=0.003 and thmu_rec<=0.003){// and vrtx.zPositionFit()<917 and vrtx.zPositionFit()>907){

// if(allmod2 and chi<20 and stub1==12 and thmu_rec>0.0002 and thmu_rec<=0.002 and the_rec<0.030 and the_rec>=0.003 and thmu_rec<=Elastic2+0.0002 and thmu_rec>=Elastic2-0.0002 and p_mu.Angle(p_e)>0.005 and abs(acoplanarity_v)<=0.4){
 if(chi>0 and stub1<=15 and thmu_rec>=0.0002 and the_rec<0.020 and the_rec>=0.005 and thmu_rec<=Elastic2+0.0002 and thmu_rec>=Elastic2-0.0002 and p_mu.Angle(p_e)>0.005 and abs(acoplanarity_v)<=0.4){// and zpos<(911.2+3) and zpos>(911.2-3)){
 if(yes2>=2){

//cout << "Event " << i << endl;

if(vrtx.zPositionFit()<(911.2+3) and vrtx.zPositionFit()>(911.2-3)){}
else{cout << "Event " << i << endl;}

std::array<double,6> t_e_hits;t_e_hits.fill({-99.});
 std::array<double,6> t_mu_hits;t_mu_hits.fill({-99.});

 for(int h=0; h<t_mu.hits().size(); h++){
  if(t_mu.hits().at(h).stationID()==1) {t_mu_hits.at(t_mu.hits().at(h).moduleID())=t_mu.hits().at(h).position();}
 }

 for(int h=0; h<t_e.hits().size(); h++){
  if(t_e.hits().at(h).stationID()==1) { t_e_hits.at(t_e.hits().at(h).moduleID())=t_e.hits().at(h).position();}
 }


if( the_rec<=0.005 ){uint64_t nID = static_cast<uint64_t>(t_e.hits().at(0).superID())*3564 + static_cast<uint64_t>(t_e.hits().at(0).bx()); nID_v_0.push_back(nID);}

if( the_rec>0.005 and the_rec<=0.01 ){uint64_t nID = static_cast<uint64_t>(t_e.hits().at(0).superID())*3564 + static_cast<uint64_t>(t_e.hits().at(0).bx()); nID_v_1.push_back(nID);}

if( the_rec>0.010 and the_rec<=0.015 ){uint64_t nID = static_cast<uint64_t>(t_e.hits().at(0).superID())*3564 + static_cast<uint64_t>(t_e.hits().at(0).bx()); nID_v_2.push_back(nID);}


				}//yes2>=2
			}//aco e chi
		}//chi!=0
	}//mu_in
yes2=0;
}; //end of ySlopefunction


   TFile *fout = TFile::Open(Form("/mnt/raid10/DATA/espedica/fairmu/reco/save_samples/%s_list_selected_events_12stubs_%s_%dhit_%s.root",version.c_str(),bend.c_str(),nhits,mc.c_str()), "RECREATE");

   fout->WriteObject(&nID_v_0, "nID_v_the5");
   fout->WriteObject(&nID_v_1, "nID_v_5the10");
   fout->WriteObject(&nID_v_2, "nID_v_10the15");
   fout->Close();

return 0;
}
