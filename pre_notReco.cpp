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

void pre_notReco(){

        //TFile *inputfile = new TFile("/mnt/raid10/DATA/espedica/fairmu/Mesmer_new_1M_1hit_bend.root");


//TChain * cbmsim = new TChain("cbmsim");
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/gen_digi_reco/WiP_v0140_commit_258faf6b_0_infmrad_1M.root");
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/bkg/mesmer_bkg_2hit_1M.root");


TChain * cbmsim = new TChain("cbmsim");
TChain * cbmsim_g = new TChain("cbmsim");
TString gen_filename;
 TTree *t2;
          cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_commit_2f4e96f4_MCsignal_bestConfig_0hit.root");
          cbmsim_g->Add("/mnt/raid10/DATA/espedica/fairmu/gen_digi/commit_2f4e96f4_MCsignal_SIM-DIGI.root");

          cbmsim->AddFriend(cbmsim_g);


/*TChain * gen_filename = new TChain("cbmsim");

//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_commit_b892679b_0_5mrad_2hitFirstModules.root");
//gen_filename->Add("/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_b892679b_0_5mrad.root",10000);
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_commit_b892679b_5_10mrad_2hitFirstModules.root",10000);
gen_filename->Add("/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_b892679b_5_10mrad.root",10000);
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_commit_b892679b_10_15mrad_2hitFirstModules.root",10000);
gen_filename->Add("/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_b892679b_10_15mrad.root",10000);
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_commit_b892679b_15_20mrad_2hitFirstModules.root",10000);
gen_filename->Add("/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_b892679b_15_20mrad.root",10000);
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_commit_b892679b_20_25mrad_2hitFirstModules.root",10000);
gen_filename->Add("/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_b892679b_20_25mrad.root",10000);
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/reco/WiP_v0140_commit_b892679b_25_32mrad_2hitFirstModules.root",10000);
gen_filename->Add("/mnt/raid10/DATA/espedica/fairmu/gen_digi/WiP_v0140_commit_b892679b_25_32mrad.root",10000);

cout << "entries cbmsim " << cbmsim->GetEntries() << endl;
cout << "entries gen_filename " << gen_filename->GetEntries() << endl;

cbmsim->AddFriend(gen_filename);
*/



//cbmsim->Add("/home/espedica/fair_install/provaFairRoot/share/MUonE/macros/prova_0_32.root");
/*cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_0-5mrad_1M_2hitFirstModules_NOoutchi2_reassign.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_5-10mrad_1M_2hitFirstModules_NOoutchi2_reassign.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_10-15mrad_1M_2hitFirstModules_NOoutchi2_reassign.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_15-20mrad_1M_2hitFirstModules_NOoutchi2_reassign.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_20-25mrad_1M_2hitFirstModules_NOoutchi2_reassign.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_25-32mrad_1M_2hitFirstModules_NOoutchi2_reassign.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/theta_32-inf_mrad_1M_2hitFirstModules_NOoutchi2_reassign.root");
*/

/*
//   TFile *f1 = new TFile("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_0-inf_mrad_2hitFirstModules_NOoutchi2_MCcorrections_reco.root");
   TFile *f1 = new TFile("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/theta_0-inf_mrad_2hitFirstModules_NOoutchi2_MCcorrections_reco_nofirst.root");
   TFile *f2 = new TFile("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/old_forRDMC_2hitFirstModules.root");

   TTree *cbmsim = (TTree*)f1->Get("cbmsim");
   TTree *t2 = (TTree*)f2->Get("cbmsim");

cout << "entries cbmsim " << cbmsim->GetEntries() << endl;
cout << "entries t2 " << t2->GetEntries() << endl;


   cbmsim->AddFriend(t2);
   cout << t2->GetEntries() << endl;
   cout << cbmsim->GetEntries() << endl;
   cbmsim->Print();
*/

/*
TChain * cbmsim = new TChain("cbmsim");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/single_sample_2hitFirstModules_NOoutchi2_MCcorrections.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/single_sample_2hitFirstModules_NOoutchi2_MCcorrections1.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/single_sample_2hitFirstModules_NOoutchi2_MCcorrections2.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/single_sample_2hitFirstModules_NOoutchi2_MCcorrections3.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/single_sample_2hitFirstModules_NOoutchi2_MCcorrections4.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/single_sample_2hitFirstModules_NOoutchi2_MCcorrections5.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/efficiency_NLO/single_sample_2hitFirstModules_NOoutchi2_MCcorrections6.root");
*/



        TClonesArray *MCTrack = 0;
        MuE::Event *MesmerEvent = 0;
        MUonERecoOutput *ReconstructionOutput = 0;

        cbmsim->SetBranchAddress("MCTrack", &MCTrack);
        cbmsim->SetBranchAddress("MesmerEvent", &MesmerEvent);
        cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);

                auto pos_on_track = [](double q, double m, double z){return (q + m*z);};

double signal=0.;

double reco=0.;
double reco_1=0.; double reco1_1=0.; double more_reco_1=0.; double reco0_1=0.;double reco3_1=0.;
double reco_2=0.; double reco1_2=0.; double more_reco_2=0.; double reco0_2=0.;double reco3_2=0.;
double reco_3=0.; double reco1_3=0.; double more_reco_3=0.; double reco0_3=0.;double reco3_3=0.;
double reco_4=0.; double reco1_4=0.; double more_reco_4=0.; double reco0_4=0.;double reco3_4=0.;
double reco_5=0.; double reco1_5=0.; double more_reco_5=0.; double reco0_5=0.;double reco3_5=0.;
double reco_6=0.; double reco1_6=0.; double more_reco_6=0.; double reco0_6=0.;double reco3_6=0.;

double bin[7]={0.};

double error=0; double error1=0.;double error2=0.;double error3=0.;double error4=0.;double error5=0.;double error6=0.;

double reco_v=0.; double more_reco_v=0.; double reco0_v=0.;
double gen=0.;

int yes2=0; int yes_v=0;
int code_mu=-99; int code_e=-99; int code_mu_in=-99;
int TrackIdreco=-99;
double z_fix=912.7;

double sumW2[6]={0.};
double sum_wgt[6]={0.};

double all=0.;

TH1::SetDefaultSumw2(kTRUE);
   const Int_t NBINS = 7;
   Double_t edges[NBINS + 1] = {0.0, 0.005, 0.010, 0.015, 0.020, 0.025, 0.032, 0.100};


   const Int_t NBINS2 = 12;
   Double_t edges2[NBINS2 + 1] = {0.,0.0002,0.0004,0.0006,0.0008,0.001,0.00125,0.0015,0.00175,0.002,0.003,0.004,0.005};


        std::vector<TH2D*> h_2d(NBINS),h_2d_pre(NBINS);

                for(int m=0; m<NBINS; m++)
        {	string name="h2D"+to_string(m);
                string title="theta mu vs theta E with all cuts"+to_string(m);
                h_2d.at(m)=new TH2D(name.c_str(),title.c_str(), 16, 0.,0.032,NBINS2,edges2);

		string name_pre="h2D_pre"+to_string(m);
                string title_pre="theta mu vs theta E with all cuts"+to_string(m);
                h_2d_pre.at(m)=new TH2D(name_pre.c_str(),title_pre.c_str(), 16, 0.,0.032,NBINS2,edges2);}


   TH1D* theta_e_all = new TH1D("theta_e", "Electron scattering reco angles from MESMER",10,0.,0.035);//NBINS,edges);
   TH1D* theta_mu_all = new TH1D("theta_mu", "Muon scattering reco angles from MESMER",20,0.,0.005);//NBINS2,edges2);

	std::vector<TH1D*>  theta_e(NBINS),theta_mu(NBINS),h_opening(NBINS),d_aco(NBINS), theta_e_pre(NBINS),theta_mu_pre(NBINS),h_opening_pre(NBINS),d_aco_pre(NBINS);


        for(int m=0; m<NBINS; m++)
        {	string name="theta_e"+to_string(m);
                string title="Electron scattering reco angles from MESMER bin"+to_string(m);
                theta_e.at(m)=new TH1D(name.c_str(),title.c_str(),10,0.,0.035);

		string name_mu="theta_mu"+to_string(m);
                string title_mu="Muon scattering reco angles from MESMER bin"+to_string(m);
                theta_mu.at(m)=new TH1D(name_mu.c_str(),title_mu.c_str(),20,0.,0.005);

                string name_op="h_opening"+to_string(m);
                string title_op="Opening angle reco events from MESMER"+to_string(m);
		h_opening.at(m)=new TH1D(name_op.c_str(),title_op.c_str(),35,0.,0.035);

                string name_aco="d_aco_MC"+to_string(m);
                string title_aco="Acoplanarity"+to_string(m);
		d_aco.at(m)=new TH1D(name_aco.c_str(),title_aco.c_str(),600,-3.2,3.2);



		string name_pre="theta_e_pre"+to_string(m);
                string title_pre="Electron scattering reco angles from MESMER  pre-cuts bin"+to_string(m);
                theta_e_pre.at(m)=new TH1D(name_pre.c_str(),title_pre.c_str(),10,0.,0.035);

                string name_mu_pre="theta_mu_pre"+to_string(m);
                string title_mu_pre="Muon scattering reco angles from MESMER  pre-cuts bin"+to_string(m);
                theta_mu_pre.at(m)=new TH1D(name_mu_pre.c_str(),title_mu_pre.c_str(),20,0.,0.005);

                string name_op_pre="h_opening_pre"+to_string(m);
                string title_op_pre="Opening angle reco events from MESMER pre-cuts"+to_string(m);
                h_opening_pre.at(m)=new TH1D(name_op_pre.c_str(),title_op_pre.c_str(),35,0.,0.035);

                string name_aco_pre="d_aco_MC_pre"+to_string(m);
                string title_aco_pre="Acoplanarity pre-cuts"+to_string(m);
                d_aco_pre.at(m)=new TH1D(name_aco_pre.c_str(),title_aco_pre.c_str(),600,-3.2,3.2);

        }


   TH1D* theta_e_gen = new TH1D("theta_e_gen", "Electron scattering generated angles from MESMER",NBINS,edges);
   TH1D* theta_mu_gen = new TH1D("theta_mu_gen", "Muon scattering generated angles from MESMER",20,0.,0.005);
   TH1D* th_mu_ris = new TH1D("th_mu_ris", "muon resolution Emu75,85GeV",300,-0.0006,0.0006);

TH1D *h_res=new TH1D("res", "(the_rec-the_true) 0<theta_e<5 GeV",100,-0.01,0.01);
TH1D *h_res1=new TH1D("res1", "(the_rec-the_true) 5<theta_e<10 GeV",50,-0.0025,0.0025);
TH1D *h_res2=new TH1D("res2", "(the_rec-the_true) 10<theta_e<15 GeV",40,-0.01,0.01);
TH1D *h_res3=new TH1D("res3", "(the_rec-the_true) 15<theta_e<20 GeV",40,-0.01,0.01);
TH1D *h_res4=new TH1D("res4", "(the_rec-the_true) 20<theta_e<25 GeV",40,-0.01,0.01);
TH1D *h_res5=new TH1D("res5", "(the_rec-the_true) 25<theta_e<32 GeV",40,-0.01,0.01);

   TH1D* th_in=new TH1D("th_in","Incoming muon theta",100,0.,0.01);


//        double r_wnorm[6]={20.786765103274643,33.313091221576336,42.396733790329876,50.584815206143652,61.828110400735824,106.88513370134392};
        double r_wnorm[7]={20.786765103274643,33.313091221576336,42.396733790329876,50.584815206143652,61.828110400735824,106.88513370134392,1019.4187737471764};

for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
		cbmsim->GetEntry(i);
		if(i%1000 == 0) cout<<"Entry "<<i<<endl;

        TVector3 p_muin_MC;
        TVector3 p_mu_MC;
        TVector3 p_e_MC;
	double the_gen, thmu_gen;
double emu=0.;


double wnorm=99.;
int index=99;

if(the_gen>=0 and the_gen<=0.005){index=0; wnorm=r_wnorm[0];}
if(the_gen>0.005 and the_gen<=0.010){index=1; wnorm=r_wnorm[1];}
if(the_gen>0.010 and the_gen<=0.015){index=2; wnorm=r_wnorm[2];}
if(the_gen>0.015 and the_gen<=0.020){index=3; wnorm=r_wnorm[3];}
if(the_gen>0.020 and the_gen<=0.025){index=4; wnorm=r_wnorm[4];}
if(the_gen>0.025 and the_gen<=0.032){index=5; wnorm=r_wnorm[5];}
if(the_gen>0.032){index=6; wnorm=r_wnorm[6];}



	for(int n = 0; n < MCTrack->GetEntries(); n++) {
	 const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
	 if(MCTr->interactionID()==0 and MCTr->pdgCode()==-13) {code_mu_in=n; p_muin_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); p_muin_MC.Unit();}
         if(MCTr->interactionID()==45 and MCTr->pdgCode()==11) {code_e=n; p_e_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); p_e_MC.Unit(); the_gen=p_muin_MC.Angle(p_e_MC);}
         if(MCTr->interactionID()==45 and MCTr->pdgCode()==-13) {code_mu=n; p_mu_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); p_mu_MC.Unit(); thmu_gen=p_muin_MC.Angle(p_mu_MC); emu=MCTr->energy();}
	}

 if(code_mu_in!=-99 and code_mu!=-99 and code_e!=-99){

vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
MUonERecoOutputVertex vrtx = ReconstructionOutput->bestVertex();

all+=MesmerEvent->wgt_full;

double chi=vrtx.chi2perDegreeOfFreedom();

int yes_mu=0;
int yes_e=0;
double th_inx,th_iny,x0_in,y0_in;
double chi2_muin;
double stubs_muin;

double th_muin=0.;

int sec0=0;
int sec1=0;
    for(int j=0; j<tracks.size();j++)
    {
        if(tracks.at(j).sector()==0) sec0++;
        if(tracks.at(j).sector()==1) sec1++;
        }

//std::array<double,2> th_rec;

//std::array<TVector3,2> p;
TVector3 p_e,p_mu;
TVector3 p_muin;

double the_rec=-99.;
double thmu_rec=-99.;
        std::vector<MUonERecoOutputHit> hits_mu=vrtx.outgoingMuon().hits();
        std::vector<double> pos_mu; pos_mu.resize(6);
        for(int p=0;p<hits_mu.size();p++)pos_mu.push_back(hits_mu.at(p).position());
        std::vector<MUonERecoOutputHit> hits_e=vrtx.outgoingElectron().hits();
        std::vector<double> pos_e; pos_e.resize(6);
        for(int p=0;p<hits_e.size();p++)pos_e.push_back(hits_e.at(p).position());

MUonERecoOutputTrack t_mu;
MUonERecoOutputTrack t_e;

for(int j=0; j<tracks.size();j++)
{

if(tracks.at(j).processIDofLinkedTrack()==45 and tracks.at(j).sector()==1) TrackIdreco=tracks.at(j).linkedTrackID();
        std::vector<MUonERecoOutputHit> hits_=tracks.at(j).hits();


if(code_mu_in==tracks.at(j).linkedTrackID() and tracks.at(j).sector()==0){
        th_inx=tracks.at(j).xSlope();
        th_iny=tracks.at(j).ySlope();
        x0_in=tracks.at(j).x0();
        y0_in=tracks.at(j).y0();
        chi2_muin=tracks.at(j).chi2perDegreeOfFreedom();
        stubs_muin=hits_.size();
        p_muin.SetXYZ(th_inx,th_iny,1.0);
        p_muin=p_muin.Unit();
	th_muin=p_muin.Theta();
                        }
if(tracks.at(j).processIDofLinkedTrack()==45 and tracks.at(j).sector()==1)
{
/*         TVector3 p1(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0); p1=p1.Unit();th_rec.at(yes2)=p1.Angle(p_muin); p.at(yes2)=p1;
         yes2++;
                 if(code_e==tracks.at(j).linkedTrackID()) {yes_e++;}
                 if(code_mu==tracks.at(j).linkedTrackID()) {yes_mu++;}*/
	std::vector<double> pos; pos.resize(6);
        for(int p=0;p<hits_.size();p++)pos.push_back(hits_.at(p).position());

         if(std::equal(pos.begin(),pos.end(),pos_mu.begin())){p_mu.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0); p_mu.Unit(); thmu_rec=p_mu.Angle(p_muin);t_mu=tracks.at(j);}
        else if(std::equal(pos.begin(),pos.end(),pos_e.begin())){p_e.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0); p_e.Unit(); the_rec=p_e.Angle(p_muin);t_e=tracks.at(j);}

         yes2++;

	}
}


//estrapolo posizione negli ultimi due moduli della seconda stazione
/*double posxIN=pos_on_track(x0_in,th_inx,899.92180);
double posyIN=pos_on_track(y0_in,th_iny,903.76930);*/

// posizione locale negli ultimi due moduli della seconda stazione

double posxIN=99.;//pos_on_track(x0_in,th_inx,z_fix);
double posyIN=99.;//pos_on_track(y0_in,th_iny,z_fix);

//h_xy->Fill(posxIN,posyIN);
std::vector<MUonERecoOutputHit> stubs=ReconstructionOutput->reconstructedHits();

int stub0 = 0;
int stub1 = 0;
         std::array<int,6> module_st1{0};

for(int s=0; s<stubs.size(); s++){
if(stubs.at(s).stationID()==0){stub0++; if(stubs.at(s).moduleID()==4){posxIN=stubs.at(s).position();} else if(stubs.at(s).moduleID()==5){posyIN=stubs.at(s).position();}   }
if(stubs.at(s).stationID()==1){stub1++;module_st1.at(stubs.at(s).moduleID())=1;}
}

if(stubs_muin==6 and abs(posxIN)<=1.5 and abs(posyIN)<=1.5 and chi2_muin<=2 and stub0==6) th_in->Fill(th_muin,MesmerEvent->wgt_full);

if(sec0==1 and stubs_muin==6 and abs(posxIN)<=1.5 and abs(posyIN)<=1.5 and chi2_muin<=2 and stub0==6 and th_muin<0.004){//and stub1<=15){


bin[index]+=MesmerEvent->wgt_full;
signal+=MesmerEvent->wgt_full;

bool allmod=std::all_of(std::begin(module_st1), std::end(module_st1), [](int i){return 1;});
 std::array<double,6> t_e_hits{-99.};
 std::array<double,6> t_mu_hits{-99.};

if(chi!=0 and the_rec!=-99 and thmu_rec!=-99){

//if(th_rec.at(0)>=th_rec.at(1)){the_rec=th_rec.at(0); p_e=p.at(0); thmu_rec=th_rec.at(1); p_mu=p.at(1);}
//else{the_rec=th_rec.at(1); p_e=p.at(1); thmu_rec=th_rec.at(0); p_mu=p.at(0);}

                                                double dotProduct_v = p_mu.Dot(p_e);
                                                TVector3 crossProduct_v = p_mu.Cross(p_e);
                                                double T_v = p_muin.Dot(crossProduct_v);
                                                TVector3 im_v= p_muin.Cross(p_mu);
                                                TVector3 ie_v= p_muin.Cross(p_e);
                                                T_v = T_v>0? 1:-1;
                                                double acoplanarity_v= T_v*(TMath::Pi()- acos( ((im_v).Dot(ie_v))/(im_v.Mag()*ie_v.Mag()) ));

 d_aco_pre[index]->Fill(acoplanarity_v,MesmerEvent->wgt_full);
 h_2d_pre[index]->Fill(the_rec,thmu_rec,MesmerEvent->wgt_full);
theta_mu_pre[index]->Fill(thmu_rec,MesmerEvent->wgt_full);
theta_e_pre[index]->Fill(the_rec,MesmerEvent->wgt_full);
h_opening_pre[index]->Fill(p_mu.Angle(p_e),MesmerEvent->wgt_full);


// if(abs(acoplanarity_v)<=1 and chi<20 and thmu_rec>0.0002 and stub1<=15){ // and vrtx.zKinematicFit()<915. and vrtx.zKinematicFit()>907.){//and the_rec<=0.032 and stub1<=15){// and vrtx.electronThea()>=0.0005 and the_rec<=0.02){//the_rec<=0.032){


double Elastic=0.5109989461*0.001*((1+(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*cos(the_rec)*cos(the_rec))/(1-(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*cos(the_rec)*cos(the_rec)));
double Elastic2=asin( (sin(the_rec)*sqrt(Elastic*Elastic-0.5109989461*0.001*0.5109989461*0.001))/sqrt( (160+0.5109989461*0.001-Elastic)*(160+0.5109989461*0.001-Elastic)-105.6583745 *0.001*105.6583745 *0.001 ) );

bool allmod=std::all_of(std::begin(module_st1), std::end(module_st1), [](int i){return 1;});
 std::array<double,6> t_e_hits{-99.};
 std::array<double,6> t_mu_hits{-99.};

 for(int h=0; h<t_mu.hits().size(); h++){
  if(t_mu.hits().at(h).stationID()==1) t_mu_hits.at(t_mu.hits().at(h).moduleID())=t_mu.hits().at(h).position();
 } 

 for(int h=0; h<t_e.hits().size(); h++){
  if(t_e.hits().at(h).stationID()==1) t_e_hits.at(t_e.hits().at(h).moduleID())=t_e.hits().at(h).position();
 }
  if(allmod and stub1<=20 and thmu_rec>0.0005 and thmu_rec<=0.002 and the_rec<0.010 and the_rec>=0.003 and thmu_rec<=Elastic2+0.0002 and thmu_rec>=Elastic2-0.0002 and p_mu.Angle(p_e)>0.005 and chi<20 and abs(acoplanarity_v)<=0.4 and abs(t_e_hits.at(2)-t_mu_hits.at(2))<0.17){

cout << "event " << i << endl;

//first

sumW2[index]+=(MesmerEvent->wgt_full*MesmerEvent->wgt_full);
sum_wgt[index]+=MesmerEvent->wgt_full;

 reco+=MesmerEvent->wgt_full; error+=MesmerEvent->wgt_full*MesmerEvent->wgt_full;

 d_aco[index]->Fill(acoplanarity_v,MesmerEvent->wgt_full);
 h_2d[index]->Fill(the_rec,thmu_rec,MesmerEvent->wgt_full);
theta_mu_all->Fill(thmu_rec,MesmerEvent->wgt_full);
theta_e_all->Fill(the_rec,MesmerEvent->wgt_full);

theta_mu[index]->Fill(thmu_rec,MesmerEvent->wgt_full);
theta_e[index]->Fill(the_rec,MesmerEvent->wgt_full);
h_opening[index]->Fill(p_mu.Angle(p_e),MesmerEvent->wgt_full);

theta_e_gen->Fill(the_gen,MesmerEvent->wgt_full);


//bin[index]+=MesmerEvent->wgt_full;

if(the_rec>=0.0 and the_rec<0.005){h_res->Fill(the_rec-the_gen,MesmerEvent->wgt_full);}

//second
if(the_rec>=0.005 and the_rec<0.01){h_res1->Fill(the_rec-the_gen,MesmerEvent->wgt_full);}

//third
if(the_rec>=0.01 and the_rec<0.015){h_res2->Fill(the_rec-the_gen,MesmerEvent->wgt_full);}

//fourth
if(the_rec>=0.015 and the_rec<0.02){h_res3->Fill(the_rec-the_gen,MesmerEvent->wgt_full);}

//fifth
if(the_rec>=0.02 and the_rec<0.025){h_res4->Fill(the_rec-the_gen,MesmerEvent->wgt_full);}

//sixth
if(the_rec>=0.025 and the_rec<=0.032){h_res5->Fill(the_rec-the_gen,MesmerEvent->wgt_full);}

			}//aco e chi
		}//chi!=0
	}//chiusura mu_in
    }// chiusura if code_x!0-99
code_e=-99;code_mu=-99;code_mu_in=-99;
TrackIdreco=-99;
yes2=0;yes_v=0;
} //end of general for

cout << "gen VS reco " << gen << " - " << reco << endl;

double Z=6.;
double sigma=1337.95331109692;

cout << "count_preM.Integral()" << all << endl;
cout << "countM.Integral() " << signal << endl;
cout << "count_afterM.Integral() " << reco << endl;
cout << "Events that passes selection over the ones that pass the fiducial " << reco/signal*100 << "%" << endl;
cout << "Events that passes selection over all " << reco/all*100 << "%" << endl;
cout << "Taking into account the CS = "<<sigma<<" micron and Z carbon: " << reco/all*sigma*Z << endl;



cout << "bin0 " << bin[0] << endl;
cout << "bin1 " << bin[1] << endl;
cout << "bin2 " << bin[2] << endl;
cout << "bin3 " << bin[3] << endl;
cout << "bin4 " << bin[4] << endl;
cout << "bin5 " << bin[5] << endl;
cout << "bin6 " << bin[6] << endl;

cout << "ALL" << endl;
cout << "Su " << signal << " eventi di segnale, " << reco << " +- " << sqrt(error) << " sono ricostruiti, con un rapporto del " << reco/signal*100 << "%"<< endl;

cout <<endl;
cout << "Theta range: [0,5] mrad" << endl;
cout << "Su " << signal << " eventi di segnale, " << reco_1 << " +- " << sqrt(error1) << " sono ricostruiti, con un rapporto del " << reco_1/signal*100 << "%"<< endl;
cout << "Su " << signal << " eventi di segnale, " << reco3_1 << " sono ricostruiti con 3 tracce, con un rapporto del " << (reco3_1/signal)*100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con piu tracce (stesso id), " << more_reco_1 << " sono ricostruiti, con un rapporto del " << more_reco_1/signal*100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con 0 tracce di segnale reco, " << reco0_1 << ", con un rapporto del " << reco0_1/signal*100 << "%"<< endl;
//cout << "Su " << signal << " eventi di con 1 sola traccia di segnale reco, " << reco1_1 << ", con un rapporto del " << reco1_1/signal*100 << "%"<< endl;
cout <<endl;
cout << "Theta range: [5,10] mrad" << endl;
cout << "Su " << signal << " eventi di segnale, " << reco_2 << " +- " << sqrt(error2) << " sono ricostruiti, con un rapporto del " << reco_2/signal *100 << "%"<< endl;
cout << "Su " << signal << " eventi di segnale, " << reco3_2 << " sono ricostruiti con 3 tracce, con un rapporto del " << (reco3_2/signal)*100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con piu tracce (stesso id), " << more_reco_2 << " sono ricostruiti, con un rapporto del " <<  more_reco_3/signal *100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con 0 tracce di segnale reco, " << reco0_2 << ", con un rapporto del " << reco0_2/signal *100 << "%"<< endl;
//cout << "Su " << signal << " eventi di con 1 sola traccia di segnale reco, " << reco1_2 << ", con un rapporto del " << reco1_2/signal *100 << "%"<< endl;
cout <<endl;
cout << "Theta range: [10,15] mrad" << endl;
cout << "Su " << signal << " eventi di segnale, " << reco_3 << " +- " << sqrt(error3) << " sono ricostruiti, con un rapporto del " << reco_3/signal *100 << "%"<< endl;
cout << "Su " << signal << " eventi di segnale, " << reco3_3 << " sono ricostruiti con 3 tracce, con un rapporto del " << (reco3_3/signal)*100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con piu tracce (stesso id), " << more_reco_3 << " sono ricostruiti, con un rapporto del " <<  more_reco_3/signal *100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con 0 tracce di segnale reco, " << reco0_3 << ", con un rapporto del " << reco0_3/signal *100 << "%"<< endl;
//cout << "Su " << signal << " eventi di con 1 sola traccia di segnale reco, " << reco1_3 << ", con un rapporto del " << reco1_3/signal *100 << "%"<< endl;
cout <<endl;
cout << "Theta range: [15,20] mrad" << endl;
cout << "Su " << signal << " eventi di segnale, " << reco_4 << " +- " << sqrt(error4) << " sono ricostruiti, con un rapporto del " << reco_4/signal *100 << "%"<< endl;
cout << "Su " << signal << " eventi di segnale, " << reco3_4 << " sono ricostruiti con 3 tracce, con un rapporto del " << (reco3_4/signal)*100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con piu tracce (stesso id), " << more_reco_4 << " sono ricostruiti, con un rapporto del " <<  more_reco_4/signal *100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con 0 tracce di segnale reco, " << reco0_4 << ", con un rapporto del " << reco0_4/signal *100 << "%"<< endl;
//cout << "Su " << signal << " eventi di con 1 sola traccia di segnale reco, " << reco1_4 << ", con un rapporto del " << reco1_4/signal *100 << "%"<< endl;
cout <<endl;
cout << "Theta range: [20,25] mrad" << endl;
cout << "Su " << signal << " eventi di segnale, " << reco_5 << " +- " << sqrt(error5) << " sono ricostruiti, con un rapporto del " << reco_5/signal *100 << "%"<< endl;
cout << "Su " << signal << " eventi di segnale, " << reco3_5 << " sono ricostruiti con 3 tracce, con un rapporto del " << (reco3_5/signal)*100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con piu tracce (stesso id), " << more_reco_5 << " sono ricostruiti, con un rapporto del " <<  more_reco_5/signal *100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con 0 tracce di segnale reco, " << reco0_5 << ", con un rapporto del " << reco0_5/signal *100 << "%"<< endl;
//cout << "Su " << signal << " eventi di con 1 sola traccia di segnale reco, " << reco1_5 << ", con un rapporto del " << reco1_5/signal *100 << "%"<< endl;
cout <<endl;
cout << "Theta range: [25,32] mrad" << endl;
cout << "Su " << signal << " eventi di segnale, " << reco_6 << " +- " << sqrt(error6) << " sono ricostruiti, con un rapporto del " << reco_6/signal *100 << "%"<< endl;
cout << "Su " << signal << " eventi di segnale, " << reco3_6 << " sono ricostruiti con 3 tracce, con un rapporto del " << (reco3_6/signal)*100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con piu tracce (stesso id), " << more_reco_6 << " sono ricostruiti, con un rapporto del " <<  more_reco_6/signal *100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con 0 tracce di segnale reco, " << reco0_6<< ", con un rapporto del " << reco0_6/signal *100 << "%"<< endl;
//cout << "Su " << signal << " eventi di con 1 sola traccia di segnale reco, " << reco1_6 << ", con un rapporto del " << reco1_6/signal *100 << "%"<< endl;

/*

TCanvas r("r","r",700,700);
r.Divide(2,3);
r.cd(1);
h_res->Draw("hist");
r.cd(2);
h_res1->Draw("hist");
r.cd(3);
h_res2->Draw("hist");
r.cd(4);
h_res3->Draw("hist");
r.cd(5);
h_res4->Draw("hist");
r.cd(6);
h_res5->Draw("hist");
//r.SaveAs("comparison_RDMC/res_bend.pdf");

TCanvas a("a","a",700,700);
th_in->Draw("E");
//a.SaveAs("comparison_RDMC/th_in_MC.pdf");


TCanvas z("z","z",200,300);
z.Divide(2,3);
for(int m=0;m<NBINS;m++){
Int_t nx = h_2d.at(m)->GetNbinsX();
Int_t ny = h_2d.at(m)->GetNbinsY();
for (Int_t i=1; i<nx+1; i++) {
for (Int_t j=1; j<ny+1; j++) {
if (h_2d.at(m)->GetBinContent(i,j)<=20) h_2d.at(m)->SetBinContent(i,j,0);}}
h_2d.at(m)->SaveAs(Form("comparison_RDMC/2D_MC_%d_new_pre2_wnorm.root",static_cast<char>(m)));
h_2d_pre.at(m)->SaveAs(Form("comparison_RDMC/preCuts_2D_MC_%d_new_pre2_wnorm.root",static_cast<char>(m)));
z.cd(m+1);
h_2d.at(m)->Draw("COLZ");
}
//z.SaveAs("comparison_RDMC/h_2D_MC.pdf");
*/
TCanvas c("c","c",700,700);
theta_mu_all->Draw("E");
//c.SaveAs("comparison_RDMC/theta_mu_all.pdf");
theta_mu_all->SaveAs("comparison_RDMC/theta_mu_MC_all.root");

TCanvas d("d","d",700,700);
theta_e_all->Draw("E");
//d.SaveAs("comparison_RDMC/theta_e_all.pdf");
theta_e_all->SaveAs("comparison_RDMC/theta_e_MC_all.root");
/*
TCanvas e("e","e",700,700);
theta_mu_gen->Draw("E");
theta_mu_gen->SaveAs("comparison_RDMC/theta_mu_gen_MC_all.root");

TCanvas f("f","f",700,700);
theta_e_gen->Draw("E");
theta_e_gen->SaveAs("comparison_RDMC/theta_e_gen_MC_all.root");

 for(int m=0;m<NBINS;m++){
theta_e[m]->SaveAs(Form("comparison_RDMC/theta_e_MC_%d_new_pre2_wnorm.root",static_cast<char>(m)));
theta_mu[m]->SaveAs(Form("comparison_RDMC/theta_mu_MC_%d_new_pre2_wnorm.root",static_cast<char>(m)));
h_opening[m]->SaveAs(Form("comparison_RDMC/opening_MC_%d_new_pre2_wnorm.root",static_cast<char>(m)));
d_aco[m]->SaveAs(Form("comparison_RDMC/d_aco_MC_%d_new_pre2_wnorm.root",static_cast<char>(m))); 


theta_e_pre[m]->SaveAs(Form("comparison_RDMC/preCuts_theta_e_MC_%d_new_pre2_wnorm.root",static_cast<char>(m)));
theta_mu_pre[m]->SaveAs(Form("comparison_RDMC/preCuts_theta_mu_MC_%d_new_pre2_wnorm.root",static_cast<char>(m)));
h_opening_pre[m]->SaveAs(Form("comparison_RDMC/preCuts_opening_MC_%d_new_pre2_wnorm.root",static_cast<char>(m)));
d_aco_pre[m]->SaveAs(Form("comparison_RDMC/preCuts_d_aco_MC_%d_new_pre2_wnorm.root",static_cast<char>(m))); 

}

for(int b=0; b< theta_e_gen->GetNbinsX(); b++){
cout << "theta_e_gen->GetBinContent(b) " << theta_e_gen->Integral() << endl;
cout << "theta_e_gen sum_wgt " << sum_wgt[b] << endl;
cout << b << " sqrt(sumw2) " << sqrt(sumW2[b]) << " GetBinError " << theta_e_gen->GetBinError(b)<< endl;
}

for(int b=0; b<7; b++){
cout << "theta_e[b]->GetBinContent(b) " << theta_e[b]->Integral() << endl;
cout << "theta_e[b] " << sum_wgt[b] << endl;
cout << b << " sqrt(sumw2) " << sqrt(sumW2[b]) << " GetBinError " << theta_e[b]->GetBinError(b)<< endl;
}
*/


}


