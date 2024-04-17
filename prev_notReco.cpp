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

void prev_notReco(){

        //TFile *inputfile = new TFile("/mnt/raid10/DATA/espedica/fairmu/theta_theta_1M_2hitFirstModules_NOoutchi2_reassign.root");
TFile *inputfile = new TFile("/mnt/raid10/DATA/espedica/fairmu/theta_0-inf_mrad_1M_2hitFirstModules_NOoutchi2_reassign.root");
//Mesmer_sample_1M.root");
        TTree* cbmsim = (TTree*) inputfile->Get("cbmsim");

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

double error=0; double error1=0.;double error2=0.;double error3=0.;double error4=0.;double error5=0.;double error6=0.;

double reco_v=0.; double more_reco_v=0.; double reco0_v=0.;
double gen=0.;

int yes2=0; int yes_v=0;
int code_mu=-99; int code_e=-99; int code_mu_in=-99;
int TrackIdreco=-99;
double z_fix=912.7;

TH2D* h_2d=new TH2D("h2D","theta mu vs theta E with all cuts", 300,0.,0.035,50,0.,0.005);
TH2D* h_xy=new TH2D("h_xy","Beam Spot Z=target",100,-5.,5.,100,-5.,5.);

TH1::SetDefaultSumw2(kTRUE);
   const Int_t NBINS = 6;
   Double_t edges[NBINS + 1] = {0.0, 0.005, 0.010, 0.015, 0.020, 0.025, 0.032};
   TH1D* d_eff = new TH1D("d_eff_MC", "Efficiency as a function of the electron's angle",NBINS,edges);
   TH1D* theta_e = new TH1D("theta_e", "Electron scattering reco angles from MESMER",10,0.,0.035);
   TH1D* theta_mu = new TH1D("theta_mu", "Muon scattering reco angles from MESMER",20,0.,0.005);

   TH1D* theta_e_gen = new TH1D("theta_e_gen", "Electron scattering generated angles from MESMER",10,0.,0.035);
   TH1D* theta_mu_gen = new TH1D("theta_mu_gen", "Muon scattering generated angles from MESMER",20,0.,0.005);
   TH1D* th_mu_ris = new TH1D("th_mu_ris", "muon resolution Emu75,85GeV",300,-0.0006,0.0006);

TH1D *h_res=new TH1D("res", "(the_rec-the_true) 0<theta_e<5 GeV",100,-0.01,0.01);
TH1D *h_res1=new TH1D("res1", "(the_rec-the_true) 5<theta_e<10 GeV",50,-0.0025,0.0025);
TH1D *h_res2=new TH1D("res2", "(the_rec-the_true) 10<theta_e<15 GeV",40,-0.01,0.01);
TH1D *h_res3=new TH1D("res3", "(the_rec-the_true) 15<theta_e<20 GeV",40,-0.01,0.01);
TH1D *h_res4=new TH1D("res4", "(the_rec-the_true) 20<theta_e<25 GeV",40,-0.01,0.01);
TH1D *h_res5=new TH1D("res5", "(the_rec-the_true) 25<theta_e<32 GeV",40,-0.01,0.01);

for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
		cbmsim->GetEntry(i);
		if(i%1000 == 0) cout<<"Entry "<<i<<endl;

        TVector3 p_muin_MC;
        TVector3 p_mu_MC;
        TVector3 p_e_MC;
	double the_gen, thmu_gen;
double emu=0.;

	for(int n = 0; n < MCTrack->GetEntries(); n++) {
	 const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
	 if(MCTr->interactionID()==0 and MCTr->pdgCode()==-13) {code_mu_in=n; p_muin_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); p_muin_MC.Unit();}
         if(MCTr->interactionID()==45 and MCTr->pdgCode()==11) {code_e=n; p_e_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); p_e_MC.Unit(); the_gen=p_muin_MC.Angle(p_e_MC);}
         if(MCTr->interactionID()==45 and MCTr->pdgCode()==-13) {code_mu=n; p_mu_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); p_mu_MC.Unit(); thmu_gen=p_muin_MC.Angle(p_mu_MC); emu=MCTr->energy();}
	}

 if(code_mu_in!=-99 and code_mu!=-99 and code_e!=-99){

vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
MUonERecoOutputVertex vrtx = ReconstructionOutput->bestVertex();

double chi=vrtx.chi2perDegreeOfFreedom();

int yes_mu=0;
int yes_e=0;
double th_inx,th_iny,x0_in,y0_in;
double chi2_muin;
double stubs_muin;
for(int j=0; j<tracks.size();j++)
{

if(tracks.at(j).processIDofLinkedTrack()==45 and tracks.at(j).sector()==1) TrackIdreco=tracks.at(j).linkedTrackID();
if(code_mu_in==tracks.at(j).linkedTrackID() and tracks.at(j).sector()==0){
        th_inx=tracks.at(j).xSlope();
        th_iny=tracks.at(j).ySlope();
        x0_in=tracks.at(j).x0();
        y0_in=tracks.at(j).y0();
        chi2_muin=tracks.at(j).chi2perDegreeOfFreedom();
        std::vector<MUonERecoOutputHit> hits_=tracks.at(j).hits();
        stubs_muin=hits_.size();
                        }
if(tracks.at(j).processIDofLinkedTrack()==45 and tracks.size()==3 and tracks.at(j).sector()==1)
{yes2++;
                 if(code_e==tracks.at(j).linkedTrackID()) yes_e++;
                 if(code_mu==tracks.at(j).linkedTrackID()) yes_mu++;
	}
}

double posxIN=pos_on_track(x0_in,th_inx,z_fix);
double posyIN=pos_on_track(y0_in,th_iny,z_fix);

//h_xy->Fill(posxIN,posyIN);

if(stubs_muin==6 and abs(posxIN)<=1.5 and abs(posyIN)<=1.5 and chi2_muin<2){

signal+=MesmerEvent->wgt_full;


                                                double dotProduct_MC = p_mu_MC.Dot(p_e_MC);
                                                TVector3 crossProduct_MC = p_mu_MC.Cross(p_e_MC);
                                                double T_MC = p_muin_MC.Dot(crossProduct_MC);
                                                TVector3 im_MC= p_muin_MC.Cross(p_mu_MC);
                                                TVector3 ie_MC= p_muin_MC.Cross(p_e_MC);
                                                T_MC = T_MC>0? 1:-1;
                                                double acoplanarity_MC= T_MC*(TMath::Pi()- acos( ((im_MC).Dot(ie_MC))/(im_MC.Mag()*ie_MC.Mag()) ));

if(abs(acoplanarity_MC)<=1 and thmu_gen>0.0002){
gen+=MesmerEvent->wgt_full;
theta_mu_gen->Fill(thmu_gen,MesmerEvent->wgt_full);
theta_e_gen->Fill(the_gen,MesmerEvent->wgt_full);
 }

if(chi!=0){
 MUonERecoOutputTrack mu_in = vrtx.incomingMuon();
 MUonERecoOutputTrack mu_out = vrtx.outgoingMuon();
 MUonERecoOutputTrack e_out = vrtx.outgoingElectron();
        TVector3 p_muin(mu_in.xSlope(),mu_in.ySlope(),1.0);
        TVector3 p_mu(mu_out.xSlope(),mu_out.ySlope(),1.0);
        TVector3 p_e(e_out.xSlope(),e_out.ySlope(),1.0);

                                                double dotProduct_v = p_mu.Dot(p_e);
                                                TVector3 crossProduct_v = p_mu.Cross(p_e);
                                                double T_v = p_muin.Dot(crossProduct_v);
                                                TVector3 im_v= p_muin.Cross(p_mu);
                                                TVector3 ie_v= p_muin.Cross(p_e);
                                                T_v = T_v>0? 1:-1;
                                                double acoplanarity_v= T_v*(TMath::Pi()- acos( ((im_v).Dot(ie_v))/(im_v.Mag()*ie_v.Mag()) ));


 if(abs(acoplanarity_v)<=1 and chi<20 and vrtx.muonTheta()>0.0002){
//first


 reco+=MesmerEvent->wgt_full; error+=MesmerEvent->wgt_full*MesmerEvent->wgt_full;

 d_eff->Fill(vrtx.electronTheta(),MesmerEvent->wgt_full);
 h_2d->Fill(vrtx.electronTheta(),vrtx.muonTheta(),MesmerEvent->wgt_full);
theta_mu->Fill(vrtx.muonTheta(),MesmerEvent->wgt_full);
theta_e->Fill(vrtx.electronTheta(),MesmerEvent->wgt_full);

if(vrtx.electronTheta()>0.0 and vrtx.electronTheta()<=0.005){
h_res->Fill(vrtx.electronTheta()-the_gen,MesmerEvent->wgt_full);
//cout << "[0-5]mrad -> Difference reco (" << vrtx.electronTheta()*1000 << ") - true ("<< the_gen*1000<< "): " << vrtx.electronTheta()-the_gen << endl;

if(yes_e>=1 and yes_mu>=1){reco_1+=MesmerEvent->wgt_full; error1+=MesmerEvent->wgt_full*MesmerEvent->wgt_full;}

if(yes_e==1 and yes_mu==1 and tracks.size()==3) {reco3_1+=MesmerEvent->wgt_full;}

if((yes_e>1 or yes_mu>1) and yes_e!=0 and yes_mu!=0){more_reco_1+=MesmerEvent->wgt_full;}

if(yes2<2 and TrackIdreco==-99){reco0_1+=MesmerEvent->wgt_full;}

if(yes2<2 and TrackIdreco!=-99){reco1_1+=MesmerEvent->wgt_full;}
}

//second
if(vrtx.electronTheta()>0.005 and vrtx.electronTheta()<=0.01){
h_res1->Fill(vrtx.electronTheta()-the_gen,MesmerEvent->wgt_full);

//cout << "[5-10]mrad -> Difference reco (" << vrtx.electronTheta()*1000 << ") - true ("<< the_gen*1000<< "): " << vrtx.electronTheta()-the_gen << endl;

if(yes_e>=1 and yes_mu>=1){reco_2+=MesmerEvent->wgt_full;error2+=MesmerEvent->wgt_full*MesmerEvent->wgt_full;}

if(yes_e==1 and yes_mu==1 and tracks.size()==3) {reco3_2+=MesmerEvent->wgt_full;}

if((yes_e>1 or yes_mu>1) and yes_e!=0 and yes_mu!=0){more_reco_2+=MesmerEvent->wgt_full;}

if(yes2<2 and TrackIdreco==-99){reco0_2+=MesmerEvent->wgt_full;}

if(yes2<2 and TrackIdreco!=-99){reco1_2+=MesmerEvent->wgt_full;}
}

//third
if(vrtx.electronTheta()>0.01 and vrtx.electronTheta()<=0.015){
h_res2->Fill(vrtx.electronTheta()-the_gen,MesmerEvent->wgt_full);

if(yes_e>=1 and yes_mu>=1){reco_3+=MesmerEvent->wgt_full; error3+=MesmerEvent->wgt_full*MesmerEvent->wgt_full;}

if(yes_e==1 and yes_mu==1 and tracks.size()==3) {reco3_3+=MesmerEvent->wgt_full;}

if((yes_e>1 or yes_mu>1) and yes_e!=0 and yes_mu!=0){more_reco_3+=MesmerEvent->wgt_full;}

if(yes2<2 and TrackIdreco==-99){reco0_3+=MesmerEvent->wgt_full;}

if(yes2<2 and TrackIdreco!=-99){reco1_3+=MesmerEvent->wgt_full;}
}

//fourth
if(vrtx.electronTheta()>0.015 and vrtx.electronTheta()<=0.02){
h_res3->Fill(vrtx.electronTheta()-the_gen,MesmerEvent->wgt_full);
if(yes_e>=1 and yes_mu>=1){reco_4+=MesmerEvent->wgt_full; error4+=MesmerEvent->wgt_full*MesmerEvent->wgt_full;}

if(yes_e==1 and yes_mu==1 and tracks.size()==3) {reco3_4+=MesmerEvent->wgt_full;}

if((yes_e>1 or yes_mu>1) and yes_e!=0 and yes_mu!=0){more_reco_4+=MesmerEvent->wgt_full;}

if(yes2<2 and TrackIdreco==-99){reco0_4+=MesmerEvent->wgt_full;}

if(yes2<2 and TrackIdreco!=-99){reco1_4+=MesmerEvent->wgt_full;}
}

//fifth
if(vrtx.electronTheta()>0.02 and vrtx.electronTheta()<=0.025){
h_res4->Fill(vrtx.electronTheta()-the_gen,MesmerEvent->wgt_full);
if(yes_e>=1 and yes_mu>=1){reco_5+=MesmerEvent->wgt_full; error5+=MesmerEvent->wgt_full*MesmerEvent->wgt_full;}

if(yes_e==1 and yes_mu==1 and tracks.size()==3) {reco3_5+=MesmerEvent->wgt_full;}

if((yes_e>1 or yes_mu>1) and yes_e!=0 and yes_mu!=0){more_reco_5+=MesmerEvent->wgt_full;}

if(yes2<2 and TrackIdreco==-99){reco0_5+=MesmerEvent->wgt_full;}

if(yes2<2 and TrackIdreco!=-99){reco1_5+=MesmerEvent->wgt_full;}
}

//sixth
if(vrtx.electronTheta()>0.025 and vrtx.electronTheta()<=0.032){
h_res5->Fill(vrtx.electronTheta()-the_gen,MesmerEvent->wgt_full);
if(yes_e>=1 and yes_mu>=1){reco_6+=MesmerEvent->wgt_full; error6+=MesmerEvent->wgt_full*MesmerEvent->wgt_full;}

if(yes_e==1 and yes_mu==1 and tracks.size()==3) {reco3_6+=MesmerEvent->wgt_full;}

if((yes_e>1 or yes_mu>1) and yes_e!=0 and yes_mu!=0){more_reco_6+=MesmerEvent->wgt_full;}

if(yes2<2 and TrackIdreco==-99){reco0_6+=MesmerEvent->wgt_full;}

if(yes2<2 and TrackIdreco!=-99){reco1_6+=MesmerEvent->wgt_full;}
}

			}//aco e chi
		}//chi!=0
	}//chiusura mu_in
    }// chiusura if code_x!0-99
code_e=-99;code_mu=-99;code_mu_in=-99;
TrackIdreco=-99;
yes2=0;yes_v=0;
} //end of general for

cout << "gen VS reco " << gen << " - " << reco << endl;

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
r.SaveAs("res_bend.pdf");

TCanvas a("a","a",700,700);
d_eff->Draw("E");
a.SaveAs("d_eff_MC.pdf");
d_eff->SaveAs("d_eff_MC_prev.root");

TCanvas b("b","b",700,700);
h_2d->Draw();
h_2d->SaveAs("2D_MC_prev.root");

TCanvas c("c","c",700,700);
theta_mu->Draw("E");
theta_mu->SaveAs("theta_mu_MC_prev.root");

TCanvas d("d","d",700,700);
theta_e->Draw("E");
theta_e->SaveAs("theta_e_MC_prev.root");

TCanvas e("e","e",700,700);
theta_mu_gen->Draw("E");
theta_mu_gen->SaveAs("theta_mu_gen_MC_prev.root");

TCanvas f("f","f",700,700);
theta_e_gen->Draw("E");
theta_e_gen->SaveAs("theta_e_gen_MC_prev.root");


}
