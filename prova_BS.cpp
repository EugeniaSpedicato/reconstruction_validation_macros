#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TSystemDirectory.h"
#include <TStyle.h>

using namespace std;

void prova_BS(){


TChain * cbmsim = new TChain("cbmsim");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/gen_digi_reco/WiP_v0140_commit_258faf6b_0_infmrad.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/bkg/mesmer_bkg_2hit.root");

        MUonERecoOutput *ReconstructionOutput = 0;
        MuE::Event *MesmerEvent = 0;
        TClonesArray *MCTrack = 0;

        cbmsim->SetBranchAddress("MCTrack", &MCTrack);
        cbmsim->SetBranchAddress("MesmerEvent", &MesmerEvent);
        cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);


TH1D *h_T1=new TH1D("T1","Theta of muon_out wrt muon_in from sig (GeV)",100,0.,0.0005);
TH1D *h_T2=new TH1D("T2","Theta muon_in - theta muon_out from sig (GeV)",200,-0.0005,0.0005);
TH1D *h_int= new TH1D("int","minimum bias interaction per particle generated",50,0,50);
TH1D *h_int_rec= new TH1D("intr","minimum bias interaction per particle reconstructed",50,0,50);

TH1D *h_them_all=new TH1D("t2_all","Generated angle of all the gen muons from sig wrt incoming muon (rad)",1143,0,0.04);
TH1D *h_thep_all=new TH1D("t3_all","Generated angle of all the gen electron from sig wrt incoming muon (rad)",1143,0,0.04);

TH1D *h_thel_PP_gen=new TH1D("t2_PP_gen","Generated angle of the reco muon from PP wrt incoming muon (rad)",1143,0,0.04);
TH1D *h_thep_PP_gen=new TH1D("t3_PP_gen","Generated angle of the reco electron from PP wrt incoming muon (rad)",1143,0,0.04);
TH1D *h_thel_PP_rec=new TH1D("t2_PP_PP_rec","Reconstructed angle of the muon from PP wrt incoming muon (rad)",1143,0,0.04);
TH1D *h_thep_PP_rec=new TH1D("t3_rec","Reconstructed angle of the electron from PP wrt incoming muon (rad)",1143,0,0.04);

TH1D *h_thel_gen=new TH1D("t2_gen","Generated angle of the reco muon from sig wrt incoming muon (rad)",1143,0,0.04);
TH1D *h_thel_rec=new TH1D("t2_rec","Reconstructed angle of the muon from sig wrt incoming muon (rad)",1143,0,0.04);

TH1D *energy_r=new TH1D("Er","Energy of the reconstructed electron (GeV)",100,0,10);
TH1D *energy_r_PP=new TH1D("Er","Energy of the reconstructed electron (GeV)",100,0,10);
TH1D *energy_g=new TH1D("Eg","Energy of all the electrons (GeV)",100,0,10);

TH1D *chi_e=new TH1D("chie","Chi2 per DOF muon or electron tracks (#trk==3)",120,0,60);
TH1D *chi_m=new TH1D("chim","Chi2 per DOF muon track (#trk==3",120,0,60);
TH1D *perc=new TH1D("perc","Quality muon and electron tracks",100,0,100);

TH1D *vrtx_chi_sig=new TH1D("chipp","Chi2 per DOF of the kinematic vrtx for sig",40,0,200);
TH1D *vrtx_chi_after02_sig=new TH1D("chiapp","Chi2 per DOF of the kinematic vrtx for sig after th_mu>0.2mrad cut",40,0,200);

TH2D *th_mu_e_sig=new TH2D("th_mu_em_sig","Reconstructed angle mu and mu/el signal with th_mu<35mrad",350,0,0.035,125,0,0.005);
TH2D *th_mu_e_chi_sig=new TH2D("th_mu_em_chi_sig","Angle mu and mu/el signal with th_mu<35mrad with Cv<100",350,0,0.035,125,0,0.005);
TH2D *th_mu_e_aco_sig=new TH2D("th_mu_em_aco_sig","Angle mu and mu/el signal with th_mu<35mrad with abs(aco)<1",350,0,0.035,125,0,0.005);
TH2D *th_mu_e_both_sig=new TH2D("th_mu_em_both_sig","Angle mu and mu/el signal with th_mu<35mrad, Cv<100 and abs(aco)<1",350,0,0.035,125,0,0.005);


TH1D *vrtx_chi_PP=new TH1D("chipp_PP","Chi2 per DOF of the kinematic vrtx for PP",40,0,200);
TH1D *vrtx_chi_after02_PP=new TH1D("chiapp_PP","Chi2 per DOF of the kinematic vrtx for PP after th_mu>0.2mrad cut",40,0,200);

TH2D *th_mu_e_PP=new TH2D("th_mu_em_PP","Reconstructed angle mu and mu/el PP with th_mu<35mrad",350,0,0.035,125,0,0.005);
TH2D *th_mu_e_chi_PP=new TH2D("th_mu_em_chi_PP","Angle mu and mu/el PP with th_mu<35mrad with Cv<100",350,0,0.035,125,0,0.005);
TH2D *th_mu_e_aco_PP=new TH2D("th_mu_em_aco_PP","Angle mu and mu/el PP with th_mu<35mrad with abs(aco)<1",350,0,0.035,125,0,0.005);
TH2D *th_mu_e_both_PP=new TH2D("th_mu_em_both_PP","Angle mu and mu/el PP with th_mu<35mrad, Cv<100 and abs(aco)<1",350,0,0.035,125,0,0.005);


TH1D * Z_sig=new TH1D("Zsig" , "Z position with adaptive fitter for sig", 400,930,1150);
TH1D * Z_PP=new TH1D("ZPP" , "Z position with adaptive fitter for PP", 400,930,1150);

TH1D* h_aco_sig=new TH1D("aco_sig","Acoplanarity of muone+electron from sig mesmer",600,-3.14,3.14);
TH1D* h_aco_02_sig=new TH1D("aco02_sig","Acoplanarity of muone+electron from sig mesmer after th_mu>0.2mrad",600,-3.14,3.14);
TH1D* h_aco_PP=new TH1D("aco_PP","Acoplanarity of muone+electron from PP mesmer",600,-3.14,3.14);
TH1D* h_aco_02_PP=new TH1D("aco02_PP","Acoplanarity of muone+electron from PP mesmer after th_mu>0.2mrad",600,-3.14,3.14);

TH2D *th_mu_e_gen_sig=new TH2D("th_mu_em_gen_sig","Generated angle mu and mu/el signal with th_mu<35mrad",350,0,0.035,125,0,0.005);
TH2D *th_mu_e_gen_PP=new TH2D("th_mu_em_gen","Generated angle mu and mu/el bkg with th_mu<35mrad",350,0,0.035,125,0,0.005);

double danger_sig=0;
double danger_02_sig=0;
double danger_PP=0;
double danger_02_PP=0;
double danger_ghost=0;
double danger_cut_sig=0;
double danger_cut_PP=0;
double all_PP=0;
double all_sig=0;

for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
		cbmsim->GetEntry(i);
		if(i%1000 == 0) cout<<"Entry "<<i<<endl;

double th_muin_rec=0; double thmu_rec=0; double thmuin_gen=0;
double thel_gen=0; double thmu_gen=0;
double thel_rec=0;
double thep_gen=0;
double thep_rec=0;
double yes=0;
	TVector3 pep,pep_rec;
	TVector3 pel,pel_rec;
	TVector3 pmu_in,pmu;
double Em,Eel,Eep;
double Z_ep;
int truth=99;
int code_ep=99;
int code_el=99;
int code_mu=99;

        for(int n = 0; n < MCTrack->GetEntries(); n++) {
         const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
         if(MCTr->interactionID()==0 and MCTr->pdgCode()==-13){pmu_in.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz());pmu_in=pmu_in.Unit();thmuin_gen=pmu_in.Theta();}

         if(MCTr->interactionID()==45){
         if(MCTr->pdgCode()==-11) {pep.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); code_ep=n;pep=pep.Unit();thep_gen=pmu_in.Angle(pep);Eep=MCTr->energy();}
         if(MCTr->pdgCode()==11){pel.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); code_el=n;pel=pel.Unit();thel_gen=pmu_in.Angle(pel);Eel=MCTr->energy();}
         if(MCTr->pdgCode()==-13){pmu.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); code_mu=n;pmu=pmu.Unit();thmu_gen=pmu_in.Angle(pmu);}
        }
	}
if(code_el!=99 and code_mu!=99) truth=1;
if(code_el!=99 and code_ep!=99 and code_mu!=99) truth=0;


if(truth==0 or truth==1)// and Z_ep<1037 and Z_ep>1031)// and thel_gen<=0.035)// and Z_ep<1037 and Z_ep>1031)
{

 if(truth==0)all_PP+=MesmerEvent->wgt_full;
 if(truth==1)all_sig+=MesmerEvent->wgt_full;

int yes_m=0; int yes_el=0; int yes_ep=0; int yes_mu2=0;

vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
TVector3 th_muin;
TVector3 th_mu;  vector<TVector3> p_mu_v;

double acoplanarity;
vector<TVector3> pel_rec_v;
vector<TVector3> pep_rec_v;

vector<double> chi_min_m,chi_min_ep,chi_min_el,chi_min_mu, th_el_rec_vec,th_ep_rec_vec, th_muin_rec_vec,thmu_rec_vec;
chi_min_m.reserve(5);chi_min_mu.reserve(5);chi_min_ep.reserve(5),chi_min_el.reserve(5);th_el_rec_vec.reserve(5);th_muin_rec_vec.reserve(5);thmu_rec_vec.reserve(5);

double th_ep;

int st0_m=0;

for(int j=0; j<tracks.size();j++)
{
 if(tracks.at(j).sector()==0) st0_m++;
}


for(int j=0; j<tracks.size();j++)
{

//cout << "int_ID " << tracks.at(j).processIDofLinkedTrack() << " sector " << tracks.at(j).sector() << " link " << tracks.at(j).linkedTrackID() << endl;

double th_inx,th_iny;
if(tracks.at(j).processIDofLinkedTrack()==0 and tracks.at(j).sector()==0){
th_inx=tracks.at(j).xSlope();
th_iny=tracks.at(j).ySlope();
th_muin.SetXYZ(th_inx,th_iny,1.0);
th_muin=th_muin.Unit();
th_muin_rec=th_muin.Theta();
}

if(tracks.at(j).processIDofLinkedTrack()==45 and tracks.at(j).linkedTrackID()==code_mu and tracks.at(j).sector()==1){
yes_mu2+=MesmerEvent->wgt_full;
th_inx=tracks.at(j).xSlope();
th_iny=tracks.at(j).ySlope();
TVector3 p;
p.SetXYZ(th_inx,th_iny,1.0);
p=p.Unit();
p_mu_v.push_back(p);
chi_min_mu.push_back(tracks.at(j).chi2perDegreeOfFreedom());
thmu_rec_vec.push_back(th_muin.Angle(p));
 }

if(tracks.at(j).processIDofLinkedTrack()==45 and tracks.at(j).linkedTrackID()==code_el and tracks.at(j).sector()==1)
{yes_el+=MesmerEvent->wgt_full;
 TVector3 p;
 p.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.);
 p=p.Unit();
 pel_rec_v.push_back(p);
 th_el_rec_vec.push_back(th_muin.Angle(p));
 chi_min_el.push_back(tracks.at(j).chi2perDegreeOfFreedom());
}

if(tracks.at(j).processIDofLinkedTrack()==45 and (tracks.at(j).linkedTrackID()==code_ep or tracks.at(j).linkedTrackID()==code_el)){
 yes_ep+=MesmerEvent->wgt_full;
 TVector3 p;
 p.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.);
 p=p.Unit();
 pep_rec_v.push_back(p);
 th_ep_rec_vec.push_back(th_muin.Angle(p));
 chi_min_ep.push_back(tracks.at(j).chi2perDegreeOfFreedom());}


}

if( thmu_rec_vec.size()!=0){ auto it = min_element(chi_min_mu.begin(),chi_min_mu.end());
                                 thmu_rec = thmu_rec_vec.at(std::distance(chi_min_mu.begin(), it));
                                 th_mu = p_mu_v.at(std::distance(chi_min_mu.begin(), it));}


if( th_el_rec_vec.size()!=0 and truth==1){
 auto it = min_element(chi_min_el.begin(),chi_min_el.end()); thel_rec = th_el_rec_vec.at(std::distance(chi_min_el.begin(), it));
                                pel_rec=pel_rec_v.at(std::distance(chi_min_el.begin(), it));
                                energy_r->Fill(Eel,MesmerEvent->wgt_full);
                                h_thel_gen->Fill(thel_gen,MesmerEvent->wgt_full);
                                h_thel_rec->Fill(thel_rec,MesmerEvent->wgt_full);
                                if(yes_mu2!=0){ double dotProduct = th_mu.Dot(pel_rec);
                                                TVector3 crossProduct = th_mu.Cross(pel_rec);
                                                double T = th_muin.Dot(crossProduct);
                                                TVector3 im= th_muin.Cross(th_mu);
                                                TVector3 ie= th_muin.Cross(pel_rec);
                                                T = T>0? 1:-1;
                                                acoplanarity= T*(TMath::Pi()- acos( ((im).Dot(ie))/(im.Mag()*ie.Mag()) ));

                                        }
		}



if( th_el_rec_vec.size()!=0 and truth==0){
 auto it = min_element(chi_min_el.begin(),chi_min_el.end()); thel_rec = th_el_rec_vec.at(std::distance(chi_min_el.begin(), it));
                                pel_rec=pel_rec_v.at(std::distance(chi_min_el.begin(), it));
                                energy_r_PP->Fill(Eel,MesmerEvent->wgt_full);
                                h_thel_PP_gen->Fill(thel_gen,MesmerEvent->wgt_full);
                                h_thel_PP_rec->Fill(thel_rec,MesmerEvent->wgt_full);
                                if(yes_mu2!=0){ double dotProduct = th_mu.Dot(pel_rec);
                                                TVector3 crossProduct = th_mu.Cross(pel_rec);
                                                double T = th_muin.Dot(crossProduct);
                                                TVector3 im= th_muin.Cross(th_mu);
                                                TVector3 ie= th_muin.Cross(pel_rec);
                                                T = T>0? 1:-1;
                                                acoplanarity= T*(TMath::Pi()- acos( ((im).Dot(ie))/(im.Mag()*ie.Mag()) ));
                                        }
                }


if( th_ep_rec_vec.size()!=0 and truth==0){
 auto it = min_element(chi_min_ep.begin(),chi_min_ep.end()); thep_rec = th_ep_rec_vec.at(std::distance(chi_min_ep.begin(), it));
                                pep_rec=pep_rec_v.at(std::distance(chi_min_ep.begin(), it));
                                energy_r_PP->Fill(Eep,MesmerEvent->wgt_full);
                                h_thep_PP_gen->Fill(thep_gen,MesmerEvent->wgt_full);
                                h_thep_PP_rec->Fill(thep_rec,MesmerEvent->wgt_full);
                                if(yes_mu2!=0){ double dotProduct = th_mu.Dot(pep_rec);
                                                TVector3 crossProduct = th_mu.Cross(pep_rec);
                                                double T = th_muin.Dot(crossProduct);
                                                TVector3 im= th_muin.Cross(th_mu);
                                                TVector3 ie= th_muin.Cross(pep_rec);
                                                T = T>0? 1:-1;
                                                acoplanarity= T*(TMath::Pi()- acos( ((im).Dot(ie))/(im.Mag()*ie.Mag()) ));
                                        }
                }


MUonERecoOutputVertex vrtx = ReconstructionOutput->bestVertex();
double chi;

	 chi=vrtx.chi2perDegreeOfFreedom();
         if(vrtx.chi2perDegreeOfFreedom()<50 and truth==1) Z_sig->Fill(vrtx.zPositionFit(),MesmerEvent->wgt_full);
         if(vrtx.chi2perDegreeOfFreedom()<50 and truth==0) Z_PP->Fill(vrtx.zPositionFit(),MesmerEvent->wgt_full);


if(yes_el==1 and yes_mu2==1 and tracks.size()==3 and st0_m==1 and chi>0 and thel_rec<0.035 and abs(acoplanarity)>=0 and truth==1){ //tracks.size()==4 and st0_m>1 and thel_rec<0.035
 danger_sig+=MesmerEvent->wgt_full;

 vrtx_chi_sig->Fill(chi,MesmerEvent->wgt_full);
 h_aco_sig->Fill(acoplanarity,MesmerEvent->wgt_full);
 th_mu_e_sig->Fill(thel_rec,thmu_rec,MesmerEvent->wgt_full);
 th_mu_e_gen_sig->Fill(thel_rec,thmu_rec,MesmerEvent->wgt_full);

 if(thmu_rec>0.0002){danger_02_sig+=MesmerEvent->wgt_full; vrtx_chi_after02_sig->Fill(chi,MesmerEvent->wgt_full); h_aco_02_sig->Fill(acoplanarity,MesmerEvent->wgt_full);}

 if(chi<20) {th_mu_e_chi_sig->Fill(thel_rec,thmu_rec,MesmerEvent->wgt_full);}
 if(abs(acoplanarity)<1) {th_mu_e_aco_sig->Fill(thel_rec,thmu_rec,MesmerEvent->wgt_full);}
 if(chi<20 and abs(acoplanarity)<1) {th_mu_e_both_sig->Fill(thel_rec,thmu_rec,MesmerEvent->wgt_full);

  if(chi<20 and abs(acoplanarity)<1) danger_cut_sig+=MesmerEvent->wgt_full;}
}
else if( yes_el==1 and yes_mu2==1 and tracks.size()==3 and st0_m==1 and chi>0 and thel_rec<0.035 and abs(acoplanarity)>=0 and truth==0){ //tracks.size()==4 and st0_m>1 and thel_rec<0.035
 danger_PP+=MesmerEvent->wgt_full;

 vrtx_chi_PP->Fill(chi,MesmerEvent->wgt_full);
 h_aco_PP->Fill(acoplanarity,MesmerEvent->wgt_full);
 th_mu_e_PP->Fill(thel_rec,thmu_rec,MesmerEvent->wgt_full);
 th_mu_e_gen_PP->Fill(thel_rec,thmu_rec,MesmerEvent->wgt_full);

 if(thmu_rec>0.0002){danger_02_PP+=MesmerEvent->wgt_full; vrtx_chi_after02_PP->Fill(chi,MesmerEvent->wgt_full); h_aco_02_PP->Fill(acoplanarity,MesmerEvent->wgt_full);}

 if(chi<20) {th_mu_e_chi_PP->Fill(thel_rec,thmu_rec,MesmerEvent->wgt_full);}
 if(abs(acoplanarity)<1) {th_mu_e_aco_PP->Fill(thel_rec,thmu_rec,MesmerEvent->wgt_full);}
 if(chi<20 and abs(acoplanarity)<1) {th_mu_e_both_PP->Fill(thel_rec,thmu_rec,MesmerEvent->wgt_full);

  if(chi<20 and abs(acoplanarity)<1) danger_cut_PP+=MesmerEvent->wgt_full;}
}
else if( yes_ep==1 and yes_mu2==1 and tracks.size()==3 and st0_m==1 and chi>0 and thep_rec<0.035 and abs(acoplanarity)>=0 and truth==0){ //tracks.size()==4 and st0_m>1 and thel_rec<0.035
 danger_PP+=MesmerEvent->wgt_full;

 vrtx_chi_PP->Fill(chi,MesmerEvent->wgt_full);
 h_aco_PP->Fill(acoplanarity,MesmerEvent->wgt_full);
 th_mu_e_PP->Fill(thep_rec,thmu_rec,MesmerEvent->wgt_full);

 if(thmu_rec>0.0002){danger_02_PP+=MesmerEvent->wgt_full; vrtx_chi_after02_PP->Fill(chi,MesmerEvent->wgt_full); h_aco_02_PP->Fill(acoplanarity,MesmerEvent->wgt_full);}

 if(chi<20) {th_mu_e_chi_PP->Fill(thep_rec,thmu_rec,MesmerEvent->wgt_full);}
 if(abs(acoplanarity)<1) {th_mu_e_aco_PP->Fill(thep_rec,thmu_rec,MesmerEvent->wgt_full);}
 if(chi<20 and abs(acoplanarity)<1) {th_mu_e_both_PP->Fill(thep_rec,thmu_rec,MesmerEvent->wgt_full);

  if(chi<20 and abs(acoplanarity)<1) danger_cut_PP+=MesmerEvent->wgt_full;}
}
/*else if((yes_el==2 and yes_mu2==0 and tracks.size()==3 and st0_m==1 and chi>0 and th_el_rec_vec.at(1)<0.035 and th_el_rec_vec.at(0)<0.035 and abs(acoplanarity)>=0) or (yes_el==0 and yes_mu2==2 and tracks.size()==3 and st0_m==1 and chi>0) or (yes_m==0 and chi>0 and yes_el==0 and tracks.size()==3 and st0_m==1 and thmu_rec_vec.size()==2 and abs(acoplanarity)>=0 and thmu_rec_vec.at(0)<0.35 and thmu_rec_vec.at(1)<0.35) or (yes_ep==2 and yes_mu2==0 and tracks.size()==3 and st0_m==1 and chi>0 and th_ep_rec_vec.at(1)<0.035 and th_ep_rec_vec.at(0)<0.035 and abs(acoplanarity)>=0))
	{danger_ghost+=MesmerEvent->wgt_full;}*/
//else if(yes_el==0 and thmu_rec_vec.size()==1 and tracks.size()==3 and st0_m==1 and chi>0) {other+=MesmerEvent->wgt_full;}

 }//end yes==2

} //end of general for


cout << "Su " << all_sig << " di muoni, hanno una sig pericolosa " << danger_sig << " con una percentuale di " << 100*(danger_sig/all_sig) << "%" << endl;
cout << "Su " << all_sig << " di muoni, hanno una sig pericolosa con th_mu>0.2mrad " << danger_02_sig << " con una percentuale di " << 100*(danger_02_sig/all_sig) << "%" << endl;
cout <<"Su " << all_sig << " di muoni, hanno una sig pericolosa con tutti i tagli " << danger_cut_sig << " con una percentuale di " << 100*(danger_cut_sig/all_sig) << "%" << endl;

cout << "Su " << all_PP << " di muoni, hanno una PP pericolosa " << danger_PP << " con una percentuale di " << 100*(danger_PP/all_PP) << "%" << endl;
cout << "Su " << all_PP << " di muoni, hanno una PP pericolosa con th_mu>0.2mrad " << danger_02_PP << " con una percentuale di " << 100*(danger_02_sig/all_PP) << "%" << endl;
cout <<"Su " << all_PP << " di muoni, hanno una PP pericolosa con tutti i tagli " << danger_cut_PP << " con una percentuale di " << 100*(danger_cut_PP/all_PP) << "%" << endl;

cout << "Su " << all_PP+all_sig << " di muoni, hanno 2 tracce con 1 ghost " << danger_ghost << " con una percentuale di " << 100*(danger_ghost/all_PP+all_sig) << "%" << endl;

/*
h_part1
h_part2
h_part3
h_part_more
h_T1
h_T2
energy_r
h_thel_gen
h_thel_rec
energy_r_PP
h_thel_PP_gen
h_thel_PP_rec
energy_r_PP
h_thep_PP_gen
h_thep_PP_rec
Z_sig
Z_PP
vrtx_chi_sig
h_aco_sig
vrtx_chi_after02_sig
h_aco_02_sig
 vrtx_chi_PP
h_aco_PP
th_mu_e_PP
vrtx_chi_after02_PP
h_aco_02_PP
vrtx_chi_PP
h_aco_PP
th_mu_e_PP
vrtx_chi_after02_PP
h_aco_02_PP
*/

TF1 *Elastic = new TF1("Elastic","0.5109989461*0.001*((1+(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*cos(x)*cos(x))/(1-(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*cos(x)*cos(x)))",0,0.035); 
TF1 *Elastic2 = new TF1("Elastic2","asin( (sin(x)*sqrt(Elastic(x)*Elastic(x)-0.5109989461*0.001*0.5109989461*0.001))/sqrt( (160+0.5109989461*0.001-Elastic(x))*(160+0.5109989461*0.001-Elastic(x))-105.6583745 *0.001*105.6583745 *0.001 ) )",0,0.035); 
TLine *t = new TLine(0.,0.0002,0.035,0.0002);
t->SetLineWidth(1);
t->SetLineColor(kRed);
Elastic2->SetLineColor(6);
Elastic2->SetLineWidth(1);

TCanvas c51("c51","c51",700,700);
c51.Divide(2,3);
c51.cd(1);
th_mu_e_gen_sig->SetMarkerColor(kBlue);
th_mu_e_gen_sig->SetMarkerStyle(43);
th_mu_e_gen_sig->SetMarkerSize(1);
th_mu_e_gen_sig->Draw();
th_mu_e_gen_PP->SetMarkerColor(kGreen);
th_mu_e_gen_PP->SetMarkerStyle(43);
th_mu_e_gen_PP->SetMarkerSize(1);
th_mu_e_gen_PP->Draw("same");
Elastic2->Draw("same");
t->Draw("same");
c51.cd(2);
th_mu_e_sig->SetMarkerColor(kBlue);
th_mu_e_sig->SetMarkerStyle(43);
th_mu_e_sig->SetMarkerSize(1);
th_mu_e_sig->Draw();
th_mu_e_PP->SetMarkerColor(kGreen);
th_mu_e_PP->SetMarkerStyle(43);
th_mu_e_PP->SetMarkerSize(1);
th_mu_e_PP->Draw("same");
Elastic2->Draw("same");
t->Draw("same");
c51.cd(3);
th_mu_e_aco_sig->SetMarkerColor(kBlue);
th_mu_e_aco_sig->SetMarkerStyle(43);
th_mu_e_aco_sig->SetMarkerSize(1);
th_mu_e_aco_sig->Draw();
th_mu_e_aco_PP->SetMarkerColor(kGreen);
th_mu_e_aco_PP->SetMarkerStyle(43);
th_mu_e_aco_PP->SetMarkerSize(1);
th_mu_e_aco_PP->Draw("same");
Elastic2->Draw("same");
t->Draw("same");
c51.cd(4);
th_mu_e_chi_sig->SetMarkerColor(kBlue);
th_mu_e_chi_sig->SetMarkerStyle(43);
th_mu_e_chi_sig->SetMarkerSize(1);
th_mu_e_chi_sig->Draw();
th_mu_e_chi_PP->SetMarkerColor(kGreen);
th_mu_e_chi_PP->SetMarkerStyle(43);
th_mu_e_chi_PP->SetMarkerSize(1);
th_mu_e_chi_PP->Draw("same");
Elastic2->Draw("same");
t->Draw("same");
c51.cd(5);
th_mu_e_both_sig->SetMarkerColor(kBlue);
th_mu_e_both_sig->SetMarkerStyle(43);
th_mu_e_both_sig->SetMarkerSize(1);
th_mu_e_both_sig->Draw();
th_mu_e_both_PP->SetMarkerColor(kGreen);
th_mu_e_both_PP->SetMarkerStyle(43);
th_mu_e_both_PP->SetMarkerSize(1);
th_mu_e_both_PP->Draw("same");
Elastic2->Draw("same");
t->Draw("same");
c51.SaveAs("bkg_mesmer/2DsigPP.pdf");


TCanvas ac("ac","ac",700,700);
ac.Divide(2,2);
ac.cd(1);
h_aco_PP->SetLineColor(kGreen);
h_aco_sig->Draw("hist");
h_aco_PP->Draw("hist same");
ac.cd(2);
h_aco_02_PP->SetLineColor(kGreen);
h_aco_02_sig->Draw("hist");
h_aco_02_PP->Draw("hist same");
ac.cd(3);
h_aco_PP->SetLineColor(kGreen);
h_aco_sig->Draw("hist");
h_aco_PP->Draw("hist same");
ac.cd(4);
h_aco_02_PP->SetLineColor(kGreen);
h_aco_02_sig->Draw("hist");
h_aco_02_PP->Draw("hist same");

ac.SaveAs("bkg_mesmer/h_aco_sigPP_min.pdf");

TCanvas z("z","z",700,700);
Z_sig->SetLineColor(kBlue);
Z_PP->SetLineColor(kGreen);
Z_sig->Draw("hist");
Z_PP->Draw("hist same");
z.SaveAs("bkg_mesmer/zz.pdf");
/*
TCanvas b2("b2","b2",700,700);
b2.Divide(2,2);
b2.cd(1);
h_thep_all->Draw("hist");
h_thel_gen->SetLineColor(kRed);
h_thel_gen->Draw("hist same");
h_thel_rec->SetLineColor(kRed);
h_thel_rec->Draw("hist same");
gPad->SetLogy();
b2.cd(2);
 TH1F *rec_t = (TH1F*)h_thel_rec->Clone("rec_t");
 rec_t->Sumw2();
 rec_t->Divide(h_thep_all);
 rec_t->Draw("ep");
b2.cd(3);
h_T1->Draw("hist");
b2.cd(4);
h_T2->Draw("hist");
b2.SaveAs("theta_signal.pdf");

TCanvas b3("b3","b3",700,700);
b3.Divide(1,2);
b3.cd(1);
h_int->Draw("hist");
b3.cd(2);
h_int_rec->SetLineColor(kRed);
h_int_rec->Draw("hist");
gPad->SetLogy();
b3.SaveAs("int_signal.pdf");

TCanvas b4("b4","b4",700,700);
b4.Divide(1,3);
b4.cd(1);
energy_g->Draw("hist");
energy_r->SetLineColor(kRed);
energy_r->Draw("hist same");
gPad->SetLogy();
b4.cd(2);
 TH1F *rec = (TH1F*)energy_r->Clone("rec");
 rec->Sumw2();
 rec->Divide(energy_g);
 rec->Draw("ep");
b4.cd(3);
chi_m->Draw("hist");
chi_e->SetLineColor(kRed);
chi_e->Draw("hist same");
gPad->SetLogy();
b4.cd(4);
perc->Draw("hist");
b4.SaveAs("energy_signal.pdf");

TCanvas b5("b5","b5",700,700);
b5.Divide(2,2);
b5.cd(1);
th_mu_e->Draw("COLZ");
b5.cd(2);
th_gst->Draw("COLZ");
b5.cd(3);
Z_sig->Draw("hist");
b5.cd(4);
Z_sig_gen->Draw("hist");
Z_sig_rec->SetLineColor(kRed);
Z_sig_rec->Draw("hist same");
b5.SaveAs("th_mu_e_signal.pdf");

TCanvas b6("b6","b6",700,700);
b6.Divide(2,3);
b6.cd(1);
mult->Draw("hist");
b6.cd(2);
h_part1->Draw("hist");
b6.cd(3);
h_part2->Draw("hist");
b6.cd(4);
h_part3->Draw("hist");
b6.cd(5);
h_part_more->Draw("hist");
b6.SaveAs("mult_signal.pdf");

TCanvas b7("b7","b7",700,700);
b7.Divide(1,2);
b7.cd(1);
vrtx_chi->Draw("hist");
gPad->SetLogy();
b7.cd(2);
vrtx_chi_after02->Draw("hist");
gPad->SetLogy();
b7.SaveAs("vrtx_chi_sig_signal.pdf");
vrtx_chi->SaveAs("vrtx_chi_sig_signal.root");
vrtx_chi_after02->SaveAs("vrtx_chi_sig_after_signal.root");
*/
}
