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

void RealDataAnalyzer(){


TChain * cbmsim = new TChain("cbmsim");
//cbmsim->Add("TRPP_minbias_offset/TRPP_minbias_1M_firstSample.root");
//cbmsim->Add("TRPP_minbias_offset/TRPP_minbias_1M_secondSample.root");
//cbmsim->Add("TRPP_minbias_offset/TRPP_minbias_1M_thirdSample.root");

cbmsim->Add("TRPP_minbias_1M_nobend_firstSample.root");

        MUonERecoOutput *ReconstructionOutput = 0;
        TClonesArray *MCTrack = 0;

        cbmsim->SetBranchAddress("MCTrack", &MCTrack);
        cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);


TH1D *h_T1=new TH1D("T1","Theta of muon_out wrt muon_in from PP (GeV)",100,0.,0.0005);
TH1D *h_T2=new TH1D("T2","Theta muon_in - theta muon_out from PP (GeV)",200,-0.0005,0.0005);
TH1D *h_int= new TH1D("int","minimum bias interaction per particle generated",50,0,50);
TH1D *h_int_rec= new TH1D("intr","minimum bias interaction per particle reconstructed",50,0,50);

TH1D *h_them_all=new TH1D("t2_all","Generated angle of all the gen electrons from PP wrt incoming muon (rad)",2857,0,0.1);
TH1D *h_thep_all=new TH1D("t3_all","Generated angle of all the gen positron from PP wrt incoming muon (rad)",2857,0,0.1);

TH1D *h_them_gen=new TH1D("t2_gen","Generated angle of the reco electron from PP wrt incoming muon (rad)",2857,0,0.1);
TH1D *h_thep_gen=new TH1D("t3_gen","Generated angle of the reco positron from PP wrt incoming muon (rad)",2857,0,0.1);
TH1D *h_them_rec=new TH1D("t2_rec","Reconstructed ngle of the electron from PP wrt incoming muon (rad)",2857,0,0.1);
TH1D *h_thep_rec=new TH1D("t3_rec","Reconstructed angle of the positron from PP wrt incoming muon (rad)",2857,0,0.1);


TH1D *energy_r=new TH1D("Er","Energy of the reconstructed electron/positron (GeV)",100,0,10);
TH1D *energy_g=new TH1D("Eg","Energy of all the electron/positron (GeV)",100,0,10);

TH1D *chi_e=new TH1D("chie","Chi2 per DOF electron or positron tracks (#trk==3)",120,0,60);
TH1D *chi_m=new TH1D("chim","Chi2 per DOF muon track (#trk==3",120,0,60);
TH1D *perc=new TH1D("perc","Quality electron and positron tracks",100,0,100);

TH1D *vrtx_chi=new TH1D("chipp","Chi2 per DOF of the kinematic vrtx for PP",40,0,200);
TH1D *vrtx_chi_after02=new TH1D("chiapp","Chi2 per DOF of the kinematic vrtx for sig after th_mu>0.2mrad cut",40,0,200);

TH2D *th_mu_e=new TH2D("th_mu_em","Reconstructed ngle mu and mu/el from sig in acceptance (th_mu<35mrad)",350,0,0.035,125,0,0.005);
TH2D *th_mu_e_chi=new TH2D("th_mu_em_chi","Angle mu and mu/el from sig in acceptance (th_mu<35mrad) with Cv<100",350,0,0.035,125,0,0.005);
TH2D *th_mu_e_aco=new TH2D("th_mu_em_aco","Angle mu and mu/el from sig in acceptance (th_mu<35mrad) with abs(aco)<1",350,0,0.035,125,0,0.005);
TH2D *th_mu_e_both=new TH2D("th_mu_em_both","Angle mu and mu/el from sig in acceptance (th_mu<35mrad) with Cv<100 and abs(aco)<1",350,0,0.035,125,0,0.005);

TH2D *th_gst=new TH2D("thgst","Angle two particles where we have 1 ghost 0r 1 particle !=sig",350,0,0.035,125,0,0.005);
TH2D *th_gst_chi=new TH2D("thgst_chi","Angle mu and mu/el from sig in acceptance (th_mu<35mrad) with Cv<100",350,0,0.035,125,0,0.005);
TH2D *th_gst_aco=new TH2D("thgst_aco","Angle mu and mu/el from sig in acceptance (th_mu<35mrad) with abs(aco)<1",350,0,0.035,125,0,0.005);
TH2D *th_gst_both=new TH2D("th_gst_both","Angle mu and mu/el from sig in acceptance (th_mu<35mrad) with Cv<100 and abs(aco)<1",350,0,0.035,125,0,0.005);


TH1D * Z_PP=new TH1D("Zsig" , "Z position with adaptive fitter for sig", 400,930,1150);
TH1D * Z_PP_gen=new TH1D("Zsig_gen" , "Generated Z position sig particles in acceptance (th_mu<35mrad)", 400,930,1150);
TH1D * Z_PP_gen_chi=new TH1D("Zsig_gen_chi" , "Generated Z position sig particles in acceptance (th_mu<35mrad) with Cv<100", 400,930,1150);
TH1D * Z_PP_gen_aco=new TH1D("Zsig_gen_aco" , "Generated Z position sig particles in acceptance (th_mu<35mrad) with abs(aco)<1", 400,930,1150);
TH1D * Z_PP_gen_both=new TH1D("Zsig_gen_both" , "Generated Z position sig particles in acceptance (th_mu<35mrad) with Cv<100 and abs(aco)<1", 400,930,1150);

TH1D * Z_PP_gen02=new TH1D("Zsig_gen02" , "Generated Z position sig particles in acceptance (th_mu<35mrad) and th_mu>0.2mrad", 400,930,1150);
TH1D * Z_PP_gen_chi02=new TH1D("Zsig_gen_chi02" , "Generated Z position sig particles in acceptance (th_mu<35mrad) with Cv<100 and th_mu>0.2mrad", 400,930,1150);
TH1D * Z_PP_gen_aco02=new TH1D("Zsig_gen_aco02" , "Generated Z position sig particles in acceptance (th_mu<35mrad) with abs(aco)<1 and th_mu>0.2mrad", 400,930,1150);
TH1D * Z_PP_gen_both02=new TH1D("Zsig_gen_both02" , "Generated Z position sig particles in acceptance (th_mu<35mrad) with Cv<100 and abs(aco)<1 and th_mu>0.2mrad", 400,930,1150);


TH1D * Z_pp=new TH1D("Zpp" , "Z position with adaptive fitter for PP", 400,930,1150);
TH1D * Z_pp_rec=new TH1D("Zpp_rec" , "Generated Z position PP particles when #reco==3", 400,930,1150);

TH1D* mult=new TH1D("mul","multiplicity of reconstructed tracks in second station when PP happens", 20,0,20);
TH1D *h_part1=new TH1D("p1","reconstructed particle ID second station's multiplicity=1", 50,0,50);
TH1D *h_part2=new TH1D("p2","reconstructed particle ID second station's multiplicity=2", 50,0,50);
TH1D *h_part3=new TH1D("p3","reconstructed particle ID second station's multiplicity=3", 50,0,50);
TH1D *h_part_more=new TH1D("pm","reconstructed particle ID multiplicity>3", 50,0,50);

TH1D* h_aco=new TH1D("aco","Acoplanarity of muone+electron/positron from PP",600,-3.14,3.14);
TH1D* h_aco_02=new TH1D("aco02","Acoplanarity of muone+electron/positron from PP after th_mu>0.2mrad",600,-3.14,3.14);


double danger=0;
double danger_02=0;
double danger_ghost=0;
double danger_ee=0;
double other=0;
double other0=0;
double danger_ghost02=0;
double other02=0;
/*vector<double> chi_v;
chi_v.reserve(500);
vector<double> chi02_v;
chi02_v.reserve(250);
*/

for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
		cbmsim->GetEntry(i);
		if(i%1000 == 0) cout<<"Entry "<<i<<endl;

double thin1_rec=0; double thin2_rec=0; double thin_gen=0;
double thep_gen=0; double them_gen=0;
double thep_rec=0; double them_rec=0;
double yes=0;
	TVector3 pem,pem_rec;
	TVector3 pep,pep_rec;
	TVector3 pmu_in;
double Em,Ep;
double Z_ep;
double code_em, code_ep;

	for(int n = 0; n < MCTrack->GetEntries(); n++) {
	 const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
         if(MCTr->interactionID()==0 and MCTr->pdgCode()==-13) {pmu_in.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); pmu_in=pmu_in.Unit();}
	 if(MCTr->interactionID()!=0) h_int->Fill(MCTr->interactionID());
	 if(MCTr->interactionID()==5){
cout << "INIZIO" << endl;
	  if(MCTr->pdgCode()==-11) {yes++; Em=MCTr->energy();
				   pem.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); code_em=n; pem=pem.Unit();
				   them_gen=pmu_in.Angle(pem);}// h_them_gen->Fill(them_gen); cout << "them_gen " << them_gen << endl;}
          if(MCTr->pdgCode()==11) {yes++; Ep=MCTr->energy();
				   pep.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); code_ep=n; pep=pep.Unit();
				   thep_gen=pmu_in.Angle(pep);Z_ep=MCTr->startZ(); cout << "Z_ep " << Z_ep << endl;}// cout << "thep_gen " <<thep_gen << endl;}

	 }
         //if(MCTr->interactionID()==9){cout << "pdgCode "<< MCTr->pdgCode() << " and mum " << MCTr->motherID() << endl;}


	}


if(yes==2)// and Z_ep<1037 and Z_ep>1031)
{
energy_g->Fill(Em);
energy_g->Fill(Ep);
//Z_pp_gen->Fill(Z_ep);
h_them_all->Fill(them_gen);
h_thep_all->Fill(thep_gen);
cout << "___________" << endl;
cout << "code_em " << code_em << " code_ep " << code_ep << endl;
int yes_m=0; int yes_p=0; int yes_mu2=0;

vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
TVector3 thin1;
TVector3 thin2; vector<TVector3> thin2_v;
double acoplanarity;
vector<TVector3> pep_rec_v;vector<TVector3> pem_rec_v;

vector<double> chi_min_m,chi_min_p,chi_min_mu, thep_rec_vec, them_rec_vec,thmu_rec_vec;
chi_min_m.reserve(5);chi_min_mu.reserve(5);chi_min_p.reserve(5);thep_rec_vec.reserve(5);them_rec_vec.reserve(5);thmu_rec_vec.reserve(5);

double th_9;

int st0_m=0;

for(int j=0; j<tracks.size();j++)
{
 if(tracks.at(j).sector()==0) st0_m++;
}

if(st0_m==1) mult->Fill(tracks.size()-1);

for(int j=0; j<tracks.size();j++)
{

cout << "int_ID " << tracks.at(j).processIDofLinkedTrack() << " sector " << tracks.at(j).sector() << " link " << tracks.at(j).linkedTrackID() << endl;

double th_inx,th_iny;
if(tracks.at(j).processIDofLinkedTrack()==0 and tracks.at(j).sector()==0){
th_inx=tracks.at(j).xSlope();
th_iny=tracks.at(j).ySlope();
thin1.SetXYZ(th_inx,th_iny,1.0);
thin1=thin1.Unit();
thin1_rec=thin1.Theta();
}

if(tracks.at(j).processIDofLinkedTrack()==0 and tracks.at(j).sector()==1){
yes_mu2++;
th_inx=tracks.at(j).xSlope();
th_iny=tracks.at(j).ySlope();
TVector3 p;
p.SetXYZ(th_inx,th_iny,1.0);
p=p.Unit();
thin2_v.push_back(p);
chi_min_mu.push_back(tracks.at(j).chi2perDegreeOfFreedom());
thmu_rec_vec.push_back(thin1.Angle(p));
 }

if(tracks.at(j).processIDofLinkedTrack()==5 and tracks.at(j).linkedTrackID()==code_em and tracks.at(j).sector()==1)
{yes_m++;
 TVector3 p;
 p.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.);
 p=p.Unit();
 pem_rec_v.push_back(p);
 them_rec_vec.push_back(thin1.Angle(p)); //h_them_rec->Fill(them_rec);
 chi_min_m.push_back(tracks.at(j).chi2perDegreeOfFreedom());
 perc->Fill(tracks.at(j).percentageOfHitsSharedWithLinkedTrack());
}

if(tracks.at(j).processIDofLinkedTrack()==5 and tracks.at(j).linkedTrackID()==code_ep and tracks.at(j).sector()==1)
{yes_p++;
 TVector3 p;
 p.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.);
 p=p.Unit();			//h_thep_gen->Fill(thep_gen);
 pep_rec_v.push_back(p);
 thep_rec_vec.push_back(thin1.Angle(p)); //h_thep_rec->Fill(thep_rec);
 chi_min_p.push_back(tracks.at(j).chi2perDegreeOfFreedom());
 perc->Fill(tracks.at(j).percentageOfHitsSharedWithLinkedTrack());
}

if(tracks.at(j).processIDofLinkedTrack()==9){TVector3 p; p.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.);p=p.Unit();th_9=thin1.Angle(p);}

if(tracks.at(j).processIDofLinkedTrack()!=0) h_int_rec->Fill(tracks.at(j).processIDofLinkedTrack());

if(tracks.at(j).sector()==1 and tracks.size()==2 and st0_m==1) {h_part1->Fill(tracks.at(j).processIDofLinkedTrack());}
if(tracks.at(j).sector()==1 and tracks.size()==3 and st0_m==1) {h_part2->Fill(tracks.at(j).processIDofLinkedTrack());}
if(tracks.at(j).sector()==1 and tracks.size()==4 and st0_m==1) {h_part3->Fill(tracks.at(j).processIDofLinkedTrack());}
if(tracks.at(j).sector()==1 and tracks.size()>4 and st0_m==1)  {h_part_more->Fill(tracks.at(j).processIDofLinkedTrack());}
}

if( thmu_rec_vec.size()!=0){ auto it = min_element(chi_min_mu.begin(),chi_min_mu.end()); 
                                thin2_rec = thmu_rec_vec.at(std::distance(chi_min_mu.begin(), it));
                                thin2 = thin2_v.at(std::distance(chi_min_mu.begin(), it));}

double diff=thin1_rec-thin2.Theta();
h_T1->Fill(thin2_rec);
h_T2->Fill(diff);

if( them_rec_vec.size()!=0){ auto it = min_element(chi_min_m.begin(),chi_min_m.end()); them_rec = them_rec_vec.at(std::distance(chi_min_m.begin(), it));
                                pem_rec=pem_rec_v.at(std::distance(chi_min_m.begin(), it));
                                energy_r->Fill(Em);
                                h_them_gen->Fill(them_gen);
                                h_them_rec->Fill(them_rec);
                                if(yes_mu2!=0){ double dotProduct = thin2.Dot(pem_rec);
                                                TVector3 crossProduct = thin2.Cross(pem_rec);
                                                double T = thin1.Dot(crossProduct);
                                                TVector3 im= thin1.Cross(thin2);
                                                TVector3 ie= thin1.Cross(pem_rec);
                                                T = T>0? 1:(-1);
                                                acoplanarity= T*(TMath::Pi()- acos( ((im).Dot(ie))/(im.Mag()*ie.Mag()) ));

                                        }
				}

if( thep_rec_vec.size()!=0){ auto it = min_element(chi_min_p.begin(),chi_min_p.end()); thep_rec = thep_rec_vec.at(std::distance(chi_min_p.begin(), it));
                                pep_rec=pep_rec_v.at(std::distance(chi_min_p.begin(), it));
                                energy_r->Fill(Ep);
                                h_thep_gen->Fill(thep_gen);
                                h_thep_rec->Fill(thep_rec);
                                if(yes_mu2!=0){ double dotProduct = thin2.Dot(pep_rec);
                                                TVector3 crossProduct = thin2.Cross(pep_rec);
                                                double T = thin1.Dot(crossProduct);
						TVector3 im= thin1.Cross(thin2);
						TVector3 ie= thin1.Cross(pep_rec);
						T = T>0? 1:-1;
						acoplanarity= T*(TMath::Pi()- acos( ((im).Dot(ie))/(im.Mag()*ie.Mag()) ));
					}
				}


vector<MUonERecoOutputVertex> vrtx = ReconstructionOutput->reconstructedVertices();
double chi;

  for(int j=0; j<vrtx.size();j++)
        {
	 if(j==0) {chi=vrtx.at(j).chi2perDegreeOfFreedom();}
         if(vrtx.at(j).chi2perDegreeOfFreedom()<50) Z_pp->Fill(vrtx.at(j).z());
        }


if(yes_p>=1 or yes_m>=1) Z_pp_rec->Fill(Z_ep);

if((yes_p==1 and yes_m==0 and tracks.size()==3 and st0_m==1 and chi>0 and thep_rec<0.035 and abs(acoplanarity)>=0) or (yes_p==0 and yes_m==1 and tracks.size()==3 and st0_m==1 and chi>0 and them_rec<0.035 and abs(acoplanarity)>=0)){ //tracks.size()==4 and st0_m>1

vrtx_chi->Fill(chi);//chi_v.push_back(chi);
h_aco->Fill(acoplanarity);
Z_PP_gen->Fill(Z_ep);

 if(yes_p==1) th_mu_e->Fill(thep_rec,thin2_rec);
 if(yes_m==1) th_mu_e->Fill(them_rec,thin2_rec);

 if(chi<100) {
		if(yes_p==1) th_mu_e_chi->Fill(thep_rec,thin2_rec);
		if(yes_m==1) th_mu_e_chi->Fill(them_rec,thin2_rec);
		Z_PP_gen_chi->Fill(Z_ep);}
 if(abs(acoplanarity)<1) {
                if(yes_p==1) th_mu_e_aco->Fill(thep_rec,thin2_rec);
                if(yes_m==1) th_mu_e_aco->Fill(them_rec,thin2_rec);
		Z_PP_gen_aco->Fill(Z_ep);}
 if(chi<100 and abs(acoplanarity)<1) {
                if(yes_p==1) th_mu_e_both->Fill(thep_rec,thin2_rec);
                if(yes_m==1) th_mu_e_both->Fill(them_rec,thin2_rec);
		Z_PP_gen_both->Fill(Z_ep);}
 danger++;

 if(thin2_rec>0.0002){danger_02++; vrtx_chi_after02->Fill(chi); h_aco_02->Fill(acoplanarity); //chi02_v.push_back(chi);
                         Z_PP_gen02->Fill(Z_ep);
                         if(chi<100) Z_PP_gen_chi02->Fill(Z_ep);
                         if(abs(acoplanarity)<1) Z_PP_gen_aco02->Fill(Z_ep);
                         if(chi<100 and abs(acoplanarity)<1) Z_PP_gen_both02->Fill(Z_ep);
}

// vector<MUonERecoOutputAdaptiveFitterVertex> vrtx = ReconstructionOutput->adaptiveFitterVertices();
// Z_pp_rec->Fill(Z_ep);
  for(int j=0; j<tracks.size();j++)
   {
    if(tracks.at(j).processIDofLinkedTrack()==5 and (tracks.at(j).linkedTrackID()==code_ep or tracks.at(j).linkedTrackID()==code_em) and tracks.at(j).sector()==1){
	chi_e->Fill(tracks.at(j).chi2perDegreeOfFreedom()); }
    if(tracks.at(j).processIDofLinkedTrack()==5 and tracks.at(j).linkedTrackID()==code_em and tracks.at(j).sector()==1){chi_m->Fill(tracks.at(j).chi2perDegreeOfFreedom());}
	}
  }
else if((yes_p==1 and yes_m==1 and tracks.size()==3 and st0_m==1 and chi>0 and them_rec<0.035 and thep_rec<0.035 and abs(acoplanarity)>=0)){danger_ee++; vrtx_chi->Fill(chi);}
else if((yes_p==2 and yes_m==0 and tracks.size()==3 and st0_m==1 and chi>0 and thep_rec<0.035 and abs(acoplanarity)>=0) or (yes_p==0 and yes_m==2 and tracks.size()==3 and st0_m==1 and chi>0 and them_rec<0.035) or (yes_m==0 and chi>0 and yes_p==0 and tracks.size()==3 and st0_m==1 and thmu_rec_vec.size()==2  and abs(acoplanarity)>=0 and them_rec<0.035 ))
	{
	 danger_ghost++;
	 if(yes_m==2){ if(them_rec_vec.at(1)>them_rec_vec.at(0))
			{th_mu_e->Fill(them_rec_vec.at(1),them_rec_vec.at(0)); 
 if(chi<100) {
              	th_mu_e_chi->Fill(them_rec_vec.at(1),them_rec_vec.at(0));
                Z_PP_gen_chi->Fill(Z_ep);}
 if(abs(acoplanarity)<1) {
                th_mu_e_aco->Fill(them_rec_vec.at(1),them_rec_vec.at(0));
                Z_PP_gen_aco->Fill(Z_ep);}
 if(chi<100 and abs(acoplanarity)<1) {
		th_mu_e_both->Fill(them_rec_vec.at(1),them_rec_vec.at(0));
                Z_PP_gen_both->Fill(Z_ep);}

			vrtx_chi->Fill(chi);h_aco->Fill(acoplanarity);//chi_v.push_back(chi);
			 if(them_rec_vec.at(0)>0.0002){danger_ghost02++; th_gst->Fill(them_rec_vec.at(1),them_rec_vec.at(0)); 
							vrtx_chi_after02->Fill(chi);h_aco_02->Fill(acoplanarity);
                         Z_PP_gen02->Fill(Z_ep);
                         if(chi<100) Z_PP_gen_chi02->Fill(Z_ep);
                         if(abs(acoplanarity)<1) Z_PP_gen_aco02->Fill(Z_ep);
                         if(chi<100 and abs(acoplanarity)<1) Z_PP_gen_both02->Fill(Z_ep);
				} //chi02_v.push_back(chi);}
			}
			else{th_mu_e->Fill(them_rec_vec.at(0),them_rec_vec.at(1)); 
 if(chi<100) {
              	th_mu_e_chi->Fill(them_rec_vec.at(0),them_rec_vec.at(1));
                Z_PP_gen_chi->Fill(Z_ep);}
 if(abs(acoplanarity)<1) {
                th_mu_e_aco->Fill(them_rec_vec.at(0),them_rec_vec.at(1));
                Z_PP_gen_aco->Fill(Z_ep);}
 if(chi<100 and abs(acoplanarity)<1) {
                th_mu_e_both->Fill(them_rec_vec.at(0),them_rec_vec.at(1));
                Z_PP_gen_both->Fill(Z_ep);}

vrtx_chi->Fill(chi);h_aco->Fill(acoplanarity);//chi_v.push_back(chi);
				if(them_rec_vec.at(1)>0.0002){danger_ghost02++; th_gst->Fill(them_rec_vec.at(0),them_rec_vec.at(1)); 
								vrtx_chi_after02->Fill(chi);h_aco_02->Fill(acoplanarity);
                         Z_PP_gen02->Fill(Z_ep);
                         if(chi<100) Z_PP_gen_chi02->Fill(Z_ep);
                         if(abs(acoplanarity)<1) Z_PP_gen_aco02->Fill(Z_ep);
                         if(chi<100 and abs(acoplanarity)<1) Z_PP_gen_both02->Fill(Z_ep);
				} //chi02_v.push_back(chi);}
			    }
		     }
         if(yes_p==2){ if(thep_rec_vec.at(1)>thep_rec_vec.at(0))
			{th_mu_e->Fill(thep_rec_vec.at(1),thep_rec_vec.at(0)); 
 if(chi<100) {
              	th_mu_e_chi->Fill(thep_rec_vec.at(1),thep_rec_vec.at(0));
                Z_PP_gen_chi->Fill(Z_ep);}
 if(abs(acoplanarity)<1) {
                th_mu_e_aco->Fill(thep_rec_vec.at(1),thep_rec_vec.at(0));
                Z_PP_gen_aco->Fill(Z_ep);}
 if(chi<100 and abs(acoplanarity)<1) {
                th_mu_e_both->Fill(thep_rec_vec.at(1),thep_rec_vec.at(0));
                Z_PP_gen_both->Fill(Z_ep);}

vrtx_chi->Fill(chi);h_aco->Fill(acoplanarity);//chi_v.push_back(chi);
			 if(thep_rec_vec.at(0)>0.0002){danger_ghost02++; th_gst->Fill(thep_rec_vec.at(1),thep_rec_vec.at(0)); 
							vrtx_chi_after02->Fill(chi);h_aco_02->Fill(acoplanarity);
                         Z_PP_gen02->Fill(Z_ep);
                         if(chi<100) Z_PP_gen_chi02->Fill(Z_ep);
                         if(abs(acoplanarity)<1) Z_PP_gen_aco02->Fill(Z_ep);
                         if(chi<100 and abs(acoplanarity)<1) Z_PP_gen_both02->Fill(Z_ep);
				} //chi02_v.push_back(chi);}
			}
			else{th_mu_e->Fill(thep_rec_vec.at(0),thep_rec_vec.at(1));
 if(chi<100) {
                th_mu_e_chi->Fill(thep_rec_vec.at(0),thep_rec_vec.at(1));
                Z_PP_gen_chi->Fill(Z_ep);}
 if(abs(acoplanarity)<1) {
                th_mu_e_aco->Fill(thep_rec_vec.at(0),thep_rec_vec.at(1));
                Z_PP_gen_aco->Fill(Z_ep);}
 if(chi<100 and abs(acoplanarity)<1) {
                th_mu_e_both->Fill(thep_rec_vec.at(0),thep_rec_vec.at(1));
                Z_PP_gen_both->Fill(Z_ep);}

 vrtx_chi->Fill(chi);//chi_v.push_back(chi);
			 	if(thep_rec_vec.at(1)>0.0002){danger_ghost02++; th_gst->Fill(thep_rec_vec.at(0),thep_rec_vec.at(1));
								 vrtx_chi_after02->Fill(chi);h_aco_02->Fill(acoplanarity);
                         Z_PP_gen02->Fill(Z_ep);
                         if(chi<100) Z_PP_gen_chi02->Fill(Z_ep);
                         if(abs(acoplanarity)<1) Z_PP_gen_aco02->Fill(Z_ep);
                         if(chi<100 and abs(acoplanarity)<1) Z_PP_gen_both02->Fill(Z_ep);
				} //chi02_v.push_back(chi);}
			}
		     }
         if(yes_p==0 and yes_m==0 and thmu_rec_vec.size()==2){if(thmu_rec_vec.at(1)>thmu_rec_vec.at(0))
			{th_mu_e->Fill(thmu_rec_vec.at(1),thmu_rec_vec.at(0));
 if(chi<100) {
                th_mu_e_chi->Fill(thmu_rec_vec.at(1),thmu_rec_vec.at(0));
                Z_PP_gen_chi->Fill(Z_ep);}
 if(abs(acoplanarity)<1) {
                th_mu_e_aco->Fill(thmu_rec_vec.at(1),thmu_rec_vec.at(0));
                Z_PP_gen_aco->Fill(Z_ep);}
 if(chi<100 and abs(acoplanarity)<1) {
                th_mu_e_both->Fill(thmu_rec_vec.at(1),thmu_rec_vec.at(0));
                Z_PP_gen_both->Fill(Z_ep);}


 vrtx_chi->Fill(chi);h_aco->Fill(acoplanarity);//chi_v.push_back(chi);
			 if(thmu_rec_vec.at(0)>0.0002){danger_ghost02++; th_gst->Fill(thmu_rec_vec.at(1),thmu_rec_vec.at(0));
							vrtx_chi_after02->Fill(chi);h_aco_02->Fill(acoplanarity);
                         Z_PP_gen02->Fill(Z_ep);
                         if(chi<100) Z_PP_gen_chi02->Fill(Z_ep);
                         if(abs(acoplanarity)<1) Z_PP_gen_aco02->Fill(Z_ep);
                         if(chi<100 and abs(acoplanarity)<1) Z_PP_gen_both02->Fill(Z_ep);
				} //chi02_v.push_back(chi);}
			}
			else{th_mu_e->Fill(thmu_rec_vec.at(0),thmu_rec_vec.at(1));
if(chi<100) {
                th_mu_e_chi->Fill(thmu_rec_vec.at(0),thmu_rec_vec.at(1));
                Z_PP_gen_chi->Fill(Z_ep);}
 if(abs(acoplanarity)<1) {
                th_mu_e_aco->Fill(thmu_rec_vec.at(0),thmu_rec_vec.at(1));
                Z_PP_gen_aco->Fill(Z_ep);}
 if(chi<100 and abs(acoplanarity)<1) {
                th_mu_e_both->Fill(thmu_rec_vec.at(0),thmu_rec_vec.at(1));
                Z_PP_gen_both->Fill(Z_ep);}

 vrtx_chi->Fill(chi);h_aco->Fill(acoplanarity);//chi_v.push_back(chi);
				if(thmu_rec_vec.at(1)>0.0002){danger_ghost02++;th_gst->Fill(thmu_rec_vec.at(0),thmu_rec_vec.at(1));
								 vrtx_chi_after02->Fill(chi);h_aco_02->Fill(acoplanarity);
                         Z_PP_gen02->Fill(Z_ep);
                         if(chi<100) Z_PP_gen_chi02->Fill(Z_ep);
                         if(abs(acoplanarity)<1) Z_PP_gen_aco02->Fill(Z_ep);
                         if(chi<100 and abs(acoplanarity)<1) Z_PP_gen_both02->Fill(Z_ep);
				} //chi02_v.push_back(chi);}
			     }
			}
 	}
else if(yes_p==0 and yes_m==0 and thmu_rec_vec.size()==1 and tracks.size()==3 and st0_m==1 and chi>0 and th_9<0.035) {other++;
		if(th_9>thmu_rec_vec.at(0)){if(thmu_rec_vec.at(0)>0.0002)other02++;}
		else{if(th_9>0.0002){other02++;}}
	}
else if(yes_p==0 and yes_m==0 and thmu_rec_vec.size()==0 and tracks.size()==3 and st0_m==1) other0++;

}


} //end of general for

/*for(int c=0;c<chi_v.size();c++){
 vrtx_chi->Fill(chi_v.at(c),1/chi_v.size());
	}

for(int c=0;c<chi02_v.size();c++){
 vrtx_chi_after02->Fill(chi02_v.at(c),1/chi02_v.size());
        }*/

cout << "Su " << cbmsim->GetEntries() << " di muoni, hanno una PP pericolosa " << danger << " con una percentuale di " << 100*(danger/cbmsim->GetEntries()) << "%" << endl;
cout << "Su " << cbmsim->GetEntries() << " di muoni, hanno una PP pericolosa con th_mu>0.2mrad " << danger_02 << " con una percentuale di " << 100*(danger_02/cbmsim->GetEntries()) << "%" << endl;
cout << "Su " << cbmsim->GetEntries() << " di muoni, hanno 2 tracce con 1 ghost " << danger_ghost << " con una percentuale di " << 100*(danger_ghost/cbmsim->GetEntries()) << "%" << endl;
cout << "Di cui " << danger_ghost02 << " hanno un angolo >0.2mrad con una percentuale di " << 100*(danger_ghost02/cbmsim->GetEntries()) << "%" << endl;
cout << "Su " << cbmsim->GetEntries() << " di muoni, hanno 2 tracce elettroniche " << danger_ee << " con una percentuale di " << 100*(danger_ee/cbmsim->GetEntries()) << "%" << endl;
cout << "Su " << cbmsim->GetEntries() << " di muoni, hanno 1 mu+1 altro " << other << " con una percentuale di " << 100*(other/cbmsim->GetEntries()) << "%" << endl;
cout << "Di cui " << other02 << " hanno un angolo >0.2mrad con una percentuale di " << 100*(other02/cbmsim->GetEntries()) << "%" << endl;
cout << "Su " << cbmsim->GetEntries() << " di muoni, hanno 2 altri (no mu/e) " << other0 << " con una percentuale di " << 100*(other0/cbmsim->GetEntries()) << "%" << endl;


TF1 *Elastic = new TF1("Elastic","0.5109989461*0.001*((1+(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*cos(x)*cos(x))/(1-(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*cos(x)*cos(x)))",0,0.030); 
TF1 *Elastic2 = new TF1("Elastic2","asin( (sin(x)*sqrt(Elastic(x)*Elastic(x)-0.5109989461*0.001*0.5109989461*0.001))/sqrt( (160+0.5109989461*0.001-Elastic(x))*(160+0.5109989461*0.001-Elastic(x))-105.6583745 *0.001*105.6583745 *0.001 ) )",0,0.030);   
TLine *t = new TLine(0.,0.0002,0.035,0.0002);
t->SetLineWidth(2);
t->SetLineColor(kRed);

TCanvas b51("b51","b51",700,700);
b51.Divide(2,2);
b51.cd(1);

th_mu_e->SetMarkerColor(kGreen);
Elastic->SetLineColor(6);
th_mu_e->SetMarkerStyle(43);
th_mu_e->SetMarkerSize(1);
th_mu_e->Draw();
th_gst->SetMarkerColor(kGreen);
th_gst->SetMarkerStyle(43);
th_gst->SetMarkerSize(1);
th_gst->Draw("same");
Elastic2->Draw("same");
t->Draw("same");
b51.cd(2);
Elastic->SetLineColor(6);
th_mu_e_chi->SetMarkerColor(kGreen);
th_mu_e_chi->SetMarkerStyle(43);
th_mu_e_chi->SetMarkerSize(1);
th_mu_e_chi->Draw();
th_gst_chi->SetMarkerColor(kGreen);
th_gst_chi->SetMarkerStyle(43);
th_gst_chi->SetMarkerSize(1);
th_gst_chi->Draw("same");
Elastic2->Draw("same");
t->Draw("same");
b51.cd(3);
Elastic->SetLineColor(6);
th_mu_e_aco->SetMarkerColor(kGreen);
th_mu_e_aco->SetMarkerStyle(43);
th_mu_e_aco->SetMarkerSize(1);
th_mu_e_aco->Draw();
th_gst_aco->SetMarkerColor(kGreen);
th_gst_aco->SetMarkerStyle(43);
th_gst_aco->SetMarkerSize(1);
th_gst_aco->Draw("same");
Elastic2->Draw("same");
t->Draw("same");
b51.cd(4);
Elastic->SetLineColor(6);
th_mu_e_both->SetMarkerColor(kGreen);
th_mu_e_both->SetMarkerStyle(43);
th_mu_e_both->SetMarkerSize(1);
th_mu_e_both->Draw();
th_gst_both->SetMarkerColor(kGreen);
th_gst_both->SetMarkerStyle(43);
th_gst_both->SetMarkerSize(1);
th_gst_both->Draw("same");
Elastic2->Draw("same");
t->Draw("same");
b51.SaveAs("2Dplots_bkg.pdf");

TCanvas b5("b5","b5",700,700);
b5.Divide(2,2);
b5.cd(1);
Z_PP_gen->Draw("hist");
gPad->SetLogy();
b5.cd(2);
Z_PP_gen_chi->Draw("hist");
gPad->SetLogy();
b5.cd(3);
Z_PP_gen_aco->Draw("hist");
gPad->SetLogy();
b5.cd(4);
Z_PP_gen_both->Draw("hist");
gPad->SetLogy();
//Z_PP_gen->SaveAs("Z_PP_gen.root");
b5.SaveAs("Z_PP_cut.pdf");

TCanvas b52("b52","b52",700,700);
b52.Divide(2,2);
b52.cd(1);
Z_PP_gen02->Draw("hist");
gPad->SetLogy();
b52.cd(2);
Z_PP_gen_chi02->Draw("hist");
gPad->SetLogy();
b52.cd(3);
Z_PP_gen_aco02->Draw("hist");
gPad->SetLogy();
b52.cd(4);
Z_PP_gen_both02->Draw("hist");
gPad->SetLogy();
//Z_PP_gen_both02->SaveAs("Z_PP_gen_both02.root");
b52.SaveAs("Z_PP_cut_02.pdf");


/*
TCanvas b7("b7","b7",700,700);
b7.Divide(2,2);
b7.cd(1);
gPad->SetLogy();
int nBins1 = vrtx_chi->GetNbinsX()+1;
TH1F *h1 = (TH1F*)(vrtx_chi->Clone("h1"));
h1->Scale(1/h1->Integral(0,nBins1));
h1->SetLineColor(kRed);
h1->Draw("hist");
gPad->SetLogy();
b7.cd(2);
int nBins2 = vrtx_chi_after02->GetNbinsX()+1;
TH1F *h2 = (TH1F*)(vrtx_chi_after02->Clone("h2"));
h2->Scale(1/h2->Integral(0,nBins2));
h2->SetLineColor(kRed);
h2->Draw("hist");
gPad->SetLogy();
b7.cd(3);
    int nx    = h1->GetNbinsX()+1;
    double *xbins= new double[nx+1];
    for (int i=0;i<nx;i++)
        xbins[i]=h1->GetBinLowEdge(i+1);
    xbins[nx]=xbins[nx-1]+h1->GetBinWidth(nx);
    //book a temporary histogram having extra bins for overflows
    TH1F *htmp = new TH1F(h1->GetName(), h1->GetTitle(), nx, xbins);
    htmp->Sumw2();
    //fill the new histogram including the overflows
    for (int i=1; i<=nx; i++) {
        htmp->SetBinContent(htmp->FindBin(htmp->GetBinCenter(i)),h1->GetBinContent(i));
        htmp->SetBinError(htmp->FindBin(htmp->GetBinCenter(i)),h1->GetBinError(i));
    }
    htmp->SetBinContent(htmp->FindBin(h1->GetBinLowEdge(1)-1), h1->GetBinContent(0));
    htmp->SetBinError(htmp->FindBin(h1->GetBinLowEdge(1)-1), h1->GetBinError(0));
    // Restore the number of entries
    htmp->SetEntries(h1->GetEffectiveEntries());

htmp->Draw("hist");
gPad->SetLogy();

b7.cd(4);
    int nx1    = h2->GetNbinsX()+1;
    double *xbins1= new double[nx1+1];
    for (int i=0;i<nx1;i++)
        xbins1[i]=h1->GetBinLowEdge(i+1);
    xbins1[nx1]=xbins1[nx1-1]+h1->GetBinWidth(nx1);
    //book a temporary histogram having extra bins for overflows
    TH1F *htmp2 = new TH1F(h2->GetName(), h2->GetTitle(), nx1, xbins1);
    htmp2->Sumw2();
    //fill the new histogram including the overflows
    for (int i=1; i<=nx1; i++) {
        htmp2->SetBinContent(htmp2->FindBin(htmp2->GetBinCenter(i)),h2->GetBinContent(i));
        htmp2->SetBinError(htmp2->FindBin(htmp2->GetBinCenter(i)),h2->GetBinError(i));
    }
    htmp2->SetBinContent(htmp2->FindBin(h2->GetBinLowEdge(1)-1), h2->GetBinContent(0));
    htmp2->SetBinError(htmp2->FindBin(h2->GetBinLowEdge(1)-1), h2->GetBinError(0));
    // Restore the number of entries
    htmp2->SetEntries(h2->GetEffectiveEntries());
htmp2->SetLineColor(kRed);
htmp2->Draw("hist");
gPad->SetLogy();

b7.SaveAs("vrtx_chi_sig_PP.pdf");
h1->SaveAs("vrtx_chi_PP_norm.root");
htmp->SaveAs("vrtx_chi_PP_norm_ovrflw.root");
h2->SaveAs("vrtx_chi_PP_norm_after02.root");
htmp2->SaveAs("vrtx_chi_PP_norm_after02_ovrflw.root");

TCanvas ac("ac","ac",700,700);
ac.Divide(1,2);
ac.cd(1);
TH1F *h_acoN = (TH1F*)(h_aco->Clone("h_acoN"));
h_acoN->Scale(1/h_acoN->Integral());
h_acoN->Draw("hist");
ac.cd(2);
TH1F *h_acoN_02 = (TH1F*)(h_aco_02->Clone("h_acoN_02"));
h_acoN_02->Scale(1/h_acoN_02->Integral());
h_acoN_02->Draw("hist");

ac.SaveAs("h_aco_PP.pdf");
h_acoN->SaveAs("h_aco_PP.root");
h_acoN_02->SaveAs("h_aco_PP_after02.root");

TCanvas b2("b2","b2",700,700);
b2.Divide(2,2);
b2.cd(1);
h_them_all->Draw("hist");
h_them_gen->SetLineColor(kOrange);
h_them_gen->Draw("hist same");
h_them_rec->SetLineColor(kRed);
h_them_rec->Draw("hist same");
gPad->SetLogy();
b2.cd(2);
h_thep_all->Draw("hist");
h_thep_gen->SetLineColor(kOrange);
h_thep_gen->Draw("hist same");
h_thep_rec->SetLineColor(kRed);
h_thep_rec->Draw("hist same");
gPad->SetLogy();
b2.cd(3);
h_T1->Draw("hist");
b2.cd(4);
h_T2->Draw("hist");
b2.SaveAs("theta.pdf");

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
b4.SaveAs("energy.pdf");


TCanvas b3("b3","b3",700,700);
b3.Divide(1,2);
b3.cd(1);
h_int->Draw("hist");
b3.cd(2);
h_int_rec->SetLineColor(kRed);
h_int_rec->Draw("hist");
gPad->SetLogy();
b3.SaveAs("int.pdf");



b5.Divide(2,2);
b5.cd(1);
th_mu_e->Draw("COLZ");
b5.cd(2);
th_gst->Draw("COLZ");
b5.cd(3);
Z_pp->Draw("hist");
b5.cd(4);
Z_pp_gen->Draw("hist");
Z_pp_rec->SetLineColor(kRed);
Z_pp_rec->Draw("hist same");
b5.SaveAs("th_mu_e.pdf");

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
b6.SaveAs("mult.pdf");

TCanvas b7("b7","b7",700,700);
vrtx_chi->Draw("hist");
gPad->SetLogy();
b7.SaveAs("vrtx_chi_PP.pdf");
*/
}
