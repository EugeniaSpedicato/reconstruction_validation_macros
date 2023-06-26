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
cbmsim->Add("TRPP_minbias_1M_firstSample.root");
cbmsim->Add("TRPP_minbias_1M_secondSample.root");


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

TH1D *vrtx_chi=new TH1D("chie","Chi2 per DOF of the kinematic vrtx for PP",500,0,6000);

TH2D *th_mu_e=new TH2D("th_mu_em","All events: angle muon and electron/positron from PP",700,0,0.07,50,0,0.005);
TH2D *th_gst=new TH2D("th_gst","Remaining events: angle muon and electron/positron from PP",700,0,0.07,50,0,0.005);


TH1D * Z_pp=new TH1D("Zpp" , "Z position with adaptive fitter for PP", 400,930,1150);
TH1D * Z_pp_gen=new TH1D("Zpp_gen" , "Generated Z position PP particles", 400,930,1150);
TH1D * Z_pp_rec=new TH1D("Zpp_rec" , "Generated Z position PP particles when #reco==3", 400,930,1150);

TH1D* mult=new TH1D("mul","multiplicity of reconstructed tracks in second station when PP happens", 20,0,20);
TH1D *h_part1=new TH1D("p1","reconstructed particle ID second station's multiplicity=1", 50,0,50);
TH1D *h_part2=new TH1D("p2","reconstructed particle ID second station's multiplicity=2", 50,0,50);
TH1D *h_part3=new TH1D("p3","reconstructed particle ID second station's multiplicity=3", 50,0,50);
TH1D *h_part_more=new TH1D("pm","reconstructed particle ID multiplicity>3", 50,0,50);

double danger=0;
double danger_02=0;
double danger_ghost=0;
double danger_ee=0;
double other=0;
double other0=0;
double danger_ghost02=0;
double other02=0;
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
				   them_gen=acos(pmu_in.Dot(pem));}// h_them_gen->Fill(them_gen); cout << "them_gen " << them_gen << endl;}
          if(MCTr->pdgCode()==11) {yes++; Ep=MCTr->energy();
				   pep.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); code_ep=n; pep=pep.Unit();
				   thep_gen=acos(pmu_in.Dot(pep));Z_ep=MCTr->startZ();}// cout << "thep_gen " <<thep_gen << endl;}

	 }
         //if(MCTr->interactionID()==9){cout << "pdgCode "<< MCTr->pdgCode() << " and mum " << MCTr->motherID() << endl;}


	}


if(yes==2)//and Z_ep<1037 and Z_ep>1031)
{
energy_g->Fill(Em);
energy_g->Fill(Ep);
Z_pp_gen->Fill(Z_ep);
h_them_all->Fill(them_gen);
h_thep_all->Fill(thep_gen);
cout << "___________" << endl;
cout << "code_em " << code_em << " code_ep " << code_ep << endl;
int yes_m=0; int yes_p=0; int yes_mu2=0;

vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
TVector3 thin1;
TVector3 thin2;

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
thin2.SetXYZ(th_inx,th_iny,1.0);
thin2=thin2.Unit();
chi_min_mu.push_back(tracks.at(j).chi2perDegreeOfFreedom());
thmu_rec_vec.push_back(acos(thin1.Dot(thin2)));
 }

if(tracks.at(j).processIDofLinkedTrack()==5 and tracks.at(j).linkedTrackID()==code_em and tracks.at(j).sector()==1)
{yes_m++;
 pem_rec.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.);
 pem_rec=pem_rec.Unit();			//h_them_gen->Fill(them_gen);
 them_rec_vec.push_back(acos(thin1.Dot(pem_rec))); //h_them_rec->Fill(them_rec);
 chi_min_m.push_back(tracks.at(j).chi2perDegreeOfFreedom());
 perc->Fill(tracks.at(j).percentageOfHitsSharedWithLinkedTrack());
}

if(tracks.at(j).processIDofLinkedTrack()==5 and tracks.at(j).linkedTrackID()==code_ep and tracks.at(j).sector()==1)
{yes_p++;
 pep_rec.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.);
 pep_rec=pep_rec.Unit();			//h_thep_gen->Fill(thep_gen);
 thep_rec_vec.push_back(acos(thin1.Dot(pep_rec))); //h_thep_rec->Fill(thep_rec);
 chi_min_p.push_back(tracks.at(j).chi2perDegreeOfFreedom());
 perc->Fill(tracks.at(j).percentageOfHitsSharedWithLinkedTrack());
}

if(tracks.at(j).processIDofLinkedTrack()==9){TVector3 p; p.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.);p=p.Unit();th_9=acos(thin1.Dot(p));}

if(tracks.at(j).processIDofLinkedTrack()!=0) h_int_rec->Fill(tracks.at(j).processIDofLinkedTrack());

if(tracks.at(j).sector()==1 and tracks.size()==2 and st0_m==1) {h_part1->Fill(tracks.at(j).processIDofLinkedTrack());}
if(tracks.at(j).sector()==1 and tracks.size()==3 and st0_m==1) {h_part2->Fill(tracks.at(j).processIDofLinkedTrack());}
if(tracks.at(j).sector()==1 and tracks.size()==4 and st0_m==1) {h_part3->Fill(tracks.at(j).processIDofLinkedTrack());}
if(tracks.at(j).sector()==1 and tracks.size()>4 and st0_m==1)  {h_part_more->Fill(tracks.at(j).processIDofLinkedTrack());}
}

if( thmu_rec_vec.size()!=0){ auto it = min_element(chi_min_mu.begin(),chi_min_mu.end()); thin2_rec = thmu_rec_vec.at(std::distance(chi_min_mu.begin(), it));}

double diff=thin1_rec-thin2.Theta();
h_T1->Fill(thin2_rec);
h_T2->Fill(diff);

if( them_rec_vec.size()!=0){ auto it = min_element(chi_min_m.begin(),chi_min_m.end()); them_rec = them_rec_vec.at(std::distance(chi_min_m.begin(), it));
				energy_r->Fill(Em);
				h_them_gen->Fill(them_gen);
				h_them_rec->Fill(them_rec);}
if( thep_rec_vec.size()!=0){ auto it = min_element(chi_min_p.begin(),chi_min_p.end()); thep_rec = thep_rec_vec.at(std::distance(chi_min_p.begin(), it));
                                energy_r->Fill(Ep);
				h_thep_gen->Fill(thep_gen);
				h_thep_rec->Fill(thep_rec);}



vector<MUonERecoOutputVertex> vrtx = ReconstructionOutput->reconstructedVertices();
double chi;
  for(int j=0; j<vrtx.size();j++)
        {
	 if(j==0) {vrtx_chi->Fill(vrtx.at(j).chi2perDegreeOfFreedom()); chi=vrtx.at(j).chi2perDegreeOfFreedom();}
         if(vrtx.at(j).chi2perDegreeOfFreedom()<50) Z_pp->Fill(vrtx.at(j).z());
        }


if(yes_p>=1 or yes_m>=1) Z_pp_rec->Fill(Z_ep);

if((yes_p==1 and yes_m==0 and tracks.size()==3 and st0_m==1 and chi<100 and thep_rec<0.035) or (yes_p==0 and yes_m==1 and tracks.size()==3 and st0_m==1 and chi<100 and them_rec<0.035)){ //tracks.size()==4 and st0_m>1

 if(yes_p==1) th_mu_e->Fill(thep_rec,thin2_rec);
 if(yes_m==1) th_mu_e->Fill(them_rec,thin2_rec);

 danger++;
 if(thin2_rec>0.0002){danger_02++;
 if(yes_p==1) th_gst->Fill(thep_rec,thin2_rec);
 if(yes_m==1) th_gst->Fill(them_rec,thin2_rec);}

// vector<MUonERecoOutputAdaptiveFitterVertex> vrtx = ReconstructionOutput->adaptiveFitterVertices();
// Z_pp_rec->Fill(Z_ep);
  for(int j=0; j<tracks.size();j++)
   {
    if(tracks.at(j).processIDofLinkedTrack()==5 and (tracks.at(j).linkedTrackID()==code_ep or tracks.at(j).linkedTrackID()==code_em) and tracks.at(j).sector()==1){
	chi_e->Fill(tracks.at(j).chi2perDegreeOfFreedom()); }
    if(tracks.at(j).processIDofLinkedTrack()==5 and tracks.at(j).linkedTrackID()==code_em and tracks.at(j).sector()==1){chi_m->Fill(tracks.at(j).chi2perDegreeOfFreedom());}
	}
  }
else if((yes_p==1 and yes_m==1 and tracks.size()==3 and st0_m==1 and chi<100 and them_rec<0.035 and thep_rec<0.035)) danger_ee++;
else if((yes_p==2 and yes_m==0 and tracks.size()==3 and st0_m==1 and chi<100 and thep_rec<0.035) or (yes_p==0 and yes_m==2 and tracks.size()==3 and st0_m==1 and chi<100 and them_rec<0.035) or (yes_m==0 and chi<100 and yes_p==0 and tracks.size()==3 and st0_m==1 and thmu_rec_vec.size()==2))
	{
	 danger_ghost++;
	 if(yes_m==2){ if(them_rec_vec.at(1)>them_rec_vec.at(0))
			{th_mu_e->Fill(them_rec_vec.at(1),them_rec_vec.at(0));
			 if(them_rec_vec.at(0)>0.0002){danger_ghost02++; th_gst->Fill(them_rec_vec.at(1),them_rec_vec.at(0));}
			}
			else{th_mu_e->Fill(them_rec_vec.at(0),them_rec_vec.at(1));
				if(them_rec_vec.at(1)>0.0002){danger_ghost02++; th_gst->Fill(them_rec_vec.at(0),them_rec_vec.at(1));}
			    }
		     }
         if(yes_p==2){ if(thep_rec_vec.at(1)>thep_rec_vec.at(0))
			{
			 if(thep_rec_vec.at(0)>0.0002){danger_ghost02++; th_gst->Fill(thep_rec_vec.at(1),thep_rec_vec.at(0));}
			}
			else{
			 	if(thep_rec_vec.at(1)>0.0002){danger_ghost02++; th_gst->Fill(thep_rec_vec.at(0),thep_rec_vec.at(1));}
			}
		     }
         if(yes_p==0 and yes_m==0 and thmu_rec_vec.size()==2){if(thmu_rec_vec.at(1)>thmu_rec_vec.at(0))
			{th_mu_e->Fill(thmu_rec_vec.at(1),thmu_rec_vec.at(0));
			 if(thmu_rec_vec.at(0)>0.0002){danger_ghost02++; th_gst->Fill(thmu_rec_vec.at(1),thmu_rec_vec.at(0));}
			}
			else{th_mu_e->Fill(thmu_rec_vec.at(0),thmu_rec_vec.at(1));
				if(thmu_rec_vec.at(1)>0.0002){danger_ghost02++;th_gst->Fill(thmu_rec_vec.at(0),thmu_rec_vec.at(1));}
			     }
			}
 	}
else if(yes_p==0 and yes_m==0 and thmu_rec_vec.size()==1 and tracks.size()==3 and st0_m==1 and chi<100 and th_9<0.035) {other++;
		if(th_9>thmu_rec_vec.at(0)){if(thmu_rec_vec.at(0)>0.0002)other02++;}
		else{if(th_9>0.0002){other02++;}}
	}
else if(yes_p==0 and yes_m==0 and thmu_rec_vec.size()==0 and tracks.size()==3 and st0_m==1) other0++;

}


} //end of general for


cout << "Su " << cbmsim->GetEntries() << " di muoni, hanno una PP pericolosa " << danger << " con una percentuale di " << 100*(danger/cbmsim->GetEntries()) << "%" << endl;
cout << "Su " << cbmsim->GetEntries() << " di muoni, hanno una PP pericolosa con th_mu>0.2mrad " << danger_02 << " con una percentuale di " << 100*(danger_02/cbmsim->GetEntries()) << "%" << endl;
cout << "Su " << cbmsim->GetEntries() << " di muoni, hanno 2 tracce con 1 ghost " << danger_ghost << " con una percentuale di " << 100*(danger_ghost/cbmsim->GetEntries()) << "%" << endl;
cout << "Di cui " << danger_ghost02 << " hanno un angolo >0.2mrad con una percentuale di " << 100*(danger_ghost02/cbmsim->GetEntries()) << "%" << endl;
cout << "Su " << cbmsim->GetEntries() << " di muoni, hanno 2 tracce elettroniche " << danger_ee << " con una percentuale di " << 100*(danger_ee/cbmsim->GetEntries()) << "%" << endl;
cout << "Su " << cbmsim->GetEntries() << " di muoni, hanno 1 mu+1 altro " << other << " con una percentuale di " << 100*(other/cbmsim->GetEntries()) << "%" << endl;
cout << "Di cui " << other02 << " hanno un angolo >0.2mrad con una percentuale di " << 100*(other02/cbmsim->GetEntries()) << "%" << endl;
cout << "Su " << cbmsim->GetEntries() << " di muoni, hanno 2 altri (no mu/e) " << other0 << " con una percentuale di " << 100*(other0/cbmsim->GetEntries()) << "%" << endl;
/*
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
*/

TCanvas b5("b5","b5",700,700);
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
/*
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
