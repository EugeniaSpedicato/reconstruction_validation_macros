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
cbmsim->Add("TRPP_minbias_offset/TRPP_minbias_1M_firstSample.root");
cbmsim->Add("TRPP_minbias_offset/TRPP_minbias_1M_secondSample.root");
cbmsim->Add("TRPP_minbias_offset/TRPP_minbias_1M_thirdSample.root");

        MUonERecoOutput *ReconstructionOutput = 0;
        TClonesArray *MCTrack = 0;

        cbmsim->SetBranchAddress("MCTrack", &MCTrack);
        cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);


TH1D *h_T1=new TH1D("T1","Theta of muon_out wrt muon_in from sig (GeV)",100,0.,0.0005);
TH1D *h_T2=new TH1D("T2","Theta muon_in - theta muon_out from sig (GeV)",200,-0.0005,0.0005);
TH1D *h_int= new TH1D("int","minimum bias interaction per particle generated",50,0,50);
TH1D *h_int_rec= new TH1D("intr","minimum bias interaction per particle reconstructed",50,0,50);

TH1D *h_them_all=new TH1D("t2_all","Generated angle of all the gen muons from sig wrt incoming muon (rad)",1143,0,0.04);
TH1D *h_thep_all=new TH1D("t3_all","Generated angle of all the gen electron from sig wrt incoming muon (rad)",1143,0,0.04);

TH1D *h_them_gen=new TH1D("t2_gen","Generated angle of the reco muon from sig wrt incoming muon (rad)",1143,0,0.04);
TH1D *h_thep_gen=new TH1D("t3_gen","Generated angle of the reco electron from sig wrt incoming muon (rad)",1143,0,0.04);
TH1D *h_them_rec=new TH1D("t2_rec","Reconstructed ngle of the muon from sig wrt incoming muon (rad)",1143,0,0.04);
TH1D *h_thep_rec=new TH1D("t3_rec","Reconstructed angle of the electron from sig wrt incoming muon (rad)",1143,0,0.04);


TH1D *energy_r=new TH1D("Er","Energy of the reconstructed electron (GeV)",100,0,10);
TH1D *energy_g=new TH1D("Eg","Energy of all the electrons (GeV)",100,0,10);

TH1D *chi_e=new TH1D("chie","Chi2 per DOF muon or electron tracks (#trk==3)",120,0,60);
TH1D *chi_m=new TH1D("chim","Chi2 per DOF muon track (#trk==3",120,0,60);
TH1D *perc=new TH1D("perc","Quality muon and electron tracks",100,0,100);

TH1D *vrtx_chi=new TH1D("chipp","Chi2 per DOF of the kinematic vrtx for PP",40,0,200);
TH1D *vrtx_chi_after02=new TH1D("chiapp","Chi2 per DOF of the kinematic vrtx for sig after th_mu>0.2mrad cut",40,0,200);

TH2D *th_mu_e=new TH2D("th_mu_em","Angle muon and muon/electron from sig",700,0,0.07,50,0,0.005);
TH2D *th_gst=new TH2D("th_mu_em","Angle two particles where we have 1 ghost 0r 1 particle !=sig",700,0,0.07,50,0,0.005);


TH1D * Z_sig=new TH1D("Zsig" , "Z position with adaptive fitter for sig", 400,930,1150);
TH1D * Z_sig_gen=new TH1D("Zsig_gen" , "Generated Z position sig particles", 400,930,1150);
TH1D * Z_sig_rec=new TH1D("Zsig_rec" , "Generated Z position sig particles when #reco==3", 400,930,1150);

TH1D* mult=new TH1D("mul","multiplicity of reconstructed tracks in second station when sig hasigens", 20,0,20);
TH1D *h_part1=new TH1D("p1","reconstructed particle ID second station's multiplicity=1", 50,0,50);
TH1D *h_part2=new TH1D("p2","reconstructed particle ID second station's multiplicity=2", 50,0,50);
TH1D *h_part3=new TH1D("p3","reconstructed particle ID second station's multiplicity=3", 50,0,50);
TH1D *h_part_more=new TH1D("pm","reconstructed particle ID multiplicity>3", 50,0,50);

TH1D* h_aco=new TH1D("aco","Acoplanarity of muone+electron from sig minbias",600,-3.14,3.14);
TH1D* h_aco_02=new TH1D("aco02","Acoplanarity of muone+electron from sig minbias after th_mu>0.2mrad",600,-3.14,3.14);

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
int yes5=0;
int yes9=0;
	for(int n = 0; n < MCTrack->GetEntries(); n++) {
	 const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));

         if(MCTr->interactionID()==0 and MCTr->pdgCode()==-13) {pmu_in.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); pmu_in=pmu_in.Unit();}
	 if(MCTr->interactionID()==5) yes5=1;//h_int->Fill(MCTr->interactionID());
	 if(MCTr->interactionID()==9){yes9=1;
cout << "INIZIO" << endl;
          if(MCTr->pdgCode()==11) {yes++; Ep=MCTr->energy();
				   pep.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); code_ep=n; pep=pep.Unit();
				   thep_gen=pmu_in.Angle(pep);Z_ep=MCTr->startZ();}// cout << "thep_gen " <<thep_gen << endl;}

	 }

	}

if(yes5==1)h_int->Fill(5);
if(yes9==1)h_int->Fill(9);

if(yes==1 and Z_ep<1037 and Z_ep>1031)// and thep_gen<=0.035)// and Z_ep<1037 and Z_ep>1031)
{
energy_g->Fill(Ep);
Z_sig_gen->Fill(Z_ep);
h_thep_all->Fill(thep_gen);

cout << "___________" << endl;
cout << "code_em " << code_em << " code_ep " << code_ep << endl;
int yes_m=0; int yes_p=0; int yes_mu2=0;

vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
TVector3 thin1;
TVector3 thin2;  vector<TVector3> thin2_v;

double acoplanarity;
vector<TVector3> pep_rec_v;

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


if(tracks.at(j).processIDofLinkedTrack()==9 and tracks.at(j).linkedTrackID()==code_ep and tracks.at(j).sector()==1)
{yes_p++;
 TVector3 p;
 p.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.);
 p=p.Unit();                    //h_thep_gen->Fill(thep_gen);
 pep_rec_v.push_back(p);
 thep_rec_vec.push_back(thin1.Angle(p)); //h_thep_rec->Fill(thep_rec);
 chi_min_p.push_back(tracks.at(j).chi2perDegreeOfFreedom());
 perc->Fill(tracks.at(j).percentageOfHitsSharedWithLinkedTrack());
}

if(tracks.at(j).processIDofLinkedTrack()==5){TVector3 p; p.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.);p=p.Unit();th_9=thin1.Angle(p);}

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

if( thep_rec_vec.size()!=0){
 auto it = min_element(chi_min_p.begin(),chi_min_p.end()); thep_rec = thep_rec_vec.at(std::distance(chi_min_p.begin(), it));
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

if(thep_rec_vec.size()!=0){

                                 energy_r->Fill(Ep);
                                h_thep_gen->Fill(thep_gen);
                                h_thep_rec->Fill(thep_rec);
}
vector<MUonERecoOutputVertex> vrtx = ReconstructionOutput->reconstructedVertices();
double chi;
  for(int j=0; j<vrtx.size();j++)
        {
	 if(j==0) {chi=vrtx.at(j).chi2perDegreeOfFreedom();}
         if(vrtx.at(j).chi2perDegreeOfFreedom()<50) Z_sig->Fill(vrtx.at(j).z());
        }


if(yes_p>=1) Z_sig_rec->Fill(Z_ep);

if(yes_p==1 and yes_mu2==1 and tracks.size()==3 and st0_m==1 and chi<100 and thep_rec<0.035 and abs(acoplanarity)<1){ //tracks.size()==4 and st0_m>1 and thep_rec<0.035
 vrtx_chi->Fill(chi);
 h_aco->Fill(acoplanarity);

 danger++;
 if(thin2_rec>0.0002){danger_02++; vrtx_chi_after02->Fill(chi); h_aco_02->Fill(acoplanarity);}

 if(yes_p==1) th_mu_e->Fill(thep_rec,thin2_rec);

// vector<MUonERecoOutputAdaptiveFitterVertex> vrtx = ReconstructionOutput->adaptiveFitterVertices();
// Z_sig_rec->Fill(Z_ep);
  for(int j=0; j<tracks.size();j++)
   {
    if(tracks.at(j).processIDofLinkedTrack()==9 and (tracks.at(j).linkedTrackID()==code_ep or tracks.at(j).linkedTrackID()==code_em) and tracks.at(j).sector()==1){
	chi_e->Fill(tracks.at(j).chi2perDegreeOfFreedom()); }
	}
  }
else if((yes_p==2 and yes_mu2==0 and tracks.size()==3 and st0_m==1 and chi<100 and thep_rec<0.035 and abs(acoplanarity)<1) or (yes_p==0 and yes_mu2==2 and tracks.size()==3 and st0_m==1 and chi<100) or (yes_m==0 and chi<100 and yes_p==0 and tracks.size()==3 and st0_m==1 and thmu_rec_vec.size()==2 and abs(acoplanarity)<1))
	{
	 danger_ghost++;
         if(yes_p==2){ if(thep_rec_vec.at(1)>thep_rec_vec.at(0))
			{th_gst->Fill(thep_rec_vec.at(1),thep_rec_vec.at(0));  vrtx_chi->Fill(chi); h_aco->Fill(acoplanarity);
				if(thep_rec_vec.at(0)>0.0002){danger_ghost02++;vrtx_chi_after02->Fill(chi); h_aco_02->Fill(acoplanarity);}
			 }
			else{th_gst->Fill(thep_rec_vec.at(0),thep_rec_vec.at(1));  vrtx_chi->Fill(chi);  h_aco->Fill(acoplanarity);
				if(thep_rec_vec.at(1)>0.0002){danger_ghost02++;vrtx_chi_after02->Fill(chi); h_aco_02->Fill(acoplanarity);}
			 }
			}
         if(yes_p==0 and thmu_rec_vec.size()==2){if(thmu_rec_vec.at(1)>thmu_rec_vec.at(0))
			{th_gst->Fill(thmu_rec_vec.at(1),thmu_rec_vec.at(0)); vrtx_chi->Fill(chi);  h_aco->Fill(acoplanarity);
				if(thmu_rec_vec.at(0)>0.0002){danger_ghost02++;vrtx_chi_after02->Fill(chi); h_aco_02->Fill(acoplanarity);}
			 }
			else{th_gst->Fill(thmu_rec_vec.at(0),thmu_rec_vec.at(1)); vrtx_chi->Fill(chi); h_aco->Fill(acoplanarity);
				if(thmu_rec_vec.at(1)>0.0002){danger_ghost02++;vrtx_chi_after02->Fill(chi); h_aco_02->Fill(acoplanarity);}
			 }
			}
	 	}
else if(yes_p==0 and thmu_rec_vec.size()==1 and tracks.size()==3 and st0_m==1 and chi<100) {other++;
		if(th_9>thmu_rec_vec.at(0)){th_gst->Fill(th_9,thmu_rec_vec.at(0)); if(thmu_rec_vec.at(0)>0.0002)other02++;}
		else{th_gst->Fill(thmu_rec_vec.at(0),th_9);if(th_9>0.0002)other02++;}
		}

 }//end yes==2

} //end of general for


cout << "Su " << cbmsim->GetEntries() << " di muoni, hanno una sig pericolosa " << danger << " con una percentuale di " << 100*(danger/cbmsim->GetEntries()) << "%" << endl;
cout << "Su " << cbmsim->GetEntries() << " di muoni, hanno una sig pericolosa con th_mu>0.2mrad " << danger_02 << " con una percentuale di " << 100*(danger_02/cbmsim->GetEntries()) << "%" << endl;
cout << "Su " << cbmsim->GetEntries() << " di muoni, hanno 2 tracce con 1 ghost " << danger_ghost << " con una percentuale di " << 100*(danger_ghost/cbmsim->GetEntries()) << "%" << endl;
cout << "Di cui " << danger_ghost02 << " hanno un angolo >0.2mrad con una percentuale di " << 100*(danger_ghost02/cbmsim->GetEntries()) << "%" << endl;
cout << "Su " << cbmsim->GetEntries() << " di muoni, hanno 1 mu+1 altro " << other << " con una percentuale di " << 100*(other/cbmsim->GetEntries()) << "%" << endl;
cout << "Di cui " << other02 << " hanno un angolo >0.2mrad con una percentuale di " << 100*(other02/cbmsim->GetEntries()) << "%" << endl;
cout << "Su " << cbmsim->GetEntries() << " di muoni, hanno 2 altri (no mu/e) " << other0 << " con una percentuale di " << 100*(other0/cbmsim->GetEntries()) << "%" << endl;

/*
TCanvas b7("b7","b7",700,700);
b7.Divide(2,2);
b7.cd(1);
gPad->SetLogy();
int nBins1 = vrtx_chi->GetNbinsX()+1;
TH1F *h1_sm = (TH1F*)(vrtx_chi->Clone("h1_sm"));
h1_sm->Scale(1/h1_sm->Integral(0,nBins1));
h1_sm->SetLineColor(kRed);
h1_sm->Draw("hist");
gPad->SetLogy();
b7.cd(2);
int nBins2 = vrtx_chi_after02->GetNbinsX()+1;
TH1F *h2_sm = (TH1F*)(vrtx_chi_after02->Clone("h2_sm"));
h2_sm->Scale(1/h2_sm->Integral(0,nBins2));
h2_sm->SetLineColor(kRed);
h2_sm->Draw("hist");
gPad->SetLogy();
b7.cd(3);
    int nx    = h1_sm->GetNbinsX()+1;
    double *xbins= new double[nx+1];
    for (int i=0;i<nx;i++)
        xbins[i]=h1_sm->GetBinLowEdge(i+1);
    xbins[nx]=xbins[nx-1]+h1_sm->GetBinWidth(nx);
    //book a temporary histogram having extra bins for overflows
    TH1F *htmp_sm = new TH1F("h1_sm_o", h1_sm->GetTitle(), nx, xbins);
    htmp_sm->Sumw2();
    //fill the new histogram including the overflows
    for (int i=1; i<=nx; i++) {
        htmp_sm->SetBinContent(htmp_sm->FindBin(htmp_sm->GetBinCenter(i)),h1_sm->GetBinContent(i));
        htmp_sm->SetBinError(htmp_sm->FindBin(htmp_sm->GetBinCenter(i)),h1_sm->GetBinError(i));
    }
    htmp_sm->SetBinContent(htmp_sm->FindBin(h1_sm->GetBinLowEdge(1)-1), h1_sm->GetBinContent(0));
    htmp_sm->SetBinError(htmp_sm->FindBin(h1_sm->GetBinLowEdge(1)-1), h1_sm->GetBinError(0));
    // Restore the number of entries
    htmp_sm->SetEntries(h1_sm->GetEffectiveEntries());
htmp_sm->Draw("hist");
gPad->SetLogy();

b7.cd(4);
    int nx1    = h2_sm->GetNbinsX()+1;
    double *xbins1= new double[nx1+1];
    for (int i=0;i<nx1;i++)
        xbins1[i]=h1_sm->GetBinLowEdge(i+1);
    xbins1[nx1]=xbins1[nx1-1]+h1_sm->GetBinWidth(nx1);
    //book a temporary histogram having extra bins for overflows
    TH1F *htmp_sm2 = new TH1F("htmp_sm2_o", h2_sm->GetTitle(), nx1, xbins1);
    htmp_sm2->Sumw2();
    //fill the new histogram including the overflows
    for (int i=1; i<=nx1; i++) {
        htmp_sm2->SetBinContent(htmp_sm2->FindBin(htmp_sm2->GetBinCenter(i)),h2_sm->GetBinContent(i));
        htmp_sm2->SetBinError(htmp_sm2->FindBin(htmp_sm2->GetBinCenter(i)),h2_sm->GetBinError(i));
    }
    htmp_sm2->SetBinContent(htmp_sm2->FindBin(h2_sm->GetBinLowEdge(1)-1), h2_sm->GetBinContent(0));
    htmp_sm2->SetBinError(htmp_sm2->FindBin(h2_sm->GetBinLowEdge(1)-1), h2_sm->GetBinError(0));
    // Restore the number of entries
    htmp_sm2->SetEntries(h2_sm->GetEffectiveEntries());
htmp_sm2->SetLineColor(kRed);
htmp_sm2->Draw("hist");
gPad->SetLogy();

b7.SaveAs("vrtx_chi_sig_minb.pdf");
h1_sm->SaveAs("vrtx_chi_sigminb_norm.root");
htmp_sm->SaveAs("vrtx_chi_sigminb_norm_ovrflw.root");
h2_sm->SaveAs("vrtx_chi_sigminb_norm_after02.root");
htmp_sm2->SaveAs("vrtx_chi_sigminb_norm_after02_ovrflw.root");


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

ac.SaveAs("h_aco_sig_min.pdf");
h_acoN->SaveAs("h_aco_sig_min.root");
h_acoN_02->SaveAs("h_aco_sig_min_after02.root");

TCanvas b2("b2","b2",700,700);
b2.Divide(2,2);
b2.cd(1);
h_thep_all->Draw("hist");
h_thep_gen->SetLineColor(kOrange);
h_thep_gen->Draw("hist same");
h_thep_rec->SetLineColor(kRed);
h_thep_rec->Draw("hist same");
gPad->SetLogy();
b2.cd(2);
 TH1F *rec_t = (TH1F*)h_thep_rec->Clone("rec_t");
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
