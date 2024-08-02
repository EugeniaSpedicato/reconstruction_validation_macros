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
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TSystemDirectory.h"
#include <TStyle.h>

using namespace std;

void analisi_sim(){

TChain * cbmsim = new TChain("cbmsim");

          cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/bkg/mesmer_bkg_2hit_1M_1.root");
          cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/bkg/mesmer_bkg_2hit_1M_2.root");


   TClonesArray *MCTrack = 0;
   MUonERecoOutput *ReconstructionOutput = 0;
   MuE::Event *MesmerEvent = 0;
   cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);
   cbmsim->SetBranchAddress("MesmerEvent", &MesmerEvent);
   cbmsim->SetBranchAddress("MCTrack", &MCTrack);

   TFile* fout = TFile::Open("bkg.root","RECREATE");
   TH1D* h_DEmm = new TH1D("hDEmm","muon energy loss ",200,0.,200.);
   TH1D* h_Eee = new TH1D("hEee","e+e- energy sum ",200,0.,200.);
   TH1D* h_Mee = new TH1D("hMee","e+e- invariant ",200,0.,1.);
   TH1D* h_Mmm = new TH1D("hMmm","muon invariant ",200,0.,1.);
   TH1D* h_Eemin = new TH1D("hEemin","e- energy ",100,0.,100.);
   TH1D* h_Eeplus = new TH1D("hEeplus","e+ energy ",100,0.,100.); 

   TH1D* h_zFit = new TH1D("hzFit","vertex position", 500, 900., 950.);   
   TH1D* h_chi2Fit_T = new TH1D("hchi2FitT","vertex chi2 Target", 500, 0., 50.);   
   TH2D* h_t2_t1_T_1 = new TH2D("h_t2_t1_T_1"," scattering angles, cut 1 ", 500, 0.,0.050, 500, 0.,0.050);   
   TH2D* h_t2_t1_T_1u = new TH2D("h_t2_t1_T_1u"," scattering angles, cut 1, usual ", 500, 0.,0.050, 50, 0.,0.005);   
   TH2D* h_t2_t1_T_2 = new TH2D("h_t2_t1_T_2"," scattering angles, cut 2 ", 500, 0.,0.050, 500, 0.,0.050);   
   TH2D* h_s_d_T = new TH2D("hsdT"," scattering angles, sum and difference of R and L ", 500, 0.,0.050, 500, 0.,0.050);
   TH2D* h_s_d_T_1 = new TH2D("hsdT1"," scattering angles, sum and difference of R and L, after cut ", 500, 0.,0.050, 500, 0.,0.050);
   TH1D* h_chi2Fit_Si = new TH1D("hchi2FitSi","vertex chi2 Silicon", 500, 0., 50.);   
   TH1D* h_ns0_T = new TH1D("hns0T","nstub in s0, Target", 30, 0., 30.);   
   TH1D* h_ns1_T = new TH1D("hns1T","nstub in s1, Target", 30, 0., 30.);   
   TH1D* h_nt1_T = new TH1D("hnt1T","n of tracksin s1, before cuts", 30, 0., 30.);   
   TH1D* h_ns10_T = new TH1D("hns10T","nstub in 0 s1, Target", 10, 0., 10.);   
   TH1D* h_ns11_T = new TH1D("hns11T","nstub in 1 s1, Target", 10, 0., 10.);   
   TH1D* h_ns12_T = new TH1D("hns12T","nstub in 2 s1, Target", 10, 0., 10.);   
   TH1D* h_ns13_T = new TH1D("hns13T","nstub in 3 s1, Target", 10, 0., 10.);   
   TH1D* h_ns14_T = new TH1D("hns14T","nstub in 4 s1, Target", 10, 0., 10.);   
   TH1D* h_ns15_T = new TH1D("hns15T","nstub in 5 s1, Target", 10, 0., 10.);   
   TH1D* h_times_T = new TH1D("htimesT","ratio s1, Target", 10, 0., 1E-5);   
   TH1D* h_chi2track1_T =  new TH1D("hchi2track1","track 1 chi2", 500, 0., 50.);
   TH1D* h_chi2track2_T =  new TH1D("hchi2track2","track 2 chi2", 500, 0., 50.);

   TH1D* h_ns10_T_1 = new TH1D("hns10T1","nstub in 0 s1, after cut 1", 10, 0., 10.);   
   TH1D* h_nt1_T_1 = new TH1D("hnt1T1","n of tracks in s1, after cut", 30, 0., 30.);   
   TH1D* h_ns1_T_1 = new TH1D("hns1T1","nstub in s1, Target, after cut 1", 30, 0., 30.);   
   TH1D* h_ns11_T_1 = new TH1D("hns11T1","nstub in 1 s1, after cut 1", 10, 0., 10.);   
   TH1D* h_ns12_T_1 = new TH1D("hns12T1","nstub in 2 s1, after cut 1", 10, 0., 10.);   
   TH1D* h_ns13_T_1 = new TH1D("hns13T1","nstub in 3 s1, after cut 1", 10, 0., 10.);   
   TH1D* h_ns14_T_1 = new TH1D("hns14T1","nstub in 4 s1, after cut 1", 10, 0., 10.);   
   TH1D* h_ns15_T_1 = new TH1D("hns15T1","nstub in 5 s1, after cut 1", 10, 0., 10.);   
   TH1D* h_times_T_1 = new TH1D("htimesT1","ratio s1, after cut 1", 10, 0., 1E-5);   
   TH1D* h_chi2Fit_T_1 = new TH1D("hchi2FitT1","vertex chi2, after cut 1", 500, 0., 50.);   
   TH1D* h_ns1_Si = new TH1D("hns1Si","nstub in s1, Silicon", 30, 0., 30.);   
   TH1D* h_th_in = new TH1D("hth_in","theta mu_in golden", 100, 0.,0.005*TMath::Pi());   
   TH1D* h_x_in = new TH1D("hx_in","x mu in",  200, -5.,+5.);
   TH1D* h_y_in = new TH1D("hy_in","y mu in",  200, -5.,+5.);
   TH2D* h_x_th_in = new TH2D("hx_th_in","x vs theta mu in",  200, -5.,+5., 100, 0.,0.005*TMath::Pi());   
   TH2D* h_y_th_in = new TH2D("hy_th_in","y vs theta mu in",  200, -5.,+5., 100, 0.,0.005*TMath::Pi());   
   TH2D* h_t2_t1_T = new TH2D("h_t2_t1_T"," scattering angles ", 500, 0.,0.050, 500, 0.,0.050);   
   TH2D* h_t2_t1_Si = new TH2D("h_t2_t1_Si"," scattering angles ", 500, 0.,0.050, 500, 0.,0.050);   
   TH2D* h_t2_t1_T_02 = new TH2D("h_t2_t1_T_02"," scattering angles ", 500, 0.,0.050, 500, 0.,0.050);   
   TH2D* h_t2_t1_Si_02 = new TH2D("h_t2_t1_Si_02"," scattering angles ", 500, 0.,0.050, 500, 0.,0.050);   
   TH1D* h_acopl = new TH1D("hacopl","acoplanarity", 200, -5, 5.);   

   TH1D* h_th_out = new TH1D("h_th_out","theta mu out", 100, 0.,0.01*TMath::Pi());   
   TH1D* h_dth = new TH1D("h_dth","delta theta mu", 200, -0.001,0.001);   
   TH1D* h_phi_in = new TH1D("h_phi_in","phi mu in", 2000,-TMath::Pi(),TMath::Pi());   
   TH1D* h_phi_out = new TH1D("h_phi_out","phi mu out", 2000,-TMath::Pi(),TMath::Pi());   
   TH1D* h_dphi = new TH1D("h_dphi","delta phi mu", 2000,-TMath::Pi(),TMath::Pi());   


cout << " Entries " <<  cbmsim->GetEntries() << endl;
//for(Long64_t i = 0; i < 1E+6; i++) {
for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
		cbmsim->GetEntry(i);
		if(i%1000 == 0) cout<<"Entry "<<i<<endl;


int nStubs_0 = 0;
int nStubs_00 = 0; 
int nStubs_01 = 0; 
int nStubs_02 = 0;
int nStubs_03 = 0;
int nStubs_04 = 0;
int nStubs_05 = 0;


int nStubs_1 = 0;
int nStubs_10 = 0; 
int nStubs_11 = 0; 
int nStubs_12 = 0;
int nStubs_13 = 0;
int nStubs_14 = 0;
int nStubs_15 = 0;


vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
vector<MUonERecoOutputVertex> vrtx = ReconstructionOutput->reconstructedVertices();

int n_0=0;
//number of reco tracks in sector0 
for(int t=0; t<tracks.size(); t++){ if(tracks.at(t).sector()==0) n_0++;}
/*
TVector3 p_out0_MC;
TVector3 p_out1_MC;
TVector3 p_out2_MC;
TVector3 p_out3_MC;
TVector3 p_out4_MC;
TLorentzVector v0;
TLorentzVector v1;
TLorentzVector v2;
TLorentzVector v3;
TLorentzVector v4;
TLorentzVector vee;
TLorentzVector vmm;
TLorentzVector vtot;


double mass_mu = 105.6583745E-3;
double mass_e  = 0.51099895069E-3;

int code_0 = -999;
int code_1 = -999;
int code_2 = -999;
int code_3 = -999;
int code_4 = -999;

for(int n = 0; n < MCTrack->GetEntries(); n++) 
{
    const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));

    if(MCTr->pdgCode()==-13 and MCTr->interactionID() == 0) 
    {
	    code_0 = n;  
	    v0.SetPxPyPzE(MCTr->px(),MCTr->py(),MCTr->pz(), MCTr->energy());
    }
    if(MCTr->pdgCode()==-13 and MCTr->interactionID() == 45) 
    { 
	    code_1 = n;  
	    p_out1_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); 
	    v1.SetPxPyPzE(MCTr->px(),MCTr->py(),MCTr->pz(), MCTr->energy());
	    p_out1_MC.Unit();
    } 
    if(MCTr->pdgCode()==-11 and MCTr->interactionID() == 45) 
    { 
	    code_2 = n;  
	    p_out2_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); 
	    v2.SetPxPyPzE(MCTr->px(),MCTr->py(),MCTr->pz(), MCTr->energy());
	    p_out2_MC.Unit();
    }
    if(MCTr->pdgCode()==11 and MCTr->interactionID() == 45) 
    { 
	    code_3 = n;  
	    p_out3_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); 
	    v3.SetPxPyPzE(MCTr->px(),MCTr->py(),MCTr->pz(), MCTr->energy());
	    p_out3_MC.Unit();
    } 
    if(MCTr->pdgCode()==22) 
    { 
	    code_4 = n;  
	    p_out4_MC.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); 
	    v4.SetPxPyPzE(MCTr->px(),MCTr->py(),MCTr->pz(), MCTr->energy());
	    p_out4_MC.Unit();
    } 
}

    h_DEmm->Fill((v0-v1).E(), MesmerEvent->wgt_full);
    h_Mmm->Fill((v0-v1).Mag(), MesmerEvent->wgt_full);
    h_Eee->Fill(v2.E()+v3.E(), MesmerEvent->wgt_full);
    h_Eemin->Fill(v2.E(), MesmerEvent->wgt_full);
    h_Eeplus->Fill(v3.E(), MesmerEvent->wgt_full);
    h_Mee->Fill((v2+v3).Mag(), MesmerEvent->wgt_full);


    if(code_0 != -999 and code_1 != -999) vmm = v0 - v1;
    if(code_2 != -999 and code_3 != -999) vee = v2 + v3;
    vtot = v0-v1-v2-v3;
    cout << " px mu in , out  " << v0.Px() << "\t " << v1.Px() << endl;
    cout << " py mu in , out  " << v0.Py() << "\t " << v1.Py() << endl;
    cout << " pz mu in , out  " << v0.Pz() << "\t " << v1.Pz() << endl;
    cout << " E  mu in , out  " << v0.E() << "\t " << v1.E() << endl;
    cout << " E mu mu   " << vmm.E() << endl;
    cout << " px mu mu   " << vmm.Px() << endl;
    cout << " py mu mu   " << vmm.Py() << endl;
    cout << " pz mu mu   " << vmm.Pz() << endl;
    cout << " inv e e      " << vee.Mag() << endl;
    cout << " E e e      " << vee.E() << endl;
    cout << " px e e     " << vee.Px() << endl;
    cout << " py e e     " << vee.Py() << endl;
    cout << " pz e e     " << vee.Pz() << endl;
    cout <<  " vtot  E   " << vtot.E() << endl;
    cout << "  vtot  px  " << vtot.Px() << endl;
    cout << "  vtot  py  " << vtot.Py() << endl;
    cout << "  vtot  pz  " << vtot.Pz() << endl;

    cout << "invariant mass of the TLorentz sum mu mu   " << vmm.Mag() << endl;
    cout << "invariant mass of the TLorentz sum " << vee.Mag() << endl;
cout << " end of searching for the associaton " << endl;
*/

std::vector<MUonERecoOutputHit> stubs=ReconstructionOutput->reconstructedHits();
// number of stubs per sector
  for (auto & hit : stubs) {
           if (hit.stationID()==0 && hit.moduleID()==0) nStubs_00++;
      else if (hit.stationID()==0 && hit.moduleID()==1) nStubs_01++;
      else if (hit.stationID()==0 && hit.moduleID()==2) nStubs_02++;
      else if (hit.stationID()==0 && hit.moduleID()==3) nStubs_03++;
      else if (hit.stationID()==0 && hit.moduleID()==4) nStubs_04++;
      else if (hit.stationID()==0 && hit.moduleID()==5) nStubs_05++;
	}
  for (auto & hit : stubs) {
           if (hit.stationID()==1 && hit.moduleID()==0) nStubs_10++;
      else if (hit.stationID()==1 && hit.moduleID()==1) nStubs_11++;
      else if (hit.stationID()==1 && hit.moduleID()==2) nStubs_12++;
      else if (hit.stationID()==1 && hit.moduleID()==3) nStubs_13++;
      else if (hit.stationID()==1 && hit.moduleID()==4) nStubs_14++;
      else if (hit.stationID()==1 && hit.moduleID()==5) nStubs_15++;
      if (hit.stationID()==1) nStubs_1++;
	}
if (nStubs_00 == 1 && nStubs_01 == 1 && nStubs_02 == 1 && nStubs_03 == 1 && nStubs_04 == 1 && nStubs_05 == 1) 
	nStubs_0 = 6; // golden muon
	else
	nStubs_0 = nStubs_00 + nStubs_01 + nStubs_02 + nStubs_03 + nStubs_04 + nStubs_05;
TVector3 p_muin;
TVector3 p_out1;
TVector3 p_out2;
double th1=-99.;
int link1=-99;

double z0 = 810. + 18.0218;
double z1 = 810. + 21.8693;
double zt = 911.2 + 1.5;

/*
for(int j=0; j<vrtx.size();j++)
{
if(vrtx.at(j).stationIndex()==1)
{
        cout << "vertex j = " << j << endl;
	cout << "There is a reconstructed vrtx with chi2 = " << vrtx.at(j).chi2perDegreeOfFreedom() << endl;
        cout << "The vrtx has a z position = " << vrtx.at(j).zPositionFit() << endl;
        cout << "The vrtx has a z knematic = " << vrtx.at(j).zKinematicFit() << endl;
	}
}
*/

MUonERecoOutputVertex best_vrtx = ReconstructionOutput->bestVertex();

if(n_0==1 and best_vrtx.chi2perDegreeOfFreedom() !=0) { // single track in S0 and at least 1 vertex
  h_zFit->Fill(best_vrtx.zPositionFit(),MesmerEvent->wgt_full); // z vertex
  h_ns0_T->Fill(nStubs_0,MesmerEvent->wgt_full); // number of stubs in station 1 
  h_chi2Fit_T->Fill( best_vrtx.chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);
  h_ns1_T->Fill(nStubs_1,MesmerEvent->wgt_full); // number of stubs in station 1 
  h_nt1_T->Fill(tracks.size() - 1,MesmerEvent->wgt_full);      // number of tracks in station 1
  h_ns10_T->Fill(nStubs_10,MesmerEvent->wgt_full);
  h_ns11_T->Fill(nStubs_11,MesmerEvent->wgt_full);
  h_ns12_T->Fill(nStubs_12,MesmerEvent->wgt_full);
  h_ns13_T->Fill(nStubs_13,MesmerEvent->wgt_full);
  h_ns14_T->Fill(nStubs_14,MesmerEvent->wgt_full);
  h_ns15_T->Fill(nStubs_15,MesmerEvent->wgt_full);

   for(int t=0; t<tracks.size(); t++) 
   { 
      if(tracks.at(t).sector()==0) 
      {
         p_muin.SetXYZ(tracks.at(t).xSlope(),tracks.at(t).ySlope(),1.0); 
         p_muin=p_muin.Unit(); 
         double theta_in = p_muin.Theta();
         h_th_in->Fill(theta_in,MesmerEvent->wgt_full); // incoming muon direction 
         std::vector<MUonERecoOutputHit> hits_track=tracks.at(t).hits();
         for(auto & hit : hits_track) 
	 {
           if(hit.moduleID()==0) 
	   { 
              h_x_in->Fill(hit.positionPerpendicular(),MesmerEvent->wgt_full); // incoming muon x entering point
              h_x_th_in->Fill(hit.positionPerpendicular(),tracks.at(t).xSlope(),MesmerEvent->wgt_full);
           }
	   else if(hit.moduleID()==1) 
	   {
               h_y_in->Fill(hit.positionPerpendicular(),MesmerEvent->wgt_full); 
               h_y_th_in->Fill(hit.positionPerpendicular(),tracks.at(t).ySlope(),MesmerEvent->wgt_full);
           } 
        }
      }
  }
   //
   // define the scattering angles 1 and 2 and the Left and Right angles
   //
  MUonERecoOutputTrack t1_out = best_vrtx.outgoingMuon();
  MUonERecoOutputTrack t2_out = best_vrtx.outgoingElectron();
  p_out1.SetXYZ(t1_out.xSlope(),t1_out.ySlope(),1.0); 
  p_out1 = p_out1.Unit(); 
  p_out2.SetXYZ(t2_out.xSlope(),t2_out.ySlope(),1.0); 
  p_out2 = p_out2.Unit(); 
  double theta_1 = TMath::ACos(p_muin*p_out1); 
  double theta_2 = TMath::ACos(p_muin*p_out2); 
  double theta_R = 0.;
  double theta_L = 0.;
  h_chi2track1_T->Fill(t1_out.chi2perDegreeOfFreedom());
  h_chi2track2_T->Fill(t2_out.chi2perDegreeOfFreedom());
  if(TMath::Cos(p_out1.Phi()) > TMath::Cos(p_out2.Phi())) 
  {
	  theta_R = theta_1;
	  theta_L = theta_2;
  }
  else
  {
	  theta_R = theta_2;
	  theta_L = theta_1;
  }
  if(theta_1 > theta_2) 
  {
	  double dummy = theta_1;
	  theta_1 = theta_2;
	  theta_2 = dummy;
  } 
  //
  // define the acoplanarity
  //
  double tripleProduct = p_muin.Dot(p_out1.Cross(p_out2));
  double acoplanarity = TMath::Pi() - (p_muin.Cross(p_out1)).Angle(p_muin.Cross(p_out2));
  if(tripleProduct < 0) acoplanarity *= -1;
  h_acopl->Fill(acoplanarity,MesmerEvent->wgt_full); 
  //
  // Histograms before cuts are applied
  //
  h_t2_t1_T->Fill(theta_L,theta_R,MesmerEvent->wgt_full);
  h_s_d_T->Fill(theta_L+theta_R,theta_R-theta_L,MesmerEvent->wgt_full); // correlations?
  h_times_T->Fill(theta_L*theta_R,MesmerEvent->wgt_full);
  //
  // Associated hit
  //
   std::vector<MUonERecoOutputHit> hits_t1_out = t1_out.hits();
   std::vector<MUonERecoOutputHit> hits_t2_out = t2_out.hits();
  //
  // Selection cuts 
  //
  if(
		  abs(acoplanarity) <= 1 and theta_L >= 0.2E-3 and theta_R >= 0.2E-3 
		  and 910. <= best_vrtx.zPositionFit() and best_vrtx.zPositionFit() <= 914. 
		  and  best_vrtx.chi2perDegreeOfFreedom() < 20.
		  and  nStubs_1 - hits_t1_out.size() - hits_t2_out.size() <= 2 // it can be refined 
		  and  theta_L * theta_R * 1E+6 > 4.
   ) 
  {
     h_t2_t1_T_1->Fill(theta_L,theta_R,MesmerEvent->wgt_full);
     h_t2_t1_T_1u->Fill(theta_2,theta_1,MesmerEvent->wgt_full);
     h_ns1_T_1->Fill(nStubs_1,MesmerEvent->wgt_full); // number of stubs in station 1 
     h_nt1_T_1->Fill(tracks.size()-1,MesmerEvent->wgt_full); // tracks in S1 after cuts
     h_ns10_T_1->Fill(nStubs_10,MesmerEvent->wgt_full);
     h_ns11_T_1->Fill(nStubs_11,MesmerEvent->wgt_full);
     h_ns12_T_1->Fill(nStubs_12,MesmerEvent->wgt_full);
     h_ns13_T_1->Fill(nStubs_13,MesmerEvent->wgt_full);
     h_ns14_T_1->Fill(nStubs_14,MesmerEvent->wgt_full);
     h_ns15_T_1->Fill(nStubs_15,MesmerEvent->wgt_full);
     h_times_T_1->Fill(theta_L*theta_R,MesmerEvent->wgt_full);
     h_chi2Fit_T_1->Fill( best_vrtx.chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);
     h_s_d_T_1->Fill(theta_L+theta_R,theta_R-theta_L,MesmerEvent->wgt_full); 
     if(theta_R * theta_L * 1E+6 > 4. and  theta_R * theta_L * 1E+6 < 16.)  h_t2_t1_T_2->Fill(theta_L,theta_R,MesmerEvent->wgt_full);
  else if(best_vrtx.zPositionFit() > 920.) 
  {
	  h_chi2Fit_Si->Fill( best_vrtx.chi2perDegreeOfFreedom(),MesmerEvent->wgt_full);
	  h_ns1_Si->Fill(nStubs_1,MesmerEvent->wgt_full);
	  h_t2_t1_Si_02->Fill(theta_L,theta_R,MesmerEvent->wgt_full);
  }
  } 
  }
}

fout->Write();
fout->Close();
}//end of file
