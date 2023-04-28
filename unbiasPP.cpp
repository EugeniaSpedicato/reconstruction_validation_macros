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

void RealDataAnalyzer(){


TChain * cbmsim = new TChain("cbmsim");
cbmsim->Add("test.root");

        TClonesArray *MCTrack = 0;

        cbmsim->SetBranchAddress("MCTrack", &MCTrack);


TH1D *h_E1=new TH1D("E1","Energy of muon_in-muon_44 from PP (GeV)",300,-30,30);
TH1D *h_E2=new TH1D("E2","Energy of electron+positron from PP (GeV)",300,-30,30);
TH1D *h_E3=new TH1D("E3","Energy E_e+E_p+E_44-E_in (GeV)",100,-10,10);

for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
		cbmsim->GetEntry(i);
		if(i%1000 == 0) cout<<"Entry "<<i<<endl;

double thep_gen=0; double them_gen=0; double the44_gen=0;
double Zep=0;double Zem=0;
int code_ep=0; int code_em=0; int code_44=0; int code_44_2=0;
double yes=0;
int yes_44=0;

double e44,z44,x44,y44;
double e_e,e_p,ein;

	TVector3 pep,pem, pmu_in, pmu_44;
                cout << "Entries: " << MCTrack->GetEntries() <<endl;
	for(int n = 0; n < MCTrack->GetEntries(); n++) {
	 const MUonETrack *MCTr = static_cast<const MUonETrack*>(MCTrack->At(n));
                cout << n <<  ") IntID : " <<MCTr->interactionID() << " and pdf " << MCTr->pdgCode() << endl;


         if(MCTr->interactionID()==0 and MCTr->pdgCode()==-13) {pmu_in.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); ein=MCTr->energy()-MCTr->mass();
									cout << "MCTr->energy() " << MCTr->energy() << endl;}
					//double mass=sqrt(ein*ein-pmu_in.Mag()*pmu_in.Mag()); cout << "mass mu "<< mass << endl;}
         if(MCTr->interactionID()==44 and MCTr->pdgCode()==-13) { yes_44++;pmu_44.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); code_44=n;
								 e44=MCTr->energy()-MCTr->mass();}
                                                                //double mass=sqrt(e44*e44-pmu_44.Mag()*pmu_44.Mag()); cout << "mass mu 44 "<< mass << endl;}

	 if(MCTr->interactionID()==5){
	 if(MCTr->pdgCode()==11) {yes++; pem.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); code_em=n; Zep=MCTr->startZ(); e_e=MCTr->energy()-MCTr->mass();}
	cout << "mum " << MCTr->motherID() << endl;
                                                                //double mass=sqrt(e_e*e_e-pem.Mag()*pem.Mag()); cout << "mass e "<< mass << endl;}

         if(MCTr->pdgCode()==-11) {yes++;pep.SetXYZ(MCTr->px(),MCTr->py(),MCTr->pz()); code_ep=n; Zem=MCTr->startZ(); e_p=MCTr->energy()-MCTr->mass();}
                                                                //double mass=sqrt(e_p*e_p-pep.Mag()*pep.Mag()); cout << "mass e "<< mass << endl;}

	 }

	}


if(yes==2 and yes_44==1)
{
// h_E1->Fill(ein-e44);
double Eres=e44+e_e+e_p-ein;
h_E1->Fill(ein-e44);
h_E2->Fill(e_e+e_p);
h_E3->Fill(Eres);
cout << "Eres " << Eres << endl;
}


	TVector3 pem_dir,pep_dir, pmu_in_dir, pmu_44_dir;
	pep_dir=pep.Unit();
	pem_dir=pem.Unit();
	pmu_in_dir=pmu_in.Unit();
	pmu_44_dir=pmu_44.Unit();

double them_sdr,thep_sdr;
		them_gen=acos(pmu_in_dir.Dot(pem_dir));
                thep_gen=acos(pmu_in_dir.Dot(pep_dir));
                the44_gen=acos(pmu_in_dir.Dot(pmu_44_dir));



} //end of general for

TCanvas b2("b2","b2",700,700);
b2.Divide(1,3);
b2.cd(1);
h_E1->Draw("hist");
gPad->SetLogy();
b2.cd(2);
//h_E2->SetLineColor(kGreen);
h_E2->Draw("hist");
gPad->SetLogy();
b2.cd(3);
h_E3->Draw("hist");
b2.SaveAs("energy.pdf");


}
