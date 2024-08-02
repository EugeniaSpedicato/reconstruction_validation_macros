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

void pre_dati_notReco(){

TChain * cbmsim = new TChain("cbmsim");

cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_0hit.root");
//dataReconstruction_skim_run6_bestConfig_2.root");

/*cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_2hit.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_2hit_1.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_2hit_2.root");
*/
/*cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_skim_run_6_v2_2hits_nochi2_MCcorr_gioskim.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_skim_run_6_v2_2hits_nochi2_MCcorr_gioskim_1.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_skim_run_6_v2_2hits_nochi2_MCcorr_gioskim_2.root");*/

/*cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_skim_run_6_v2_2hits_nochi2_MCcorr.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_skim_run_6_v2_2hits_nochi2_MCcorr_1.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/skim_GA_divided/dataReconstruction_skim_run_6_v2_2hits_nochi2_MCcorr_2.root");*/


        TClonesArray *MCTrack = 0;
        MuE::Event *MesmerEvent = 0;
        MUonERecoOutput *ReconstructionOutput = 0;

        cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);

                auto pos_on_track = [](double q, double m, double z){return (q + m*z);};

double signal=0.;

double reco_1=0.; double reco1_1=0.; double more_reco_1=0.; double reco0_1=0.;double reco3_1=0.;
double reco_2=0.; double reco1_2=0.; double more_reco_2=0.; double reco0_2=0.;double reco3_2=0.;
double reco_3=0.; double reco1_3=0.; double more_reco_3=0.; double reco0_3=0.;double reco3_3=0.;
double reco_4=0.; double reco1_4=0.; double more_reco_4=0.; double reco0_4=0.;double reco3_4=0.;
double reco_5=0.; double reco1_5=0.; double more_reco_5=0.; double reco0_5=0.;double reco3_5=0.;
double reco_6=0.; double reco1_6=0.; double more_reco_6=0.; double reco0_6=0.;double reco3_6=0.;

double error1=0.;double error2=0.;double error3=0.;double error4=0.;double error5=0.;double error6=0.;

double reco_v=0.; double more_reco_v=0.; double reco0_v=0.;
int yes_v=0;
int code_mu=-99; int code_e=-99; int code_mu_in=-99;
int TrackIdreco=-99;
double z_fix=912.7;

TH2D* h_xy=new TH2D("h_xy","Beam Spot Z=target",100,-5.,5.,100,-5.,5.);


   const Int_t NBINS2 = 9;
   Double_t edges2[NBINS2 + 1] = {0.,0.0002,0.0004,0.0006,0.0008,0.001,0.002,0.003,0.004,0.005};

TH2D* h_2d=new TH2D("h2D","theta mu vs theta E with all cuts", 16, 0.,0.032,NBINS2,edges2);

TH1::SetDefaultSumw2(kTRUE);

   const Int_t NBINS = 6;
   Double_t edges[NBINS + 1] = {0.0, 0.005, 0.010, 0.015, 0.020, 0.025, 0.032};
   TH1D* d_aco = new TH1D("d_aco_real", "Acoplanarity",600,-3.2,3.2);
   TH1D* theta_e = new TH1D("theta_e", "Electron scattering angles from MESMER",10,0.,0.035);
   TH1D* theta_mu = new TH1D("theta_mu", "Muon scattering angles from MESMER",20,0.,0.005);
   TH1D* th_in=new TH1D("th_in","Incoming muon theta",100,0.,0.01);
TH1D *h_opening=new TH1D("h_opening", "Opening angle reco events from MESMER",35,0.,0.035);

for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
		cbmsim->GetEntry(i);
//		if(i%1000 == 0) cout<<"Entry "<<i<<endl;

vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
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
int yes2=0;

    for(int j=0; j<tracks.size();j++)
    {
        if(tracks.at(j).sector()==0) sec0++;
        if(tracks.at(j).sector()==1) sec1++;
        }

//std::array<double,2> th_rec;
double the_rec=-99.;
double thmu_rec=-99.;

//std::array<TVector3,2> p;
TVector3 p_e,p_mu;

        std::vector<MUonERecoOutputHit> hits_mu=vrtx.outgoingMuon().hits();
        std::vector<double> pos_mu; pos_mu.resize(6);
	for(int p=0;p<hits_mu.size();p++)pos_mu.push_back(hits_mu.at(p).position());
	std::vector<MUonERecoOutputHit> hits_e=vrtx.outgoingElectron().hits();
        std::vector<double> pos_e; pos_e.resize(6);
        for(int p=0;p<hits_e.size();p++)pos_e.push_back(hits_e.at(p).position());

TVector3 p_muin;


MUonERecoOutputTrack t_mu;
MUonERecoOutputTrack t_e;

for(int j=0; j<tracks.size();j++)
{
	std::vector<MUonERecoOutputHit> hits_=tracks.at(j).hits();

        if(tracks.at(j).sector()==0 and sec0==1){

	stubs_muin=hits_.size();
        th_inx=tracks.at(j).xSlope();
        th_iny=tracks.at(j).ySlope();
        x0_in=tracks.at(j).x0();
        y0_in=tracks.at(j).y0();
        chi2_muin=tracks.at(j).chi2perDegreeOfFreedom();
	p_muin.SetXYZ(th_inx,th_iny,1.0);
        p_muin=p_muin.Unit();
        th_muin=p_muin.Theta();
                        }
	//if(tracks.at(j).sector()==1 and sec1==2)
	if(tracks.at(j).sector()==1)
	{
        std::vector<double> pos; pos.resize(6);
        for(int p=0;p<hits_.size();p++)pos.push_back(hits_.at(p).position());

	 if(std::equal(pos.begin(),pos.end(),pos_mu.begin())){p_mu.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0); p_mu.Unit(); thmu_rec=p_mu.Angle(p_muin);t_mu=tracks.at(j);}
	else if(std::equal(pos.begin(),pos.end(),pos_e.begin())){p_e.SetXYZ(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0); p_e.Unit(); the_rec=p_e.Angle(p_muin);t_e=tracks.at(j);}

	 //TVector3 p1(tracks.at(j).xSlope(),tracks.at(j).ySlope(),1.0); p1=p1.Unit();th_rec.at(yes2)=p1.Angle(p_muin); p.at(yes2)=p1;
	 yes2++;
	}
}

//estrapolo posizione negli ultimi due moduli della seconda stazione
/*double posxIN=pos_on_track(x0_in,th_inx,899.92180);
double posyIN=pos_on_track(y0_in,th_iny,903.76930);*/

// posizione locale negli ultimi due moduli della seconda stazione
double posxIN=99.;
double posyIN=99.;

//posxIN=pos_on_track(x0_in,th_inx,z_fix);
//posyIN=pos_on_track(y0_in,th_iny,z_fix);

std::vector<MUonERecoOutputHit> stubs=ReconstructionOutput->reconstructedHits();
         std::array<int,6> module_st1{0};
        std::array<std::vector<float>,6> localX{{std::vector<float>(0.0)}};
int stub0 = 0;
int stub1 = 0;
for(int s=0; s<stubs.size(); s++){
if(stubs.at(s).stationID()==0){stub0++; if(stubs.at(s).moduleID()==4){posxIN=stubs.at(s).position(); } else if(stubs.at(s).moduleID()==5){posyIN=stubs.at(s).position();}   }
if(stubs.at(s).stationID()==1){stub1++;module_st1.at(stubs.at(s).moduleID())+=1; int link= stubs.at(s).moduleID(); localX.at(link).push_back(stubs.at(s).position());}
}


 if(sec0==1 and stubs_muin==6 and abs(posxIN)<=1.5 and abs(posyIN)<=1.5 and chi2_muin<=2 and stub0==6)th_in->Fill(th_muin);

if(sec0==1 and stubs_muin==6 and abs(posxIN)<=1.5 and abs(posyIN)<=1.5 and chi2_muin<=2 and stub0==6 and th_muin<0.004){//and stub1<=15){

 signal++;

bool allmod=std::all_of(std::begin(module_st1), std::end(module_st1), [](int i){return i==2;});
cout <<"---------------"<<endl;
double u=-99.;
  if(allmod and stub1==12){
			double real_pos1= localX.at(2).at(1);
		        double real_pos0= localX.at(2).at(0);
			u=real_pos1-real_pos0;
			if(abs(real_pos1-real_pos0)<0.1 and chi==0 ) cout << "allmod and stub1==12 event and chi==0 " << i << endl;
}

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



//  if(abs(acoplanarity_v)<=1 and chi<20 and thmu_rec>0.0002 and stub1<=15){// and vrtx.zKinematicFit()<915. and vrtx.zKinematicFit()>907.){//the_rec>=0.0005 and the_rec<=0.02){//the_rec<=0.032){
double Elastic=0.5109989461*0.001*((1+(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*cos(the_rec)*cos(the_rec))/(1-(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*(sqrt(160*160-(105.6583745 *0.001*105.6583745 *0.001))/(160+0.5109989461*0.001))*cos(the_rec)*cos(the_rec)));
double Elastic2=asin( (sin(the_rec)*sqrt(Elastic*Elastic-0.5109989461*0.001*0.5109989461*0.001))/sqrt( (160+0.5109989461*0.001-Elastic)*(160+0.5109989461*0.001-Elastic)-105.6583745 *0.001*105.6583745 *0.001 ) );

//bool allmod=std::all_of(std::begin(module_st1), std::end(module_st1), [](int i){return i==2;});
 std::array<double,6> t_e_hits{-99.};
 std::array<double,6> t_mu_hits{-99.};

 for(int h=0; h<t_mu.hits().size(); h++){
  if(t_mu.hits().at(h).stationID()==1) t_mu_hits.at(t_mu.hits().at(h).moduleID())=t_mu.hits().at(h).position();
 } 

 for(int h=0; h<t_e.hits().size(); h++){
  if(t_e.hits().at(h).stationID()==1) t_e_hits.at(t_e.hits().at(h).moduleID())=t_e.hits().at(h).position();
 }

 // if(allmod and stub1==12){
  if(stub1<=20 and thmu_rec>0.0005 and thmu_rec<=0.002 and the_rec<0.020 and the_rec>=0.003 and thmu_rec<=Elastic2+0.0002 and thmu_rec>=Elastic2-0.0002){


//first
 if(yes2>=2){

cout <<"selected event " << i <<endl;

d_aco->Fill(acoplanarity_v);

 h_2d->Fill(the_rec,thmu_rec);
 theta_mu->Fill(thmu_rec);
 theta_e->Fill(the_rec);
 h_opening->Fill(p_mu.Angle(p_e));

if(thmu_rec>=the_rec-0.00005 and thmu_rec<=the_rec+0.00005){cout<<2*0.511*thmu_rec<<endl; cout <<"thmu_rec " << thmu_rec << " the_rec " << the_rec << endl;}
}

if(the_rec>0.0 and the_rec<=0.005){
if(yes2>=2){reco_1++;}

if(yes2==2 and tracks.size()==3) {reco3_1++;}

}

//second
if(the_rec>0.005 and the_rec<=0.01){
if(yes2>=2){reco_2++;}

if(yes2==2 and tracks.size()==3) {reco3_2++;}

}

//third
if(the_rec>0.01 and the_rec<=0.015){
if(yes2>=2){reco_3++;}

if(yes2==2 and tracks.size()==3) {reco3_3++;}

}

//fourth
if(the_rec>0.015 and the_rec<=0.02){
if(yes2>=2){reco_4++;}

if(yes2==2 and tracks.size()==3) {reco3_4++;}

}

//fifth
if(the_rec>0.02 and the_rec<=0.025){
if(yes2>=2){reco_5++;}

if(yes2==2 and tracks.size()==3) {reco3_5++;}

}

//sixth
if(the_rec>0.025 and the_rec<=0.032){
if(yes2>=2){reco_6++;}

if(yes2==2 and tracks.size()==3) {reco3_6++;}

}

			}//aco e chi
else{

	if(allmod and stub1==12 and abs(u)<0.2){ cout << "not selected event " << i << endl;}
}
		}//chi!=0
	}//mu_in
yes2=0;
cout <<"------end---------"<<endl;
} //end of general for


cout <<endl;
cout << "Theta range: [0,5] mrad" << endl;
cout << "Su " << signal << " golden muons, " << reco_1 << " +- " << sqrt(reco_1) << " sono ricostruiti, con un rapporto del " << reco_1/signal*100 << "%"<< endl;
cout << "Su " << signal << " golden muons, " << reco3_1 << " sono ricostruiti con 3 tracce, con un rapporto del " << (reco3_1/signal)*100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con piu tracce (stesso id), " << more_reco_1 << " sono ricostruiti, con un rapporto del " << more_reco_1/signal*100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con 0 tracce di segnale reco, " << reco0_1 << ", con un rapporto del " << reco0_1/signal*100 << "%"<< endl;
//cout << "Su " << signal << " eventi di con 1 sola traccia di segnale reco, " << reco1_1 << ", con un rapporto del " << reco1_1/signal*100 << "%"<< endl;
cout <<endl;
cout << "Theta range: [5,10] mrad" << endl;
cout << "Su " << signal << " golden muons, " << reco_2 << " +- " << sqrt(reco_2) << " sono ricostruiti, con un rapporto del " << reco_2/signal *100 << "%"<< endl;
cout << "Su " << signal << " golden muons, " << reco3_2 << " sono ricostruiti con 3 tracce, con un rapporto del " << (reco3_2/signal)*100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con piu tracce (stesso id), " << more_reco_2 << " sono ricostruiti, con un rapporto del " <<  more_reco_3/signal *100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con 0 tracce di segnale reco, " << reco0_2 << ", con un rapporto del " << reco0_2/signal *100 << "%"<< endl;
//cout << "Su " << signal << " eventi di con 1 sola traccia di segnale reco, " << reco1_2 << ", con un rapporto del " << reco1_2/signal *100 << "%"<< endl;
cout <<endl;
cout << "Theta range: [10,15] mrad" << endl;
cout << "Su " << signal << " golden muons, " << reco_3 << " +- " << sqrt(reco_3) << " sono ricostruiti, con un rapporto del " << reco_3/signal *100 << "%"<< endl;
cout << "Su " << signal << " golden muons, " << reco3_3 << " sono ricostruiti con 3 tracce, con un rapporto del " << (reco3_3/signal)*100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con piu tracce (stesso id), " << more_reco_3 << " sono ricostruiti, con un rapporto del " <<  more_reco_3/signal *100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con 0 tracce di segnale reco, " << reco0_3 << ", con un rapporto del " << reco0_3/signal *100 << "%"<< endl;
//cout << "Su " << signal << " eventi di con 1 sola traccia di segnale reco, " << reco1_3 << ", con un rapporto del " << reco1_3/signal *100 << "%"<< endl;
cout <<endl;
cout << "Theta range: [15,20] mrad" << endl;
cout << "Su " << signal << " golden muons, " << reco_4 << " +- " << sqrt(reco_4) << " sono ricostruiti, con un rapporto del " << reco_4/signal *100 << "%"<< endl;
cout << "Su " << signal << " golden muons, " << reco3_4 << " sono ricostruiti con 3 tracce, con un rapporto del " << (reco3_4/signal)*100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con piu tracce (stesso id), " << more_reco_4 << " sono ricostruiti, con un rapporto del " <<  more_reco_4/signal *100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con 0 tracce di segnale reco, " << reco0_4 << ", con un rapporto del " << reco0_4/signal *100 << "%"<< endl;
//cout << "Su " << signal << " eventi di con 1 sola traccia di segnale reco, " << reco1_4 << ", con un rapporto del " << reco1_4/signal *100 << "%"<< endl;
cout <<endl;
cout << "Theta range: [20,25] mrad" << endl;
cout << "Su " << signal << " golden muons, " << reco_5 << " +- " << sqrt(reco_5) << " sono ricostruiti, con un rapporto del " << reco_5/signal *100 << "%"<< endl;
cout << "Su " << signal << " golden muons, " << reco3_5 << " sono ricostruiti con 3 tracce, con un rapporto del " << (reco3_5/signal)*100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con piu tracce (stesso id), " << more_reco_5 << " sono ricostruiti, con un rapporto del " <<  more_reco_5/signal *100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con 0 tracce di segnale reco, " << reco0_5 << ", con un rapporto del " << reco0_5/signal *100 << "%"<< endl;
//cout << "Su " << signal << " eventi di con 1 sola traccia di segnale reco, " << reco1_5 << ", con un rapporto del " << reco1_5/signal *100 << "%"<< endl;
cout <<endl;
cout << "Theta range: [25,32] mrad" << endl;
cout << "Su " << signal << " golden muons, " << reco_6 << " +- " << sqrt(reco_6) << " sono ricostruiti, con un rapporto del " << reco_6/signal *100 << "%"<< endl;
cout << "Su " << signal << " golden muons, " << reco3_6 << " sono ricostruiti con 3 tracce, con un rapporto del " << (reco3_6/signal)*100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con piu tracce (stesso id), " << more_reco_6 << " sono ricostruiti, con un rapporto del " <<  more_reco_6/signal *100 << "%"<< endl;
//cout << "Su " << signal << " eventi di segnale con 0 tracce di segnale reco, " << reco0_6<< ", con un rapporto del " << reco0_6/signal *100 << "%"<< endl;
//cout << "Su " << signal << " eventi di con 1 sola traccia di segnale reco, " << reco1_6 << ", con un rapporto del " << reco1_6/signal *100 << "%"<< endl;
/*
Int_t nx = h_2d->GetNbinsX();
Int_t ny = h_2d->GetNbinsY();
for (Int_t i=1; i<nx+1; i++) {
for (Int_t j=1; j<ny+1; j++) {
if (h_2d->GetBinContent(i,j)<=20) h_2d->SetBinContent(i,j,0);}}


//TCanvas a("a","a",700,700);
//th_in->Draw("E");
//a.SaveAs("comparison_RDMC/th_in_RD.pdf");


d_aco->SaveAs("comparison_RDMC/d_aco_RD_pre2.root");

//TCanvas b("b","b",700,700);
//h_2d->Draw();
h_2d->SaveAs("comparison_RDMC/2D_RD_pre2.root");

//TCanvas c("c","c",700,700);
//theta_mu->Draw("E");
theta_mu->SaveAs("comparison_RDMC/theta_mu_RD_pre2.root");

//TCanvas d("d","d",700,700);
//theta_e->Draw("E");
theta_e->SaveAs("comparison_RDMC/theta_e_RD_pre2.root");

h_opening->SaveAs("comparison_RDMC/opening_RD_pre2.root");
*/
cout << "muon" << endl;
for(int b=1; b< h_opening->GetNbinsX(); b++){
cout << "h_opening[0]->GetBinContent(b) " << h_opening->GetBinContent(b) << endl;
cout << b << " sqrt(N) " << sqrt(h_opening->GetBinContent(b)) << " GetBinError " << h_opening->GetBinError(b)<< endl;
}
cout << "muon" << endl;
for(int b=1; b< theta_mu->GetNbinsX(); b++){
cout << "theta_mu_all->GetBinContent(b) " << theta_mu->GetBinContent(b) << endl;
cout << b << " sqrt(N) " << sqrt(theta_mu->GetBinContent(b)) << " GetBinError " << theta_mu->GetBinError(b)<< endl;
}

}
