
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

int Uhole_study(string version,int nhits,string bend){



TChain * cbmsim = new TChain("cbmsim");

if(version=="default" and nhits==2){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_1.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_2.root");
}
else if(version=="default" and nhits==0 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_0.root");
//cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_0hit.root");
//cbmsim->Add("/mnt/raid00/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_0hit_1.root");
}
else if(version=="default" and nhits==0 and bend=="noBend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_noBend_0hit.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_noBend_0hit_1.root");
}
else if(version=="ricMis"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_commit_b2ed7c3b.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_commit_b2ed7c3b_1.root");
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_commit_b2ed7c3b_2.root");
}
else if(version=="chi2out50" and nhits==0 and bend=="Bend"){
cbmsim->Add("/mnt/raid10/DATA/espedica/fairmu/real_data_skim_run6/dataReconstruction_skim_run6_bestConfig_0hit_chi2out50.root");
}
else if(version=="wip11" and nhits==0 and bend=="Bend"){

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


        TClonesArray *MCTrack = 0;
        MuE::Event *MesmerEvent = 0;
        MUonERecoOutput *ReconstructionOutput = 0;

        cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);


auto trackatZ = [](double q, double m,double z) {return q + (z ) * m;};

double tilt[6]={0.233,-0.233,0.,0.,0.233,-0.233};
//double tilt[6]={0.233,0.233,0.,0.,0.233,0.233};//,0.233,-0.233,0.,0.,0.233,-0.233};
double posZ[6]={18.0218,21.8693,55.3635,56.6205,89.9218,93.7693};//, 119.2218, 123.0693,156.4535, 157.8205,191.1218, 194.9693};
double rot[6]={0  *TMath::DegToRad(),90 *TMath::DegToRad(),135*TMath::DegToRad(),45 *TMath::DegToRad(),0  *TMath::DegToRad(),90 *TMath::DegToRad()};
double Zsens[6]={-0.09,+0.09,+0.2,-0.2,-0.09,+0.09};
double alpha[6]={0  *TMath::DegToRad(), 90 *TMath::DegToRad(), 135*TMath::DegToRad(), 45 *TMath::DegToRad(),0  *TMath::DegToRad(), 90 *TMath::DegToRad()};


//double offsetX[6]={-0.041151187339928,0.,0.07044990681274964,0.07366436386828577,-0.02216311771211892,0., //prima staz
double offsetX[6]={0.05309770641837935,0.,0.01831731454809257,0.1400357237114471,-0.0564074004379715,0.};//seconda staz
//double offsetY[6]={0.,-0.01962944331210647,-0.07050890136986443, 0.07323347608588825,0.,-0.01011926023475611,//prima staz
double offsetY[6]={0.,0.2324750543793093,-0.01838517522060464,0.1391336460004981,0.,-0.04404067629216706};
//double offset_alpha[6]={-0.1164515010010689*TMath::DegToRad(),-0.03411354091431071*TMath::DegToRad(),0.1369298349690934*TMath::DegToRad(),-0.4577137342911528*TMath::DegToRad(),0.007725616313160143*TMath::DegToRad(),-0.1283674677074805*TMath::DegToRad(),//prima staz
double offset_alpha[6]={0.08712552021046044*TMath::DegToRad(),0.005738391422644699*TMath::DegToRad(),0.1284617932104747*TMath::DegToRad(),0.2515280983199862*TMath::DegToRad(),-0.179018408807327*TMath::DegToRad()};//seconda staz


auto trackFit = [&] (std::vector<MUonERecoOutputHit> stubs){

        double tmpA[4*(4+1)/2] = {0.};
        ROOT::Math::SVector<double, 4> B;
        for(int j = 0; j<stubs.size(); j++){

                int linkID= stubs.at(j).moduleID();
                if(stubs.at(j).stationID()==0 or linkID==2 or linkID==3) continue;
                double cos_term = cos(rot[linkID]+offset_alpha[linkID]);
                double sin_term = sin(rot[linkID]+offset_alpha[linkID]);
//BEND
                double z_term   = stubs.at(j).z();
                double ui	= stubs.at(j).positionPerpendicular();
                double du2	= 0.002*0.002;

                tmpA[0] += cos_term*cos_term/du2;
                tmpA[1] += cos_term*sin_term/du2;
                tmpA[2] += sin_term*sin_term/du2;
                tmpA[3] += cos_term*cos_term*z_term/du2;
                tmpA[4] += cos_term*sin_term*z_term/du2;
                tmpA[5] += cos_term*cos_term*z_term*z_term/du2;
                tmpA[6] += cos_term*sin_term*z_term/du2;
                tmpA[7] += sin_term*sin_term*z_term/du2;
                tmpA[8] += cos_term*sin_term*z_term*z_term/du2;
                tmpA[9] += sin_term*sin_term*z_term*z_term/du2;

                B(0)    += (ui*cos_term + offsetX[linkID]*cos_term*cos_term + offsetY[linkID]*sin_term*cos_term)/du2;//x0
                B(1)    += (ui*sin_term + offsetX[linkID]*sin_term*cos_term + offsetY[linkID]*sin_term*sin_term)/du2;//y0
                B(2)    += (ui*z_term*cos_term + offsetX[linkID]*z_term*cos_term*cos_term + offsetY[linkID]*z_term*sin_term*cos_term)/du2;//xSlope
                B(3)    += (ui*z_term*sin_term + offsetX[linkID]*z_term*sin_term*cos_term + offsetY[linkID]*z_term*sin_term*sin_term)/du2;//ySlope
        }
       //note MatRepSym, hence not all terms are needed to define the matrix (it's symmetric)
        ROOT::Math::SMatrix<double, 4, 4, ROOT::Math::MatRepSym<double, 4>> A(tmpA, 4*(4+1)/2);

        //I'm inverting A, so that it becomes the covariance matrix:
        if(!A.Invert()) {
                cout<<"ERROR!! FAILED TO INVERT THE DERIVATIVES MATRIX!! CANNOT CREATE THE COVARIANCE MATRIX"<<endl; //failed to invert
        }

	ROOT::Math::SVector<Double_t, 4> c_svd = A*B;

        double *bestFitPars = new double[4];
        for(int i = 0; i < 4; i++) bestFitPars[i] = c_svd[i];

        return bestFitPars;};




auto trackPos = [&](int linkID,double globalz, double *bestFitPars) {
        double xtrack = bestFitPars[2]*globalz + bestFitPars[0];
        double ytrack = bestFitPars[3]*globalz + bestFitPars[1];
        return (xtrack - offsetX[linkID])*cos(rot[linkID]+offset_alpha[linkID]) +
               (ytrack - offsetY[linkID])*sin(rot[linkID]+offset_alpha[linkID]);
};



std::vector<uint64_t> nID_v;


for(Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
                cbmsim->GetEntry(i);
                if(i%1000 == 0) cout<<"Entry "<<i<<endl;

vector<MUonERecoOutputTrack> tracks = ReconstructionOutput->reconstructedTracks();
MUonERecoOutputVertex vrtx = ReconstructionOutput->bestVertex();
std::vector<MUonERecoOutputHit> stubs=ReconstructionOutput->reconstructedHits();

	 double chi=vrtx.chi2perDegreeOfFreedom();

        std::array<uint16_t,6> nstubs_per_mod_BX_st1{0};
        std::array<std::vector<double>,6> localX{{std::vector<double>(0.0)}};
        std::array<std::vector<double>,6> localZ{{std::vector<double>(0.0)}};
        int n_stubs=0.;

        std::vector<MUonERecoOutputHit> all_hits;
        for(auto&& stub : stubs)
        {
	if(stub.stationID()==1){
                int link= stub.moduleID();
        localX.at(link).push_back(stub.positionPerpendicular());
	localZ.at(link).push_back(stub.z());
        nstubs_per_mod_BX_st1.at(link)+=1;
	                    n_stubs+=1;
	all_hits.push_back(stub);}
	}

int  size1=0;
        std::array<std::vector<double>,6> posT{{std::vector<double>(0.0)}};
        std::array<std::vector<double>,6> posTZ{{std::vector<double>(0.0)}};
        std::vector<MUonERecoOutputHit> hits1;

int sec0=0;
int stubs_muin=0.;

         for (auto&& track : tracks) {
        std::vector<MUonERecoOutputHit> hits_=track.hits();

  if(track.sector()==0){
        stubs_muin=hits_.size();
	sec0++;
                        }
  if(track.sector()==1)
        {
          for(int p=0;p<hits_.size();p++){posT.at(hits_.at(p).moduleID()).push_back(hits_.at(p).positionPerpendicular());posTZ.at(hits_.at(p).moduleID()).push_back(hits_.at(p).z());}
        size1=hits_.size();
	hits1=hits_;
	}
 }

if(sec0==1 and stubs_muin==6 and n_stubs==12 and std::all_of(std::begin(nstubs_per_mod_BX_st1), std::end(nstubs_per_mod_BX_st1), [](int i){return i==2;}) ){


double u=99.;
double v=99.;
for(int m=0;m<6;m++){
      	double real_pos1= localX.at(m).at(1);
        double real_pos0= localX.at(m).at(0);

                            if(m==2){ //h_hist_distance2->Fill(real_pos1-real_pos0);
                             //h_hist_distance_zoom2->Fill(real_pos1-real_pos0);
				u=real_pos1-real_pos0;}

                            if(m==3){ //h_hist_distance3->Fill(real_pos1-real_pos0);
                             //h_hist_distance_zoom3->Fill(real_pos1-real_pos0);
				v=real_pos1-real_pos0;}

		} //end for


if(chi==0 and size1==6){
double tilt[6]={0.233,0.233,0.,0.,0.233,0.233};
        std::array<std::vector<double>,6> not_used{{std::vector<double>(0.0)}};
        std::array<std::vector<double>,6> not_usedZ{{std::vector<double>(0.0)}};


int aiut=0;
for(int p=0;p<6;p++){
if(localX.at(p).at(0)!=posT.at(p).at(0)){not_used.at(p).push_back(localX.at(p).at(0));not_usedZ.at(p).push_back(localZ.at(p).at(0));}
else if(localX.at(p).at(1)!=posT.at(p).at(0)){not_used.at(p).push_back(localX.at(p).at(1));not_usedZ.at(p).push_back(localZ.at(p).at(1));}
else{aiut=1;}
}

if(aiut==0 and (abs(u) <0.17 or abs(v)<0.17)){
        std::vector<MUonERecoOutputHit> not_used_hits;
         for (int s=0; s<all_hits.size(); s++) {

if(all_hits.at(s).positionPerpendicular()==not_used.at(all_hits.at(s).moduleID()).at(0) ){
not_used_hits.push_back(all_hits.at(s));
		}
        }
     auto paramsT = trackFit(hits1);
     auto params = trackFit(not_used_hits);

TVector3 pmuin;
TVector3 pmu;
TVector3 pe;
         for (auto&& track : tracks) {
        std::vector<MUonERecoOutputHit> hits_=track.hits();
  if(track.sector()==0)
         {pmuin.SetXYZ(track.xSlope(),track.ySlope(),1.0);}
	 }

  double mx=params[2];
  double my=params[3];
  double mxT=paramsT[2];
  double myT=paramsT[3];

pmu.SetXYZ(mxT,myT,1.0);
pe.SetXYZ(mx,my,1.0);


if( (pe.Angle(pmuin)>pmu.Angle(pmuin) and pmu.Angle(pmuin)>0.0003) or (pe.Angle(pmuin)<pmu.Angle(pmuin) and pe.Angle(pmuin)>0.0003) ){

         for (auto&& track : tracks) {
        std::vector<MUonERecoOutputHit> hits_=track.hits();
 	if(track.sector()==1)
         {
	for(int h=0;h<hits_.size();h++){
if(h==0){uint64_t nID = static_cast<uint64_t>(hits_.at(h).superID())*3564 + static_cast<uint64_t>(hits_.at(h).bx()); nID_v.push_back(nID);}
			}
		 }
}

	}//0.2mrad
  }//aiut==0
}

	}//end ==12
}//end while


cout <<"nID_v size " << nID_v.size() << endl;
for( int f=0; f<nID_v.size(); f++){
cout << nID_v.at(f) <<", ";
}

}

