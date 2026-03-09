#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TGraph.h"
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

const double zT = (667.3 - 2.7);

//const double zT = (784.6 - 3.9);
const double Dz = 2.E-2; // cm
const double z0 = 0.;   // common reference point 
const double thetaMS = 13.6 * 1E-3 / 100. * sqrt(4./19.3) * (1. + 0.038 * log(4./19.3)); // good for muons

void readMCT0(string type){

   long int nTotEntries = 0;
    long int nTotRecoAtS2 = 0;
    long int nTotRecoAtS2E = 0;
    long int nElasticS2 = 0;
    long int nElasticS2E = 0;


       string name=Form("/mnt/raid10/DATA/espedica/fairmu/TB2025/MC/ideal_%s_noECAL_WiP_v1_1_x.root",type.c_str());
	TChain * cbmsim = new TChain("cbmsim");
	cbmsim->Add(name.c_str());

       TString fileNameSelected = "selected_T0_" + type;
       TFile *fout = TFile::Open(fileNameSelected, "RECREATE");
       TTree *tout = cbmsim->CloneTree(0);


        MUonERecoOutputAnalysis *ReconstructionOutput = 0;
	MuE::Event *MesmerEvent = 0;
        cbmsim->SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);
        cbmsim->SetBranchAddress("MesmerEvent", &MesmerEvent);

        cbmsim->SetBranchStatus("ReconstructionOutput*", 1);
        cbmsim->SetBranchStatus("MesmerEvent*", 1);


        long int nEntries = cbmsim->GetEntries();
        nTotEntries += nEntries;

    for(Long64_t i = 0; i < nEntries; i++) {

    cbmsim->GetEntry(i);

    if( !ReconstructionOutput->isReconstructed() ) continue;

    MUonERecoOutputVertexAnalysis vrtx = ReconstructionOutput->bestVertex();

    if(vrtx.stationIndex() != 1 ) continue;

    double wgt_full=MesmerEvent->wgt_full;

    nTotRecoAtS2+=wgt_full;
    vector<MUonERecoOutputTrackAnalysis> tracks = ReconstructionOutput->reconstructedTracks();
    vector<MUonERecoOutputHitAnalysis> hits = ReconstructionOutput->reconstructedHits();
    MUonERecoOutputTrackAnalysis e_out = vrtx.outgoingElectron();
    MUonERecoOutputTrackAnalysis mu_out = vrtx.outgoingMuon();

    double clusterEnergy = ReconstructionOutput->reconstructedCalorimeterCluster().clusterEnergy();
    double electronDeltaXfromCaloCentroid =  ReconstructionOutput->reconstructedCalorimeterCluster().xCentroid();
    double electronDeltaYfromCaloCentroid =  ReconstructionOutput->reconstructedCalorimeterCluster().yCentroid();


    if(clusterEnergy > 1.) nTotRecoAtS2E+=wgt_full; 

    // selection of the elastic events
    MUonERecoOutputTrackAnalysis mu_in_v = vrtx.incomingMuon();
    MUonERecoOutputTrackAnalysis mu_out_v = vrtx.outgoingMuon();
    MUonERecoOutputTrackAnalysis e_out_v = vrtx.outgoingElectron();

    std::array<int,3> nTr{0};
    std::vector<MUonERecoOutputTrackAnalysis> track_0;
    std::vector<MUonERecoOutputTrackAnalysis> track_1;
    std::vector<MUonERecoOutputTrackAnalysis> track_vrtx1;
    std::vector<MUonERecoOutputTrackAnalysis> track_vrtx2;
    std::vector<TVector3> p_v0;
    std::vector<TVector3> p_v1;
    std::vector<TVector3> p_v2;


    for (auto&& track : tracks) {
        if(track.sector()==0) {
           nTr.at(0)++;
           track_0.push_back(track);
           TVector3 p(track.xSlope(),track.ySlope(),1.0);
           p.Unit();
           p_v0.push_back(p);
        }
        else if(track.sector()==1) {
           nTr.at(1)++;
           track_1.push_back(track);
           TVector3 p(track.xSlope(),track.ySlope(),1.0);
           p.Unit();
           p_v1.push_back(p);
        }
        else if(track.sector()==2) {
           nTr.at(2)++;track_vrtx2.push_back(track);
           TVector3 p(track.xSlope(),track.ySlope(),1.0);
           p.Unit();
           p_v2.push_back(p);
        }
   }
   double theta1 = p_v0.at(0).Angle(p_v1.at(0));
   double theta2 = p_v0.at(0).Angle(p_v1.at(0));
   cout << " theta1 " << theta1 <<"  theta2  " << theta2 << endl;
/*
    for (auto&& track : tracks) {
        if(track.sector()==0) { 
		nTr.at(0)++; 
		track_0.push_back(track); 
	}
        else if(track.sector()==1){ 
		nTr.at(1)++; 
		track_1.push_back(track);
	        if(track.index() == mu_out_v.index() or track.index() == e_out_v.index()) { track_vrtx1.push_back(track);}              
	}
        else if(track.sector()==2) {
		nTr.at(2)++;
	        if(track.index() == mu_out_v.index() or track.index()==e_out_v.index()) { track_vrtx2.push_back(track);}
	}
    }
*/

    MUonERecoOutputHitAnalysis muin_4;
    MUonERecoOutputHitAnalysis muin_5;

    if(nTr.at(0)==1){
    std::vector<Short_t> ids_=track_0.at(0).hitIds();
    for (auto&& stub : hits) {
        for(auto&& h : ids_){
         if(h==stub.index()){
          if(stub.moduleID()==4){muin_4 = stub;}
          else if(stub.moduleID()==5){muin_5 = stub;}
          }
         }
        }
     }


        int stub0=0;int stub1=0;int stub2=0;

        for (auto&& hit : hits) {
        if(hit.stationID()==0){stub0++;}
        else if(hit.stationID()==1){stub1++;}
        else if(hit.stationID()==2){stub2++;}
        }

        double chi2 = vrtx.chi2();
//	double the_rec  = vrtx.electronTheta();
//	double thmu_rec = vrtx.muonTheta();
        double acoplanarity_v=vrtx.modifiedAcoplanarity();
        double zpos =vrtx.zPositionFit();

        double IP_out =  ReconstructionOutput->bestVertex().distanceBetweenOutgoingTracksAtTargetZ();

	// selection cuts
	if(chi2 > 20 ) continue;
       	if(IP_out > 0.2) continue;
        if( stub1 > 15 ) continue;
	if(abs(acoplanarity_v) > 0.3) continue;
        if(abs(zpos - zT) > 3.) continue;
//        if( the_rec0 > 20.E-3) continue;
//       	if( thmu_rec < 0.2E-3) continue;
        if( theta1 > 20.E-3 or theta2 < 0.2E-3) continue;
       	if( theta2 > 20.E-3 or theta1 < 0.2E-3) continue;

        if( abs(muin_4.position()) > 1.4985 || abs(muin_5.position()) > 1.4985) continue;


	// saving the events
	tout->Fill();
    nElasticS2+=wgt_full;
 
    if(clusterEnergy > 1) {
	nElasticS2E+=wgt_full;
    }

    }
      fout->Write();
      fout->Close();

     cout << "     nTotEntries                  " << nTotEntries   << endl; 
     cout << "     nTotRecoAtS2                 " << nTotRecoAtS2  << endl; 
     cout << "     nTotRecoAtS2E E>1 Gev        " << nTotRecoAtS2E << endl; 
     cout << "     nElasticS2                   " << nElasticS2    <<endl;
     cout << "     nElastiS2E                   " << nElasticS2E   << endl;

/*
     // Create a ROOT output file
      TString filename;
      filename.Form("histT0_%d_%d_%d.root", nrun, nFirstFile, nFirstFile+nOfFiles);
      TFile *fhist = new TFile(filename, "RECREATE");
      fhist->Close();
*/
}

