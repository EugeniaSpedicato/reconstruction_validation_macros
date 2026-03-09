#include <TChain.h>
auto pos_on_track = [](double q, double m, double z){return (q + m*z);};

void MC_effRecoEcal_tar0(void) {
TChain chain("cbmsim");  // replace "myTree" with your TTree name

chain.Add("/eos/user/e/ehess/production_Mesmer/FairMUonE_WiP_v1.1.x/Condor_files/MCsignalJob_target0_1.root");
chain.Add("/eos/user/e/ehess/production_Mesmer/FairMUonE_WiP_v1.1.x/Condor_files/MCsignalJob_target0_3.root");
chain.Add("/eos/user/e/ehess/production_Mesmer/FairMUonE_WiP_v1.1.x/Condor_files/MCsignalJob_target0_6.root");

TH1D *hIPdist = new TH1D("hIPdist","",100,0,1);
TH1D *hCopdist = new TH1D("hCopdist","",200,-1,1);

TH1D *hE = new TH1D("hE","all",100,0,200);
TH1D *hEreco = new TH1D("hEreco","reconstructed",100,0,200);
TH1D *hEreco12 = new TH1D("hEreco12","reco 1 in 2 out",100,0,200);
TH1D *hElastic = new TH1D("hElastic","all cuts",100,0,200);
TH1D *hAcc= new TH1D("hAcc","IP cut",100,0,200);
TH1D *hIP= new TH1D("hIP","IP cut",100,0,200);
TH1D *hT1T2= new TH1D("hT1T2","angles window",100,0,200);
TH1D *hCop= new TH1D("hCop","acoplanarity",100,0,200);


TH1D *hr = new TH1D("hr","",100,0,200);
TH1D *hr12 = new TH1D("hr12","",100,0,200);
TH1D *hrel = new TH1D("hrel","",100,0,200);
TH1D *hrIP= new TH1D("hrIP","",100,0,200);
TH1D *hrIPAcc = new TH1D("hrIPAcc","",100,0,200);
TH1D *hrAcc= new TH1D("hrAcc","",100,0,200);
TH1D *hrT1T2= new TH1D("hrT1T2","",100,0,200);
TH1D *hrCop= new TH1D("hrCop","",100,0,200);

TH2D *hT1T2E = new TH2D("hT1T2E","",100,0,200,100,-10E-3,10E-3);
TH2D *hT1vsT2R12 = new TH2D("hT1vsT2R12","",100,0,10e-3,100,0,10E-3);
TH2D *hT1vsT2Acc = new TH2D("hT1vsT2Acc","",100,0,10e-3,100,0,10E-3);
TH2D *hT1vsT2IP = new TH2D("hT1vsT2IP","",100,0,10e-3,100,0,10E-3);
TH2D *hT1T2ESum = new TH2D("hT1T2ESum","",100,0,200,100,-10E-3,10E-3);
TH2D *hT1T2ED = new TH2D("hT1T2ED","",200,-10E-3,10E-3,400,-20E-3,20E-3);
TH2D *hT1vsT2Elas = new TH2D("hT1vsT2Elas","",200,-10E-3,10E-3,400,-20E-3,20E-3);
TH1D *hEbeam = new TH1D("hEbeam","",300,70,100);
TH2D *hT1T2Emin = new TH2D("hT1T2Emin","",200,-10E-3,10E-3,200,-10E-3,10E-3);

MUonERecoOutputAnalysis *ReconstructionOutput = 0;
chain.SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);
        MuE::Event *MesmerEvent = 0;
chain.SetBranchAddress("MesmerEvent", &MesmerEvent);

Long64_t nentries = chain.GetEntries();
cout << " entries " << nentries <<endl;
for (Long64_t i = 0; i < nentries; ++i) {
    if (i % 10000 == 0) {
       cout << " event " << i << endl;
    }
    chain.GetEntry(i);
    double wgt_full=MesmerEvent->wgt_full;
//    double totalEnergy = ReconstructionOutput->reconstructedCalorimeterCluster().totalEnergy();
//    double electronDeltaXfromCaloCentroid =  ReconstructionOutput->reconstructedCalorimeterCluster().xCentroid();
//    double electronDeltaYfromCaloCentroid =  ReconstructionOutput->reconstructedCalorimeterCluster().yCentroid();
    double clusterEnergy = ReconstructionOutput->reconstructedCalorimeterCluster().clusterEnergy();
    bool isReco =  ReconstructionOutput->isReconstructed();
    hE->Fill(clusterEnergy,wgt_full);
    if(isReco) hEreco->Fill(clusterEnergy,wgt_full);

    if(clusterEnergy<30.) continue;

    vector<MUonERecoOutputHitAnalysis> hits = ReconstructionOutput->reconstructedHits();
    vector<MUonERecoOutputTrackAnalysis> tracks = ReconstructionOutput->reconstructedTracks();
    std::array<int,3> nTr{0};
    std::vector<MUonERecoOutputTrackAnalysis> track_0;
    std::vector<MUonERecoOutputTrackAnalysis> track_1;
    std::vector<MUonERecoOutputTrackAnalysis> track_vrtx1;
    std::vector<MUonERecoOutputTrackAnalysis> track_vrtx2;
    std::vector<TVector3> p_v0;
    std::vector<TVector3> p_v1;
    std::vector<TVector3> p_v2;

   //
   // counts number of tracks per station
   //

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
           nTr.at(2)++;
	   track_vrtx2.push_back(track);
           TVector3 p(track.xSlope(),track.ySlope(),1.0); 
           p.Unit(); 
           p_v2.push_back(p);
        }
   }

    if( nTr.at(0) != 1 ||  nTr.at(1) != 2) continue;

    hEreco12->Fill(clusterEnergy,wgt_full);
   //
   // elastic selection variables
   //
   double dotProduct_v = p_v1.at(0).Dot(p_v1.at(1));
   TVector3 crossProduct_v = p_v1.at(0).Cross(p_v1.at(1));
   double T_v = p_v0.at(0).Dot(crossProduct_v);
   TVector3 i1_v = p_v0.at(0).Cross(p_v1.at(0));
   TVector3 i2_v = p_v0.at(0).Cross(p_v1.at(1));
   T_v = T_v > 0 ? 1 : -1;
   double acoplanarity_v = T_v * (TMath::Pi() - acos( ((i1_v).Dot(i2_v))/(i1_v.Mag()*i2_v.Mag()) ));
   
    //
    // Fiducial cut on the incoming muon
    //

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

   //
   // lowest chi2 tracks scateering angles
   //
   double theta1 = p_v0.at(0).Angle(p_v1.at(0));
   double theta2 = p_v0.at(0).Angle(p_v1.at(1));
    hT1vsT2R12->Fill(theta1,theta2,wgt_full);
    hT1vsT2R12->Fill(theta2,theta1,wgt_full);
    /*
   if(i%2 == 0) 
    hT1vsT2R12->Fill(theta1,theta2);
   else
    hT1vsT2R12->Fill(theta2,theta1);
    */
     //
     // number of stubs per station
     //

        int stub0=0;int stub1=0;int stub2=0;

        for (auto&& hit : hits)
        {
           if(hit.stationID()==0){stub0++;}
           else if(hit.stationID()==1){stub1++;}
           else if(hit.stationID()==2){stub2++;}
        }

      //
      // impact point
      //


//       double IP_out =  ReconstructionOutput->bestVertex().distanceBetweenOutgoingTracksAtTargetZ();

        if( abs(muin_4.position()) > 1.4985 || abs(muin_5.position()) > 1.4985) continue;
        if( stub1 > 15 ) continue;

double IPx=pos_on_track(track_1.at(1).x0(),track_1.at(1).xSlope(),(667.3-2.7)) - pos_on_track(track_1.at(0).x0(),track_1.at(0).xSlope(),(667.3-2.7));
double IPy=pos_on_track(track_1.at(1).y0(),track_1.at(1).ySlope(),(667.3-2.7)) - pos_on_track(track_1.at(0).y0(),track_1.at(0).ySlope(),(667.3-2.7));
              double IP_out=sqrt(IPx*IPx + IPy*IPy);

	hAcc->Fill(clusterEnergy,wgt_full);
	hIPdist->Fill(IP_out,wgt_full);
	hCopdist->Fill(acoplanarity_v,wgt_full);
          hT1vsT2Acc->Fill(theta1,theta2,wgt_full);
          hT1vsT2Acc->Fill(theta2,theta1,wgt_full);
	  /*
        if(i%2 == 0) 
          hT1vsT2Acc->Fill(theta1,theta2);
        else
          hT1vsT2Acc->Fill(theta2,theta1);
          */

        if( theta1 > 20.E-3 or theta1 < 0.2E-3) continue;
       	if( theta2 > 20.E-3 or theta2 < 0.2E-3) continue;
	hT1T2->Fill(clusterEnergy,wgt_full);


	if(IP_out > 0.2) continue;
	hIP->Fill(clusterEnergy,wgt_full);
	 hT1vsT2IP->Fill(theta1,theta2,wgt_full);
	  hT1vsT2IP->Fill(theta2,theta1,wgt_full);
	  /*
        if(i%2 == 0) 
          hT1vsT2IP->Fill(theta1,theta2);
        else
          hT1vsT2IP->Fill(theta2,theta1);
         */


	if(abs(acoplanarity_v) > 0.3) continue; 
	hCop->Fill(clusterEnergy,wgt_full);
	hT1vsT2Elas->Fill(theta1,theta2,wgt_full);
	hT1vsT2Elas->Fill(theta2,theta1,wgt_full);


  //
   // elastic event candidate: one track in with two outgoing tracks
   //
       hElastic->Fill(clusterEnergy,wgt_full);
//       if(theta1+theta2 < 4.E-3) continue;
       if(abs(theta1-theta2) < 0.5E-3) hEbeam->Fill(clusterEnergy,wgt_full);
       hT1T2E->Fill(clusterEnergy,(theta1-theta2),wgt_full);   
       hT1T2ESum->Fill(clusterEnergy,(theta1+theta2),wgt_full);   
       hT1T2ED->Fill(theta1-theta2,theta1+theta2,wgt_full);   
       if(abs(theta1-theta2) < 0.5E-3) hT1T2Emin->Fill(theta2,theta1,wgt_full);       
  }
   hE->Sumw2();
   hEreco->Sumw2();
   hElastic->Sumw2();
   hIP->Sumw2();
   hT1T2->Sumw2();
   hCop->Sumw2();
   hr->Divide(hEreco,hE);
   hr->SetMarkerStyle(20);
   hr->SetMarkerSize(0.5);
   hr12->Divide(hEreco12,hE);
   hr12->SetMarkerStyle(20);
   hr12->SetMarkerSize(0.5);
   hrAcc->Divide(hAcc,hE);
   hrAcc->SetMarkerStyle(20);
   hrAcc->SetMarkerSize(0.5);
   hrIP->Divide(hIP,hE);
   hrIP->SetMarkerStyle(20);
   hrIP->SetMarkerSize(0.5);
   hrIPAcc->Divide(hIP,hAcc);
   hrIPAcc->SetMarkerStyle(20);
   hrIPAcc->SetMarkerSize(0.5);
 
   hrT1T2->Divide(hT1T2,hE);
   hrT1T2->SetMarkerStyle(20);
   hrT1T2->SetMarkerSize(0.5);
   hrCop->Divide(hCop,hE);
   hrCop->SetMarkerStyle(20);
   hrCop->SetMarkerSize(0.5);


   hrel->Divide(hElastic,hE);
   hrel->SetMarkerStyle(20);
   hrel->SetMarkerSize(0.5);



      TString filename = "effRecoT1ECAL30GeV_tar0.root";

      TFile *fhist = new TFile(filename, "RECREATE");
      gROOT->GetList()->Write();
      fhist->Close();
}
