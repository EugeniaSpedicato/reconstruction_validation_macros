#include <TChain.h>
#include <TMath.h>
void buildHistograms(void) {

TH1D *hE = new TH1D("hE","",200,0.,200.);
TH1D *hEtot = new TH1D("hEtot","",200,0.,200.);
TH1D *hEtot0 = new TH1D("hEtot0","",200,0.,200.);
TH1D *hEtot1 = new TH1D("hEtot1","",200,0.,200.);
TH1D *hEtot2 = new TH1D("hEtot2","",200,0.,200.);
TH1D *hthe     = new TH1D("hthe", "", 300, 0, 30.E-3);
TH1D *htheE0  = new TH1D("htheE0", "", 300, 0, 30.E-3);
TH1D *htheE1  = new TH1D("htheE1", "", 300, 0, 30.E-3);
TH1D *htheE2  = new TH1D("htheE2", "", 300, 0, 30.E-3);
TH1D *htheE10  = new TH1D("htheE10", "", 300, 0, 30.E-3);
TH1D *htheE105  = new TH1D("htheE105", "", 300, 0, 30.E-3);
TH1D *hthmu    = new TH1D("hthmu", "",50, 0, 5E-3);
TH1D *hthmuE0 = new TH1D("hthmuE0", "",50, 0, 5E-3);
TH1D *hthmuE1 = new TH1D("hthmuE1", "",50, 0, 5E-3);
TH1D *hthmuE2 = new TH1D("hthmuE2", "",50, 0, 5E-3);
TH1D *hthmuE10 = new TH1D("hthmuE10", "",50, 0, 5E-3);
TH1D *hthmuE105 = new TH1D("hthmuE105", "",50, 0, 5E-3);
TH2D *hthethmu = new TH2D("hthethmu", "", 300, 0, 30.E-3, 50, 0, 5E-3);
TH2D *hthethmuE0 = new TH2D("hthethmuE0", "", 300, 0, 30.E-3, 50, 0, 5E-3);
TH2D *hthethmuE1 = new TH2D("hthethmuE1", "", 300, 0, 30.E-3, 50, 0, 5E-3);
TH2D *hthethmuE2 = new TH2D("hthethmuE2", "", 300, 0, 30.E-3, 50, 0, 5E-3);
TH2D *hthethmuE10 = new TH2D("hthethmuE10", "", 300, 0, 30.E-3, 50, 0, 5E-3);
TH2D *hthethmuE10Swap = new TH2D("hthethmuE10Swap", "", 300, 0, 30.E-3, 50, 0, 5E-3);
TH2D *hthethmuE105 = new TH2D("hthethmuE105", "", 300, 0, 30.E-3, 50, 0, 5E-3);
TH2D *hthethmuE105Swap = new TH2D("hthethmuE105Swap", "", 300, 0, 30.E-3, 50, 0, 5E-3);

TH1D *hDXele  = new TH1D("hDXele","",200,-10E-2,+10E-2);
TH1D *hDXele0 = new TH1D("hDXele0","",200,-10E-2,+10E-2);
TH1D *hDXele1 = new TH1D("hDXele1","",200,-10E-2,+10E-2);
TH1D *hDXele2 = new TH1D("hDXele2","",200,-10E-2,+10E-2);
TH1D *hDXele10 = new TH1D("hDXele10","",500,-2.,+3.);
TH1D *hDXele105 = new TH1D("hDXele105","",500,-2.,+3.);
TH1D *hDYele  = new TH1D("hDYele","",200,-10E-2,+10E-2);
TH1D *hDYele0 = new TH1D("hDYele0","",200,-10E-2,+10E-2);
TH1D *hDYele1 = new TH1D("hDYele1","",200,-10E-2,+10E-2);
TH1D *hDYele2 = new TH1D("hDYele2","",200,-10E-2,+10E-2);
TH1D *hDYele10 = new TH1D("hDYele10","",400,-2.,+2.);
TH1D *hDYele105 = new TH1D("hDYele105","",400,-2.,+2.);

TH1D *hDXmu  = new TH1D("hDXmu","",200,-10E-2,+10E-2);
TH1D *hDXmu0 = new TH1D("hDXmu0","",200,-10E-2,+10E-2);
TH1D *hDXmu1 = new TH1D("hDXmu1","",200,-10E-2,+10E-2);
TH1D *hDXmu2 = new TH1D("hDXmu2","",200,-10E-2,+10E-2);
TH1D *hDXmu10 = new TH1D("hDXmu10","",500,-2.,+3.);
TH1D *hDXmu105 = new TH1D("hDXmu105","",500,-2.,+3.);
TH1D *hDYmu  = new TH1D("hDYmu","",200,-10E-2,+10E-2);
TH1D *hDYmu0 = new TH1D("hDYmu0","",200,-10E-2,+10E-2);
TH1D *hDYmu1 = new TH1D("hDYmu1","",200,-10E-2,+10E-2);
TH1D *hDYmu2 = new TH1D("hDYmu2","",200,-10E-2,+10E-2);
TH1D *hDYmu10 = new TH1D("hDYmu10","",400,-2.,+2.);
TH1D *hDYmu105 = new TH1D("hDYmu105","",400,-2.,+2.);

TH2D *htheD = new TH2D("htheD","", 300, 0, 30.E-3, 300, 0., 3.);
TH2D *hthmuD = new TH2D("hthmuD","", 300, 0, 30.E-3, 300, 0., 3.);
TH2D *htDeDmu = new TH2D("htDeDmu","", 300, 0, 3, 300, 0., 3.);

TH2D *hwrong = new TH2D("hwrong","", 300, 0, 30.E-3, 50, 0., 5.E-3);
TH2D *hthethmuE105SwapWrong = new TH2D("hthethmuE105SwapWrong", "", 300, 0, 30.E-3, 50, 0, 5E-3);

TH2D *hRL1 = new TH2D("hRL1", "", 300, 0, 30.E-3, 300, 0, 30E-3);
TH2D *hRL10 = new TH2D("hRL10", "", 100, 0, 10.E-3, 100, 0, 10E-3);
TH2D *hRL10noE = new TH2D("hRL10noE", "", 100, 0, 10.E-3, 100, 0, 10E-3);
TH2D *hRL105 = new TH2D("hRL105", "", 100, 0, 10.E-3, 100, 0, 10E-3);
TH2D *hRL105noE = new TH2D("hRL105noE", "", 100, 0, 10.E-3, 100, 0, 10E-3);

TChain chain("cbmsim");  // replace "myTree" with your TTree name
TString name = " /home/espedica/WIP_0_17_6_fair_install/instFairRoot/share/MUonE/macros/";
name.Append("selected_T0_*.root");
cout << " chaining files " << endl; 
chain.Add(name.Data());
cout << " chaining completed " << endl;

MUonERecoOutputAnalysis *ReconstructionOutput = 0;
chain.SetBranchAddress("ReconstructionOutput", &ReconstructionOutput);
        MuE::Event *MesmerEvent = 0;
        chain.SetBranchAddress("MesmerEvent", &MesmerEvent);

long int nSwap10 = 0;
long int nTot10 = 0;
long int nSwap105 = 0;
long int nSwap105Wrong = 0;
long int nTot105 = 0;
long int nEleIsMuon105 = 0;
long int nMuIsNotMuon105 = 0;

Long64_t nentries = chain.GetEntries();

for (Long64_t i = 0; i < nentries; ++i) {
//for (Long64_t i = 0; i < 100; ++i) {
    if (i % 10000 == 0) {
       cout << " event " << i << endl;
    }

    double wgt_full=MesmerEvent->wgt_full;

    chain.GetEntry(i);
    double totalEnergy = ReconstructionOutput->reconstructedCalorimeterCluster().totalEnergy();
//    double electronDeltaXfromCaloCentroid =  ReconstructionOutput->reconstructedCalorimeterCluster().xCentroid();
//    double electronDeltaYfromCaloCentroid =  ReconstructionOutput->reconstructedCalorimeterCluster().yCentroid();
    double clusterEnergy = ReconstructionOutput->reconstructedCalorimeterCluster().clusterEnergy();
    double the  = ReconstructionOutput->bestVertex().electronTheta();
    bool   eleIsMuon =  ReconstructionOutput->bestVertex().outgoingElectron().isMuon();
    bool   muIsMuon =  ReconstructionOutput->bestVertex().outgoingMuon().isMuon();
    double thmu = ReconstructionOutput->bestVertex().muonTheta();
    double dXele = ReconstructionOutput->bestVertex().outgoingElectronDeltaXfromCaloCentroid();
    double dYele = ReconstructionOutput->bestVertex().outgoingElectronDeltaYfromCaloCentroid();
    double dXmu = ReconstructionOutput->bestVertex().outgoingMuonDeltaXfromCaloCentroid();
    double dYmu = ReconstructionOutput->bestVertex().outgoingMuonDeltaYfromCaloCentroid();
    
    vector<MUonERecoOutputTrackAnalysis> tracks = ReconstructionOutput->reconstructedTracks();

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
   double theta2 = p_v0.at(0).Angle(p_v1.at(1));
    
    
    /*
    TVectorD p_in(3);
    p_in[0] = ReconstructionOutput->bestVertex().incomingMuon().xSlope();
    p_in[1] = ReconstructionOutput->bestVertex().incomingMuon().ySlope();
    p_in[2] = 1.;
    double p_in_norm = sqrt(p_in.Norm2Sqr());
    p_in *= 1.0 / p_in_norm;

    TVectorD p1_out(3);
    p1_out[0] = ReconstructionOutput->bestVertex().outgoingMuon().xSlope();
    p1_out[1] = ReconstructionOutput->bestVertex().outgoingMuon().ySlope();
    p1_out[2] = 1.;
 
    TVectorD p2_out(3);
    p2_out[0] = ReconstructionOutput->bestVertex().outgoingElectron().xSlope();
    p2_out[1] = ReconstructionOutput->bestVertex().outgoingElectron().ySlope();
    p2_out[2] = 1.;

    double p1Dotpin = p1_out * p_in;
    double p2Dotpin = p2_out * p_in;
    //
    // orthogonal outgoing vectors with respect to the incoming muon
    //
    for(int i =0; i < 3; i++) {
	    p1_out[i] =  p1_out[i] - p1Dotpin * p_in[i];
	    p2_out[i] =  p2_out[i] - p2Dotpin * p_in[i];
    }
    //
    //check they sum to zero
    //
    double vsum = 0.;
    for(int i =0; i < 3; i++) {
	    vsum = p1_out[i] + p2_out[i];
    }

    TVectorD xvers(3);
    xvers[0] = 1.; xvers[1] =0.; xvers[2] =0.;
    double side = p1_out * xvers;
    if (side >0) nSideRight++;
    else nSideLeft++;
    //
    // check if the electron track is flagged as muon by the muon PID
    //
*/

   if(i%2 == 0) {
	   the = theta1;
	   thmu = theta2;
   } else {
	   the = theta2;
	   thmu = theta1;
   }

   hE->Fill(clusterEnergy,wgt_full);
    hthethmu->Fill(the,thmu,wgt_full);
    hthe->Fill(the,wgt_full);
    hthmu->Fill(thmu,wgt_full);
    hDXele->Fill(dXele,wgt_full);
    hDYele->Fill(dXele,wgt_full);
    hDXmu->Fill(dXmu,wgt_full);
    hDYmu->Fill(dXmu,wgt_full);


    hEtot->Fill(totalEnergy,wgt_full);
    if(clusterEnergy == 0) {
        hEtot0->Fill(totalEnergy,wgt_full);
    }

    if(clusterEnergy > 0. ) {
        hthethmuE0->Fill(the,thmu,wgt_full);
        htheE0->Fill(the,wgt_full);
        hthmuE0->Fill(thmu,wgt_full);
        hDXele0->Fill(dXele,wgt_full);
        hDYele0->Fill(dYele,wgt_full);
        hDXmu0->Fill(dXmu,wgt_full);
        hDYmu0->Fill(dYmu,wgt_full);
    }


     if(clusterEnergy > 1. ) {
        hRL1->Fill(the,thmu,wgt_full);
        hthethmuE1->Fill(the,thmu,wgt_full);
        htheE1->Fill(the,wgt_full);
        hthmuE1->Fill(thmu,wgt_full);
        hEtot1->Fill(totalEnergy,wgt_full);
        hDXele1->Fill(dXele,wgt_full);
        hDYele1->Fill(dYele,wgt_full);
        hDXmu1->Fill(dXmu,wgt_full);
        hDYmu1->Fill(dYmu,wgt_full);
    }   
      if(clusterEnergy > 2. ) {
        hthethmuE2->Fill(the,thmu,wgt_full);
        htheE2->Fill(the,wgt_full);
        hthmuE2->Fill(thmu,wgt_full);
        hEtot2->Fill(totalEnergy,wgt_full);
        hDXele2->Fill(dXele,wgt_full);
        hDYele2->Fill(dYele,wgt_full);
        hDXmu2->Fill(dXmu,wgt_full);
        hDYmu2->Fill(dYmu,wgt_full);
    }   

     if(the < 10.E-3 and thmu < 10.E-3) {
       hRL10noE->Fill(the,thmu,wgt_full);
     }
     if(clusterEnergy > 10. and the < 10.E-3 and thmu < 10.E-3) {
       hRL10->Fill(the,thmu,wgt_full);
        nTot10+=wgt_full;
        double De = sqrt(dXele*dXele + dYele*dYele);
        double Dmu = sqrt(dXmu*dXmu + dYmu*dYmu);
        if(De > Dmu)  { 
           nSwap10+=wgt_full; 
           hthethmuE10Swap->Fill(the,thmu,wgt_full);
        } 
        htheD->Fill(the,De,wgt_full);
        hthmuD->Fill(thmu,Dmu,wgt_full);
        hthethmuE10->Fill(the,thmu,wgt_full);
        htheE10->Fill(the,wgt_full);
        hthmuE10->Fill(thmu,wgt_full);
        hDXele10->Fill(dXele,wgt_full);
        hDYele10->Fill(dYele,wgt_full);
        hDXmu10->Fill(dXmu,wgt_full);
        hDYmu10->Fill(dYmu,wgt_full);
    }   

     if(the < 5.E-3 and thmu < 5.E-3) {
                 hRL105noE->Fill(the,thmu,wgt_full);
     }

     if(clusterEnergy > 10. and the < 5E-3 and thmu < 5E-3) {
         hRL105->Fill(the,thmu,wgt_full);
         double De = sqrt(dXele*dXele + dYele*dYele);
         double Dmu = sqrt(dXmu*dXmu + dYmu*dYmu);
    if(eleIsMuon or !muIsMuon) {
	    hwrong->Fill(the,thmu,wgt_full);
            if(eleIsMuon) nEleIsMuon105+=wgt_full;
            if(!muIsMuon) nMuIsNotMuon105+=wgt_full;
            if(De > Dmu)  { 
                nSwap105Wrong+=wgt_full;
                hthethmuE105SwapWrong->Fill(the,thmu,wgt_full);
            } 
	    continue;
    }

       nTot105+=wgt_full;
        if(De > Dmu)  { 
           nSwap105+=wgt_full;
           hthethmuE105Swap->Fill(the,thmu,wgt_full);
                      } 
        htDeDmu->Fill(De,Dmu,wgt_full);
        hthethmuE105->Fill(the,thmu,wgt_full);
        htheE105->Fill(the,wgt_full);
        hthmuE105->Fill(thmu,wgt_full);
        hDXele105->Fill(dXele,wgt_full);
        hDYele105->Fill(dYele,wgt_full);
        hDXmu105->Fill(dXmu,wgt_full);
        hDYmu105->Fill(dYmu,wgt_full);
    }
  }

      TString filename = "out_histograms.root";
      TFile *fhist = new TFile(filename, "RECREATE");
      gROOT->GetList()->Write();
      fhist->Close();
      cout << "number of entries " << nentries  << endl;
      cout << "nEleIsMuon105     " << nEleIsMuon105 << endl; 
      cout << "nMuIsNotMuon105   " << nMuIsNotMuon105 << endl; 
      cout << " nTot10           " << nTot10 << endl;
      cout << " nSwap10          " << nSwap10 << endl;
      cout << " nTot105          " << nTot105 << endl;
      cout << " nSwap105         " << nSwap105 << endl;
      cout << " nSwap105Wrong    " << nSwap105Wrong << endl;

}
