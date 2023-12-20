#include <iostream>
#include <TROOT.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>


using namespace std;





int main(int argc, char *argv[]) {

	TString filename = argv[1];
	TFile *infileRawData = new TFile(filename);

        TString outfileName = argv[2];
        TFile* outfileRawData = new TFile(outfileName, "UPDATE");

	outfileRawData->cd();
	TList* tl = (TList*)infileRawData->Get("BranchList");
	tl->Write("BranchList", 1);


	outfileRawData->Close();
	infileRawData->Close();

	cout<<"\nCreated outputfile: "<<outfileName<<endl;

return 0;


}
