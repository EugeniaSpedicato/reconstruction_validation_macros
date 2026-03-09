void sync_inside_run(int run, std::array<std::string,10> nfiles){


std::array<std::time_t,nfiles.size()> time;
std::array<std::array<TEfficiency*,nfiles.size()>,18> eff;
std::array<TLegend*,18> leg;
for(int mod=0; mod<18; mod++){leg.at(mod)= new TLegend(0.65, 0.10, 0.90, 0.55);}
int id=0;

for(auto&& nfile : nfiles){

        // Trova la posizione del trattino
        size_t pos = nfile.find('-');
        if (pos != std::string::npos && pos + 1 < nfile.size()) {
            std::string number = nfile.substr(pos + 1);  // prendi la parte dopo '-'
            std::cout << "Input: " << nfile << "  →  Numero estratto: " << number << std::endl;
		time.at(id)=std::stol(number);
        } else {
            std::cout << "Formato non valido: " << nfile << std::endl;
        }

	for(int mod=0; mod<18; mod++){
	TFile *f0=TFile::Open(Form("/home/espedica/trackeranalysis/results/analysis_run%i_muedaq%s.root",run,nfile.c_str()));
	TDirectoryFile* Synchronisation_Link=(TDirectoryFile*)f0->Get("Synchronisation_Link");
	Synchronisation_Link->cd();
	TCanvas *cEff_0 = (TCanvas*)gDirectory->Get(Form("cEff_%i",mod));
	eff.at(mod).at(id)=(TEfficiency*)cEff_0->FindObject(Form("eff_IN TIME_%i",mod));
	leg.at(mod)->AddEntry(eff.at(mod).at(id),Form("Run %i",run),"l");
	f0->Close();
	}
	id++;
}


TMultiGraph *mgx = new TMultiGraph();
std::array<TGraphAsymmErrors*,time.size()> gr_eff;

TCanvas c_gr("c_gr","c_gr",3000,2000);

//for(auto&& run : time){
//int i=0;l
for (int i = 0; i < time.size(); i++) {
 gr_eff.at(i)=new TGraphAsymmErrors(18);
	for (int m = 0; m < 18; m++) {
	 gr_eff.at(i)->AddPoint(m,eff.at(m).at(i)->GetEfficiency(13));
	 gr_eff.at(i)->SetPointError(m,0.5,0.5,eff.at(m).at(i)->GetEfficiencyErrorLow(13),eff.at(m).at(i)->GetEfficiencyErrorUp(13));
	}
 if (i < 9)  gr_eff.at(i)->SetMarkerColor(1+i);
 else gr_eff.at(i)->SetMarkerColor(i + 30);
 gr_eff.at(i)->SetMarkerSize(2);
 gr_eff.at(i)->SetMarkerStyle(20+i);
 gr_eff.at(i)->SetTitle(Form("Run %i",time.at(i)));
 mgx->Add(gr_eff.at(i));//,"APE");
// i++;
}
    mgx->SetTitle("Efficiency per mod");
    mgx->Draw("AP");
    mgx->SetMinimum(0.7);
    gPad->BuildLegend(0.55,0.7,0.85,0.9);
    c_gr.SaveAs(Form("/home/espedica/trackeranalysis/results/Alltime_run%i.pdf",run));


TCanvas sync("sync","sync",6000,3000);
sync.Divide(6,3);

for (int m = 0; m < 18; m++) {
sync.cd(m+1);
	for (int i = 0; i < eff.at(m).size(); i++) {


	    if (!eff.at(m).at(i)) continue;  // sicurezza

		eff.at(m).at(i)->SetTitle(Form("eff mod %i",m));

	    if (i < 9) eff.at(m).at(i)->SetLineColor(i + 1);
	    else eff.at(m).at(i)->SetLineColor(i + 30);

	    if (i == 0) {
	        eff.at(m).at(i)->Draw("");
	        gPad->Update();
	        eff.at(m).at(i)->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.7, 1.1);
	    } else {
	        eff.at(m).at(i)->Draw("same");
	    }
 leg.at(m)->Draw();
	}
 }

sync.SaveAs(Form("/home/espedica/trackeranalysis/results/AlltimePerMod_run%i.pdf",run));




TMultiGraph *mgx_1 = new TMultiGraph();
std::array<TGraphAsymmErrors*,18> gr_eff_1;


TLine* l1= new TLine(17.,0.7,17. ,1.01);
l1->SetLineColor(kGray);
l1->SetLineWidth(1);
l1->SetLineStyle(3);
TLine* l2= new TLine(20.,0.7,20. ,1.01);
l2->SetLineColor(kGray);
l2->SetLineWidth(1);
l2->SetLineStyle(3);
TLine* l3= new TLine(23.,0.7,23. ,1.01);
l3->SetLineColor(kGray);
l3->SetLineWidth(1);
l3->SetLineStyle(3);
TLine* l4= new TLine(24.,0.7,24. ,1.01);
l4->SetLineColor(kGray);
l4->SetLineWidth(1);
l4->SetLineStyle(3);
TLine* l5= new TLine(27.,0.7,27. ,1.01);
l5->SetLineColor(kGray);
l5->SetLineWidth(1);
l5->SetLineStyle(3);
TLine* l6= new TLine(28.,0.7,28. ,1.01);
l6->SetLineColor(kGray);
l6->SetLineWidth(1);
l6->SetLineStyle(3);
TLine* l7= new TLine(29.,0.7,29. ,1.01);
l7->SetLineColor(kGray);
l7->SetLineWidth(1);
l7->SetLineStyle(3);
TLine* l8= new TLine(30.,0.7,30. ,1.01);
l8->SetLineColor(kGray);
l8->SetLineWidth(1);
l8->SetLineStyle(3);


TCanvas c_gr_1("c_gr_1","c_gr_1",3500,2000);
c_gr_1.Divide(6,3);
for (int m = 0; m < 18; m++) {
c_gr_1.cd(m+1);
 gr_eff_1.at(m)=new TGraphAsymmErrors(time.size());

	for (int i = 0; i < time.size(); i++) {
//std::cout << "mod " << m << " run " << i << " eff.at(m).at(i)->GetEfficiency(13) " << eff.at(m).at(i)->GetEfficiency(13) << std::endl;
	std::cout << "mod " << m << "time " << time.at(i) << " eff " << eff.at(m).at(i)->GetEfficiency(13) << std::endl;
         gr_eff_1.at(m)->SetPoint(i,time.at(i),eff.at(m).at(i)->GetEfficiency(13));
//         gr_eff_1.at(m)->SetPointError(time.at(i),0.5,0.5,eff.at(m).at(i)->GetEfficiencyErrorLow(13),eff.at(m).at(i)->GetEfficiencyErrorUp(13));
        }

 if (m < 9)  gr_eff_1.at(m)->SetMarkerColor(1+m);
 else gr_eff_1.at(m)->SetMarkerColor(m + 30);

 gr_eff_1.at(m)->SetMarkerSize(2);
 gr_eff_1.at(m)->SetMarkerStyle(20+m);
 gr_eff_1.at(m)->SetMinimum(0.7);
// gr_eff_1.at(m)->GetXaxis()->SetLimits(1752895463,1752895463);
 gr_eff_1.at(m)->SetTitle(Form("Efficiency during time Mod %i",m));
 gr_eff_1.at(m)->GetXaxis()->SetTimeDisplay(1);
 gr_eff_1.at(m)->GetXaxis()->SetLabelOffset(0.03);
 gr_eff_1.at(m)->GetXaxis()->SetTimeFormat("#splitline{%H:%M}{%d/%m}");
 gr_eff_1.at(m)->GetXaxis()->SetTimeOffset(0,"CEST");
 gr_eff_1.at(m)->GetXaxis()->SetTitle("Time");
 gr_eff_1.at(m)->GetYaxis()->SetTitle("Efficiency");
 gr_eff_1.at(m)->Draw("ALP");
 l1->Draw("same");
 l2->Draw("same");
 l3->Draw("same");
 l4->Draw("same");
 l5->Draw("same");
 l6->Draw("same");
 l7->Draw("same");
 l8->Draw("same");
}
/* mgx_1->Add(gr_eff_1.at(m));//,"APE");
// i++;
}
    mgx_1->SetTitle("Efficiency per run");
    mgx_1->Draw("AP");
    mgx_1->SetMinimum(0.7);
    mgx_1->GetXaxis()->SetLimits(6,48);
    gPad->BuildLegend(0.1,0.3,0.25,0.8);*/

    c_gr_1.SaveAs(Form("/home/espedica/trackeranalysis/results/AllMod_run%i.pdf",run));



}
