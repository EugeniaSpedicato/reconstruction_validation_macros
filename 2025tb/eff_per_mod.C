void eff_per_mod() {
    TFile *f_V530=TFile::Open("/eos/experiment/mu-e/daq/2025_testData/decoded/vcthScan_highIntensity_station1_530.root");
    TFile *f_V535=TFile::Open("/eos/experiment/mu-e/daq/2025_testData/decoded/vcthScan_highIntensity_station1_540.root");
    TFile *f_V540=TFile::Open("/eos/experiment/mu-e/daq/2025_testData/decoded/vcthScan_highIntensity_station1_550.root");
    TFile *f_V545=TFile::Open("/eos/experiment/mu-e/daq/2025_testData/decoded/vcthScan_highIntensity_station1_554.root");
    TFile *f_V550=TFile::Open("/eos/experiment/mu-e/daq/2025_testData/decoded/vcthScan_highIntensity_station1_556.root");
    TFile *f_V555=TFile::Open("/eos/experiment/mu-e/daq/2025_testData/decoded/vcthScan_highIntensity_station1_560.root");
    TFile *f_V560=TFile::Open("/eos/experiment/mu-e/daq/2025_testData/decoded/vcthScan_highIntensity_station1_564.root");
    TFile *f_V565=TFile::Open("/eos/experiment/mu-e/daq/2025_testData/decoded/vcthScan_highIntensity_station1_564.root");
    TFile *f_V570=TFile::Open("/eos/experiment/mu-e/daq/2025_testData/decoded/vcthScan_highIntensity_station1_570.root");
    
    TCanvas* c_530=(TCanvas*)f_V530->Get("c_CLlim_perModule");
    std::array<TGraphAsymmErrors*,6> hmod530;
    for(int i=0; i<6; i++){
        TPad* mod530=(TPad*)c_530->FindObject(Form("c_CLlim_perModule_%i",i+1));
        hmod530.at(i)=(TGraphAsymmErrors*)mod530->FindObject(Form("g_CLlim_perModule_station0_linkID%i",i));
    }
    
    TCanvas* c_535=(TCanvas*)f_V535->Get("c_CLlim_perModule");
    std::array<TGraphAsymmErrors*,6> hmod535;
    for(int i=0; i<6; i++){
        TPad* mod535=(TPad*)c_535->FindObject(Form("c_CLlim_perModule_%i",i+1));
        hmod535.at(i)=(TGraphAsymmErrors*)mod535->FindObject(Form("g_CLlim_perModule_station0_linkID%i",i));
    }
    
    TCanvas* c_540=(TCanvas*)f_V540->Get("c_CLlim_perModule");
    std::array<TGraphAsymmErrors*,6> hmod540;
    for(int i=0; i<6; i++){
        TPad* mod540=(TPad*)c_540->FindObject(Form("c_CLlim_perModule_%i",i+1));
        hmod540.at(i)=(TGraphAsymmErrors*)mod540->FindObject(Form("g_CLlim_perModule_station0_linkID%i",i));
    }
    
    TCanvas* c_545=(TCanvas*)f_V545->Get("c_CLlim_perModule");
    std::array<TGraphAsymmErrors*,6> hmod545;
    for(int i=0; i<6; i++){
        TPad* mod545=(TPad*)c_545->FindObject(Form("c_CLlim_perModule_%i",i+1));
        hmod545.at(i)=(TGraphAsymmErrors*)mod545->FindObject(Form("g_CLlim_perModule_station0_linkID%i",i));
    }
    
    TCanvas* c_550=(TCanvas*)f_V550->Get("c_CLlim_perModule");
    std::array<TGraphAsymmErrors*,6> hmod550;
    for(int i=0; i<6; i++){
        TPad* mod550=(TPad*)c_550->FindObject(Form("c_CLlim_perModule_%i",i+1));
        hmod550.at(i)=(TGraphAsymmErrors*)mod550->FindObject(Form("g_CLlim_perModule_station0_linkID%i",i));
    }
    
    TCanvas* c_555=(TCanvas*)f_V555->Get("c_CLlim_perModule");
    std::array<TGraphAsymmErrors*,6> hmod555;
    for(int i=0; i<6; i++){
        TPad* mod555=(TPad*)c_555->FindObject(Form("c_CLlim_perModule_%i",i+1));
        hmod555.at(i)=(TGraphAsymmErrors*)mod555->FindObject(Form("g_CLlim_perModule_station0_linkID%i",i));
    }
    
    TCanvas* c_560=(TCanvas*)f_V560->Get("c_CLlim_perModule");
    std::array<TGraphAsymmErrors*,6> hmod560;
    for(int i=0; i<6; i++){
        TPad* mod560=(TPad*)c_560->FindObject(Form("c_CLlim_perModule_%i",i+1));
        hmod560.at(i)=(TGraphAsymmErrors*)mod560->FindObject(Form("g_CLlim_perModule_station0_linkID%i",i));
    }
    
    TCanvas* c_565=(TCanvas*)f_V565->Get("c_CLlim_perModule");
    std::array<TGraphAsymmErrors*,6> hmod565;
    for(int i=0; i<6; i++){
        TPad* mod565=(TPad*)c_565->FindObject(Form("c_CLlim_perModule_%i",i+1));
        hmod565.at(i)=(TGraphAsymmErrors*)mod565->FindObject(Form("g_CLlim_perModule_station0_linkID%i",i));
    }
    
    TCanvas* c_570=(TCanvas*)f_V570->Get("c_CLlim_perModule");
    std::array<TGraphAsymmErrors*,6> hmod570;
    for(int i=0; i<6; i++){
        TPad* mod570=(TPad*)c_570->FindObject(Form("c_CLlim_perModule_%i",i+1));
        hmod570.at(i)=(TGraphAsymmErrors*)mod570->FindObject(Form("g_CLlim_perModule_station0_linkID%i",i));
    }
    
    

    TCanvas a1("a1","a1",3000,2500);
    a1.Divide(3,2);
    for(int m=0; m<6; m++){
        a1.cd(m+1);
        hmod530.at(m)->SetLineColor(kRed);
        hmod535.at(m)->SetLineColor(kOrange);
        hmod540.at(m)->SetLineColor(kBlue);
        hmod545.at(m)->SetLineColor(kGreen);
        hmod550.at(m)->SetLineColor(kBlack);
        hmod555.at(m)->SetLineColor(kYellow);
        hmod560.at(m)->SetLineColor(kCyan);
        hmod565.at(m)->SetLineColor(kViolet);
        hmod570.at(m)->SetLineColor(kGray);
        //hmod530.at(m)->GetYaxis()->SetRangeUser(0.85,1.10);
        hmod530.at(m)->Draw();
        hmod535.at(m)->Draw("same");
        hmod540.at(m)->Draw("same");
        hmod545.at(m)->Draw("same");
        hmod550.at(m)->Draw("same");
        hmod555.at(m)->Draw("same");
        hmod560.at(m)->Draw("same");
        hmod565.at(m)->Draw("same");
        hmod570.at(m)->Draw("same");
        TLegend *legend = new TLegend(0.50, 0.05, 0.88, 0.28);
        legend->AddEntry(hmod530.at(m),"Vcth 530","l");
        legend->AddEntry(hmod535.at(m),"Vcth 535","l");
        legend->AddEntry(hmod540.at(m),"Vcth 540","l");
        legend->AddEntry(hmod545.at(m),"Vcth 545","l");
        legend->AddEntry(hmod550.at(m),"Vcth 550","l");
        legend->AddEntry(hmod555.at(m),"Vcth 555","l");
        legend->AddEntry(hmod560.at(m),"Vcth 560","l");
        legend->AddEntry(hmod565.at(m),"Vcth 565","l");
        legend->AddEntry(hmod570.at(m),"Vcth 570","l");
        legend->Draw();

    }
    a1.SaveAs("eff_per_mod.pdf");
    
    
    
    TCanvas a2("a2","a2",3000,3000);
    a2.Divide(3,3);
    a2.cd(1);
    TLegend *legend1 = new TLegend(0.50, 0.05, 0.88, 0.28);
    hmod530.at(0)->SetTitle("Efficiency VCth 530");
    //hmod530.at(0)->GetYaxis()->SetRange(0.85,1.1);
    for(int m=0; m<6; m++){
        hmod530.at(m)->SetLineColor(1+m);
        if(m==0)hmod530.at(m)->Draw();
        hmod530.at(m)->Draw("same");
        TLegend *legend = new TLegend(0.50, 0.05, 0.88, 0.28);
        legend1->AddEntry(hmod530.at(m),Form("Mod %i",m),"l");
    }
    legend1->Draw();
    a2.cd(2);
    hmod535.at(0)->SetTitle("Efficiency VCth 535");
    TLegend *legend2 = new TLegend(0.50, 0.05, 0.88, 0.28);
        for(int m=0; m<6; m++){
            hmod535.at(m)->SetLineColor(1+m);
            if(m==0)hmod535.at(m)->Draw();
            hmod535.at(m)->Draw("same");
            legend2->AddEntry(hmod535.at(m),Form("Mod %i",m),"l");
        }
        legend2->Draw();
    a2.cd(3);
    TLegend *legend3 = new TLegend(0.50, 0.05, 0.88, 0.28);
    hmod540.at(0)->SetTitle("Efficiency VCth 540");
                for(int m=0; m<6; m++){
                    hmod540.at(m)->SetLineColor(1+m);
                    if(m==0)hmod540.at(m)->Draw();
                    hmod540.at(m)->Draw("same");
                    legend3->AddEntry(hmod540.at(m),Form("Mod %i",m),"l");
                }
                legend3->Draw();
    a2.cd(4);
    TLegend *legend4 = new TLegend(0.50, 0.05, 0.88, 0.28);
    hmod545.at(0)->SetTitle("Efficiency VCth 545");
                    for(int m=0; m<6; m++){
                        hmod545.at(m)->SetLineColor(1+m);
                        if(m==0)hmod545.at(m)->Draw();
                        hmod545.at(m)->Draw("same");
                        legend4->AddEntry(hmod545.at(m),Form("Mod %i",m),"l");
                    }
                    legend4->Draw();
    a2.cd(5);
    TLegend *legend5 = new TLegend(0.50, 0.05, 0.88, 0.28);
    hmod550.at(0)->SetTitle("Efficiency VCth 550");
                        for(int m=0; m<6; m++){
                            hmod550.at(m)->SetLineColor(1+m);
                            if(m==0)hmod550.at(m)->Draw();
                            hmod550.at(m)->Draw("same");
                            legend5->AddEntry(hmod550.at(m),Form("Mod %i",m),"l");
                        }
                        legend5->Draw();
    a2.cd(6);
    TLegend *legend6 = new TLegend(0.50, 0.05, 0.88, 0.28);
    hmod555.at(0)->SetTitle("Efficiency VCth 555");
                        for(int m=0; m<6; m++){
                                hmod555.at(m)->SetLineColor(1+m);
                            if(m==0)hmod555.at(m)->Draw();
                            hmod555.at(m)->Draw("same");
                                legend6->AddEntry(hmod555.at(m),Form("Mod %i",m),"l");
                        }
                        legend6->Draw();
    a2.cd(7);
    TLegend *legend7 = new TLegend(0.50, 0.05, 0.88, 0.28);
    hmod560.at(0)->SetTitle("Efficiency VCth 560");
                            for(int m=0; m<6; m++){
                                    hmod560.at(m)->SetLineColor(1+m);
                                if(m==0)hmod560.at(m)->Draw();
                                    hmod560.at(m)->Draw("same");
                                    legend7->AddEntry(hmod560.at(m),Form("Mod %i",m),"l");
                            }
    legend7->Draw();
    a2.cd(8);
    TLegend *legend8 = new TLegend(0.50, 0.05, 0.88, 0.28);
    hmod565.at(0)->SetTitle("Efficiency VCth 565");
                                    for(int m=0; m<6; m++){
                                        hmod565.at(m)->SetLineColor(1+m);
                                        if(m==0)hmod565.at(m)->Draw();
                                        hmod565.at(m)->Draw("same");
                                        legend8->AddEntry(hmod565.at(m),Form("Mod %i",m),"l");
                                    }
                                    legend8->Draw();
    a2.cd(9);
    TLegend *legend9 = new TLegend(0.50, 0.05, 0.88, 0.28);
    hmod570.at(0)->SetTitle("Efficiency VCth 570");
                                        for(int m=0; m<6; m++){
                                            hmod570.at(m)->SetLineColor(1+m);
                                            if(m==0)hmod570.at(m)->Draw();
                                            hmod570.at(m)->Draw("same");
                                            legend9->AddEntry(hmod570.at(m),Form("Mod %i",m),"l");
                                        }
                                        legend9->Draw();
                                                    
                                                
    a2.SaveAs("eff_per_Vcth.pdf");

    std::array<TGraph*,6> eff_v;
    for(int m=0;m<6;m++){
        //        eff_v.at(m)=new TGraph(Form("eff_v_mod%i",m),"Efficiency vs Vcth",9,9);
        eff_v.at(m)=new TGraph(9);
    }
    TMultiGraph *mgx = new TMultiGraph();
    TCanvas a3("a3","a3",1000,1000);
    for(int m=0;m<6;m++){
//        std::cout<<"hmod530.at(m)->GetBinSize and content " << hmod530.at(m)->GetNbinsX() << " " << hmod530.at(m)->GetY()[0] << std::endl;
        eff_v.at(m)->AddPoint(530,hmod530.at(m)->GetY()[0]);
        eff_v.at(m)->AddPoint(535,hmod535.at(m)->GetY()[0]);
        eff_v.at(m)->AddPoint(540,hmod540.at(m)->GetY()[0]);
        eff_v.at(m)->AddPoint(545,hmod545.at(m)->GetY()[0]);
        eff_v.at(m)->AddPoint(550,hmod550.at(m)->GetY()[0]);
        eff_v.at(m)->AddPoint(555,hmod555.at(m)->GetY()[0]);
        eff_v.at(m)->AddPoint(560,hmod560.at(m)->GetY()[0]);
        eff_v.at(m)->AddPoint(565,hmod565.at(m)->GetY()[0]);
        eff_v.at(m)->AddPoint(570,hmod570.at(m)->GetY()[0]);
        eff_v.at(m)->SetTitle(Form("mod%i",m));
        eff_v.at(m)->SetLineColor(1+m);
        mgx->Add(eff_v.at(m),"AL*");
    }
    mgx->SetTitle("Efficiency vs Vcth");
    mgx->Draw("AL*");
    mgx->GetXaxis()->SetLimits(530.,570.);
    mgx->SetMinimum(0.8);
    gPad->BuildLegend(0.55,0.1,0.85,0.3);
    a3.SaveAs("eff_vcth.pdf");

return ;

}
