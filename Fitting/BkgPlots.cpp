
// Fits all gamma peaks from Th232
void FitThPeaks(TH1D *dHisto, bool bSavePlots = false)
{
    // Pb212 238.6
    double dFitMinPb    =   230.;
    double dFitMaxPb    =   247.;

    // 511
    double dFitMinE     =   495;
    double dFitMaxE     =   525;

    // Tl208 583.2
    double dFitMinTl2   =   575.;
    double dFitMaxTl2   =   590.;

    // Tl208 860.56
    double dFitMinTl3   =   854.;
    double dFitMaxTl3   =   868.;

    // Ac228 - 911
    double dFitMinAc2   =   900;
    double dFitMaxAc2   =   920;

    // Ac228 - 968
    double dFitMinAc1   =   958;
    double dFitMaxAc1   =   978;

    // Double Escape - 1593
    double dFitMinDE    =   1570;
    double dFitMaxDE    =   1605;

    // Single Escape - 2104
    double dFitMinSE    =   2094;
    double dFitMaxSE    =   2115;

    // Tl208 - 2615
    double dFitMinTl1   =   2600;
    double dFitMaxTl1   =   2630;


    // Fit functions

    TF1 *fFitPb = new TF1("fFitPb", "gaus(0) + pol1(3)", dFitMinPb, dFitMaxPb);
    fFitPb->SetParameters(10000., 238.6, 2.8, 0., 0.);

    TF1 *fFitE = new TF1("fFitE", "gaus(0) + pol1(3)", dFitMinE, dFitMaxE);
    // fFitE->SetParameters(1000., 510.77., 2.8, 1000., 511, 2.8, 0., 0.);
    fFitE->SetParameters(1000, 510.77, 3., 0., 0.);

    TF1 *fFitTl2 = new TF1("fFitTl2","gaus(0) + pol1(3)", dFitMinTl2, dFitMaxTl2);
    fFitTl2->SetParameters(1000., 583.2, 2.8, 0., 0.);

    TF1 *fFitTl3 = new TF1("fFitTl3","gaus(0) + pol1(3)", dFitMinTl3, dFitMaxTl3);
    fFitTl3->SetParameters(1000., 860.56, 2.8, 0., 0.);

    TF1 *fFitAc2 = new TF1("fFitAc2","gaus(0) + pol1(3)", dFitMinAc2, dFitMaxAc2);
    fFitAc2->SetParameters(1000., 911.2, 2.8, 0., 0.);

    TF1 *fFitAc1 = new TF1("fFitAc1","gaus(0) + gaus(3) + pol1(6)", dFitMinAc1, dFitMaxAc1); // Fit 968 with 2 gaussians
    fFitAc1->SetParameters(50., 964.77, 2.8, 1000., 968.98, 2.8, 0., 0.);


    TF1 *fFitDE = new TF1("fFitDE","gaus(0) + gaus(3) + pol1(6)", dFitMinDE, dFitMaxDE);
    // fFitDE->SetParameters(5000., 1588.19, 2.8, 10000., 1593., 2.8, 0., 0.);
    fFitDE->SetParameters(500, 1580, 3., 5000, 1588, 3., 0., 0.);

    TF1 *fFitSE = new TF1("fFitSE","gaus(0) + pol1(3)", dFitMinSE, dFitMaxSE);
    fFitSE->SetParameters(1000., 2104., 3, 0, 0);

    TF1 *fFitTl1 = new TF1("fFitTl1","gaus(0) + pol1(3)", dFitMinTl1, dFitMaxTl1);
    fFitTl1->SetParameters(1000, 2615, 2.8, 0, 0);



    TCanvas *cth = new TCanvas(Form("%s",dHisto->GetName()), Form("%s",dHisto->GetName()), 1600, 1200);
    cth->Divide(3,3);

    // Fit and draw all peaks
    // TCanvas *cpb = new TCanvas("cpb", "cpb", 1100, 750);
    cth->cd(1);
    TH1D *hpb   = (TH1D*)dHisto->Clone("hpb");
    hpb->SetTitle("Pb212 (238.632 keV)");
    hpb->SetAxisRange(dFitMinPb, dFitMaxPb);
    hpb->Fit("fFitPb","R");


    // TCanvas *ce = new TCanvas("ce", "ce", 1100, 750);
    cth->cd(2);
    TH1D *he   = (TH1D*)dHisto->Clone("he");
    he->SetTitle("Tl208 and electron (510.77 and 511 keV)");
    he->SetAxisRange(dFitMinE, dFitMaxE);
    he->Fit("fFitE","R");

    // TCanvas *ctl2 = new TCanvas("ctl2", "ctl2", 1100, 750);
    cth->cd(3);
    TH1D *htl2   = (TH1D*)dHisto->Clone("htl2");
    htl2->SetTitle("Tl208 (583.2 keV)");
    htl2->SetAxisRange(dFitMinTl2, dFitMaxTl2);
    htl2->Fit("fFitTl2","R");

    // TCanvas *ctl3 = new TCanvas("ctl3", "ctl3", 1100, 750);
    cth->cd(4);
    TH1D *htl3   = (TH1D*)dHisto->Clone("htl3");
    htl3->SetTitle("Tl208 (860.56 keV)");
    htl3->Rebin();
    htl3->SetAxisRange(dFitMinTl3, dFitMaxTl3);
    htl3->Fit("fFitTl3","R");

    // TCanvas *cac2 = new TCanvas("cac2", "cac2", 1100, 750);
    cth->cd(5);
    TH1D *hac2   = (TH1D*)dHisto->Clone("hac2");
    hac2->SetTitle("Ac228 (911.2 keV)");
    hac2->SetAxisRange(dFitMinAc2, dFitMaxAc2);
    hac2->Fit("fFitAc2","R");

    // TCanvas *cac1 = new TCanvas("cac1", "cac1", 1100, 750);
    cth->cd(6);
    TH1D *hac1   = (TH1D*)dHisto->Clone("hac1");
    hac1->SetTitle("Ac228 (964.77 and 968.98 keV)");
    hac1->SetAxisRange(dFitMinAc1, dFitMaxAc1);
    hac1->Fit("fFitAc1","R");

    // TCanvas *cde = new TCanvas("cde", "cde", 1100, 750);
    cth->cd(7);
    TH1D *hde   = (TH1D*)dHisto->Clone("hde");
    hde->SetTitle("Ac228 (1580.53 and 1588.19)");
    hde->Rebin(1.5);
    hde->SetAxisRange(dFitMinDE, dFitMaxDE);
    hde->Fit("fFitDE","R");

    // TCanvas *cse = new TCanvas("cse", "cse", 1100, 750);
    cth->cd(8);
    TH1D *hse   = (TH1D*)dHisto->Clone("hse");
    hse->SetTitle("Tl208 (2104 keV)");
    hse->Rebin();
    hse->SetAxisRange(dFitMinSE, dFitMaxSE);
    hse->Fit("fFitSE","R");
    // hse->Fit("gaus");

    // TCanvas *ctl1 = new TCanvas("ctl1", "ctl1", 1100, 750);
    cth->cd(9);
    TH1D *htl1  = (TH1D*)dHisto->Clone("htl1");
    htl1->SetTitle("Tl208 (2615 keV)");
    htl1->SetAxisRange(dFitMinTl1, dFitMaxTl1);
    htl1->Fit("fFitTl1","R");




    // Ratios
    TCanvas *cratioth = new TCanvas(Form("Th_Ratio%s",dHisto->GetName()), Form("Th_Ratio%s",dHisto->GetName()), 1200, 800);

    double dFitEnergy[11] = {};
    double dFitEnergyErr[11] = {};
    double dFitNorm[11] = {};
    double dFitArea[11] = {};
    double dFitAreaErr[11] = {};    


    // Mean
    dFitEnergy[0] = fFitPb->GetParameter(1);
    dFitEnergy[1] = fFitE->GetParameter(1);
    dFitEnergy[2] = fFitTl2->GetParameter(1);
    dFitEnergy[3] = fFitTl3->GetParameter(1);
    dFitEnergy[4] = fFitAc2->GetParameter(1);
    dFitEnergy[5] = fFitAc1->GetParameter(1);
    dFitEnergy[6] = fFitAc1->GetParameter(4);
    dFitEnergy[7] = fFitDE->GetParameter(1);
    dFitEnergy[8] = fFitDE->GetParameter(4);
    dFitEnergy[9] = fFitSE->GetParameter(1);
    dFitEnergy[10] = fFitTl1->GetParameter(1);

    // Total Area
    double dTotAreaTl = TMath::Abs(fFitTl2->GetParameter(0)/fFitTl2->GetParameter(2) 
            + fFitE->GetParameter(0)/fFitE->GetParameter(2)
            + fFitTl3->GetParameter(0)/fFitTl3->GetParameter(2) + fFitTl1->GetParameter(0)/fFitTl1->GetParameter(2));

    double dTotAreaAc = TMath::Abs(fFitAc2->GetParameter(0)/fFitAc2->GetParameter(2) + fFitAc1->GetParameter(0)/fFitAc1->GetParameter(2)
            + fFitAc1->GetParameter(3)/fFitAc1->GetParameter(5) + fFitDE->GetParameter(0)/fFitDE->GetParameter(2)
            + fFitDE->GetParameter(3)/fFitDE->GetParameter(5)); 

    double dTotAreaPb = TMath::Abs(fFitPb->GetParameter(0)/fFitPb->GetParameter(2));


    // Norm/Sigma, SE and DE just set to be 1
    dFitNorm[0] = TMath::Abs(fFitPb->GetParameter(0)/fFitPb->GetParameter(2));
    // dFitNorm[1] = TMath::Abs(fFitE->GetParameter(0)/fFitE->GetParameter(2));
    dFitNorm[1] = 1;
    dFitNorm[2] = TMath::Abs(fFitTl2->GetParameter(0)/fFitTl2->GetParameter(2));
    dFitNorm[3] = TMath::Abs(fFitTl3->GetParameter(0)/fFitTl3->GetParameter(2));
    dFitNorm[4] = TMath::Abs(fFitAc2->GetParameter(0)/fFitAc2->GetParameter(2));
    dFitNorm[5] = TMath::Abs(fFitAc1->GetParameter(0)/fFitAc1->GetParameter(2));
    dFitNorm[6] = TMath::Abs(fFitAc1->GetParameter(3)/fFitAc1->GetParameter(5));
    dFitNorm[7] = TMath::Abs(fFitDE->GetParameter(0)/fFitDE->GetParameter(2));
    dFitNorm[8] = TMath::Abs(fFitDE->GetParameter(3)/fFitDE->GetParameter(5));
    dFitNorm[9] = 1;
    dFitNorm[10] = TMath::Abs(fFitTl1->GetParameter(0)/fFitTl1->GetParameter(2));



    // Total Tl BR is 195.92 without e- peak, 218.52 with e-
/*
    // with e-
    dFitArea[0] = dFitNorm[0]/dTotAreaPb;
    dFitArea[1] = dFitNorm[1]/dTotAreaTl * 218.52/22.6;
    dFitArea[2] = dFitNorm[2]/dTotAreaTl * 218.52/84.5;
    dFitArea[3] = dFitNorm[3]/dTotAreaTl * 218.52/12.42;
    dFitArea[4] = dFitNorm[4]/dTotAreaAc * 49.81/25.8;
    dFitArea[5] = dFitNorm[5]/dTotAreaAc * 49.81/4.99;
    dFitArea[6] = dFitNorm[6]/dTotAreaAc * 49.81/15.8;
    dFitArea[7] = 1;
    dFitArea[8] = dFitNorm[8]/dTotAreaAc * 49.81/3.22;
    dFitArea[9] = 1;
    dFitArea[10] = dFitNorm[10]/dTotAreaTl * 218.52/99.0;
*/

    // without e-
    dFitArea[0] = dFitNorm[0]/dTotAreaPb;
    dFitArea[1] = 1;
    dFitArea[2] = dFitNorm[2]/dTotAreaTl * 195.92/84.5;
    dFitArea[3] = dFitNorm[3]/dTotAreaTl * 195.92/12.42;
    dFitArea[4] = dFitNorm[4]/dTotAreaAc * 50.41/25.8;
    dFitArea[5] = dFitNorm[5]/dTotAreaAc * 50.41/4.99;
    dFitArea[6] = dFitNorm[6]/dTotAreaAc * 50.41/15.8;
    dFitArea[7] = dFitNorm[7]/dTotAreaAc * 50.41/0.6;
    dFitArea[8] = dFitNorm[8]/dTotAreaAc * 50.41/3.22;
    dFitArea[9] = 1;
    dFitArea[10] = dFitNorm[10]/dTotAreaTl * 195.92/99.0;


    // Error: Ratio * sqrt(Tot - Area/(Tot*Area))
    dFitAreaErr[0] = dFitArea[0]*TMath::Sqrt((dTotAreaPb - dFitNorm[0])/(dTotAreaPb*dFitNorm[0]));
    // dFitAreaErr[1] = dFitArea[1]*TMath::Sqrt((dTotAreaTl - dFitNorm[1])/(dTotAreaTl*dFitNorm[1]));
    dFitAreaErr[1] = 0;
    dFitAreaErr[2] = dFitArea[2]*TMath::Sqrt((dTotAreaTl - dFitNorm[2])/(dTotAreaTl*dFitNorm[2]));
    dFitAreaErr[3] = dFitArea[3]*TMath::Sqrt((dTotAreaTl - dFitNorm[3])/(dTotAreaTl*dFitNorm[3]));
    dFitAreaErr[4] = dFitArea[4]*TMath::Sqrt((dTotAreaAc - dFitNorm[4])/(dTotAreaAc*dFitNorm[4]));
    dFitAreaErr[5] = dFitArea[5]*TMath::Sqrt((dTotAreaAc - dFitNorm[5])/(dTotAreaAc*dFitNorm[5]));
    dFitAreaErr[6] = dFitArea[6]*TMath::Sqrt((dTotAreaAc - dFitNorm[6])/(dTotAreaAc*dFitNorm[6]));
    dFitAreaErr[7] = dFitArea[7]*TMath::Sqrt((dTotAreaAc - dFitNorm[7])/(dTotAreaAc*dFitNorm[7]));
    dFitAreaErr[8] = dFitArea[8]*TMath::Sqrt((dTotAreaAc - dFitNorm[8])/(dTotAreaAc*dFitNorm[8]));
    dFitAreaErr[9] = 0;
    dFitAreaErr[10] = dFitArea[10]*TMath::Sqrt((dTotAreaTl - dFitNorm[10])/(dTotAreaTl*dFitNorm[10]));



    TGraphErrors *g1 = new TGraphErrors(11, dFitEnergy, dFitArea, dFitEnergyErr, dFitAreaErr);
    g1->SetMarkerStyle(21);
    g1->SetMarkerColor(4);
    g1->Draw("AP");
    g1->GetXaxis()->SetTitle("Energy (keV)");
    g1->GetYaxis()->SetTitle("Area/Total Area * Total Branching Ratio/Branching Ratio");

    TLine *line = new TLine();
    line->SetLineStyle(10);
    line->DrawLine(100, 1, 2800, 1);

    TLegend *leg;
    leg = new TLegend(0.72,0.75,0.9,0.9);

    double dSERatio = fFitSE->GetParameter(0)/fFitSE->GetParameter(2)/ (fFitTl1->GetParameter(0)/fFitTl1->GetParameter(2));

    leg->AddEntry((TObject*)0, Form("Single Escape/(2615 keV) = %.3f", dSERatio), "");
    // leg->AddEntry((TObject*)0, Form("Double Escape Ratio = %.3f%", dDERatio), "");

    leg->Draw();

    if(bSavePlots)
    {
        cth->SaveAs(Form("%s-Th232.png",dHisto->GetName()));
        cth->SaveAs(Form("%s-Th232.pdf",dHisto->GetName()));


    }



}

// Fits gamma peaks from Ra226 chain
void FitRaPeaks(TH1D *dHisto, bool bSavePlots = false)
{
    // Pb214 241.997
    double dFitMinPb1    =   234.;
    double dFitMaxPb1    =   250.;

    // Pb214 295.224
    double dFitMinPb2    =   288.;
    double dFitMaxPb2    =   302.;

    // Pb214 351.932
    double dFitMinPb3    =   344.;
    double dFitMaxPb3    =   360.;

    // Bi214 609.31
    double dFitMinBi1    =   601.;
    double dFitMaxBi1    =   618.;

    // Bi214 665.45
    double dFitMinBi2    =   658.;
    double dFitMaxBi2    =   672.;

    // Bi214 768.36
    double dFitMinBi3    =   760.;
    double dFitMaxBi3    =   776.;

    // Bi214 806.17
    double dFitMinBi4    =   798.;
    double dFitMaxBi4    =   814.;

    // Bi214 934.06
    double dFitMinBi5    =   926.;
    double dFitMaxBi5    =   942.;

    // Bi214 1120.29
    double dFitMinBi6    =   1110.;
    double dFitMaxBi6    =   1128.;

    // Bi214 1155.19
    double dFitMinBi7    =   1148.;
    double dFitMaxBi7    =   1157.;

    // Bi214 1238.11
    double dFitMinBi8    =   1230.;
    double dFitMaxBi8    =   1246.;

    // Bi214 1377.67
    double dFitMinBi9    =   1370.;
    double dFitMaxBi9    =   1385.;

    // Bi214 1407.98
    double dFitMinBi10   =   1394.;
    double dFitMaxBi10   =   1420.;

    // Bi214 1729.6
    double dFitMinBi11   =   1720.;
    double dFitMaxBi11   =   1738.;

    // Bi214 1764.49
    double dFitMinBi12   =   1756.;
    double dFitMaxBi12   =   1772.;

    // Bi214 1847.42
    double dFitMinBi13   =   1841.;
    double dFitMaxBi13   =   1855.;

    // Bi214 2204.21
    double dFitMinBi14   =   2196.;
    double dFitMaxBi14   =   2212.;

    // Bi214 2447.86
    double dFitMinBi15   =   2440.;
    double dFitMaxBi15   =   2455.;

    // Fit functions

    TF1 *fFitPb1 = new TF1("fFitPb1", "gaus(0) + pol1(3)", dFitMinPb1, dFitMaxPb1);
    fFitPb1->SetParameters(100., 241.997, 2.8, 0., 0.);

    TF1 *fFitPb2 = new TF1("fFitPb2", "gaus(0) + pol1(3)", dFitMinPb2, dFitMaxPb2);
    fFitPb2->SetParameters(100., 295.224, 2.8, 0., 0.);

    TF1 *fFitPb3 = new TF1("fFitPb3", "gaus(0) + pol1(3)", dFitMinPb3, dFitMaxPb3);
    fFitPb3->SetParameters(100., 351.932, 2.8, 0., 0.);

    TF1 *fFitBi1 = new TF1("fFitBi1","gaus(0) + pol1(3)", dFitMinBi1, dFitMaxBi1);
    fFitBi1->SetParameters(1000, 609.31, 2.8, 0, 0);

    TF1 *fFitBi2 = new TF1("fFitBi2","gaus(0) + pol1(3)", dFitMinBi2, dFitMaxBi2);
    fFitBi2->SetParameters(1000, 665.45, 2.8, 0, 0);

    TF1 *fFitBi3 = new TF1("fFitBi3","gaus(0) + pol1(3)", dFitMinBi3, dFitMaxBi3);
    fFitBi3->SetParameters(1000, 768.36, 2.8, 0, 0);

    TF1 *fFitBi4 = new TF1("fFitBi4","gaus(0) + pol1(3)", dFitMinBi4, dFitMaxBi4);
    fFitBi4->SetParameters(1000, 806.17, 2.8, 0, 0);

    TF1 *fFitBi5 = new TF1("fFitBi5","gaus(0) + pol1(3)", dFitMinBi5, dFitMaxBi5);
    fFitBi5->SetParameters(1000, 934.06, 2.8, 0, 0);    

    TF1 *fFitBi6 = new TF1("fFitBi6","gaus(0) + pol1(3)", dFitMinBi6, dFitMaxBi6);
    fFitBi6->SetParameters(1000, 1120.29, 2.8, 0, 0);

    TF1 *fFitBi7 = new TF1("fFitBi7","gaus(0) + pol1(3)", dFitMinBi7, dFitMaxBi7);
    fFitBi7->SetParameters(100, 1155.19, 3.5, 0, 0);

    TF1 *fFitBi8 = new TF1("fFitBi8","gaus(0) + pol1(3)", dFitMinBi8, dFitMaxBi8);
    fFitBi8->SetParameters(1000, 1238.11, 2.8, 0, 0);    

    TF1 *fFitBi9 = new TF1("fFitBi9","gaus(0) + pol1(3)", dFitMinBi9, dFitMaxBi9);
    fFitBi9->SetParameters(100, 1377.67, 2.8, 0, 0);

    TF1 *fFitBi10 = new TF1("fFitBi10","gaus(0) + gaus(3) + pol1(6)", dFitMinBi10, dFitMaxBi10);
    fFitBi10->SetParameters(40, 1401.515, 2.8, 100, 1407.98, 2.8, 0, 0);

    TF1 *fFitBi11 = new TF1("fFitBi11","gaus(0) + pol1(3)", dFitMinBi11, dFitMaxBi11);
    fFitBi11->SetParameters(1000, 1729.6, 2.8, 0, 0);

    TF1 *fFitBi12 = new TF1("fFitBi12","gaus(0) + pol1(3)", dFitMinBi12, dFitMaxBi12);
    fFitBi12->SetParameters(1000, 1764.49, 2.8, 0, 0);

    TF1 *fFitBi13 = new TF1("fFitBi13","gaus(0) + pol1(3)", dFitMinBi13, dFitMaxBi13);
    fFitBi13->SetParameters(10, 1847.42, 2.8, 0, 0);

    TF1 *fFitBi14 = new TF1("fFitBi14","gaus(0) + pol1(3)", dFitMinBi14, dFitMaxBi14);
    fFitBi14->SetParameters(1000, 2204.21, 2.8, 0, 0);

    TF1 *fFitBi15 = new TF1("fFitBi15","gaus(0) + pol1(3)", dFitMinBi15, dFitMaxBi15);
    fFitBi15->SetParameters(10, 2447.86, 2.8, 0, 0);

            

    TCanvas *cra1 = new TCanvas(Form("%s1",dHisto->GetName()), Form("%s1",dHisto->GetName()), 1600, 1200);
    cra1->Divide(3,3);

    TCanvas *cra2 = new TCanvas(Form("%s2",dHisto->GetName()), Form("%s2",dHisto->GetName()), 1600, 1200);
    cra2->Divide(3,3);


    // Fit and draw all peaks
    // TCanvas *cpb = new TCanvas("cpb", "cpb", 1100, 750);
    cra1->cd(1);
    TH1D *hpb1   = (TH1D*)dHisto->Clone("hpb1");
    hpb1->SetTitle("Pb214 (241.997 keV)");
    hpb1->SetAxisRange(dFitMinPb1, dFitMaxPb1);
    hpb1->Fit("fFitPb1","R");

    cra1->cd(2);
    TH1D *hpb2   = (TH1D*)dHisto->Clone("hpb2");
    hpb2->SetTitle("Pb214 (295.224 keV)");
    // hpb2->Rebin();
    hpb2->SetAxisRange(dFitMinPb2, dFitMaxPb2);
    hpb2->Fit("fFitPb2","R");

    cra1->cd(3);
    TH1D *hpb3   = (TH1D*)dHisto->Clone("hpb1");
    hpb3->SetTitle("Pb214 (351.932 keV)");
    hpb3->SetAxisRange(dFitMinPb3, dFitMaxPb3);
    hpb3->Fit("fFitPb3","R");

    cra1->cd(4);
    TH1D *hbi1   = (TH1D*)dHisto->Clone("hbi1");
    hbi1->SetTitle("Bi214 (609.31 keV)");
    hbi1->SetAxisRange(dFitMinBi1, dFitMaxBi1);
    hbi1->Fit("fFitBi1","R");

    cra1->cd(5);
    TH1D *hbi2   = (TH1D*)dHisto->Clone("hbi2");
    hbi2->SetTitle("Bi214 (665.45 keV)");
    hbi2->SetAxisRange(dFitMinBi2, dFitMaxBi2);
    hbi2->Fit("fFitBi2","R");

    cra1->cd(6);
    TH1D *hbi3   = (TH1D*)dHisto->Clone("hbi3");
    hbi3->SetTitle("Bi214 (768.36 keV)");
    hbi3->SetAxisRange(dFitMinBi3, dFitMaxBi3);
    hbi3->Fit("fFitBi3","R");

    cra1->cd(7);
    TH1D *hbi4   = (TH1D*)dHisto->Clone("hbi4");
    hbi4->SetTitle("Bi214 (806.17 keV)");
    hbi4->SetAxisRange(dFitMinBi4, dFitMaxBi4);
    hbi4->Fit("fFitBi4","R");

    cra1->cd(8);
    TH1D *hbi5   = (TH1D*)dHisto->Clone("hbi5");
    hbi5->SetTitle("Bi214 (934.06 keV)");
    hbi5->SetAxisRange(dFitMinBi5, dFitMaxBi5);
    hbi5->Fit("fFitBi5","R");

    cra1->cd(9);
    TH1D *hbi6   = (TH1D*)dHisto->Clone("hbi6");
    hbi6->SetTitle("Bi214 (1120.29 keV)");
    hbi6->SetAxisRange(dFitMinBi6, dFitMaxBi6);
    hbi6->Fit("fFitBi6","R");

    cra2->cd(1);
    TH1D *hbi7   = (TH1D*)dHisto->Clone("hbi7");
    hbi7->SetTitle("Bi214 (1155.19 keV)");
    // hbi7->Rebin();
    hbi7->SetAxisRange(dFitMinBi7, dFitMaxBi7+4);
    hbi7->Fit("fFitBi7","R");

    cra2->cd(2);
    TH1D *hbi8   = (TH1D*)dHisto->Clone("hbi8");
    hbi8->SetTitle("Bi214 (1238.11 keV)");
    hbi8->SetAxisRange(dFitMinBi8, dFitMaxBi8);
    hbi8->Fit("fFitBi8","R");

    cra2->cd(3);
    TH1D *hbi9   = (TH1D*)dHisto->Clone("hbi9");
    hbi9->SetTitle("Bi214 (1377.67 keV)");
    hbi9->SetAxisRange(dFitMinBi9, dFitMaxBi9);
    hbi9->Fit("fFitBi9","R");

    cra2->cd(4);
    TH1D *hbi10   = (TH1D*)dHisto->Clone("hbi10");
    hbi10->SetTitle("Bi214 (1401.515 and 1407.98 keV)");
    hbi10->Rebin();
    hbi10->SetAxisRange(dFitMinBi10, dFitMaxBi10);
    hbi10->Fit("fFitBi10","R");

    cra2->cd(5);
    TH1D *hbi11   = (TH1D*)dHisto->Clone("hbi11");
    hbi11->SetTitle("Bi214 (1729.6 keV)");
    hbi11->SetAxisRange(dFitMinBi11, dFitMaxBi11);
    hbi11->Fit("fFitBi11","R");

    cra2->cd(6);
    TH1D *hbi12   = (TH1D*)dHisto->Clone("hbi12");
    hbi12->SetTitle("Bi214 (1764.49 keV)");
    hbi12->SetAxisRange(dFitMinBi12, dFitMaxBi12);
    hbi12->Fit("fFitBi12","R");

    cra2->cd(7);
    TH1D *hbi13   = (TH1D*)dHisto->Clone("hbi13");
    hbi13->SetTitle("Bi214 (1847.42 keV)");
    hbi13->SetAxisRange(dFitMinBi13, dFitMaxBi13);
    hbi13->Fit("fFitBi13","R");

    cra2->cd(8);
    TH1D *hbi14   = (TH1D*)dHisto->Clone("hbi14");
    hbi14->SetTitle("Bi214 (2204.21 keV)");
    hbi14->SetAxisRange(dFitMinBi14, dFitMaxBi14);
    hbi14->Fit("fFitBi14","R");

    cra2->cd(9);
    TH1D *hbi15   = (TH1D*)dHisto->Clone("hbi15");
    hbi15->SetTitle("Bi214 (2447.86 keV)");
    hbi15->Rebin();
    hbi15->SetAxisRange(dFitMinBi15, dFitMaxBi15);
    hbi15->Fit("fFitBi15","R");





    if(bSavePlots)
    {
        cra1->SaveAs(Form("%s-Ra226-1.png",dHisto->GetName()));
        cra1->SaveAs(Form("%s-Ra226-1.pdf",dHisto->GetName()));
        cra2->SaveAs(Form("%s-Ra226-2.png",dHisto->GetName()));
        cra2->SaveAs(Form("%s-Ra226-2.pdf",dHisto->GetName()));

    }





}


void DrawBkg(int dMult = 1, bool bSavePlots = false)
{
    gStyle->SetOptStat(0);
	gStyle->SetOptFit();	

    int dMult = 1;

    int bin = 3500;
    double binsize = 3500/bin;

    TH1D *hBkg = new TH1D("hBkg", "", bin, 0, 3500);

    TChain *qtree = new TChain("qtree");
    qtree->Add("/Users/brian/macros/CUOREZ/Bkg/ReducedBkg-ds*.root");

    TCut base_cut;
    base_cut = base_cut && "(TimeUntilSignalEvent_SameChannel > 4.0 || TimeUntilSignalEvent_SameChannel < 0)";
    base_cut = base_cut && "(TimeSinceSignalEvent_SameChannel > 3.1 || TimeSinceSignalEvent_SameChannel < 0)";
    base_cut = base_cut && "abs(BaselineSlope)<0.1";
    base_cut = base_cut && "OF_TVR < 1.75 && OF_TVL < 2.05";


    qtree->Project("hBkg", "Energy", base_cut && "Multiplicity_OFTime == 1");

    FitThPeaks(hBkg, bSavePlots);
    // FitRaPeaks(hBkg, bSavePlots);


}
