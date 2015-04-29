// Fitter

void TBkgFitter()
{
	




}


// Draws background spectrum
void TBkgModel::DrawBkg()
{

 	gStyle->SetOptStat(0);
  gStyle->SetOptFit();
 	// gStyle->SetOptTitle(0);	

  TLegend *leg = new TLegend(0.75,0.75,0.97,0.97);
  leg->AddEntry(fDataHistoTot, "Total", "l");
  leg->AddEntry(fDataHistoM1, "M1", "l");
  leg->AddEntry(fDataHistoM2, "M2", "l");
  leg->Draw();

  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 1800);
  c1->Divide(1, 2);
  c1->cd(1);
  gPad->SetLogy();
  fAdapDataHistoM1->GetXaxis()->SetTitle("Energy (keV)");
  fAdapDataHistoM1->GetYaxis()->SetTitle("Counts/bin");  
  fAdapDataHistoM1->Draw("");


  c1->cd(2); 
  gPad->SetLogy();
  fAdapDataHistoM2->GetXaxis()->SetTitle("Energy (keV)");
  fAdapDataHistoM2->GetYaxis()->SetTitle("Counts/bin");   
  fAdapDataHistoM2->Draw("");

}

TH1D *TBkgModel::CalculateResidualsAdaptive(TH1D *h1, TH1D *h2, TH1D *hResid, int binMin, int binMax, int dMult)
{

  if(binMin >= binMax)
  {
    cout << " Residuals: min bin >= max bin" << endl;
  }

  if(dMult == 1)
  {
  TH1D  *hOut       = new TH1D("hOutResidualM1", "Fit Residuals M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  }
  if(dMult == 2)
  {
  TH1D  *hOut       = new TH1D("hOutResidualM2", "Fit Residuals M2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  }
  if(dMult == 3)
  {
  TH1D  *hOut       = new TH1D("hOutResidualM2Sum", "Fit Residuals M2Sum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  }

  // Clone histograms for rebinning
  TH1D  *hCloneBkg    = (TH1D*)h1->Clone("hCloneBkg");
  TH1D  *hCloneMC   = (TH1D*)h2->Clone("hCloneMC");

  // Variables used in Residual calculations
  double dResidualX = 0, dResidualY = 0, dResidualXErr = 0, dResidualYErr = 0;

  double datam1_i = 0, modelm1_i = 0;

  // Residual plot and distribution
  for(int j = binMin ; j < binMax ; j++)
  {

    if( hCloneBkg->GetBinCenter(j) >= 3150 && hCloneBkg->GetBinCenter(j) <= 3400)continue;    
    // if( hCloneBkg->GetBinCenter(j) >= 800 && hCloneBkg->GetBinCenter(j) <= 808)continue;
    // if( hCloneBkg->GetBinCenter(j) >= 1060 && hCloneBkg->GetBinCenter(j) <= 1068)continue;

    // if( hCloneBkg->GetBinCenter(j) >= 506 && hCloneBkg->GetBinCenter(j) <= 515)continue;
    // if( hCloneBkg->GetBinCenter(j) >= 579 && hCloneBkg->GetBinCenter(j) <= 589)continue;
    // if( hCloneBkg->GetBinCenter(j) >= 605 && hCloneBkg->GetBinCenter(j) <= 615)continue;
    // if( hCloneBkg->GetBinCenter(j) >= 906 && hCloneBkg->GetBinCenter(j) <= 917)continue;
    // if( hCloneBkg->GetBinCenter(j) >= 1450 && hCloneBkg->GetBinCenter(j) <= 1475)continue;
    // if( hCloneBkg->GetBinCenter(j) >= 1755 && hCloneBkg->GetBinCenter(j) <= 1780)continue;
    // if( hCloneBkg->GetBinCenter(j) >= 2090 && hCloneBkg->GetBinCenter(j) <= 2130)continue;
    // if( hCloneBkg->GetBinCenter(j) >= 2200 && hCloneBkg->GetBinCenter(j) <= 2220)continue;
    // if( hCloneBkg->GetBinCenter(j) >= 2440 && hCloneBkg->GetBinCenter(j) <= 2450)continue;
    // if( hCloneBkg->GetBinCenter(j) >= 2600 && hCloneBkg->GetBinCenter(j) <= 2630)continue;

    dResidualX    = hCloneBkg->GetBinCenter(j);

    datam1_i = h1->GetBinContent(j)*h1->GetBinWidth(j);
    modelm1_i = h2->GetBinContent(j)*h1->GetBinWidth(j);
    
    // Re-multiply bin content by bin width (for # of counts)
    dResidualY = (datam1_i - modelm1_i)/TMath::Sqrt(datam1_i);

    hOut->SetBinContent(j, dResidualY);
    hOut->SetBinError(j, 0.01);
    hResid->Fill(dResidualY);
  }

  return hOut;
}