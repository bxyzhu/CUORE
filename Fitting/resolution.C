{
//=========Macro generated from canvas: resolution/resolution
//=========  (Wed Apr 30 10:54:57 2014) by ROOT version5.34/09
   TCanvas *resolution = new TCanvas("resolution", "resolution",1809,45,1555,1036);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   resolution->Range(-721.5324,1.050733e-05,6675.168,0.005982627);
   resolution->SetFillColor(0);
   resolution->SetBorderMode(0);
   resolution->SetBorderSize(2);
   resolution->SetFrameBorderMode(0);
   resolution->SetFrameBorderMode(0);
   
   TGraphErrors *gre = new TGraphErrors(8);
   gre->SetName("Graph");
   gre->SetTitle("Resolution vs Energy");
   gre->SetFillColor(1);
   gre->SetMarkerStyle(3);
   gre->SetPoint(0,511.4316,0.004627107);
   gre->SetPointError(0,0.1804386,0.0003601674);
   gre->SetPoint(1,911.7092,0.002379718);
   gre->SetPointError(1,0.2004599,0.0002539386);
   gre->SetPoint(2,968.9005,0.002855411);
   gre->SetPointError(2,1.074554,0.00184955);
   gre->SetPoint(3,1461.188,0.001675077);
   gre->SetPointError(3,0.06281654,4.129996e-05);
   gre->SetPoint(4,2615.364,0.001103712);
   gre->SetPointError(4,0.1334105,5.599414e-05);
   gre->SetPoint(5,3291.408,0.002031246);
   gre->SetPointError(5,0.3748868,0.0001529872);
   gre->SetPoint(6,5346.33,0.004780345);
   gre->SetPointError(6,0.8816534,0.0001682918);
   gre->SetPoint(7,5441.71,0.004720157);
   gre->SetPointError(7,0.674658,0.0001142878);
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","Resolution vs Energy",100,18.13775,5935.498);
   Graph_Graph1->SetMinimum(0.0006077194);
   Graph_Graph1->SetMaximum(0.005385415);
   Graph_Graph1->SetDirectory(0);
   Graph_Graph1->SetStats(0);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#000099");
   Graph_Graph1->SetLineColor(ci);
   Graph_Graph1->GetXaxis()->SetTitle("Energy (keV)");
   Graph_Graph1->GetXaxis()->SetRange(1,100);
   Graph_Graph1->GetXaxis()->SetLabelFont(42);
   Graph_Graph1->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetXaxis()->SetTitleFont(42);
   Graph_Graph1->GetYaxis()->SetTitle("Resolution (%)");
   Graph_Graph1->GetYaxis()->SetLabelFont(42);
   Graph_Graph1->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetYaxis()->SetTitleFont(42);
   Graph_Graph1->GetZaxis()->SetLabelFont(42);
   Graph_Graph1->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph1);
   
   
   TPaveStats *ptstats = new TPaveStats(0.7066409,0.8392857,0.9793681,0.9355159,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetTextAlign(12);
   ptstats->SetTextFont(42);
   TText *text = ptstats->AddText("#chi^{2} / ndf = 8.978 / 4");
   text = ptstats->AddText("p0       = 0.007227 #pm 0.0005765 ");
   text = ptstats->AddText("p1       = -6.322e-06 #pm 7.246e-07 ");
   text = ptstats->AddText("p2       = 1.947e-09 #pm 2.67e-10 ");
   text = ptstats->AddText("p3       = -1.596e-13 #pm 2.838e-14 ");
   ptstats->SetOptStat(0);
   ptstats->SetOptFit(111);
   ptstats->Draw();
   gre->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(gre->GetListOfFunctions());
   
   TF1 *PrevFitTMP = new TF1("PrevFitTMP","pol3",18.13775,5935.498);
   PrevFitTMP->SetFillColor(19);
   PrevFitTMP->SetFillStyle(0);
   PrevFitTMP->SetLineColor(2);
   PrevFitTMP->SetLineWidth(2);
   PrevFitTMP->SetChisquare(8.978214);
   PrevFitTMP->SetNDF(4);
   PrevFitTMP->GetXaxis()->SetLabelFont(42);
   PrevFitTMP->GetXaxis()->SetLabelSize(0.035);
   PrevFitTMP->GetXaxis()->SetTitleSize(0.035);
   PrevFitTMP->GetXaxis()->SetTitleFont(42);
   PrevFitTMP->GetYaxis()->SetLabelFont(42);
   PrevFitTMP->GetYaxis()->SetLabelSize(0.035);
   PrevFitTMP->GetYaxis()->SetTitleSize(0.035);
   PrevFitTMP->GetYaxis()->SetTitleFont(42);
   PrevFitTMP->SetParameter(0,0.007226989);
   PrevFitTMP->SetParError(0,0.0005764916);
   PrevFitTMP->SetParLimits(0,0,0);
   PrevFitTMP->SetParameter(1,-6.321793e-06);
   PrevFitTMP->SetParError(1,7.24598e-07);
   PrevFitTMP->SetParLimits(1,0,0);
   PrevFitTMP->SetParameter(2,1.946817e-09);
   PrevFitTMP->SetParError(2,2.67043e-10);
   PrevFitTMP->SetParLimits(2,0,0);
   PrevFitTMP->SetParameter(3,-1.595559e-13);
   PrevFitTMP->SetParError(3,2.837877e-14);
   PrevFitTMP->SetParLimits(3,0,0);
   gre->GetListOfFunctions()->Add(PrevFitTMP);
   gre->Draw("ap");
   
   TPaveText *pt = new TPaveText(0.3175547,0.9351662,0.6824453,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   text = pt->AddText("Resolution vs Energy");
   pt->Draw();
   resolution->Modified();
   resolution->cd();
   resolution->SetSelected(resolution);
}
