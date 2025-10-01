#include "Style_jaebeom.h"
#include "sPhenixStyle.h"

void draw_photon_efficiency() {
  const char * runtext = "#it{p+p #sqrt{s_{N}}=200 GeV}";
  TFile * f = TFile::Open(Form("data_hists/mass_data.root"),"READ");
  const char * type = "pt";
  TH1D * h = (TH1D*)f->Get( "hmbd");
  TH1D * h3 = (TH1D*)f->Get("hphoton3");
  TH1D * h4 = (TH1D*)f->Get("hphoton4");
  TH1D * h5 = (TH1D*)f->Get("hphoton5");
  
  h->Rebin(5);
  h3->Rebin(5);
  h4->Rebin(5);
  h5->Rebin(5);

  h->GetXaxis()->SetRangeUser(0,20);
  h3->GetXaxis()->SetRangeUser(0,20);
  h4->GetXaxis()->SetRangeUser(0,20);
  h5->GetXaxis()->SetRangeUser(0,20);
  
  h->SetMarkerStyle(20);
  h3->SetMarkerStyle(24);
  h4->SetMarkerStyle(24);
  h5->SetMarkerStyle(24);

  h->SetMarkerColor(kBlack);
  h3->SetMarkerColor(kRed+1);
  h4->SetMarkerColor(kBlue+1);
  h5->SetMarkerColor(kGreen+2);
  
  h->SetLineColor(kBlack);
  h3->SetLineColor(kRed+1);
  h4->SetLineColor(kBlue+1);
  h5->SetLineColor(kGreen+2);

  TCanvas * c = new TCanvas("c","",600,800);

  TPad * pcounts = new TPad("pcounts","",0,.33,1,1);
  TPad * peff = new TPad("peff","",0,0,1,.33);

  pcounts->cd();
  pcounts->SetBottomMargin(0.00);
  pcounts->SetLeftMargin(0.1); 
  gPad->SetLogy();
  gPad->SetTicks(1,1);
  gStyle->SetOptStat(0);
  for (int i = 0 ; i < h->GetNbinsX(); i++) {
    h->SetBinError(i,TMath::Sqrt(h->GetBinContent(i)));
    h3->SetBinError(i,TMath::Sqrt(h3->GetBinContent(i)));
    h4->SetBinError(i,TMath::Sqrt(h4->GetBinContent(i)));
    h5->SetBinError(i,TMath::Sqrt(h5->GetBinContent(i)));
  }

  h->GetXaxis()->SetTitle();
  h3->GetXaxis()->SetTitle();
  h4->GetXaxis()->SetTitle();
  h5->GetXaxis()->SetTitle();
  
  h->GetXaxis()->SetLabelSize(0);
  h3->GetXaxis()->SetLabelSize(0);
  h4->GetXaxis()->SetLabelSize(0);
  h5->GetXaxis()->SetLabelSize(0);

  h->GetYaxis()->SetRangeUser(1,h->GetMaximum() * 100000);

  h->Draw("e p");
  h3->Draw("e p same");
  h4->Draw("e p same");
  h5->Draw("e p same");

  TLegend * l = new TLegend(.45,.5,.88,.68);
  l->SetLineWidth(0);
  l->SetFillStyle(0);
  l->AddEntry(h,"MBD triggered events");
  l->AddEntry(h3,"MBD + Photon 3 GeV");
  l->AddEntry(h4,"MBD + Photon 4 GeV");
  l->AddEntry(h5,"MBD + Photon 5 GeV");
  l->Draw();

  drawText("#bf{#it{sPHENIX}} Internal",0.5,0.82,1,22);
  drawText("#it{p+p} #sqrt{s_{N}} = 200 GeV",0.5,0.77,1,18);
  //drawText(Form("Run %s",runtext),0.3,0.72,1,18);
  
  drawText("|v_{z}| < 30 cm",0.55,0.45,1,18);
  drawText("|y^{#gamma}| < 1",0.55,0.40,1,18);
  drawText("#mathcal{L}_{int}=20.9 nb%{-1}",0.55,0.35,1,18);
  //drawText("prob_#gamma > 0.05",0.55,0.40,1,16);
  
  peff->cd();
  peff->SetTopMargin(0);
  peff->SetBottomMargin(0.3);
  peff->SetLeftMargin(0.1); 
  gPad->SetTicks(1,1);
  gStyle->SetOptStat(0);

  TEfficiency * h3divide = new TEfficiency(*h3,*h);
  TEfficiency * h4divide = new TEfficiency(*h4,*h);
  TEfficiency * h5divide = new TEfficiency(*h5,*h);
  h3divide->SetName("photon3eff");
  h4divide->SetName("photon4eff");
  h5divide->SetName("photon5eff");
  TH1D * dummy = new TH1D("dummy",";Max p_{T,cluster} [GeV];Ratio",100,0,20);
  dummy->GetYaxis()->SetRangeUser(0,2);
  dummy->GetXaxis()->SetTitleSize(0.06);
  dummy->GetXaxis()->SetTitleOffset(1.5);
  dummy->GetXaxis()->SetLabelSize(0.08);
  dummy->GetYaxis()->SetTitleSize(0.07);
  dummy->GetYaxis()->SetLabelSize(0.08);
  dummy->GetYaxis()->SetTitleOffset(.62);

  h3divide->SetMarkerStyle(20);
  h4divide->SetMarkerStyle(20);
  h5divide->SetMarkerStyle(20);

  h3divide->SetMarkerColor(kRed+1);
  h4divide->SetMarkerColor(kBlue+1);
  h5divide->SetMarkerColor(kGreen+2);
  
  h3divide->SetLineColor(kRed+1);
  h4divide->SetLineColor(kBlue+1);
  h5divide->SetLineColor(kGreen+2);
  
  /*
  TF1 * func = new TF1("func", [=](double *x, double *par) {
      Double_t num = par[0];
      Double_t den = 1 + TMath::Exp(-par[1]*(x[0]-par[2]));
      Double_t val = num / den;
      return val;
      }, 3, 12, 3);
  func->SetParameters(1.0,2.0,3.5);
  //func->FixParameter(0,1);
  h4divide->Fit(func,"RM0");
  func->SetRange(0,15);
*/
  TLine * line = new TLine(0,1,32,1);
  //line->SetNDC(1);
  line->SetLineStyle(9);
  dummy->Draw();
  h3divide->Draw("e1 p same");
  h4divide->Draw("e1 p same");
  h5divide->Draw("e1 p same");
  //func->SetLineColor(kBlack);
  //func->Draw("same");
  line->Draw("same");
  
  TLegend * l2 = new TLegend(.15,.7,.75,.95);
  l2->SetLineWidth(0);
  l2->AddEntry(h3divide,"MBD + Photon 4 GeV");
  l2->AddEntry(h4divide,"MBD + Photon 5 GeV");
  l2->AddEntry(h5divide,"MBD + Photon 6 GeV");
  l2->Draw();
  
  c->cd();
  pcounts->Draw();
  peff->Draw();

  gStyle->SetPalette(kInvertedDarkBodyRadiator);
  gPad->SetLogz();


}
